import io, os, csv, sys, re
import numpy as np

from selftarget.indel import tokFullIndel
from selftarget.data import getWTDir
from selftarget.oligo import getSummaryFileSuffix

#KL Divergence between the two indel profiles (non-symmetric)
def KL(p1, p2, ignore_null=True, missing_count=0.5):
    
    p1_indels = set([x for x in p1 if p1[x]>0 and (x != '-' or not ignore_null)])
    p2_indels = set([x for x in p2 if p2[x]>0 and (x != '-' or not ignore_null)])
    common = p1_indels.intersection(p2_indels)
    p1_only = p1_indels.difference(p2_indels)
    p2_only = p2_indels.difference(p1_indels)
    
    p1_total = sum([p1[x] for x in p1_indels]) + missing_count*len(p2_only)
    p2_total = sum([p2[x] for x in p2_indels]) + missing_count*len(p1_only)
    
    if p1_total > 0 and p2_total > 0:
        norm1, norm2 = 1.0/p1_total, 1.0/p2_total
        score = 0.0
        for indel in common:
            score += p1[indel]*norm1*np.log2(p1[indel]*norm1/(p2[indel]*norm2))
        for indel in p1_only:
            score += p1[indel]*norm1*np.log2(p1[indel]*norm1/(missing_count*norm2))
        for indel in p2_only:
            score += missing_count*norm1*np.log2(missing_count*norm1/(p2[indel]*norm2))
    else: score = np.nan
    return score

#Percent Overlap between the two indel profiles
def percentOverlap( p1, p2, ignore_null ):
    score=0.0
    if ignore_null:
        if len(p1) == 1 and '-' in p1 or len([x for x in p2 if x != '-']) == 0:
            norm1, norm2 = 0.0,0.0
        else:
            norm1 = 100.0/sum([p1[x] for x in p1 if x != '-'])
            norm2 = 100.0/sum([p2[x] for x in p2 if x != '-'])
    else:
        norm1 = 100.0/sum([p1[x] for x in p1])
        norm2 = 100.0/sum([p2[x] for x in p2])
    for indel in p1:
        if ignore_null and indel=='-':
            continue
        if indel in p2:
            score += min(p1[indel]*norm1,p2[indel]*norm2)
    return score

#Entropy of an indel profile
def entropy(p1, ignore_null):
    score = 0.0
    if len(p1) == 0 or sum([p1[x] for x in p1]) == 0:
        return score
    if ignore_null:
        if len(p1) == 1 and '-' in p1:
            norm1 = 0.0
    if ignore_null:
        if len(p1) == 1 and '-' in p1:
            norm1 = 0.0
        else:
            norm1 = 1.0/sum([p1[x] for x in p1 if x != '-'])
    else:
        norm1 = 1.0/sum([p1[x] for x in p1])
    for indel in p1:
        if ignore_null and indel=='-':
            continue
        score += -p1[indel]*norm1*np.log2(p1[indel]*norm1)
    return score

def isAllowableOligoIndel(oligo_indel):
    itype, isize, details, muts = tokFullIndel(oligo_indel)
    #Exclude reads from oligos with any mutations in the guide or PAM sequence
    is_ok = True
    mut_locs = [x for x in muts if x[0] not in ['N','I','D']]
    if len(mut_locs) > 0: 
        if any([x[1] > -20 and x[1] < 6 for x in mut_locs]):
            is_ok = False
        if len(mut_locs) > 5:
            is_ok = False
    #Only allow oligo indels if they're size 1 or 2 insertion/deletions outside the guide or PAM sequence
    ins_del_muts = [x for x in muts if x[0] in ['I','D']]
    if len(ins_del_muts) > 0:
        if any([x[1] > 2 for x in ins_del_muts]):
            is_ok = False
    if oligo_indel[0] != '-':
        if isize > 2 or (details['L'] < 6 and details['R'] > -20):   
            is_ok = False
    return is_ok
    
#Read in profile from indel summary file
def readSummaryToProfile(filename, profile, oligoid=None, noexclude=False, remove_long_indels=False, remove_wt=True, wt_thresh=3.0):

    if not os.path.isfile(filename): return 0,0,0
    
    dirname = '/'.join(filename.split('/')[:-3])
    filename_suffix = '/'.join(filename.split('/')[-3:])
    wt_p, wt_p_wfilter = {}, {}
    if 'WT' not in dirname and dirname != '' and not noexclude and remove_wt:
        wt_filename = getWTDir(dirname) + '/' + filename_suffix
        #if wt_filename[0] == '/' and wt_filename[1:7] != 'lustre': wt_filename = wt_filename[1:]
        if not os.path.isfile(wt_filename):
            print('Warning: Could not find', wt_filename)
        else:
            readSummaryToProfile(wt_filename, wt_p, oligoid=oligoid, noexclude=True, remove_wt=False)
            _, wt_acc, _ = readSummaryToProfile(wt_filename, wt_p_wfilter, oligoid=oligoid, noexclude=False, remove_wt=False)
            if wt_acc < 10.0: return 0,0,0    #Need at least 20% acceptable reads in the wild type
                                          #(to remove oligos that are really messed up)

    total, accepted = 0,0
    f = io.open(filename)
    reader = csv.reader(f, delimiter='\t')
    if '-' not in profile:
        profile['-'] = 0
    orig_null = profile['-']
    curr_oligo_id = None
    wt_indels = []
    for toks in reader:
        if toks[0][:3] == '@@@':
            curr_oligo_id = toks[0][3:].split()[0]
            continue
        if oligoid != curr_oligo_id:
            continue
        indel = toks[0]
        oligo_indel = toks[1]
        num_reads = eval(toks[2])
        total += num_reads
        if not noexclude:
            if oligo_indel != '-':
                if not isAllowableOligoIndel(oligo_indel):
                    continue
            #Only allow indels that span the cut site and which are
            #not present in the corresponding WT sample
            if indel != '-':
                itype, isize, details, muts = tokFullIndel(indel)
                if itype != '-' and (details['L'] > 5 or details['R'] < -5):
                    continue
                if remove_long_indels and isize > 30:
                    continue
                if indel in wt_p and remove_wt: 
                    #Check the levels of the indel in the WT sample,
                    #only include it if present at at least 3 x that level (including NULLS)
                    # - will need to wait til we know total reads to do this
                    wt_indels.append((indel, num_reads))
                    continue
        if indel not in profile:
            profile[indel] = 0
        profile[indel] += num_reads
        accepted += num_reads
    for indel, num_reads in wt_indels:
        if num_reads*1.0/total > wt_p[indel]*wt_thresh/sum([wt_p[x] for x in wt_p]):
            if indel not in profile: profile[indel] = 0
            profile[indel] += num_reads
            accepted += num_reads
    f.close()
    if total == 0:
        perc_accepted = 0.0
    else:
        perc_accepted = accepted*100.0/total
    return accepted, perc_accepted, profile['-']-orig_null
    
def loadMergedProfile(oligo_id, sample_dirs=[]):
    profile, sumfilename = {}, getSummaryFileSuffix(oligo_id)
    mut_reads = 0
    for sample_dir in sample_dirs:
        acc, pacc, null = readSummaryToProfile(sample_dir + 'mapped_reads/' + sumfilename, profile, oligoid=oligo_id )
        mut_reads += (acc-null)
    return profile, mut_reads

def readNullSummaryToProfile(filename, profile, oligoid=None, unedited_only=True):
    total, accepted = 0,0
    f = io.open(filename)
    reader = csv.reader(f, delimiter='\t')
    curr_oligo_id = None
    for toks in reader:
        if toks[0][:3] == '@@@':
            curr_oligo_id = toks[0][3:]
            continue
        if oligoid != curr_oligo_id:
            continue
        indel = toks[0]
        if unedited_only and indel != '-':
            continue
        oligo_indel = toks[1]
        num_reads = eval(toks[2])
        if oligo_indel not in profile:
            profile[oligo_indel] = 0
        profile[oligo_indel] += num_reads
    f.close()

def getProfileCounts(profile):
    total = sum([profile[x] for x in profile])
    if total == 0:
        return []
    indel_total = total
    if '-' in profile:
        indel_total -= profile['-']
        null_perc = profile['-']*100.0/indel_total if indel_total != 0 else 100.0
        null_profile = (profile['-'],'-',profile['-']*100.0/total, null_perc)
    counts = [(profile[x],x, profile[x]*100.0/total, profile[x]*100.0/indel_total) for x in profile if x != '-']
    counts.sort(reverse = True)
    if '-' in profile:
        counts = [null_profile] + counts
    return counts

def printCompareProfiles( profile1, profile2 ):
    counts1 = getProfileCounts(profile1)
    counts2 = getProfileCounts(profile2)
    i = 0
    while i < len(counts1) and i < len(counts2) and i < 50:
        (cnt1,indel1,perc1a,perc1b) = counts1[i] 
        (cnt2,indel2,perc2a,perc2b) = counts2[i] 
        print('%s\t%.1f\t%d\t\t%s\t%.1f\t%d' % (indel1,perc1b,cnt1, indel2,perc2b,cnt2)) 
        i += 1
    print('...%d more\t\t\t\t...%d more' % (len(counts1)-i, len(counts2)-i))
    
def symmetricKL( profile1, profile2, ignore_null=True ):
    return 0.5*KL(profile1, profile2, ignore_null) + 0.5*KL(profile2, profile1, ignore_null)
    
def compareTopIndels( p1, p2 ):
    cnt1 = [[p1[x],x] for x in p1 if x != '-']
    cnt2 = [[p2[x],x] for x in p2 if x != '-']
    cnt1.sort(reverse=True); cnt2.sort(reverse=True)
    idx, nonmatch_idx = 0,-1
    top_indels_1 = {3:[[],0.0],5:[[],0.0],10:[[],0.0]}
    top_indels_2 = {3:[[],0.0],5:[[],0.0],10:[[],0.0]}
    while idx < len(cnt1) and idx < len(cnt2):
        i1, i2 = cnt1[idx][1], cnt2[idx][1]
        if i1 != i2 and nonmatch_idx < 0:
            nonmatch_idx = idx
        for thresh in top_indels_1:
            if idx < thresh:
                top_indels_1[thresh][0].append(i1)
                top_indels_1[thresh][1] += cnt1[idx][0]
                top_indels_2[thresh][0].append(i2)
                top_indels_2[thresh][1] += cnt2[idx][0]
        if idx >= max(top_indels_1.keys()) and nonmatch_idx >= 0:
            break
        idx += 1
    top_common, top_percs = {3:-1,5:-1,10:-1}, {3:(-1,-1),5:(-1,-1),10:(-1,-1)}
    if len(cnt1) > 0 and len(cnt2) > 0:
        for thresh in top_indels_1:
            top_common[thresh] = len(set(top_indels_1[thresh][0]).intersection(set(top_indels_2[thresh][0])))
            top_percs[thresh] = (top_indels_1[thresh][1]*100.0/sum([x[0] for x in cnt1]),top_indels_2[thresh][1]*100.0/sum([x[0] for x in cnt2]))
    return nonmatch_idx+1, top_common, top_percs

def getHighestIndel(p1):
    cnts = [(p1[x],x) for x in p1]
    cnts.sort(reverse=True)
    if (len(p1) == 1 and '-' in p1) or len(p1) == 0:
        return '-'
    if cnts[0][1] == '-' and len(cnts)>1:
        return cnts[1][1]
    else:
        return cnts[0][1]
        
def fetchRepresentativeCleanReads(mapped_profile_filename, rep_reads, oligoid=None):
    curr_oligo_id = None
    f = io.open(mapped_profile_filename)
    used_null_oligo = 'SOME'
    for toks in csv.reader(f, delimiter='\t'):
        if toks[0][:3] == '@@@':
            curr_oligo_id = toks[0][3:]
            continue
        if oligoid != curr_oligo_id:
            continue
        indel = toks[1]
        if indel not in rep_reads:
            rep_reads[indel] = toks[0]
            used_null_oligo = toks[-2]
        elif toks[-2] == '-' and used_null_oligo != '-':
            rep_reads[indel] = toks[0]
            used_null_oligo = toks[-2]
        elif toks[-2] == '-' and toks[-1] == '-':
            rep_reads[indel] = toks[0]
            used_null_oligo = toks[-2]
            
    f.close()
    
def fetchIndelSizeCounts(p1):
    inframe, outframe, size_counts, = 0,0,{'I':{},'D':{}}
    for i in range(1,21):
        size_counts['I'][i] = 0
        size_counts['D'][i] = 0
    for indel in p1:
        if indel == '-':
            continue
        itype,isize,details, muts = tokFullIndel(indel)
        net_isize = isize - details['I'] - details['D']
        if net_isize % 3 == 0:
            inframe += p1[indel]
        else:
            outframe += p1[indel]
        if net_isize not in size_counts[itype]:
            size_counts[itype][net_isize] = 0
        size_counts[itype][net_isize] += p1[indel]
    return inframe, outframe, size_counts

def makeClassProfile(p1):
    p1_class = {}
    for indel in p1:
        if indel == '-': continue
        indel_class = indel.split('_')[0]
        if indel_class not in p1_class:
            p1_class[indel_class] = 0
        p1_class[indel_class] += p1[indel]
    return p1_class

def classSymmetricKL(p1, p2):
    p1_class, p2_class = makeClassProfile(p1), makeClassProfile(p2)
    return symmetricKL(p1_class, p2_class)

def limProfile(p1, N):
    cnts = [[p1[x],x] for x in p1 if x != '-']
    cnts.sort(reverse=True)
    return {x[1]: x[0] for x in cnts[:N]}

def symmetricClassKLTopNIndels(p1, p2, N=10):
    p1_lim, p2_lim = limProfile(p1,N), limProfile(p2,N)
    return classSymmetricKL(p1_lim, p2_lim)

def symmetricClassKLTop5Indels(p1, p2):
    return symmetricClassKLTopNIndels(p1, p2, N=5)
