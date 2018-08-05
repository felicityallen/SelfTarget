import io
import os
import sys
import csv
from Bio import SeqIO
import numpy as np
import re

from selftarget.util import loadFastaReadsById
from selftarget.indel import tokFullIndel
from selftarget.oligo import loadPamLookup, loadExpOligoLookup
from selftarget.data import getPamLocFile, getExpOligoFile

def shortenLookupIds(lookup):
    new_lookup = {}
    for x in lookup:
        new_lookup[x.split('_')[0]] = lookup[x]
    return new_lookup

def addNullRead(profile, indel, seq, read_id, oligo_indel,mut):
    if indel not in profile:
        profile[indel] = {}
    if oligo_indel not in profile[indel]:
        profile[indel][oligo_indel] = []
    profile[indel][oligo_indel].append((seq, read_id, mut))

def countReads(profile, indel):
    return sum([len(profile[indel][oligo_indel]) for oligo_indel in profile[indel]])

def writeReads(profile, indel, fout, fout_sum):
    oligo_counts = [(len(profile[indel][x]),x) for x in profile[indel]]
    oligo_counts.sort(reverse=True)
    for (count, oligo_indel) in oligo_counts:
        fout_sum.write(u'%s\t%s\t%d\n' % (indel, oligo_indel, count))
        for (seq, read_id, muts) in profile[indel][oligo_indel]:
            fout.write(u'%s\t%s\t%s\t%s\t%s\n' % (seq, indel, read_id, oligo_indel, muts))

def updatePam(indel, orig_pam_loc, pam_dir):
    pam_loc = orig_pam_loc
    itype, isize, details, muts = tokFullIndel(indel)
    if itype != '-':
        if pam_dir == 'REVERSE':
            left_pos = pam_loc + 2 - (details['R']-1) + details['C']
            right_pos = pam_loc + 2 - (details['L']+1) + details['C']
        else:
            left_pos = pam_loc - 3 + (details['L']+1) + details['C']
            right_pos = pam_loc - 3 + (details['R']-1) + details['C']

        if itype == 'D':
            delsize = isize - details['I']
        else:
            delsize = -isize + details['D']

        if left_pos < pam_loc:
            pam_loc = max(pam_loc-delsize, left_pos)
    
    for (muttype,mutpos,nucl) in muts:
        if muttype == 'D':
            msize = mutpos
        if muttype == 'I':
            msize = -mutpos
        if muttype != 'S':
            continue
        if pam_dir == 'REVERSE':
            mutidx = pam_loc + 2 - mutpos
        else:
            mutidx = pam_loc - 3 + mutpos
        if mutidx < pam_loc:
            pam_loc = pam_loc-msize
        
    return pam_loc

def compileMappedNull(file_prefix, read_lookup, pam_lookup, exp_oligo_lookup):
    read_profiles, indel_seqs = {}, {}
    if not os.path.isfile(file_prefix + '_mappings.txt'):
        print('Could not find file', file_prefix + '_mappings.txt')
    else:

        #Add 5 pseudo reads for the NULL indel for all oligos (in case poorly represented in the NULL measure)
        if file_prefix.split('/')[-1] in exp_oligo_lookup:
            for (oligo_id, pam_loc, pam_dir, seq) in exp_oligo_lookup[file_prefix.split('/')[-1]]:
                read_profiles[oligo_id] = {'-': 5}
                indel_seqs[oligo_id] = {'-': seq}

        f = io.open(file_prefix + '_mappings.txt')
        rdr = csv.reader(f, delimiter='\t')
        for toks in rdr:
            oligo_id = toks[1].split('_')[0]
            read_id = oligo_id + '.' + toks[0].split()[0]
            if oligo_id not in read_profiles:
                read_profiles[oligo_id] = {}
                indel_seqs[oligo_id] = {}
            seq = read_lookup[read_id]
            indel = toks[2]+'_'+toks[3] #combine mutations with indels
            itype, isize, details, muts = tokFullIndel(indel)
            if indel == '-_-':
                indel = '-'
            if indel not in read_profiles[oligo_id]:
                read_profiles[oligo_id][indel]=0
                indel_seqs[oligo_id][indel] = seq
            read_profiles[oligo_id][indel] += 1
        f.close()

        fout = io.open(file_prefix + '_nullsummary.txt','w')
        oligo_ids = [x for x in read_profiles.keys()]
        oligo_ids.sort()
        for oligo_id in oligo_ids:

            orig_pam_loc, pam_dir = pam_lookup[oligo_id]
            fout.write(u'@@@%s\n' % oligo_id)
            indel_counts = [(read_profiles[oligo_id][x],x) for x in read_profiles[oligo_id]]
            indel_counts.sort(reverse=True)
            total_counts = sum([x[0] for x in indel_counts])	
            for (count, indel) in indel_counts:
                seq = indel_seqs[oligo_id][indel]
                perc = count*100.0/total_counts
                pam_loc = updatePam(indel, orig_pam_loc, pam_dir)
                fout.write(u'%s\t%s\t%d\t%s\t%.3f\n' % (seq,indel,pam_loc,pam_dir,perc))
        
        fout.close()

def convertToExpFile(null_sum_file,  output_file, discard_long=True):
    f = io.open(null_sum_file)
    fout = io.open(output_file,'w')
    for line in f:
        if line[:3] == '@@@':
            oligo_id = line[3:-1]
            continue
        toks = line.split('\t')
        seq, indel, pam_idx, pam_dir, perc = toks
        if discard_long and len(seq) > 120: 
            print('Too long:', oligo_id, seq)
            continue     #Leave out long templates
        if eval(perc) < 0.5:
            continue
        fout.write('>%s:%s:%s %s %s\n%s\n' % (oligo_id, indel, perc[:-1], pam_idx, pam_dir, seq))	
    fout.close()
        
if __name__ == '__main__':

    if len(sys.argv) != 5:
        print('Usage: compile_mapped_null_profile.py <subdir> <alt_pam> <pam_file> <exp_oligo_file>')
    else:
        subdir = sys.argv[1]
        alt_pam = eval(sys.argv[2])
        pam_adj = 40 if alt_pam else 20
        pam_file = sys.argv[3]
        exp_oligo_file = sys.argv[4]
        
        fasta_files = [x for x in os.listdir(subdir) if x[-6:] == '.fasta']
        for fasta_file in fasta_files:

            file_prefix = subdir + '/' + fasta_file[:-6]
            read_lookup = loadFastaReadsById(file_prefix + '.fasta')
            pam_lookup = shortenLookupIds(loadPamLookup(pam_file,adj=pam_adj))
            exp_oligo_lookup = loadExpOligoLookup(subdir, exp_oligo_file=exp_oligo_file)
            compileMappedNull(file_prefix, read_lookup, pam_lookup, exp_oligo_lookup)
            convertToExpFile(file_prefix + '_nullsummary.txt', file_prefix + '_exptargets.txt')

