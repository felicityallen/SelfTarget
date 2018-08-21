import pylab as PL
PL.switch_backend('Agg')
import io, os, csv, sys, re
import numpy as np
import pylab as PL
import Bio.Seq
from selftarget.util import getPlotDir
from selftarget.indel import tokFullIndel
from selftarget.profile import readSummaryToProfile, fetchRepresentativeCleanReads, getProfileCounts
from selftarget.oligo import getFileForOligoIdx
from selftarget.plot import saveFig

from matplotlib import rcParams
rcParams['font.family'] = 'monospace'
rcParams['font.monospace'] = ['DejaVu Sans Mono'] + rcParams['font.monospace']

FORECAST_GREEN = '#%02x%02x%02x' % (0, 103, 111)
FORECAST_ORANGE = '#%02x%02x%02x' % (244, 87, 33)
MIN_X = -30
MAX_X = 30

def padReadForIndel(read_seq, indel, pam_idx):
    itype,isize,details,muts = tokFullIndel(indel)
    red_idxs, green_idxs = set(), set()
    if itype == 'D':
        read_seq = read_seq[:pam_idx-3+details['L']+details['C']+1] + ' '*isize + read_seq[pam_idx-3+details['L']+details['C']+1:]
        green_idxs = set(range(pam_idx-3+details['L']+1, pam_idx-3+details['L']+1+details['C']))
    if itype == 'I':
        green_idxs = set(range(pam_idx-3+details['L']+1, pam_idx-3+details['L']+1+details['C']))
        red_idxs = set(range(pam_idx-3+details['L']+1+details['C'],pam_idx-3+details['L']+details['C']+1+isize))
    return read_seq, red_idxs, green_idxs
    
def plotSeqLetterwise(seq, y, pam_idx, red_idxs=set(), green_idxs=set(), default_clr='black'):
    min_xloc,max_xloc  = 0, 0
    for i,nt in enumerate(seq):
        clr = FORECAST_ORANGE if i in red_idxs else FORECAST_GREEN if i in green_idxs else default_clr
        weight = 'normal' if clr == default_clr else 'bold'
        xloc = i-pam_idx+3
        if xloc > -1:    xloc += 0.1
        if xloc > MIN_X and xloc <= MAX_X:
            PL.text(xloc,y, nt, verticalalignment='bottom', horizontalalignment='left', color=clr, fontweight=weight)
            if xloc < min_xloc:
                min_xloc = xloc
            if xloc > max_xloc:
                max_xloc = xloc
    return min_xloc,max_xloc

def plotProfiles(profiles, rep_reads, pam_idxs, reverses, labels, title='', max_lines=10):
    if len(profiles) == 0: raise Exception('Empty list of profiles')
    
    colors = [FORECAST_GREEN,'C0', 'C2', 'C2','C1','C1','C3','C3','C4','C4','C5','C5','C6']

    PL.rcParams['svg.fonttype'] = 'none'
    ocounts = [getProfileCounts(p1) for p1 in profiles]
    counts = [{indel: (cnt,indel,perc1a,perc1b) for (cnt,indel,perc1a,perc1b) in x} for x in ocounts]
    
    #Count total non-null reads for each sample (to report in labels)
    nonnull_reads = [sum([x[indel][0] for indel in x if indel != '-']) for x in counts]
    labels = ['%s(%d Reads)' % (tit,nn) for (tit,nn) in zip(labels, nonnull_reads)]

    #Fetch the indels to display as union of top N indels across profiles
    num_top = 20
    top_indels = [[y[1] for y in x[:num_top]] for x in ocounts]
    union_top_indels = set()
    for x in top_indels: union_top_indels = union_top_indels.union(set(x))
    
    for indel in union_top_indels:
        for count in counts:
            if indel not in count:
                count[indel] = (0,indel,0.0,0.0)
    union_top_indels = [x for x in union_top_indels]
    indel_toks = [tokFullIndel(indel) for indel in union_top_indels]
    max_insert = max([0] + [toks[1] for toks in indel_toks if toks[0] == 'I'])

    #Order indels by decreasing average percentage across profiles
    top_av_percs = [(np.mean([x[indel][-1] for x in counts]),indel) for indel in union_top_indels]
    top_av_percs.sort(reverse=True)
    max_indels = max_lines/len(profiles)

    #Figure out Trims
    null_reads = [x['-'] if '-' in x else [x[y[1]] for y in ocnt if y[1] in x][0] for x,ocnt in zip(rep_reads, ocounts)]
    null_reads = [Bio.Seq.reverse_complement(x) if rev else x for x,rev in zip(null_reads, reverses)]
    pam_idxs = [len(x) - pam if rev else pam for x,pam,rev in zip(null_reads, pam_idxs, reverses)]
    min_null, pam_idx = min([(len(null),pidx) for (null, pidx) in zip(null_reads, pam_idxs)])
    Ls = [x - pam_idx for x in pam_idxs]
    Rs = [L + min_null - len(null) for (L,null) in zip(Ls, null_reads)]

    #Plot
    scale_factor = 10.0/max([x[1][3] for x in ocounts])
    fig = PL.figure(figsize=(9,5*len(labels)))
    fig.patch.set_visible(False)
    ax = PL.gca()
    ax.axis('off')
    N = min(len(union_top_indels), max_indels)
    line_height = 0.8
    min_xloc,max_xloc = MIN_X, MAX_X
    PL.ylim( (0,(N+1.0)*line_height) )
    bar_ypos, bar_len = [[] for x in profiles], [[] for x in profiles]
    for i, (av_perc, indel) in enumerate(top_av_percs):
        if i > max_indels: break
        for repr, cnts, rev, L1, R1, j in zip(rep_reads, counts, reverses, Ls, Rs, range(len(Rs))):
            (cnt1,indel1,perc1a,perc1b) = cnts[indel] 
            if indel in repr:
                if R1 == 0: R1 = len(repr[indel])
                seq = Bio.Seq.reverse_complement(repr[indel])[L1:R1] if rev else repr[indel][L1:R1]
                padded_seq, red_idxs, green_idxs = padReadForIndel(seq, indel, pam_idx)
                min_xloc,max_xloc = plotSeqLetterwise(padded_seq, (N-i+(j+0.3)*1.0/len(profiles))*line_height, pam_idx, red_idxs=red_idxs, green_idxs=green_idxs)
            if indel != '-':
                bar_ypos[j].append((N-i+(j+0.4)*1.0/len(profiles))*line_height)
                bar_len[j].append(perc1b*scale_factor)
    hist_loc = max_xloc + 10
    for bar1_ypos, bar1_len, label1,clr in zip(bar_ypos, bar_len, labels,colors):
        PL.barh(bar1_ypos, bar1_len, height=0.8*line_height/len(profiles), left=hist_loc, label=label1, color=clr )
        for (ypos, blen) in zip(bar1_ypos, bar1_len):
            PL.text(hist_loc+blen+1,ypos-0.5/len(profiles)*line_height, '%.1f%%' % (blen/scale_factor))
    xlims = (min_xloc-10,MAX_X+20+(min_xloc-MIN_X))
    PL.xlim( xlims )
    for i, (av_perc, indel) in enumerate(top_av_percs):
        if i > max_indels: break
        if indel == '-':
            PL.text(xlims[0],(N-i+0.4)*line_height,'Target:', fontweight='bold')
        else:
            PL.text(xlims[0],(N-i+0.4)*line_height,indel.split('_')[0], fontweight='bold')
        PL.plot([min_xloc-10,max_xloc + 10],[(N-i)*line_height, (N-i)*line_height],'lightgrey')
    PL.plot([0,0],[0,(N+1)*line_height],'k--')
    PL.plot([min_xloc-10,hist_loc],[N*line_height,N*line_height],'k')
    PL.plot([hist_loc,hist_loc],[0,N*line_height],'k')
    PL.xticks([])
    PL.yticks([])
    if len(labels) > 1: PL.legend(loc='upper right')
    PL.text(hist_loc,(N+0.5)*line_height, title, fontweight='bold')
    PL.subplots_adjust(left=0.05,right=0.95,top=0.95, bottom=0.05)
    PL.show(block=False)
    PL.axis('off')
    saveFig('%s_%d' % (title.replace(' ','_'), len(labels)), bbox=False)
    return fig
