import io, os, csv, sys, re, random, itertools
import numpy as np
import pylab as PL
import Bio.Seq

from selftarget.profile import readSummaryToProfile, symmetricKL, classSymmetricKL, fetchIndelSizeCounts
from selftarget.oligo import getOligoIdxFromId, getFileForOligoIdx, loadOldNewMapping
from selftarget.data import getAllDataDirs,getSampleSelectors

MIN_READS = 100

updir = '/lustre/scratch117/cellgen/team227/fa9/self-target/indel_analysis'

#New Samples
new_dirs = [updir + x for x in ['/ST_June_2017/data/K562_800x_LV7A_DPI7/',
            '/ST_June_2017/data/K562_800x_LV7B_DPI7/',
            '/ST_Feb_2018/data/CAS9_12NA_1600X_DPI7/']]

#Old Samples
old_dirs = [updir + x for x in ['/ST_June_2017/data/K562_1600x_6OA_DPI7/',
            '/ST_April_2017/data/K562_800x_6OA_DPI7_Old8/',
            '/ST_April_2017/data/K562_800x_6OB_DPI7_Old11/']]

def getFileSuffix(oligo_id):
    oligo_idx = getOligoIdxFromId(oligo_id)
    subdir, sumfilename = getFileForOligoIdx(oligo_idx, ext='_mappedindelsummary.txt')
    return subdir + '/' + sumfilename

def loadProfilePair(old_id, new_id):
    p_old, p_new = {}, {}
    old_file, new_file = getFileSuffix(old_id), getFileSuffix(new_id)
    mut_reads_old, mut_reads_new = 0, 0
    for new_dir in new_dirs:
        acc, pacc, null = readSummaryToProfile(new_dir + 'mapped_reads/' + new_file, p_new, oligoid=new_id )
        mut_reads_new += (acc-null)
    for old_dir in old_dirs:
        acc, pacc, null = readSummaryToProfile(old_dir + 'mapped_reads/' + old_file, p_old, oligoid=old_id )
        mut_reads_old += (acc-null)
    return p_old, p_new, mut_reads_old, mut_reads_new

def loadSeparateProfilePairs(old_id, new_id):
    old_ps, new_ps = [], []
    old_file, new_file = getFileSuffix(old_id), getFileSuffix(new_id)
    for new_dir in new_dirs:
        p_new, mut_reads_new = {}, 0
        acc, pacc, null = readSummaryToProfile(new_dir + 'mapped_reads/' + new_file, p_new, oligoid=new_id )
        mut_reads_new += (acc-null)
        new_ps.append(p_new)
    for old_dir in old_dirs:
        p_old, mut_reads_old = {}, 0
        acc, pacc, null = readSummaryToProfile(old_dir + 'mapped_reads/' + old_file, p_old, oligoid=old_id )
        mut_reads_old += (acc-null)
        old_ps.append(p_old)
    return old_ps, new_ps

def meanSymKL(profs, kl_func = symmetricKL):
    sum_tot, n_tot = 0.0, 0
    for p1, p2 in itertools.combinations(profs, 2):
        sum_tot += kl_func(p1, p2)
        n_tot += 1
    return sum_tot/n_tot

if __name__ == '__main__':

    data = loadOldNewMapping()
    fout = io.open('old_new_kl_summaries.txt', 'w')

    fout.write(u'Old Oligo Id\tNew Oligo Id\tOld Mut Reads\tNew Mut Reads\tOld In Frame Perc\tNew In Frame Perc')
    fout.write(u'\tOld v New KL\tOld v Old KL\tNew v New KL')
    fout.write(u'\tOld v New Class KL\tOld v Old Class KL\tNew v New Class KL')
    fout.write(u'\tAlt Old Id\tAlt New Id\tAlt Old Mut Reads\tAlt New Mut Reads')
    fout.write(u'\tAlt Old v New KL\tAlt Old v Old KL\tAlt New v New KL')
    fout.write(u'\tAlt Old v New Class KL\tAlt Old v Old Class KL\tAlt New v New Class KL')
    fout.write(u'\tAlt2 Old Id\tAlt2 New Id\tAlt2 Old Mut Reads\tAlt2 New Mut Reads')
    fout.write(u'\tAlt2 Old v New KL\tAlt2 Old v Old KL\tAlt2 New v New KL')
    fout.write(u'\tAlt2 Old v New Class KL\tAlt2 Old v Old Class KL\tAlt2 New v New Class KL\n')

    id_pairs_all = [y for y in zip([x for x in data['Old']], [x for x in data['New']])]
    id_pairs = []
    for (old_id, new_id) in id_pairs_all:
        p_old, p_new, mut_reads_old, mut_reads_new = loadProfilePair(old_id, new_id)
        print(old_id, new_id, mut_reads_old, mut_reads_new)
        if mut_reads_old < MIN_READS or mut_reads_new < MIN_READS: continue
        id_pairs.append((old_id, new_id))

    alt_id_pairs = [x for x in id_pairs]; random.shuffle(alt_id_pairs)
    alt2_id_pairs = [x for x in id_pairs]; random.shuffle(alt2_id_pairs)

    print(len(id_pairs),len(alt_id_pairs),len(alt2_id_pairs))
    for (old_id, new_id), (alt_old_id, alt_new_id), (alt2_old_id, alt2_new_id) in zip(id_pairs, alt_id_pairs, alt2_id_pairs):
        print(old_id, new_id)
        p_old, p_new, mut_reads_old, mut_reads_new = loadProfilePair(old_id, new_id)
        old_ps, new_ps = loadSeparateProfilePairs(old_id, new_id)
        p_old_alt, p_new_alt, mut_reads_old_alt, mut_reads_new_alt= loadProfilePair(alt_old_id, alt_new_id)
        p_old_alt2, p_new_alt2, mut_reads_old_alt2, mut_reads_new_alt2= loadProfilePair(alt2_old_id, alt2_new_id)
        
        old_if, old_of, _ = fetchIndelSizeCounts(p_old)
        new_if, new_of, _ = fetchIndelSizeCounts(p_new)
        old_if_perc = old_if*100.0/(old_if+old_of)
        new_if_perc = new_if*100.0/(new_if+new_of)
 
        out_str = '' 
        for kl_func in [symmetricKL, classSymmetricKL]:
            str_args = (kl_func(p_old, p_new), meanSymKL(old_ps, kl_func=kl_func), meanSymKL(new_ps, kl_func=kl_func))
            out_str += u'\t%.5f\t%.5f\t%.5f' % str_args
        out_str += u'\t%s\t%s\t%d\t%d' % (alt_old_id, alt_new_id, mut_reads_old_alt, mut_reads_new_alt)
        for kl_func in [symmetricKL, classSymmetricKL]:
            str_args = (np.mean([kl_func(p_old, p_new_alt),kl_func(p_new, p_old_alt)]), kl_func(p_old, p_old_alt), kl_func(p_new, p_new_alt))
            out_str += u'\t%.5f\t%.5f\t%.5f' % str_args
        out_str += u'\t%s\t%s\t%d\t%d' % (alt2_old_id, alt2_new_id, mut_reads_old_alt2, mut_reads_new_alt2)
        for kl_func in [symmetricKL, classSymmetricKL]:
            str_args = (np.mean([kl_func(p_old, p_new_alt2),kl_func(p_new, p_old_alt2)]), kl_func(p_old, p_old_alt2), kl_func(p_new, p_new_alt2))
            out_str += u'\t%.5f\t%.5f\t%.5f' % str_args
        fout.write(u'%s\t%s\t%d\t%d\t%.3f\t%.3f\t%s\n' % (old_id, new_id, mut_reads_old, mut_reads_new, old_if_perc, new_if_perc, out_str))
        fout.flush()
    fout.close()
        
        





        


