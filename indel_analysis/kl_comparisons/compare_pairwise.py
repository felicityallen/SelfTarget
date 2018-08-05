import io, sys, os, csv
import numpy as np

from selftarget.profile import readSummaryToProfile, symmetricKL, compareTopIndels, entropy, percentOverlap
from selftarget.oligo import getOligoIdsFromFile
from selftarget.data import getDirLabel, getIndelSummaryFiles, getSubdirs

def filterLargeI(profile):
    return {x:profile[x] for x in profile if (x[0] == '-' or x[0] != 'I' or x[1] != '1')}

if len(sys.argv) != 4:
    print('compare_pairwise.py <dirname1> <dirname2> <subdir>')
else:

    remove_largeI = True

    dirname1, dirname2 = sys.argv[1], sys.argv[2]
    subdir = sys.argv[3]

    if subdir not in getSubdirs(dirname1, withpath=False):
        raise Exception('No subdir %s in %s' % (subdir, dirname1) )
    if subdir not in getSubdirs(dirname2, withpath=False):
        raise Exception('No subdir %s in %s' % (subdir, dirname2) )    
        
    out_dir = 'profile_comparison_summaries' if not remove_largeI else 'profile_comparison_summaries_nolargeI'
    if not os.path.isdir(out_dir): os.mkdir(out_dir)
    out_dir += '/%s_vs_%s' % (getDirLabel(dirname1),getDirLabel(dirname2))
    if not os.path.isdir(out_dir): os.mkdir(out_dir)
    out_file = out_dir + '/%s.txt' % subdir
    
    fout = io.open(out_file,'w')
    fout.write(u'ID\tNum Reads 1\tNum Reads 2\tNum States 1\tNum States 2\tNum null reads 1\tNum null reads 2\tKL with Null\tKL without null\tPerc Accepted Reads 1\tPerc Accepted Reads 2\t1st Nonmatch Indel\tNum Top 3 Common\tNum Top 5 Common\tNum Top 10 Common\tProfile 1 Entropy (before mods)\tProfile 2 Entropy (before mods)\tProfile 1 Entropy (after mods)\tProfile 2 Entropy (after mods)\tPerc Overlap\tP1 Perc in Top 3\tP2 Perc in Top 3\tP1 Perc in Top 5\tP2 Perc in Top 5\tP1 Perc in Top 10\tP2 Perc in Top 10\n')
    
    dir1_files = getIndelSummaryFiles(dirname1 + '/mapped_reads/' + subdir, withpath=False)
    dir2_files = getIndelSummaryFiles(dirname2 + '/mapped_reads/' + subdir, withpath=False)
    common_files =  set(dir1_files).intersection(set(dir2_files))
    for filename in common_files:
        
        filename1 = dirname1 + '/mapped_reads/' + subdir + '/' + filename
        filename2 = dirname2 + '/mapped_reads/' + subdir + '/' + filename
        
        oligo_ids1 = getOligoIdsFromFile( filename1 )
        oligo_ids2 = getOligoIdsFromFile( filename2 )
        common_oligos =  set(oligo_ids1).intersection(set(oligo_ids2))
        for oligo_id in common_oligos:

            profile1, profile2 = {}, {}
            num_reads1, perc_acc1, nonull1 = readSummaryToProfile(filename1, profile1, oligoid=oligo_id)	
            num_reads2, perc_acc2, nonull2 = readSummaryToProfile(filename2, profile2, oligoid=oligo_id)
            ns1, ns2 = len(profile1), len(profile2)

            if remove_largeI:
                profile1 = filterLargeI(profile1)
                profile2 = filterLargeI(profile2)

            ent1a, ent2a = entropy(profile1,True), entropy(profile2,True)
            poverlap = percentOverlap( profile1, profile2, True )

            score1 = symmetricKL( profile1, profile2, False )
            score2 = symmetricKL( profile1, profile2, True )

            ent1b, ent2b = entropy(profile1,True), entropy(profile2,True)	#Since comparing the profiles appends missing states to both profiles

            nonmatch_idx, top_common, top_percs = compareTopIndels(profile1, profile2)
            fout.write(u'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (oligo_id, num_reads1, num_reads2, ns1, ns2, nonull1, nonull2, score1, score2, perc_acc1, perc_acc2, nonmatch_idx, top_common[3], top_common[5], top_common[10], ent1a, ent2a, ent1b, ent2b, poverlap, top_percs[3][0],top_percs[3][1],top_percs[5][0],top_percs[5][1], top_percs[10][0], top_percs[10][1]))	
    
    fout.close()


