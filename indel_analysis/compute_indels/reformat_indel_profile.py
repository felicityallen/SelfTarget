import io
import os
import sys
import csv
from Bio import SeqIO
import numpy as np
import re

from selftarget.indel import indelOutofGuideSeedPAM
from selftarget.oligo import getFileForOligoIdx
from selftarget.util import loadFastaReadsById

def addRead(profile, indel, seq, read_id, oligo_indel,mut):
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

def reformatIndelProfile( file_prefix, read_lookup ):
    read_profiles = {}
    if not os.path.isfile(file_prefix + '_mappedindels.txt'):
        print('Could not find file', file_prefix + '_mappedindels.txt')
    else:
        f = io.open(file_prefix + '_mappedindels.txt')
        rdr = csv.reader(f, delimiter='\t')
        for toks in rdr:
            if '@@@' in toks[0]: continue
            oligo_id = toks[0].split('.')[0]
            read_id = toks[0].split('.')[1]
            if oligo_id not in read_profiles:
                read_profiles[oligo_id] = {}
            seq = read_lookup[toks[0]]
            if toks[1] == '':
                continue #unmapped read?
            nulls = [(indel, null.split(':')[1], eval(null.split(':')[2]), muts) for (null,indel, muts) in zip(toks[1].split(','), toks[2].split(','), toks[3].split(','))]

            if len(nulls) == 1:
                indel, oligo_indel, perc_oligo_indel,muts = nulls[0]
                addRead( read_profiles[oligo_id], indel, seq, read_id, oligo_indel, muts )
            else:		
                #If there is a null oligo that results in a null indel, select that
                for (indel, oligo_indel, perc_oligo_indel, muts) in nulls:
                    if indel == '-':
                        addRead( read_profiles[oligo_id], indel, seq, read_id, oligo_indel, muts )
                        break
                if indel == '-':
                    continue

                #Otherwise, if 18 closest region in gRNA or PAM of oligo is affected, exclude from consideration,
                #otherwise assign randomly according to oligo percentages in plasmids
                ok_nulls = [x for x in nulls if indelOutofGuideSeedPAM(x[1])]  
                if len(ok_nulls) == 0:
                    ok_nulls = nulls
                percs = [x[2] for x in ok_nulls]
                total_perc = sum(percs)
                idx = np.where(np.random.multinomial(1,np.array(percs)/total_perc)>0)[0][0]
                selected_null = ok_nulls[idx]
                addRead( read_profiles[oligo_id], selected_null[0], seq, read_id, selected_null[1], selected_null[3] )
        f.close()

        fout = io.open(file_prefix + '_mappedindelprofiles.txt','w')
        fout_sum = io.open(file_prefix + '_mappedindelsummary.txt','w')
        for oligo_id in read_profiles:

            fout.write(u'@@@%s\n' % oligo_id)
            fout_sum.write(u'@@@%s\n' % oligo_id)
            indel_counts = [(countReads(read_profiles[oligo_id],x),x) for x in read_profiles[oligo_id]]
            indel_counts.sort(reverse=True)	
            if '-' in read_profiles[oligo_id]:
                writeReads(read_profiles[oligo_id], '-', fout, fout_sum)
            for (count, indel) in indel_counts:
                if indel == '-':
                    continue
                writeReads(read_profiles[oligo_id], indel, fout, fout_sum)
        fout.close()
        fout_sum.close()
            
if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: reformat_indel_profile.py <file_prefix>')
    else:
        file_prefix = sys.argv[1]
        
        read_lookup = loadFastaReadsById(file_prefix + '.fasta')
        reformatIndelProfile( file_prefix, read_lookup )
