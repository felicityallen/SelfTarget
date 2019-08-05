import io, sys, os, csv
import Bio.Seq

from selftarget.data import getIndelSummaryFiles, getDirLabel, createResultDirectory, getExpOligoFile, setHighDataDir
from selftarget.indel import tokFullIndel
from selftarget.profile import readSummaryToProfile, getHighestIndel, getProfileCounts
from selftarget.oligo import getOligoIdsFromFile, getShortOligoId, loadAllOligoDetails, loadExpOligoLookup

def getSequence(oligo_det, left, right):
    pam_loc, pam_dir, seq = oligo_det
    pam_loc = eval(pam_loc)
    if pam_dir not in ['REVERSE', 'FORWARD']:
        raise Exception('Unexpected PAM direction: ' + str(pam_dir))
    if pam_dir == 'REVERSE':
        seq = Bio.Seq.reverse_complement(seq)
        pam_loc = len(seq) - pam_loc
    return seq[pam_loc-3+left:pam_loc-3+right+1]

def writeMCISummary(fout, id, p1, stats1, oligo_det, more_indels=False):
    if not more_indels: mcis = [getHighestIndel(p1)]
    else: mcis = [x[1] for x in getProfileCounts(p1) if x[1] != '-']
    for mci in mcis:
        mci_reads = p1[mci]
        total_reads = stats1[0]-stats1[2]
        itype, isize, details, muts = tokFullIndel(mci)
        pam_loc, pam_dir, seq = oligo_det

        mh_seq, altered_seq = '', ''
        if itype == 'D' and ('I' not in details or details['I']==0):
            if details['C'] > 0:
                left_c_seq = getSequence(oligo_det, details['L']+1, details['L'] + details['C'])
                right_c_seq = getSequence(oligo_det, details['R']-details['C'], details['R']-1)
                if left_c_seq == right_c_seq:
                    mh_seq = left_c_seq
            altered_seq = getSequence(oligo_det, details['L'] +1, details['R']-1)  #Note includes MH seq at both ends

        str_args = (id, mci, details['L'], details['R'],details['C'],itype,isize,mci_reads,total_reads,mh_seq,altered_seq)
        fout.write(u'%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%s\n' % str_args)

if __name__ == '__main__':

    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print('Usage: compile_most_common_indels.py <high_dir> <subdir> <more_indels>')
    else:

        #Set up output file
        high_dir = sys.argv[1]
        if high_dir != '.': setHighDataDir(high_dir)
        subdir = sys.argv[2]
        more_indels = True
        if len(sys.argv) > 3: more_indels = eval(sys.argv[3])
    
        if more_indels: out_dir = createResultDirectory(high_dir + '/more_indel_summaries',subdir)
        else: out_dir = createResultDirectory(high_dir + '/most_common_indel_summaries', subdir)
        fout = io.open(out_dir + '/' + subdir.split('/')[-1] + '.txt', 'w')
        oligo_lookup = loadExpOligoLookup(subdir)

        #For each Oligo, summarise details of its most common indel
        fout.write(u'Oligo Id\tMost Common Indel\tLeft\tRight\tCentral\tType\tSize\tMCI Reads\tTotal reads\tMicrohomology Sequence\tAltered Sequence\n')
        sum_files = getIndelSummaryFiles(subdir)
        for filename in sum_files:
            file_prefix = filename.split('/')[-1][:-23]
            oligo_details = {x[0]: x[1:] for x in oligo_lookup[file_prefix]}
            oligo_ids = getOligoIdsFromFile(filename)
            for id in oligo_ids:
    
                #Read in the profile (if it exists)	
                p1 = {}
                stats1 = readSummaryToProfile(filename, p1, oligoid=getShortOligoId(id))

                if len(p1) == 0 or p1.keys() == ['-']:
                    continue

                #Compute and summarise its MCI details
                writeMCISummary(fout, id, p1, stats1, oligo_details[id], more_indels)

        fout.close()

            

