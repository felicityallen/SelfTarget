import io, sys, os, csv
import Bio.Seq

from selftarget.data import isOldLib, getDirNameFromSubdir, getIndelSummaryFiles, getDirLabel, createResultDirectory, getExpOligoFile
from selftarget.indel import tokFullIndel
from selftarget.profile import readSummaryToProfile, getHighestIndel, getProfileCounts, fetchRepresentativeCleanReads
from selftarget.oligo import getOligoIdsFromFile, getShortOligoId, loadAllOligoDetails, loadExpOligoLookup


def loadI1to3Details(filename):
    f = io.open(filename)
    f.readline()    #Git commit line
    rdr = csv.DictReader(f, delimiter='\t')
    i1_details = {getShortOligoId(row['Oligo Id']): row for row in rdr}
    f.close
    return i1_details, rdr.fieldnames[1:-2]

def writeI1to3Summary(fout, id, p1, stats1, i1_details_row, indel_hdrs):

    i1_indels = [i1_details_row[x] for x in indel_hdrs]
    i1_indel_map = {1:{},2:{},3:{}}
    for indel, hdr in zip(i1_indels, indel_hdrs):
        ilen = len(hdr.split('_')[-1])
        if indel not in i1_indel_map[ilen]:
            i1_indel_map[ilen][indel] = 0
        i1_indel_map[ilen][indel] += 1
        
    read_count = lambda indel: p1[indel] if indel in p1 else 0 
    rpt_left = i1_details_row['Repeat Nucleotide Left']
    rpt_right = i1_details_row['Repeat Nucleotide Right']
    fout.write(u'%s' % id )
    for ilen in [1,2,3]:
        left_indel = i1_details_row['I1_%s' % (rpt_left*ilen)]
        right_indel = i1_details_row['I1_%s' % (rpt_right*ilen)]
        left_cnt = read_count(left_indel)/(1.0+(left_indel==right_indel))
        right_cnt = read_count(right_indel)/(1.0+(left_indel==right_indel))
        other_cnt = sum([read_count(x) for x in i1_indel_map[ilen]]) - left_cnt - right_cnt
        fout.write(u'\t%d\t%d\t%d' % (left_cnt, right_cnt, other_cnt))
    fout.write(u'\t%d\n' % (stats1[0]-stats1[2]))


if __name__ == '__main__':

    if len(sys.argv) != 3:
        raise Exception('Usage: compile_i1.py <highdir> <subdir>')

    #Set up output file
    highdir = sys.argv[1]
    subdir = sys.argv[2]
    out_dir = createResultDirectory(highdir + '/i1_summaries', subdir)
    fout = io.open(out_dir + '/' + subdir.split('/')[-1] + '.txt', 'w')

    #Load I1 details file
    dirname = getDirNameFromSubdir(subdir)
    if isOldLib(dirname): i1_details_file = 'exp_target_pam_old_gen_i1_indels.txt'
    else: i1_details_file = 'exp_target_pam_new_gen_i1_indels.txt'
    
    i1_loc = '.' if highdir == '.' else highdir + '/ST_June_2017/data'
    i1_details, i1_hdrs = loadI1to3Details(i1_loc + '/' + i1_details_file)

    #For each Oligo, summarise details of I1 insertions
    fout.write(u'Oligo Id')
    for i in [1,2,3]:
        fout.write(u'\tI%d_Rpt Left Reads\tI%d_Rpt Right Reads\tI%d_NonRpt Reads' % (i,i,i))
    fout.write(u'\tTotal reads\t\n')
    sum_files = getIndelSummaryFiles(subdir)
    for filename in sum_files:
        file_prefix = filename.split('/')[-1][:-23]
        oligo_ids = getOligoIdsFromFile(filename)
        for id in oligo_ids:
    
            #Read in the profile (if it exists)	
            p1 = {}
            stats1 = readSummaryToProfile(filename, p1, oligoid=getShortOligoId(id))

            if len(p1) == 0 or p1.keys() == ['-']:
                continue

            #Compute and summarise its i1
            writeI1to3Summary(fout, id, p1, stats1, i1_details[id], i1_hdrs)

    fout.close()

            

