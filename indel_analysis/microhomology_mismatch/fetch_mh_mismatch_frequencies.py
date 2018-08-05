import io, os, sys, csv

from selftarget.data import isOldLib, createResultDirectory, getDirLabel, getHighDataDir
from selftarget.oligo import getFileForOligoIdx, getOligoIdxFromId, getShortOligoId
from selftarget.profile import readSummaryToProfile

def fetchMhMismatchFrequencies(dirname, outdir='mh_mismatch_indel_frequencies'):

    if not os.path.isdir(outdir): os.makedirs(outdir)
    if isOldLib(dirname): raise Exception('Old Lib not supported')
        
    mh_exp_indels_file = getHighDataDir() + '/mh_mismatch_indels.txt'

    fout = io.open(outdir + '/' + getDirLabel(dirname) + '.txt','w')
    hdr_str = '\t'.join(['\t'.join([x + ' Indel Reads in ' + y for x in ['Orig', 'Left Mut', 'Right Mut', 'Merged Mut1', 'Merged Mut2']]) for y in ['Mut', 'Orig']])

    f = io.open(mh_exp_indels_file)
    rdr = csv.DictReader(f, delimiter='\t')
    fout.write(u'%s\t%s\tMut Non-Null Reads\tOrig Non-Null Reads\n' % ('\t'.join(rdr.fieldnames), hdr_str))
    for row in rdr:        

        #Load Indel Profiles for both the original and mutated micrhomology forms
        mut_oligo_id = row['Oligo ID'].replace('_','')
        orig_oligo_id = row['Mapped Oligo Id'].replace('_','')
        
        mut_filepath, mut_filename = getFileForOligoIdx(getOligoIdxFromId(mut_oligo_id), ext='_mappedindelsummary.txt')
        orig_filepath, orig_filename = getFileForOligoIdx(getOligoIdxFromId(orig_oligo_id), ext='_mappedindelsummary.txt')

        p_mut, p_orig = {},{}
        stats_mut = readSummaryToProfile(dirname + '/mapped_reads/' + mut_filepath + '/' + mut_filename, p_mut, oligoid=mut_oligo_id)
        stats_orig = readSummaryToProfile(dirname + '/mapped_reads/' + orig_filepath + '/' + orig_filename, p_orig, oligoid=orig_oligo_id)    
    
        indels = [row['Orig Indel'], row['Left Mut-MH Indel'], row['Right Mut-MH Indel'], row['Merge Mut 1 Indel'], row['Merge Mut 2 Indel']]
        reads = lambda indel, profile: profile[indel] if (indel in profile and indel != '') else 0
        mut_read_str = '\t'.join(['%d' % reads(indel, p_mut) for indel in indels])
        orig_read_str = '\t'.join(['%d' % reads(indel, p_orig) for indel in indels])

        str_args = ('\t'.join([row[col] for col in rdr.fieldnames]), mut_read_str, orig_read_str, stats_mut[0]-stats_mut[2], stats_orig[0]-stats_orig[2])
        fout.write(u'%s\t%s\t%s\t%d\t%d\n' % str_args)

    f.close()
    fout.close()
    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise Exception('Usage: fetch_mh_mismatch_frequencies.py dirname')

    dirname = sys.argv[1]    
    
    
    
    
    


