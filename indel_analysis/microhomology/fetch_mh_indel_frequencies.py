import io, os, sys, csv

from selftarget.data import isOldLib, createResultDirectory, getDirNameFromSubdir, getIndelSummaryFiles
from selftarget.oligo import getFileForOligoIdx, getOligoIdsFromFile, getShortOligoId
from selftarget.profile import readSummaryToProfile

def loadMhExpIndels(filename, oligo_ids):
    
    mh_indels = {}
    f = io.open(filename)
    for toks in csv.reader(f,delimiter='\t'):
        
        oligo_id = getShortOligoId(toks[0])
        if oligo_id not in oligo_ids: continue
        mhs = toks[1][:-1].split(',')
        indels = toks[2][:-1].split(',') 
        mh_indels[oligo_id] = ( mhs,indels )
        
    f.close()
    return mh_indels

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('Usage: fetch_mh_indel_frequencies.py highdir subdir')

    else:

        highdir = sys.argv[1]
        subdir = sys.argv[2]
        outdir = createResultDirectory(highdir + '/mh_indel_frequencies', subdir, with_subdir=True)
        dirname = getDirNameFromSubdir(subdir)
        if isOldLib(dirname):
            mh_exp_indels_file = 'exp_target_old_mh_indels.txt'
        else:
            mh_exp_indels_file = 'exp_target_new_mh_indels.txt'

        indel_files = getIndelSummaryFiles(subdir, withpath=False)
        for indel_file in indel_files:
    
            oligo_ids = getOligoIdsFromFile(subdir + '/' + indel_file)
            mh_loc = '.' if highdir == '.' else highdir + '/ST_June_2017/data'
            mh_indels = loadMhExpIndels(mh_loc + '/' + mh_exp_indels_file, set(oligo_ids))
    
            fout = io.open(outdir + '/' + indel_file[:-23] + '_mhindels.txt','w')
            for oligo_id in oligo_ids:
            
                profile = {}
                acc, pacc, nullr = readSummaryToProfile(subdir + '/' + indel_file, profile, oligoid=oligo_id)
            
                fout.write(u'@@@%s:%d:%d\n' % (oligo_id,acc,acc-nullr))
                mhs,indels = mh_indels[oligo_id]
            
                for (mh, indel) in zip(mhs, indels):
                    left, right, mh_len = mh.split(':')
                    if indel == 'Unmappable': continue
                    if indel in profile: nreads = profile[indel]
                    else: nreads = 0
                    fout.write(u'%s\t%s\t%s\t%s\t%d\n' % (left, right, mh_len, indel, nreads))
            fout.close()
                
    
    
    
    
    
    


