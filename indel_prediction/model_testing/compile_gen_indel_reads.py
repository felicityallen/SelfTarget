import io, sys, os, csv
import Bio.Seq

from selftarget.data import getAllDataDirs, isOldLib, getDirLabel, getHighDataDir
from selftarget.profile import readSummaryToProfile
from selftarget.oligo import getOligoIdxFromId, getFileForOligoIdx

def compileGenIndelReads(gen_indel_dir='generated_indels', out_dir = 'reads_for_gen_indels_all_samples', sample_dirs=[]):

    if not os.path.isdir(out_dir): os.mkdir(out_dir)

    for gen_file in os.listdir(gen_indel_dir):

        oligo_id = gen_file.split('_')[0]
        oligo_idx = getOligoIdxFromId(oligo_id)
        oligo_subdir, sum_filename = getFileForOligoIdx(oligo_idx, ext='_mappedindelsummary.txt')

        out_subdir = out_dir + '/' + oligo_subdir
        if not os.path.isdir(out_subdir): os.mkdir(out_subdir)

        #Read all profiles for this oligo
        profiles, mut_read_totals = [], []
        for dirname in sample_dirs:
            profiles.append({})
            filename = getHighDataDir() + '/' + dirname + '/mapped_reads/' + oligo_subdir + '/' + sum_filename
            stats = readSummaryToProfile(filename, profiles[-1], oligoid=oligo_id)
            mut_read_totals.append('%d' % (stats[0]-stats[2]))

        #Compile reads for each indel across all samples
        f = io.open(gen_indel_dir + '/' + gen_file)
        fout = io.open(out_subdir + '/%s_gen_indel_reads.txt' % oligo_id, 'w')
        fout.write(f.readline())    #Git commit
        fout.write(u'Indel\tDetails\t%s\n' % '\t'.join([getDirLabel(x) for x in sample_dirs]))
        fout.write(u'All Mutated\t[]\t%s\n' % '\t'.join(mut_read_totals))
        for toks in csv.reader(f,delimiter='\t'):
            indel, indel_details = toks[0], toks[2]
            read_str = '\t'.join(['%d' % (p1[indel] if indel in p1 else 0) for p1 in profiles])
            fout.write(u'%s\t%s\t%s\n' % (indel, indel_details, read_str))
        fout.close()
        f.close()

if __name__ == '__main__':
    sample_dirs = [ '/ST_June_2017/data/K562_800x_LV7A_DPI7/',
                    '/ST_June_2017/data/K562_800x_LV7B_DPI7/',
                    '/ST_Feb_2018/data/CAS9_12NA_1600X_DPI7/',
                    '/ST_June_2017/data/K562_1600x_6OA_DPI7/',
                    '/ST_April_2017/data/K562_800x_6OA_DPI7_Old8/',
                    '/ST_April_2017/data/K562_800x_6OB_DPI7_Old11/']
    compileGenIndelReads(sample_dirs=sample_dirs)


            

