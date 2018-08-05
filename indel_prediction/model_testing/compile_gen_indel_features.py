import io, sys, os, csv
import Bio.Seq

from predictor.features import calculateFeatures

from selftarget.oligo import loadAllOligoDetails, getShortOligoId, loadPamLookup, getOligoIdxFromId, getFileForOligoIdx
from selftarget.data import getHighDataDir

def computeFeaturesForGenIndels(gen_indel_dir = 'generated_indels', out_dir='features_for_gen_indels'):

    if not os.path.isdir(out_dir): os.mkdir(out_dir)

    #Load Oligo details
    oligo_details = loadAllOligoDetails(oligo_detail_dir=getHighDataDir() + '/ST_June_2017/data')
    oligo_details = {id.replace('_',''): row for (id,row) in oligo_details.items()}

    for gen_file in os.listdir(gen_indel_dir):
        print(gen_file)

        oligo_id = gen_file.split('_')[0]
        oligo_idx = getOligoIdxFromId(oligo_id)
        oligo_subdir, _ = getFileForOligoIdx(oligo_idx, ext='')

        out_subdir = out_dir + '/' + oligo_subdir
        if not os.path.isdir(out_subdir): os.mkdir(out_subdir)

        row = oligo_details[oligo_id]

        uncut_seq = row['Target'] if row['PAM Direction'] != 'REVERSE' else Bio.Seq.reverse_complement(row['Target'])
        cut_site = eval(row['PAM Location'])-3 if row['PAM Direction'] != 'REVERSE' else (79 - eval(row['PAM Location']) - 3)
            
        f = io.open(gen_indel_dir + '/' + gen_file )
        fout = io.open(out_subdir + '/%s_gen_indel_features.txt' % oligo_id, 'w')
        fout.write(f.readline())    #Git commit line (pass on)
        fout.write(u'###%s\t%d\t%s\n' % (uncut_seq, cut_site, row['PAM Direction']))
        first = True
        A,T,G,C = 'A','T','G','C'
        AA,AT,AC,AG,CG,CT,CA,CC = 'AA','AT','AC','AG','CG','CT','CA','CC'
        GT,GA,GG,GC,TA,TG,TC,TT = 'GT','GA','GG','GC','TA','TG','TC','TT'

        for toks in csv.reader(f,delimiter='\t'):
            indel, indel_locs = toks[0], eval(toks[2])
            for indel_loc in indel_locs:
                ins_seq = indel_loc[2] if len(indel_loc) > 2 else ''
                left = indel_loc[0] if row['PAM Direction'] != 'REVERSE' else (78 - indel_loc[1]) 
                right = indel_loc[1] if row['PAM Direction'] != 'REVERSE' else (78 - indel_loc[0]) 
                ins_seq = ins_seq if row['PAM Direction'] != 'REVERSE' else Bio.Seq.reverse_complement(ins_seq)
                indel_details = (uncut_seq, cut_site, left, right, ins_seq)
             
                features, feature_labels = calculateFeatures(indel_details)
                
                if first: fout.write(u'Indel\tLeft\tRight\tInserted Seq\t%s\n' % '\t'.join(feature_labels))
                feature_str = '\t'.join(['%d' % x for x in features])
                fout.write(u'%s\t%d\t%d\t%s\t%s\n' % (indel, left, right, ins_seq, feature_str))
                first = False
        fout.close()
        f.close()

if __name__ == '__main__':

    setHighDataDir('predicted_vs_measured_example')
    computeFeaturesForGenIndels()