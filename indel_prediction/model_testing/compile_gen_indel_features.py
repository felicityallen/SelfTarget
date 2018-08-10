import io, sys, os, csv
import Bio.Seq

from predictor.features import calculateFeaturesForGenIndelFile

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
        generated_indel_file = gen_indel_dir + '/' + gen_file    
        out_file = out_subdir + '/%s_gen_indel_features.txt' % oligo_id
        is_reverse = (row['PAM Direction'] == 'REVERSE')
        calculateFeaturesForGenIndelFile( generated_indel_file, uncut_seq, cut_site, out_file, is_reverse=is_reverse)

if __name__ == '__main__':

    setHighDataDir('predicted_vs_measured_example')
    computeFeaturesForGenIndels()