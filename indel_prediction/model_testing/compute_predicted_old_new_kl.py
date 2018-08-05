import io, os, csv, sys, re, random, itertools
import numpy as np
import pylab as PL
import Bio.Seq
import pandas as pd

from selftarget.profile import fetchRepresentativeCleanReads,readSummaryToProfile, symmetricKL, fetchIndelSizeCounts, getSummaryFileSuffix
from selftarget.oligo import getOligoIdxFromId, getFileForOligoIdx, loadOldNewMapping
from selftarget.data import getHighDataDir, setHighDataDir
from selftarget.view import plotProfiles
from selftarget.util import setPlotDir, getIndelGenExe, setIndelGenExe

from predictor.model import computePredictedProfile, readTheta, setFeaturesDir, setReadsDir

from compile_gen_indel_features import computeFeaturesForGenIndels
from compile_gen_indel_reads import compileGenIndelReads

#New Samples
new_dirs = ['/ST_June_2017/data/K562_800x_LV7A_DPI7',
            '/ST_June_2017/data/K562_800x_LV7B_DPI7',
            '/ST_Feb_2018/data/CAS9_12NA_1600X_DPI7']

#Old Samples
old_dirs = ['/ST_June_2017/data/K562_1600x_6OA_DPI7',
            '/ST_April_2017/data/K562_800x_6OA_DPI7_Old8',
            '/ST_April_2017/data/K562_800x_6OB_DPI7_Old11']

def loadProfilePair(old_id, new_id):
    p_old, p_new = {}, {}
    old_file, new_file = getSummaryFileSuffix(old_id), getSummaryFileSuffix(new_id)
    mut_reads_old, mut_reads_new = 0, 0
    for new_dir in [getHighDataDir() + '/' + x for x in new_dirs]:
        acc, pacc, null = readSummaryToProfile(new_dir + '/mapped_reads/' + new_file, p_new, oligoid=new_id )
        mut_reads_new += (acc-null)
    for old_dir in [getHighDataDir() + '/' + x for x in old_dirs]:
        acc, pacc, null = readSummaryToProfile(old_dir + '/mapped_reads/' + old_file, p_old, oligoid=old_id )
        mut_reads_old += (acc-null)
    return p_old, p_new, mut_reads_old, mut_reads_new

def loadRepReads(new_id):
    oligo_idx = getOligoIdxFromId(new_id)
    subdir, profilefilename = getFileForOligoIdx(oligo_idx, ext='_mappedindelprofiles.txt')
    profile_file = getHighDataDir() + '/' + new_dirs[0] + '/mapped_reads/' + subdir + '/' + profilefilename
    rep_reads = {}
    fetchRepresentativeCleanReads(profile_file, rep_reads, oligoid=new_id) 
    return rep_reads

def combineProfiles(p_old, p_new, mut_reads_old, mut_reads_new):
    p_combined = {indel: count for indel,count in p_old.items()}
    for indel, count in p_new.items():
        if indel not in p_combined: p_combined[indel] = 0
        p_combined[indel] += count
    return p_combined, mut_reads_old + mut_reads_new

def loadValidationPairs():
    f = io.open(getHighDataDir() + '/old_new_validation_guides.txt')
    id_pairs = [[row['Old Oligo Id'],row['New Oligo Id']] for row in csv.DictReader(f, delimiter='\t')]
    f.close()
    return id_pairs

def computeAndComparePredicted(theta_file, selected_id=None, out_dir='.'):

    features_dir = getHighDataDir() + '/gen_indels/features_for_gen_indels'
    theta, train_set, feature_columns = readTheta(theta_file)

    #Note: here old refers to conventional scaffold library, new refers to improved scaffold library
    fout = io.open(out_dir + '/old_new_kl_predicted_summaries.txt', 'w')
    fout.write(u'Old Oligo Id\tNew Oligo Id\tOld Mut Reads\tNew Mut Reads\tCombined Mut Reads\tOld In Frame Perc\tNew In Frame Perc\tCombined in Frame Perc\tPredicted In Frame Per')
    fout.write(u'\tOld v New KL\tOld v Predicted KL\tNew v Predicted KL\tCombined v Predicted KL\n')

    id_pairs = loadValidationPairs()
    for (old_id, new_id) in id_pairs:
        if old_id in train_set or new_id in train_set:
            raise Exception('Bad!!! Testing on Training data: %s %s' % (old_id, new_id))

        if selected_id is not None and selected_id != old_id:
            continue #Guide pair selected for plotting

        #Load Old and new profiles, and produce combined profile from the two
        p_old, p_new, mut_reads_old, mut_reads_new = loadProfilePair(old_id, new_id)
        p_comb, mut_reads_comb = combineProfiles(p_old, p_new, mut_reads_old, mut_reads_new)

        #Predict the profile (old and new will be the same so just do one)
        p_predict, _ = computePredictedProfile(new_id, theta, feature_columns)

        #Compute in frame percentages
        old_if, old_of, _ = fetchIndelSizeCounts(p_old)
        new_if, new_of, _ = fetchIndelSizeCounts(p_new)
        comb_if, comb_of, _ = fetchIndelSizeCounts(p_comb)
        pred_if, pred_of, _ = fetchIndelSizeCounts(p_predict)
        old_if_perc = old_if*100.0/(old_if+old_of)
        new_if_perc = new_if*100.0/(new_if+new_of)
        comb_if_perc = comb_if*100.0/(comb_if+comb_of)
        pred_if_perc = pred_if*100.0/(pred_if+pred_of)

        #Plot the comparison
        if selected_id is not None:
            rrds = loadRepReads(new_id)
            plotProfiles([p_old,p_new,p_predict], [rrds,rrds,rrds], [42,42,42], [False,False,False], ['Replicate 1','Replicate 2','Predicted'], title='%s (KL=%.2f, KL=%.2f)' % (new_id, symmetricKL(p_old, p_new), symmetricKL(p_comb,p_predict)))

        str_args = (symmetricKL(p_old, p_new), symmetricKL(p_old, p_predict), symmetricKL(p_new, p_predict), symmetricKL(p_comb, p_predict) )
        kl_str = u'\t%.5f\t%.5f\t%.5f\t%.5f' % str_args
        fout.write(u'%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f%s\n' % (old_id, new_id, mut_reads_old, mut_reads_new,mut_reads_comb, old_if_perc, new_if_perc, comb_if_perc, pred_if_perc, kl_str))
        fout.flush()
    fout.close()
        
def prepareExample(example_dir):

    setHighDataDir(example_dir)

    #Generate all possible indels
    new_gen_dir, old_gen_dir = getHighDataDir() + '/generated_indels_new', getHighDataDir() + '/generated_indels_old'
    if not os.path.isdir(new_gen_dir): os.makedirs(new_gen_dir)
    if not os.path.isdir(old_gen_dir): os.makedirs(old_gen_dir)
    cmd = getIndelGenExe() + ' ' + getHighDataDir() + '/exp_target_pam_new.fasta ' + new_gen_dir  +'/'
    print(cmd); os.system(cmd)
    cmd = getIndelGenExe() + ' ' + getHighDataDir() + '/exp_target_pam_old.fasta ' + old_gen_dir  +'/'
    print(cmd); os.system(cmd)
    
    #Compile number of reads per sample for each indel
    reads_dir = getHighDataDir() + '/reads_for_gen_indels'
    compileGenIndelReads(gen_indel_dir=new_gen_dir, out_dir = reads_dir, sample_dirs=new_dirs)
    compileGenIndelReads(gen_indel_dir=old_gen_dir, out_dir = reads_dir, sample_dirs=old_dirs)
    setReadsDir(reads_dir)
    
    #Compute features for each indel
    features_dir = getHighDataDir() + '/features_for_gen_indels'
    computeFeaturesForGenIndels(gen_indel_dir = new_gen_dir, out_dir=features_dir)
    computeFeaturesForGenIndels(gen_indel_dir = old_gen_dir, out_dir=features_dir)
    setFeaturesDir(features_dir)

if __name__ == '__main__':
 
    setIndelGenExe('C:/Users/fa9/postdoc/indelmap/build/Release/indelgen.exe')
    setPlotDir('plots')
    prepareExample('predicted_vs_measured_example')

    #Predict mutations using pre-trained model and compare to actual (for one oligo only)
    theta_file = getHighDataDir() + '/model_output_2000_0.01000000_1.835_theta.txt_cf0.txt'
    computeAndComparePredicted(theta_file, selected_id='Oligo69919', out_dir='.')


        


