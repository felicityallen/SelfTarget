import io, os, csv, sys, re, random, itertools, shutil
import numpy as np
import pylab as PL
import Bio.Seq
import pandas as pd

from selftarget.profile import fetchRepresentativeCleanReads,readSummaryToProfile, symmetricKL, fetchIndelSizeCounts, getSummaryFileSuffix
from selftarget.oligo import getOligoIdxFromId, getFileForOligoIdx, loadOldNewMapping
from selftarget.data import getHighDataDir, setHighDataDir
from selftarget.view import plotProfiles
from selftarget.util import setPlotDir, getIndelGenExe, setIndelGenExe
from selftarget.plot import setFigType

from predictor.model import computePredictedProfile, readTheta, setFeaturesDir, setReadsDir, loadOligoFeaturesAndReadCounts

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

def loadProfilesSeparately(old_id, new_id):

    p_olds, p_news, old_sep_mr, new_sep_mr = [{},{}], [{},{}], [0,0],[0,0]
    old_file, new_file = getSummaryFileSuffix(old_id), getSummaryFileSuffix(new_id)
    for new_dir in [getHighDataDir() + '/' + x for x in new_dirs]:
        idx = 0 if '800' in new_dir else 1
        acc, pacc, null = readSummaryToProfile(new_dir + '/mapped_reads/' + new_file, p_news[idx], oligoid=new_id )
        new_sep_mr[idx] += acc-null
    for old_dir in [getHighDataDir() + '/' + x for x in old_dirs]:
        idx = 0 if '800' in old_dir else 1
        acc, pacc, null = readSummaryToProfile(old_dir + '/mapped_reads/' + old_file, p_olds[idx], oligoid=old_id )
        old_sep_mr[idx] += acc-null
    return p_olds, p_news, old_sep_mr, new_sep_mr

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

def getInFramePerc(profile):
    p_if, p_of, _ = fetchIndelSizeCounts(profile)
    p_if_perc = p_if*100.0/(p_if+p_of)
    return p_if_perc

def getCombs(obj):
    return [(x,y) for (x,y) in itertools.combinations(obj,2)]
    

def computeAndComparePredicted(theta_file, selected_id=None, out_dir='.', start_count=0, end_count=10000):

    features_dir = getHighDataDir() + '/gen_indels/features_for_gen_indels'
    theta, train_set, feature_columns = readTheta(theta_file)

    new_sep_labels = 'New 2x800x', 'New 1600x'
    old_sep_labels = 'Old 2x800x', 'Old 1600x'

    #Note: here old refers to conventional scaffold library, new refers to improved scaffold library
    fout = io.open(out_dir + '/old_new_kl_predicted_summaries.txt' % (start_count, end_count), 'w')
    fout.write(u'Old Oligo Id\tNew Oligo Id\tOld Mut Reads\tNew Mut Reads\tCombined Mut Reads\t')
    fout.write(u'\t'.join('%s Mut Reads' % x.split('/')[-1] for x in new_sep_labels + old_sep_labels))
    fout.write(u'\tOld In Frame Perc\tNew In Frame Perc\tCombined in Frame Perc\tPredicted In Frame Per\t')
    fout.write(u'\t'.join('%s In Frame Perc' % x.split('/')[-1] for x in new_sep_labels + old_sep_labels))
    fout.write(u'\tOld v New KL\tOld v Predicted KL\tNew v Predicted KL\tCombined v Predicted KL\t')
    fout.write(u'\t'.join('%s vs Predicted KL' % x.split('/')[-1] for x in new_sep_labels + old_sep_labels) + '\t')
    fout.write(u'\t'.join(['%s vs %s KL' % (x.split('/')[-1],y.split('/')[-1]) for x,y in (getCombs(new_sep_labels) + getCombs(old_sep_labels))]) + '\n')

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
        feature_data = loadOligoFeaturesAndReadCounts(new_id, [])
        p_predict, _ = computePredictedProfile(feature_data, theta, feature_columns)

        #Load separate profiles too
        p_old_sep, p_new_sep, old_sep_mr, new_sep_mr =  loadProfilesSeparately(old_id, new_id)

        #Compute in frame percentages
        old_if_perc = getInFramePerc(p_old)
        new_if_perc = getInFramePerc(p_new)
        comb_if_perc = getInFramePerc(p_comb)
        pred_if_perc = getInFramePerc(p_predict)
        new_sep_if_percs = [getInFramePerc(profile) if len(profile)>1 else -1 for profile in p_new_sep]
        old_sep_if_percs = [getInFramePerc(profile) if len(profile)>1 else -1 for profile in p_old_sep]

        #Plot the comparison
        if selected_id is not None:
            rrds = loadRepReads(new_id)
            plotProfiles([p_new_sep[0],p_new_sep[1],p_predict], [rrds,rrds,rrds], [56,56,56], [False,False,False], ['Replicate 1','Replicate 2','Predicted'], title='%s (KL=%.2f, KL=%.2f)' % (new_id, symmetricKL(p_new_sep[0], p_new_sep[1]), symmetricKL(p_new,p_predict)))

        str_args = (symmetricKL(p_old, p_new), symmetricKL(p_old, p_predict), symmetricKL(p_new, p_predict), symmetricKL(p_comb, p_predict) )
        kl_str = u'\t%.5f\t%.5f\t%.5f\t%.5f\t' % str_args
        kl_str += u'\t'.join(['%.5f' % symmetricKL(p_predict,x) for x in p_new_sep + p_old_sep])
        kl_str += u'\t' + u'\t'.join(['%.5f' % symmetricKL(x,y) for (x,y) in (getCombs(p_new_sep) + getCombs(p_old_sep))])
        if_str = u'\t'.join(['%.3f' % x for x in new_sep_if_percs + old_sep_if_percs])
        mut_str = u'\t'.join(['%d' % x for x in new_sep_mr + old_sep_mr])
        fout.write(u'%s\t%s\t%d\t%d\t%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s%s\n' % (old_id, new_id, mut_reads_old, mut_reads_new,mut_reads_comb,mut_str, old_if_perc, new_if_perc, comb_if_perc, pred_if_perc, if_str, kl_str))
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
 
    setIndelGenExe('/usr/local/bin/indelgen')
    setPlotDir('/results/plots')
    setFigType('png')
	    
    shutil.copytree('/data/predicted_vs_measured_example', '/results/predicted_vs_measured_example')
    prepareExample('/results/predicted_vs_measured_example')

    #Predict mutations using pre-trained model and compare to actual (for one oligo only)
    theta_file = getHighDataDir() + '/model_output_10000_0.01000000_0.01000000_-0.607_theta.txt_cf0.txt'
    computeAndComparePredicted(theta_file, selected_id='Oligo35785', out_dir='.')


        


