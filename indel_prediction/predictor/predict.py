import io, os, csv, sys, re, random, itertools, subprocess
import numpy as np
import pylab as PL
import Bio.Seq
import pandas as pd

from selftarget.profile import fetchRepresentativeCleanReads,readSummaryToProfile, symmetricKL, fetchIndelSizeCounts, getSummaryFileSuffix
from selftarget.oligo import getOligoIdxFromId, getFileForOligoIdx, loadOldNewMapping
from selftarget.data import getHighDataDir, setHighDataDir
from selftarget.view import plotProfiles
from selftarget.plot import setFigType
from selftarget.util import setPlotDir, getIndelGenExe, setIndelGenExe

from predictor.model import computePredictedProfile, readTheta, setFeaturesDir, setReadsDir
from predictor.features import calculateFeaturesForGenIndelFile, readFeaturesData

INDELGENTARGET_EXE = os.getenv("INDELGENTARGET_EXE", "C:/Users/fa9/postdoc/indelmap/build/Release/indelgentarget.exe")


def setIndelGenTargetExeLoc(val):
    global INDELGENTARGET_EXE
    INDELGENTARGET_EXE = val

def fetchRepReads(genindels_file):
    f = io.open(genindels_file)
    rep_reads = {toks[0]:toks[-1] for toks in csv.reader(f, delimiter='\t')}
    f.close()
    return rep_reads

def predictMutations(theta_file, target_seq, pam_idx):

    theta, train_set, theta_feature_columns = readTheta(theta_file)

    #generate indels
    tmp_genindels_file = 'tmp_genindels_%s_%d.txt' % (target_seq, random.randint(0,100000))
    cmd = INDELGENTARGET_EXE + ' %s %d %s' % (target_seq, pam_idx, tmp_genindels_file)
    print(cmd); subprocess.check_call(cmd.split())
    rep_reads = fetchRepReads(tmp_genindels_file)

    #compute features for all generated indels
    tmp_features_file = 'tmp_features_%s_%d.txt' % (target_seq, random.randint(0,100000))
    calculateFeaturesForGenIndelFile( tmp_genindels_file, target_seq, pam_idx-3, tmp_features_file)
    os.remove(tmp_genindels_file)
    feature_data, feature_columns = readFeaturesData(tmp_features_file)
    os.remove(tmp_features_file)

    if len(set(theta_feature_columns).union(set(feature_columns))) != len(theta_feature_columns):
        raise Exception('Stored feature names associated with model thetas do not match those computed')

    #Predict the profile
    p_predict, _ = computePredictedProfile(feature_data, theta, theta_feature_columns)
    return p_predict, rep_reads


def plot_predictions(theta_file, target_seq, pam_idx, title='Title'):

    if pam_idx < 0 or pam_idx >= (len(target_seq)-3):
        raise Exception('PAM idx out of range')

    if sum([x in ['A','T','G','C'] for x in target_seq]) != len(target_seq):
        raise Exception('Sequence must be composed of A,T,G,or C only')

    if len(target_seq) < 20 or pam_idx < 13 or pam_idx > len(target_seq)-7:
        raise Exception('Sequence too short or PAM too close to edge of sequence (must have at least 10nt either side of cut site)')

    if target_seq[pam_idx+1:pam_idx+3] != 'GG':
        raise Exception('Non NGG PAM (check correct index of PAM)')


    profile, rep_reads = predictMutations(theta_file, target_seq, pam_idx)
    profile['-'] = 1000
    rep_reads['-'] = target_seq
    setFigType('png')
    return plotProfiles([profile], [rep_reads], [pam_idx], [False], ['Predicted'], title=title)

if __name__ == '__main__':
    
    theta_file =  'model_output_2000_0.01000000_1.835_theta.txt_cf0.txt'
    target_seq = 'GGCCAGCGGAGCATGCATGCAGGGAAGATGAGAGTGATGTAGCAGCGATCGATGCTAGCACG'
    pam_idx = 20
    plot_predictions(theta_file, target_seq, pam_idx)
    import pdb; pdb.set_trace()
    
