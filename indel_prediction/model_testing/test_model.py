import pandas as pd
import numpy as np
import random

from mpi4py import MPI

import io, os, sys, csv, time
from multiprocessing import Process, Pipe

from scipy.stats import pearsonr, spearmanr

from sklearn.model_selection import KFold
from sklearn import metrics

from selftarget.data import getSampleSelectors, getAllDataDirs
from selftarget.oligo import loadOldNewMapping, partitionGuides, getFileForOligoIdx, getOligoIdxFromId
from selftarget.profile import getProfileCounts

from predictor.model import writeTheta, readTheta, printAndFlush, trainModelParallel, testModelParallel, recordPredictions

comm = MPI.COMM_WORLD
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()

NUM_OLIGO = -1
FOLD = 2
OUT_PREFIX = 'model_output'

def getModelDevGuideSet(guideset_file):
    f = io.open(guideset_file)
    guideset = [line[:-1] for line in f]
    f.cloes()
    return np.array(guideset)

def loadFeatureLabels(oligo_id):
    data = loadOligoFeaturesAndReadCounts(oligo_id, [], FEATURES_DIR)
    return [x for x in data.columns if x not in ['Oligo ID','Indel','Frac Sample Reads','Left','Right','Inserted Seq']]

def runAnalysis(guideset_file = 'model_development_guideset.txt'):

    guideset = getFullModelDevGuideSet(guideset_file)
    sample_names = ['ST_Feb_2018_CAS9_12NA_1600X_DPI7', 'ST_June_2017_K562_800x_LV7A_DPI7', 'ST_June_2017_K562_800x_LV7B_DPI7']
    
    feature_columns= loadFeatureLabels([x for x in guideset][0])
    if NUM_OLIGO != -1:
        guideset = random.sample([x for x in guideset],NUM_OLIGO)

    kf = KFold(n_splits=2)
    for i,(train_idx, test_idx) in enumerate(kf.split(guideset)):
        printAndFlush('Cross Validation Fold %d' % (i+1))
        train_set, test_set = np.array(guideset)[train_idx], np.array(guideset)[test_idx]

        outfile = OUT_THETA_FILE + '_cf%d.txt' % i

        theta0 = None
        tmp_file = 'tmp_%s_%d.txt' % (OUT_THETA_FILE, i)
        if os.path.isfile(tmp_file):
            printAndFlush('Loading from previous tmp file')
            theta0, rec_train_set, feature_columns = readTheta(tmp_file)
            test_set = [x for x in ([y for y in train_set] + [y for y in test_set]) if x not in rec_train_set][:int(NUM_OLIGO/2)]
            train_set = rec_train_set

        printAndFlush('Training')
        theta = trainModelParallel(train_set, sample_names, feature_columns, theta0, cv_idx=i)
        testModelParallel( theta, train_set, sample_names, feature_columns )    #Check final training result with lambda=0
        writeTheta(OUT_THETA_FILE + '_cf%d.txt' % i, feature_columns, theta, train_set)
        recordPredictions(OUT_PROFILE_DIR + '_train_%d' % i, theta, train_set, feature_columns )
        
        printAndFlush('Testing')
        testModelParallel( theta, test_set, sample_names, feature_columns )
        recordPredictions(OUT_PROFILE_DIR + '_test_%d' % i, theta, test_set, feature_columns )

if __name__ == '__main__':
    if len(sys.argv) > 1: NUM_OLIGO = eval(sys.argv[1])
    if len(sys.argv) > 3: REG_CONST = eval(sys.argv[3])
    if len(sys.argv) > 4: OUT_PREFIX = sys.argv[4]
    else:
        rand_val = np.random.normal(loc=0.0, scale=1.0)
        rand_val = comm.bcast(rand_val, root=0)
        OUT_PREFIX = 'model_output_%d_%.8f_%.3f' % (NUM_OLIGO, REG_CONST, rand_val )
    OUT_PROFILE_DIR = OUT_PREFIX + '_predictions'
    OUT_THETA_FILE = OUT_PREFIX + '_theta.txt'
    runAnalysis()
   

