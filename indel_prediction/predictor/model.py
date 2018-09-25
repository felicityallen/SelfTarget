import pandas as pd
import numpy as np
import random

# from mpi4py import MPI

import io, os, sys, csv, time
from multiprocessing import Process, Pipe

from scipy.stats import pearsonr, spearmanr
from scipy.optimize import minimize

from selftarget.data import getSampleSelectors, getAllDataDirs
from selftarget.oligo import loadOldNewMapping, partitionGuides, getFileForOligoIdx, getOligoIdxFromId
from selftarget.profile import getProfileCounts

from predictor.features import readFeaturesData

# comm = MPI.COMM_WORLD
mpi_rank = 0
mpi_size = 1

REG_CONST = 0.01
I1_REG_CONST = 0.01
GEN_INDEL_FEATURES_DIR = '/features_ext_for_gen_indels'
GEN_INDEL_LOC = '../../indel_analysis/gen_indels'
FEATURES_DIR = GEN_INDEL_LOC + GEN_INDEL_FEATURES_DIR
READS_DIR = GEN_INDEL_LOC + '/reads_for_gen_indels_all_samples'
OUT_THETA_FILE = 'model_thetas.txt'

def setOutputThetaFile(filename):
    global OUT_THETA_FILE
    OUT_THETA_FILE = filename

def setFeaturesDir(dirname):
    global FEATURES_DIR
    FEATURES_DIR = dirname

def setReadsDir(dirname):
    global READS_DIR
    READS_DIR = dirname

def setRegConst(val):
    global REG_CONST
    REG_CONST = val

def setI1RegConst(val):
    global I1_REG_CONST
    I1_REG_CONST = val

def writeTheta(out_file, feature_columns, theta, train_set):
    fout = io.open(out_file, 'w')
    fout.write(u'%s\n' % ','.join([x for x in train_set]))
    lines = [u'%s\t%s\n' % (x,y) for (x,y) in zip(feature_columns, theta)]
    fout.write(''.join(lines))
    fout.close()

def readTheta(theta_file):
    f = io.open(theta_file)
    train_set = f.readline()[:-1].split(',')
    feature_columns, theta = [], []
    for toks in csv.reader(f, delimiter='\t'):
        feature_columns.append(toks[0])
        theta.append(eval(toks[1]))
    return theta, train_set, feature_columns

def printAndFlush(msg, master_only=True):
    if not master_only or mpi_rank == 0:
        print(msg)
        sys.stdout.flush()

def getCutSite(features_file):
    f = io.open(features_file); f.readline()
    cut_site = eval(f.readline().split('\t')[1])
    return cut_site

def loadOligoFeaturesAndReadCounts(oligo_id, sample_names):

    oligo_idx = getOligoIdxFromId(oligo_id)
    oligo_subdir, _ = getFileForOligoIdx(oligo_idx, ext='')

    features_file = FEATURES_DIR + '/' + oligo_subdir + '/%s_gen_indel_features.txt' % oligo_id
    reads_file = READS_DIR + '/' + oligo_subdir + '/%s_gen_indel_reads.txt' % oligo_id

    cut_site = getCutSite(features_file)
    indel_feature_data, feature_cols = readFeaturesData(features_file)
    
    if len(sample_names) > 0:
        read_data =  pd.read_csv(reads_file, skiprows=1, sep='\t')
        read_data['Sum Sample Reads'] = read_data[sample_names].sum(axis=1) + 0.5
        read_data = read_data.loc[read_data['Indel']!='All Mutated']
        total_mut_reads = read_data['Sum Sample Reads'].sum()
        if total_mut_reads == 0: raise Exception('No Mutated Reads in %s' % reads_file)
        read_data['Frac Sample Reads'] = read_data['Sum Sample Reads']/total_mut_reads
        merged_data = pd.merge(indel_feature_data, read_data[['Indel','Frac Sample Reads']], left_index=True, right_on='Indel', how='inner')
    else:
        merged_data = indel_feature_data
        merged_data['Indel'] = merged_data.index

    return merged_data

def calcThetaX(row, theta, feature_columns):
    return sum(theta*row[feature_columns])

def computeRegularisers(theta, feature_columns, reg_const, i1_reg_const):
    Q_reg = sum([i1_reg_const*val**2.0 if 'I' in name else reg_const*val**2.0 for (val, name) in zip(theta, feature_columns)])
    grad_reg = theta*np.array([i1_reg_const if 'I' in name else reg_const for name in feature_columns])
    return Q_reg, grad_reg

def computeKLObjAndGradients(theta, guideset, sample_names, feature_columns, reg_const, i1_reg_const):
    N = len(feature_columns)
    Q, jac, minQ, maxQ = 0.0, np.zeros(N), 0.0, 1000.0
    Qs = []
    for oligo_id in guideset:
        data = loadOligoFeaturesAndReadCounts(oligo_id, sample_names)
        Y = data['Frac Sample Reads']
        data['ThetaX'] = data.apply(calcThetaX, axis=1, args=(theta,feature_columns))
        sum_exp = np.exp(data['ThetaX']).sum()
        Q_reg, grad_reg =  computeRegularisers(theta, feature_columns, reg_const, i1_reg_const)
        tmpQ = (np.log(sum_exp) + sum(Y*(np.log(Y) - data['ThetaX'])) + Q_reg)
        Q += tmpQ
        Qs.append(tmpQ)
        jac += np.matmul(np.exp(data['ThetaX']),data[feature_columns].astype(int))/sum_exp - np.matmul(Y,data[feature_columns].astype(int)) + grad_reg
    return Q, jac, Qs 

def assessFit(theta, guideset, sample_names, feature_columns, cv_idx=0, reg_const=REG_CONST, i1_reg_const=I1_REG_CONST, test_only=False):
    #Send out thetas
    theta, done = comm.bcast((theta, False), root=0)
    while not done:
        #Compute objective and gradients
        Q, jac, Qs = computeKLObjAndGradients(theta, guideset, sample_names, feature_columns, reg_const, i1_reg_const)

        #Combine all
        full_guideset = comm.gather([x for x in guideset], root=0)
        flatten = lambda l: [item for sublist in l for item in sublist]
        objs_and_grads = comm.gather((Q, jac, Qs), root=0)
        if mpi_rank == 0:
            Q, jac, Qs = sum([x[0] for x in objs_and_grads]), sum([x[1] for x in objs_and_grads]), []
            for x in objs_and_grads: Qs.extend(x[2]) 
            Q, jac, Qs = Q/len(Qs), jac/len(Qs), Qs
            printAndFlush(' '.join(['Q=%.5f' % Q, 'Min=%.3f' % min(Qs), 'Max=%.3f' % max(Qs), 'Num=%d' % len(Qs), 'Lambda=%e' % reg_const, 'I1_Lambda=%e' % i1_reg_const]))
            writeTheta('tmp_%s_%d.txt' % (OUT_THETA_FILE, cv_idx), feature_columns, theta, flatten(full_guideset))
    
        Q, jac, Qs = comm.bcast((Q, jac, Qs), root=0)
        if mpi_rank == 0 or test_only: done = True
        else:
            done = False
            theta, done = comm.bcast((theta, done), root=0)
    return Q, jac, Qs

def debugIndel(theta, data, indel, feature_columns):
    indel_data = data.loc[data['Indel']==indel]
    for index,row in indel_data.iterrows():
        print(indel, [(x,theta) for (x,y,theta) in zip(feature_columns,[row[x] for x in feature_columns],theta) if y])

def computePredictedProfile(data, theta, feature_columns):
    data['expThetaX'] = np.exp(data.apply(calcThetaX, axis=1, args=(theta,feature_columns)))
    sum_exp = data['expThetaX'].sum()
    profile = {x: expthetax*1000/sum_exp for (x,expthetax) in zip(data['Indel'],data['expThetaX'])}
    counts = getProfileCounts(profile)
    return profile, counts
       
def recordProfiles( output_dir, theta, guideset, feature_columns ):
    while not os.path.isdir(output_dir):
        if mpi_rank == 0: 
            os.mkdir(output_dir)
        else: time.sleep(5)
    for oligo_id in guideset:
        profile, counts = computePredictedProfile(oligo_id, theta, feature_columns)
        idx = getOligoIdxFromId(oligo_id)
        filepath, filename = getFileForOligoIdx(idx)
        if not os.path.isdir(output_dir + '/' + filepath):
            os.mkdir(output_dir + '/' + filepath)
        fout = io.open(output_dir + '/' + filepath + '/%s_mappedindelsummary_predicted.txt' % oligo_id, 'w')
        fout.write(u'@@@%s\n' % oligo_id)
        for val, indel, perc1, perc2 in counts:
            if val >= 1:
                fout.write('%s\t-\t%d\n' % (indel, val))
        fout.close()

def trainModelParallel(guideset, sample_names, feature_columns, theta0, cv_idx=0):
    
    guidesubsets = [guideset[i:len(guideset):mpi_size] for i in range(mpi_size)]
    if theta0 is None: theta0 = np.array([np.random.normal(loc=0.0, scale=1.0) for x in feature_columns])
    args=(guidesubsets[mpi_rank], sample_names, feature_columns, cv_idx, REG_CONST, I1_REG_CONST)
    if mpi_rank == 0:
        result = minimize(assessFit, theta0, args=args, method='L-BFGS-B', jac=True, tol=1e-4)
        theta = result.x
        printAndFlush("Optimization Result: " + str(result.success))
        done = True
        theta, done = comm.bcast((theta, done), root=0)
    else:
        assessFit(theta0, guidesubsets[mpi_rank], sample_names, feature_columns, cv_idx, REG_CONST, I1_REG_CONST)
        theta, done = None, True
    theta, done = comm.bcast((theta, done), root=0)
    return theta

def testModelParallel(theta, guideset, sample_names, feature_columns):
    guidesubsets = [guideset[i:len(guideset):mpi_size] for i in range(mpi_size)]
    assessFit( theta, guidesubsets[mpi_rank], sample_names, feature_columns,reg_const=0.0,i1_reg_const=0.0,test_only=True )

def recordPredictions(output_dir, theta, guideset, feature_columns ):
    guidesubsets = [guideset[i:len(guideset):mpi_size] for i in range(mpi_size)]
    recordProfiles( output_dir, theta, guidesubsets[mpi_rank], feature_columns )



