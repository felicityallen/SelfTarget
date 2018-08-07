import io, sys, os, csv
import pylab as PL
import numpy as np
import itertools
import pandas

from selftarget.oligo import partitionGuides
from selftarget.util import getPickleDir
from selftarget.data import getAllDataDirs, getSampleSelectors, sortSampleNames, getDirLabel, getSimpleName, parseSampleName, getHighDataDir, setHighDataDir
from selftarget.plot import saveFig

from scipy.stats import pearsonr

ST_COMPARISON_RESULTS_DIR = 'kl_comparisons/kl_comparison_summaries'
MIN_READS = 20 

def getDirsFromFilename( filename ):
    idx = filename.find('_vs_')
    if idx < 0: return(None, None)
    return(filename[:idx],filename[idx+4:-4])

def loadAllData( guideset, sample_selector=lambda x: True, label='', cols=['KL without null'], allow_pickle=False):
    pickle_file = '%s/kl_analysis_%s.pickle' % (getPickleDir(), label.replace(' ','_'))
    if os.path.exists(pickle_file) and allow_pickle:
        merged_data = pandas.read_pickle(pickle_file)
    else:
        cmp_files = os.listdir(getHighDataDir() + '/' + ST_COMPARISON_RESULTS_DIR)
        merged_data = None
        for filename in cmp_files:
            dir1, dir2 = getDirsFromFilename(filename)
            if not sample_selector(dir1) or not sample_selector(dir2): continue
        
            #Load data from file
            data = pandas.read_csv(getHighDataDir() + '/' + ST_COMPARISON_RESULTS_DIR + '/' + filename, sep='\t')
            data['Mutated Reads 1'] = data['Num Reads 1'] - data['Num null reads 1']
            data['Mutated Reads 2'] = data['Num Reads 2'] - data['Num null reads 2']
            data = data.loc[data['Mutated Reads 1'] > MIN_READS]
            data = data.loc[data['Mutated Reads 2'] > MIN_READS]
            data = data.loc[data['ID'].isin(guideset)][['ID'] + cols]
            if merged_data is not None and len(data) < 0.75*len(merged_data):
                print('Skipping %s, data for insufficient guides (%d vs %d)' % (filename, len(data), len(merged_data)))
                continue

            #Merge with the other data (keep only common Oligos)
            suffix_fn = lambda x: '$' + x
            if merged_data is None:  
                merged_data, first_suffix = data, suffix_fn(filename)
            else:
                merged_data = merged_data.merge(data, how='inner', on='ID', suffixes=('',suffix_fn(filename)))
            print(len(merged_data), filename)
           
        merged_data = merged_data.rename(columns={x: (x + first_suffix) for x in cols})
        if allow_pickle:
            merged_data.to_pickle(pickle_file)
    return merged_data

def getUniqueSamples(colnames):
    samples = set()
    for colname in colnames:
        dir1, dir2 = getDirsFromFilename(colname.split('$')[-1])
        samples.add(dir1)
        samples.add(dir2)
    return [x for x in samples]

def plotHeatMap(data, col='KL without null', label=''):

    #Compute and collate medians
    sel_cols = [x for x in data.columns if col in x]
    cmp_meds = data[sel_cols].median(axis=0)
    samples = sortSampleNames(getUniqueSamples(sel_cols))
    cell_lines = ['CHO', 'E14TG2A', 'BOB','RPE1', 'HAP1','K562','eCAS9','TREX2']
    sample_idxs = [(cell_lines.index(parseSampleName(x)[0]),x) for x in getUniqueSamples(sel_cols)]
    sample_idxs.sort()
    samples = [x[1] for x in sample_idxs]

    N = len(samples)
    meds = np.zeros((N,N))
    for colname in sel_cols:
        dir1, dir2 = getDirsFromFilename(colname.split('$')[-1])
        idx1, idx2 = samples.index(dir1), samples.index(dir2)
        meds[idx1,idx2] = cmp_meds[colname]
        meds[idx2,idx1] = cmp_meds[colname]

    for i in range(N):
        print(' '.join(['%.2f' % x for x in meds[i,:]]))
    print( np.median(meds[:,:-4],axis=0))

	#Display in Heatmap
    PL.figure(figsize=(5,5))
    PL.imshow(meds, cmap='hot_r', vmin = 0.0, vmax = 3.0, interpolation='nearest')
    PL.colorbar()
    PL.xticks(range(N))
    PL.yticks(range(N))
    PL.title("Median KL") # between %d mutational profiles (for %s with >%d mutated reads)" % (col, len(data), label, MIN_READS))
    ax1 = PL.gca()
    ax1.set_yticklabels([getSimpleName(x) for x in samples], rotation='horizontal')
    ax1.set_xticklabels([getSimpleName(x) for x in samples], rotation='vertical')
    PL.subplots_adjust(left=0.25,right=0.95,top=0.95, bottom=0.25)
    PL.show(block=False) 
    saveFig('median_kl_heatmap_cell_lines')

def runAnalysis():
	
    partitions = partitionGuides(oligo_detail_dir=getHighDataDir()+ '/ST_June_2017/data')

    for part_desc in ['Real Guides']:

        selector = getSampleSelectors()['DPI7']
        guideset = partitions[part_desc]

        desc = part_desc + ' DPI7'
        data = loadAllData(guideset, sample_selector=selector, label=desc)
        plotHeatMap(data, label=desc)


if __name__ == '__main__':
    setHighDataDir('..')
    runAnalysis()
    import pdb; pdb.set_trace()
	