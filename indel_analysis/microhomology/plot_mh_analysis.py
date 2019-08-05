import io, csv, sys, os
import pylab as PL
import numpy as np
import random
import scipy.stats
from sklearn import linear_model
import pandas as pd

from scipy.stats import pearsonr, spearmanr

MIN_READS = 20
NUM_PSEUDO_READS = 5

from selftarget.util import getPlotDir, analyseResultsPerPartition, defaultLoadData
from selftarget.data import parseSampleName, getSimpleName, setHighDataDir, getHighDataDir
from selftarget.plot import saveFig

def loadData(filename, guideset=set()):
    mhlen = eval(filename.split('_')[-1][:-4])
    data = pd.read_csv(filename, sep='\t')
    if len(guideset) > 0: data = data.loc[data['Oligo ID'].isin(guideset)] 
    data['MH Dist'] = data['Right Position'] - data['Left Position'] - mhlen
    data['Percent Non-Null Reads'] = data['Indel Reads']*100.0/data['Non-Null Reads']
    data['MH Len'] = mhlen
    return data

def loadAllMHLenData(dirname, guideset=set()):
    files = [dirname + '/' + x for x in os.listdir(dirname)]
    return pd.concat([loadData(x, guideset=guideset) for x in files])
    
def loadAllMHLenAndOtherData(dirname, guideset=set()):
    data_mh =loadAllMHLenData(dirname, guideset=guideset)
    return data_mh

def passData(data, label='test'):
    return data

def plotPercScatterAnalysis(data, label='test', y_axis = 'Percent Non-Null Reads', plot_scatters=False, plot_regr_lines=False, scatter_mh_lens=[], mh_lens=[9]):
    
    plot_dir = getPlotDir()
    regr_lines = {}
    for mh_len in mh_lens:
        mh_data = data.loc[data['MH Len'] == mh_len]
        mh_rdata = mh_data.loc[(mh_data['MH Dist'] >= 0) & (mh_data['MH Dist'] < (30-mh_len)) ]
        
        regr = linear_model.LinearRegression()
        rx, ry = mh_rdata[['MH Dist']], mh_rdata[[y_axis]] #np.log(mh_rdata[[y_axis]])
        regr.fit(rx, ry)
        corr = scipy.stats.pearsonr(rx, ry)
        min_x, max_x = rx.min()[0], rx.max()[0]
        x_pts = [min_x, max_x]
        regr_lines[mh_len] = (x_pts,[regr.predict(x)[0] for x in x_pts],corr[0])
        
        if plot_scatters and mh_len in scatter_mh_lens:
            fig = PL.figure(figsize=(5,5))
            PL.plot( mh_data['MH Dist'], mh_data[y_axis], '.', alpha=0.4 )
            PL.plot(regr_lines[mh_len][0],regr_lines[mh_len][1],'dodgerblue',linewidth=3)
        
            PL.xlabel('Distance between nearest ends of\nmicrohomologous sequences',fontsize=14)
            PL.ylabel('Percent of mutated reads of corresponding\nMH-mediated deletion',fontsize=14)
            PL.tick_params(labelsize=14)
            PL.xlim((0,20))
            PL.title('Microhomology of length %d (r=%.2f)' % (mh_len,corr[0]),fontsize=14)
            PL.show(block=False)  
            saveFig('mh_scatter_len%d_%s' % (mh_len,label.split('/')[-1])) 
    
    if plot_regr_lines:
        fig = PL.figure()
        output_data = {}
        for mh_len in mh_lens:
            fit_data = regr_lines[mh_len]
            if mh_len > 15:
                continue
            lsty = '--' if mh_len < 9 else '-'
            PL.plot(fit_data[0], fit_data[1], linewidth=2, linestyle=lsty, label='MH length %d (R=%.1f)' % (mh_len, fit_data[2]))
        PL.title(label,fontsize=18)
        PL.xlabel('Distance between nearest ends of\nmicrohomologous sequences',fontsize=14)
        PL.ylabel('Percent of mutated reads of corresponding\nMH-mediated deletion',fontsize=14)
          
        PL.tick_params(labelsize=18)
        PL.legend()
        PL.ylim((0,100))
        PL.show(block=False)  
        saveFig(plot_dir + '/mh_scatter_all_len_%s' % label.split('/')[-1]) 
    return regr_lines

def plotK562PercScatterAnalysis(data, label=''):
    return plotPercScatterAnalysis(data, label=label, plot_scatters=True, plot_regr_lines=False, scatter_mh_lens=[9], mh_lens=[3,4,5,6,7,8,9,10,11,12,13,14,15])

def compareMHlines(all_result_outputs, label='', y_axis = 'Percent Non-Null Reads', data_label='RegrLines'):

    color_map = {'K562':'b','K562_1600x':'lightblue', 'BOB':'g', 'RPE1':'purple', 'TREX2':'k', 'TREX2_2A':'gray', 'HAP1':'r', 'E14TG2A':'orange','eCAS9':'c', 'WT':'pink', 'CHO':'salmon'}
    lysty_map = {2:'--',3:'--',5:'--',7:'-',10:'-.',16:':',20:':'}

    dirnames = [x[1] for x in all_result_outputs]
    lystys = [lysty_map[parseSampleName(x)[1]] for x in dirnames]
    clrs = [color_map[parseSampleName(x)[0]] for x in dirnames]

    for mh_len in [9]:
        PL.figure()
        regr_lines = [x[0][data_label][mh_len] for x in all_result_outputs]
        for dirname, regr_line,lysty,clr in zip(dirnames,regr_lines,lystys, clrs):
            PL.plot(regr_line[0], regr_line[1], label='%s (R=%.1f)' % (getSimpleName(dirname), regr_line[2]), linewidth=2, color=clr, linestyle=lysty, alpha=0.5 )
        PL.xlabel('Distance between nearest ends of microhomologous sequences',fontsize=14)
        PL.ylabel('Corresponding microhomology-mediated deletion\n as percent of total mutated reads',fontsize=14)
        PL.tick_params(labelsize=16)
        PL.legend(loc='upper right')
        PL.ylim((0,70))
        PL.xlim((0,20))
        PL.xticks(range(0,21,5))
        PL.title('Microhomology Length %d' % mh_len,fontsize=18)
        PL.show(block=False)  
        saveFig('mh_regr_lines_all_samples__%d' % mh_len)

def compareMHK562lines(all_result_outputs, label='', y_axis = 'Percent Non-Null Reads', data_label='RegrLines'):

    dirnames = [x[1] for x in all_result_outputs]
    clrs = ['silver','grey','darkgreen','green','lightgreen','royalblue','dodgerblue','skyblue','mediumpurple','orchid','red','orange','salmon']

    fig = PL.figure(figsize=(6,6))
    leg_handles = []
    mh_lens = [3,4,5,6,7,8,9,10,11,12,13,14,15]
    for mh_len, clr in zip(mh_lens,clrs):
        regr_lines = [x[0][data_label][mh_len] for x in all_result_outputs]
        mean_line = np.mean([x[:2] for x in regr_lines], axis=0)
        leg_handles.append(PL.plot(mean_line[0], mean_line[1], label='MH Len=%d  (R=%.1f)' % (mh_len,np.mean([x[2] for x in regr_lines])) , linewidth=2, color=clr )[0])
    PL.xlabel('Distance between nearest ends of\nmicrohomologous sequences',fontsize=16)
    PL.ylabel('Correspondng microhomology-mediated deletion\n as percent of total mutated reads',fontsize=16)
    PL.tick_params(labelsize=16)
    PL.legend(handles=[x for x in reversed(leg_handles)], loc='upper right')
    PL.ylim((0,80))
    PL.xlim((0,20))
    PL.xticks(range(0,21,5))
    PL.show(block=False)  
    saveFig('mh_regr_lines_K562')

def getRegrValue(mh_len, mh_dist, mean_lines):
    x1,x2 = mean_lines[mh_len][0]
    y1,y2 = mean_lines[mh_len][1]
    grad = (y2-y1)*1.0/(x2-x1)
    intercept = y1-grad*x1
    return mh_dist*grad + intercept
    
def plotGCContent(all_result_outputs, label=''):
    
    #Merge data across samples
    unique_cols = ['Oligo ID','Indel', 'GC Content', 'MH Len', 'MH Dist'] 
    datas = [x[0]['Data'][unique_cols + ['Indel Reads', 'Non-Null Reads']] for x in all_result_outputs]
    merged_data = datas[0]
    for i, data in enumerate(datas[1:]):
        merged_data = pd.merge(merged_data, data, on=unique_cols, suffixes=('','%d' % (i+2)), how='outer')
    suffix = lambda i: '%d' % (i+1) if i > 0 else ''
    merged_data['Indel Reads Sum'] = merged_data[['Indel Reads' + suffix(i) for i in range(len(datas))]].sum(axis=1)
    merged_data['Non-Null Reads Sum'] = merged_data[['Non-Null Reads' + suffix(i) for i in range(len(datas))]].sum(axis=1)

    #Compute mean regression lines across samples for each MH length
    mean_lines = {}
    for mh_len in range(2,16):
        if mh_len not in all_result_outputs[0][0]['RegrLines']: continue
        regr_lines = [x[0]['RegrLines'][mh_len][:2] for x in all_result_outputs]
        mean_lines[mh_len] = np.mean(regr_lines, axis=0)

    #Restrict to only MH dist in (0,10) and adjust for mh len-dist relationship
    for mh_len in [9]:
        compute_resid = lambda row: row['Perc Reads']# - getRegrValue(row['MH Len'],row['MH Dist'],mean_lines)
        sel_data = merged_data.loc[(merged_data['MH Len'] == mh_len) & (merged_data['MH Dist'] >= 0) & (merged_data['MH Dist'] <= 10)]
        sel_data['Perc Reads'] = sel_data['Indel Reads Sum']*100.0/sel_data['Non-Null Reads Sum']
        sel_data['Perc Reads Residual'] = sel_data.apply(compute_resid, axis=1)
        PL.figure(figsize=(4,4))
        gcs = sel_data['GC Content'].unique(); gcs.sort()
        boxdata_lk = {gc: sel_data.loc[sel_data['GC Content'] == gc]['Perc Reads Residual'] for gc in gcs}
        gcs = [gc for gc in gcs if len(boxdata_lk[gc])>20]  #Limit to GC with at least 20 data points
        boxdata = [boxdata_lk[gc] for gc in gcs]
        print([len(x) for x in boxdata])
        PL.boxplot(boxdata)
        PL.ylabel('Percent total mutated reads of MH-mediated deletion')
        PL.xlabel('GC content of microhomologous sequence')
        PL.title('Microhomology of length %d\n(at max 10 distance)' % mh_len)
        PL.xticks(range(1,len(gcs)+1),gcs)
        PL.show(block=False)
        saveFig('gc_content_mh%d' % mh_len)

def runAnalysis():

    spec = {'results_dir': getHighDataDir() + '/microhomology/mh_freqs_by_len',
            'dirname_to_result_fn': lambda x: x,
            'result_to_dirname_fn': lambda x: x,
            'py_func_load': loadAllMHLenData,
            'py_funcs_per_result':  [(plotK562PercScatterAnalysis,'RegrLines'), (passData, 'Data')],
            'py_funcs_all_results': [compareMHK562lines, plotGCContent],
            'reads_colname': 'Non-Null Reads',
            'check_output_fn': lambda x: True, 
            'id_colname': 'Oligo ID',
            'min_reads': MIN_READS,
            'partitions': ['Non-Targeting'],
            'samples': ['K562 New']
            }
    analyseResultsPerPartition( spec ) 

    spec = {'results_dir': getHighDataDir() + '/microhomology/mh_freqs_by_len',
            'dirname_to_result_fn': lambda x: x,
            'result_to_dirname_fn': lambda x: x,
            'py_func_load': loadAllMHLenData,
            'py_funcs_per_result':  [(plotPercScatterAnalysis,'RegrLines')],
            'py_funcs_all_results': [compareMHlines],
            'reads_colname': 'Non-Null Reads',
            'check_output_fn': lambda x: True, 
            'id_colname': 'Oligo ID',
            'min_reads': MIN_READS,
            'partitions': ['Non-Targeting'],
            'samples': ['DPI7']
            }
    analyseResultsPerPartition( spec ) 

if __name__ == '__main__':
    setHighDataDir('..')
    runAnalysis()
    import pdb; pdb.set_trace()