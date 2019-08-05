import io, os, sys, csv
import pandas as pd
import numpy as np
import pylab as PL
from sklearn import linear_model
from matplotlib import colors as mcolors
import Bio.Seq

from selftarget.indel import tokFullIndel
from selftarget.util import mergeSamples, getPlotDir, analyseResultsPerPartition, defaultLoadData
from selftarget.plot import plotBarSummary, saveFig
from selftarget.data import setHighDataDir, getHighDataDir

MIN_READS = 20
PLOT_COLORS = ['red','blue','green','orange','purple','cyan','lightblue','darkgreen','yellow','magenta','gray','navy']

def getRegrLine(xdata, ydata):
    regr = linear_model.LinearRegression(fit_intercept=False)
    regr.fit(xdata, ydata)
    min_x, max_x = xdata.min()[0], xdata.max()[0]
    return [min_x, max_x],[regr.predict(min_x)[0], regr.predict(max_x)[0]], regr.predict(1)

def loadData(results_file, guideset=set()):
    data = pd.read_csv(results_file, sep='\t')
    if len(guideset) > 0: data = data.loc[data['Oligo ID'].isin(guideset)] 
    data = data.loc[(data['Mut Non-Null Reads'] > MIN_READS) & (data['Orig Non-Null Reads'] > MIN_READS)]

    fout = io.open('mhmm_guides_in_k562_samples.txt', 'w')
    for oligo_id in data['Oligo ID'].unique(): fout.write(u'%s\n' % oligo_id)
    fout.close()
    return data

def passData(data,label=''):
    return data

def getMismatch(row, idx=0):
    orig_seq = row['Orig MH']
    mismatch_seq = row['Left Mut-MH'] if row['Left Mut-MH'] != orig_seq else row['Right Mut-MH']
    mismatches = [(x,y) for (x,y) in zip(orig_seq,mismatch_seq) if (x!=y)]
    if len(mismatches) == 0: return ''
    orig_nt, mm_nt = mismatches[idx]
    return orig_nt + '-' + mm_nt

def getLastMismatch(row):
    return getMismatch(row, idx=-1)

def getMhGC(row):
    return sum([1 if (x==y and x in ['G','C']) else 0  for (x,y) in zip(row['Left Mut-MH'],row['Right Mut-MH'])])

def isMatched(orow):
    grna_seq_19mer = orow['Guide'][1:]
    return (grna_seq_19mer in orow['Target']) or (Bio.Seq.reverse_complement(grna_seq_19mer) in orow['Target'])

def passData(data, label=''):
    return data

def plotMicrohomologyMismatches(all_result_outputs, label=''):
    
    mut_hdrs =  ['Left Mut', 'Right Mut','Merged Mut1', 'Merged Mut2']
    cols_to_sum = [x + ' Indel Reads in Mut' for x in mut_hdrs] + ['Orig Indel Reads in Orig', 'Mut Non-Null Reads', 'Orig Non-Null Reads'] 
    common_cols = ['Oligo ID','Mapped Oligo Id','Num Mismatches', 'Orig Indel','Orig MH','Left Mut-MH','Right Mut-MH'] #,'Merged Mut 1 MH','Merged Mut 2 MH','Orig Indel','Left Mut-MH Indel','Right Mut-MH Indel','Merge Mut 1 Indel','Merge Mut 2 Indel']
    data =  mergeSamples(all_result_outputs, cols_to_sum, merge_on=common_cols)

    getLeft = lambda indel: tokFullIndel(indel)[2]['L']
    getRight = lambda indel: tokFullIndel(indel)[2]['R']
    getMHSize = lambda indel: tokFullIndel(indel)[2]['C']

    oligo_data = pd.read_csv(getHighDataDir() + '/ST_June_2017/data/self_target_oligos_details_with_pam_details.csv', sep='\t')
    oligo_data['Guide is matched'] = oligo_data.apply(isMatched, axis=1)
    reverse_lookup = {x: y == 'REVERSE' for (x,y) in zip(oligo_data['ID'],oligo_data['PAM Direction'])}
    is_reverse = lambda x: reverse_lookup[x]

    data = pd.merge(data, oligo_data[['ID','Guide is matched']], left_on='Oligo ID', right_on='ID', how='inner')

    data['MH Size'] = data['Orig Indel'].apply(getMHSize)
    data = data.loc[(data['MH Size'] != 0) & (data['Guide is matched'])]
    data['MH Left Loc'] = data['Orig Indel'].apply(getLeft) + data['MH Size']
    data['MH Right Loc'] = data['Orig Indel'].apply(getRight) - data['MH Size']
    data['Is Reverse'] = data['Oligo ID'].apply(is_reverse)

    for hdrL,hdrR in [mut_hdrs[:2], mut_hdrs[2:]]:
        data[hdrL + ' Reads'] = data['Is Reverse']*data[hdrR + ' Indel Reads in Mut Sum'] + (1- data['Is Reverse'])*data[hdrL + ' Indel Reads in Mut Sum']
        data[hdrR + ' Reads'] = data['Is Reverse']*data[hdrL + ' Indel Reads in Mut Sum'] + (1- data['Is Reverse'])*data[hdrR + ' Indel Reads in Mut Sum']
        data[hdrL + ' Reads Ratio'] =  data[hdrL + ' Reads']*100.0/data['Mut Non-Null Reads Sum']
        data[hdrR + ' Reads Ratio'] =  data[hdrR + ' Reads']*100.0/data['Mut Non-Null Reads Sum']
    data['Orig Indel Reads Ratio'] = data['Orig Indel Reads in Orig Sum']*100.0/data['Orig Non-Null Reads Sum']
    data['All Mut Reads Ratio'] = (data[[x + ' Reads' for x in mut_hdrs]].sum(axis=1))*100.0/data['Mut Non-Null Reads Sum']
    data['MH Dist'] = data['MH Right Loc'] - data['MH Left Loc']
    data['1st Mismatch'] = data.apply(getMismatch, axis=1)
    data['Last Mismatch'] = data.apply(getLastMismatch, axis=1)
    data['MH GC Content'] = data.apply(getMhGC, axis=1)

    mh_indel_types = [('Orig Indel','Left Mut'), ('Orig Indel','Right Mut'), ('Orig Indel','All Mut'),('Left Mut','Right Mut') ]
 
    label_lookup = {'Orig Indel': 'Percent of mutated reads of\ncorresponding deletion (no mismatch)',
                    'Left Mut': 'Percent of mutated reads of corresponding\ndeletion with retained left sequence',
                    'Right Mut': 'Percent of mutated reads of corresponding\ndeletion with retained right sequence',
                    'All Mut': 'Percent of mutated reads of\ncorresponding deletion (mismatches)'
        }

    fig1 = PL.figure(figsize=(4,4))
    fig_all = PL.figure(figsize=(10,10))
    for i, (mh_typex, mh_typey) in enumerate(mh_indel_types):
        figs = [(fig_all, True), (fig1,False)] if i==2 else [(fig_all, True)]
        for fig, is_all in figs:
            PL.figure(fig.number)
            if is_all: PL.subplot(2,2,i+1)
            for nm,clr  in zip([1,2],['royalblue','orange']):
                nm_data = data.loc[data['Num Mismatches'] == nm]

                sty, lsty = 'o', '-'
                sel_data = nm_data.loc[(nm_data['MH Size'] >= 6) & (nm_data['MH Size'] <= 15)]

                rx, ry, grad = getRegrLine(sel_data[[mh_typex + ' Reads Ratio']], sel_data[[mh_typey + ' Reads Ratio']])
                PL.plot(sel_data[mh_typex + ' Reads Ratio'], sel_data[mh_typey + ' Reads Ratio'], sty, color=clr, markersize=4, label='%d mismatch%s (R=%.2f, %d gRNAs)' % (nm, 'es' if nm==2 else '', grad, len(sel_data)), alpha=0.5)

                if not is_all: print(grad, nm, mh_typex, mh_typey)
                if i != 3: PL.plot(rx, ry, lsty, color=clr, linewidth=2)

            PL.xlabel(label_lookup[mh_typex])
            PL.ylabel(label_lookup[mh_typey])
            PL.xlim((0,80))
            PL.ylim((0,80))
            PL.plot([0,80],[0,80],'k--')
            PL.legend()
            PL.show(block=False)
    saveFig('mm_mismatch_all')
    PL.figure(fig1.number)
    saveFig('mm_mismatch_one')

def runAnalysis():
	
    spec = {'results_dir':getHighDataDir() + '/microhomology_mismatch/mh_mismatch_indel_frequencies',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(passData,'Data')],
            'py_funcs_all_results': [plotMicrohomologyMismatches],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Orig Non-Null Reads',
            'min_reads': MIN_READS,
            'id_colname': 'Oligo ID',
            'partitions': ['Non-Targeting'],
            'samples': ['K562 New']
            }
    analyseResultsPerPartition( spec )

if __name__ == '__main__':
    setHighDataDir('..')
    runAnalysis()
    import pdb; pdb.set_trace()
