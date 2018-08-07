import io, os, sys, csv, random
import pandas as pd
import numpy as np
import pylab as PL

from scipy.stats import pearsonr, spearmanr
from sklearn.manifold import TSNE

from selftarget.util import getPlotDir, analyseResultsPerPartition, mergeSamples
from selftarget.plot import plotBarSummary, plotCorrelations, saveFig
from selftarget.data import parseSampleName, setHighDataDir, getHighDataDir

COLORS = ['C0','C9','C3','C1','darkgreen','C2', 'gray','C7','C8']

MIN_READS = 20

MH_DEL_LABELS = ['Large D, MH','Small D, MH']
NON_MH_DEL_LABELS = ['Large D, No MH','Small D, No MH']
DEL_LABELS = MH_DEL_LABELS + NON_MH_DEL_LABELS
INS_LABELS = ['I1','I>1']
ALL_LABELS = INS_LABELS + DEL_LABELS

def loadData(results_file, guideset=set()):
    data = pd.read_csv(results_file, sep='\t')
    if len(guideset) > 0: data = data.loc[data['Oligo Id'].isin(guideset)] 
    data['MCI Ratio'] = data['MCI Reads']*100.0/data['Total reads']
    data['All Deletions'] = data[DEL_LABELS].sum(axis=1)
    data['All Insertions'] = data[INS_LABELS].sum(axis=1)
    data['MH-Mediated Deletions'] = data[MH_DEL_LABELS].sum(axis=1)
    data['Non-MH-Mediated Deletions'] = data[NON_MH_DEL_LABELS].sum(axis=1)
    return data

def passData(data, label=''):
    return data

def perOligoMCI(data, label=''):
    po_data =  data.groupby('Oligo Id')['MCI Type','Most Common Indel'].sum()
    po_data['Oligo Id'] = po_data.index
    return po_data

def perOligoCounts(data, label=''):
    po_data =  data.groupby('Oligo Id')[ALL_LABELS + ['Total reads']].mean()
    po_data['Oligo Id'] = po_data.index
    return po_data

def computePieData(data,label=''):
    pie_labels = ALL_LABELS
    data = perOligoCounts(data)
    pie_data = {x: (data[x]*100.0/(data['Total reads'])).mean(axis=0) for x in pie_labels}
    return pie_data, pie_labels, (data['Total reads']).median()

def computePercentages(data, label='', norm = 'All'):
    data = perOligoCounts(data)
    norm_col, col_labels = 'Total reads', ALL_LABELS
    data['Norm Column'] = data[norm_col]
    for col_label in ALL_LABELS:
        data[col_label + ' Perc'] = data[col_label]*100.0/data[norm_col]
    return data[['Oligo Id', 'Norm Column'] + [x + ' Perc' for x in col_labels]], data[norm_col].median()

def plotPercCorrelations(all_result_outputs, label='', data_label='PercData'):
   scatter_fig = PL.figure()
   col_labels = [x for x in all_result_outputs[0][0]['PercData'][0].columns if 'Perc' in x]
   for i, col_label in enumerate(col_labels):
        scatter_samples = {('K562 800x LV7A DPI7','K562 800x LV7B DPI7'): i+1}
        plotCorrelations(all_result_outputs, label=label, data_label=data_label, y_label=col_label[:-5], val_str=col_label, oligo_id_str='Oligo Id', total_reads_str='Norm Column', scatter_samples=scatter_samples, sdims=(2,3), plot_scatters=True, scatter_fig=scatter_fig, plot_label='K562_cat_correlations', add_leg=(i==5))

def plotSumPie(all_result_outputs, label=''):
    mapping = {'Large D, No MH': 'D>=4,\nno MH', 'Large D, MH': 'D>=4,\nMH', 'Small D, No MH': 'D<4, no MH', 'Small D, MH': 'D<4, MH'}
    merged_data =  mergeSamples(all_result_outputs, ['Total reads'] + ALL_LABELS, data_label='perOligoCounts')
    for col in ALL_LABELS:
        merged_data[col + ' Perc'] = merged_data[col + ' Sum']*100.0/merged_data['Total reads Sum']
    merged_data.to_csv('data_dump_indel_pie.txt',sep='\t',columns=['Oligo Id'] + [col + ' Perc' for col in ALL_LABELS])
    pie_vals = [merged_data[col + ' Perc'].mean() for col in ALL_LABELS]
    PL.figure(figsize=(4,4))
    wedge_labels = [mapping[x] if x in mapping else x for x in ALL_LABELS]
    PL.pie(pie_vals, labels=wedge_labels, autopct='%.1f', labeldistance=1.05, startangle=90.0, counterclock=False, colors=COLORS)
    PL.title('Average distribution\n of mutations\n per gRNA')
    PL.show(block=False)
    saveFig('pie_chart_cats')

def plotMCIPie(all_result_outputs, label=''):
    mapping = {'Large D, No MH': 'D>=4,\nno MH', 'Large D, MH': 'D>=4,\nMH', 'Small D, No MH': 'D<4, no MH', 'Small D, MH': 'D<4, MH'}
    mci_merged_data = mergeSamples(all_result_outputs, ['MCI Type','Most Common Indel'], data_label='perOligoMCI')
    mci_common = mci_merged_data.loc[(mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 2']) & (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 3'])]
    pie_vals, pie_labels = [], []
    for mci_type in ALL_LABELS:
        pie_vals.append(len(mci_common.loc[mci_common['MCI Type'] == mci_type]))
        pie_labels.append(mci_type)
    pie_vals.append(len(mci_merged_data)-len(mci_common))
    pie_labels.append('Inconsistent\nbetween\nreplicates')

    PL.figure(figsize=(4,4))
    PL.pie(pie_vals, labels=[mapping[x] if x in mapping else x for x in pie_labels], autopct='%.1f', labeldistance=1.05, startangle=90.0, counterclock=False, colors=COLORS)
    PL.title('Most frequent\nmutation per gRNA')
    PL.show(block=False)
    saveFig('pie_chart_cats_dominant')

def plotBarSummaryPieIndels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_indel_bar_plots', stacked=True, combine_reps=True, colors=COLORS)

def runAnalysis():
	
    spec = {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(perOligoCounts,'perOligoCounts'),(perOligoMCI,'perOligoMCI'),(computePercentages,'PercData')],
            'py_funcs_all_results': [plotSumPie,plotMCIPie,plotPercCorrelations],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname':'Oligo Id',
            'partitions': ['Real Guides'],
            'samples': ['K562 New']
            }
    analyseResultsPerPartition( spec ) 

    spec = {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(computePieData, 'PieData')],
            'py_funcs_all_results': [plotBarSummaryPieIndels],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname':'Oligo Id',
            'partitions': ['Real Guides'],
            'samples': ['DPI7']
            }
    analyseResultsPerPartition( spec ) 

if __name__ == '__main__':
    setHighDataDir('..')
    runAnalysis()
    import pdb; pdb.set_trace()