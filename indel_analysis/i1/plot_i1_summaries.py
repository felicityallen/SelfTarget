import io, os, sys, csv
import pandas as pd
import numpy as np
import pylab as PL

from selftarget.util import getPlotDir, analyseResultsPerPartition, defaultLoadData, mergeSamples
from selftarget.plot import saveFig, plotBarSummary, plotVerticalHistSummary
from selftarget.oligo import getShortOligoId
from selftarget.data import getHighDataDir, setHighDataDir

    
MIN_READS = 20

def loadI1andMCIData(results_file, guideset=set()):
    data = pd.read_csv(results_file[0], sep='\t')
    mci_data = pd.read_csv(results_file[1], sep='\t')
    po_mci_data =  mci_data.groupby('Oligo Id')[['Most Common Indel', 'MCI Type']].sum()    #Note: concatenates strings if Oligos have multiple most common indels, so when matching replicates, has to be the same for the other replicates (unlikely)
    po_mci_data['Oligo Id'] = po_mci_data.index
    rd_mci_data = mci_data.groupby('Oligo Id')[['MCI Reads']].mean()
    rd_mci_data['Oligo Id'] = rd_mci_data.index
    merged_data = pd.merge(data, po_mci_data, on='Oligo Id', how='inner')
    merged_rd_data = pd.merge(merged_data, rd_mci_data, on='Oligo Id', how='inner')
    if len(guideset) > 0: merged_rd_data = merged_rd_data.loc[merged_rd_data['Oligo Id'].isin(guideset)] 

    return merged_rd_data

def loadIndelData():
    indel_data_new = pd.read_csv(getHighDataDir() + '/i1/exp_target_pam_new_gen_i1_indels.txt', sep='\t', header=1)
    indel_data_old = pd.read_csv(getHighDataDir() + '/i1/exp_target_pam_old_gen_i1_indels.txt', sep='\t', header=1)
    indel_data = pd.concat([indel_data_new, indel_data_old])[['Oligo Id','Repeat Nucleotide Left','Repeat Nucleotide Right']]
    indel_data['Short Oligo Id'] = indel_data['Oligo Id'].apply(getShortOligoId)
    return indel_data


def mergeWithIndelData(data, label=''):
    indel_data = loadIndelData()
    indel_data = indel_data[[x for x in indel_data.columns if x != 'Oligo Id']]
    merged_data = pd.merge(data, indel_data, how='inner', left_on='Oligo Id', right_on='Short Oligo Id')
    merged_data['Is Ambiguous'] = (merged_data['Repeat Nucleotide Left'] == merged_data['Repeat Nucleotide Right'])
    merged_data['I1_Rpt Left Reads - NonAmb'] = merged_data['I1_Rpt Left Reads']*(1-merged_data['Is Ambiguous'])
    merged_data['I1_Rpt Right Reads - NonAmb'] = merged_data['I1_Rpt Right Reads']*(1-merged_data['Is Ambiguous'])
    merged_data['Ambiguous Rpt Reads'] = merged_data['Is Ambiguous']*(merged_data['I1_Rpt Left Reads'] + merged_data['I1_Rpt Right Reads'])
    pie_labels = ['I1_Rpt Left Reads - NonAmb','I1_Rpt Right Reads - NonAmb','Ambiguous Rpt Reads','I1_NonRpt Reads']
    merged_data['I1 Total'] = merged_data[pie_labels].sum(axis=1)
    return merged_data

def computePieData(data,label=''):
    pie_labels = ['I1_Rpt Left Reads','I1_Rpt Right Reads','I1_NonRpt Reads']
    pie_data = {x: (data[x]*100.0/data['Total reads']).mean(axis=0) for x in pie_labels}
    return pie_data, pie_labels, data['Total reads'].median()

def computePieDataWithAmbig(data,label='', norm='I1 Total'):
    merged_data = mergeWithIndelData(data)
    pie_labels = ['I1_Rpt Left Reads - NonAmb','Ambiguous Rpt Reads','I1_Rpt Right Reads - NonAmb','I1_NonRpt Reads']
    labels = ['Repeated\nleft nucleotide', 'Ambiguous\n(Left = Right)', 'Repeated\nright nucleotide', 'Non-repeated\nnucleotide']
    pie_data = {x: (merged_data[x]*100.0/merged_data[norm]).mean(axis=0) for x in pie_labels}
    PL.figure(figsize=(3,3))
    PL.pie([pie_data[x] for x in pie_labels], labels=labels, autopct='%.1f', labeldistance=1.05, startangle=120.0, counterclock=False)
    PL.title('Single nucleotide insertions (I1)')
    PL.show(block=False)
    saveFig('ambig_pie_%s' % label)
    return pie_data, pie_labels, data['Total reads'].median()

def i1RepeatNucleotides(data, label=''):
    merged_data = mergeWithIndelData(data)
    nt_mean_percs, nts = [], ['A','T','G','C']
    for nt in nts:
        nt_data  = merged_data.loc[merged_data['Repeat Nucleotide Left'] == nt]  
        nt_mean_percs.append((nt_data['I1_Rpt Left Reads - NonAmb']*100.0/nt_data['Total reads']).mean())
    PL.figure(figsize=(3,3))
    PL.bar(range(4),nt_mean_percs)
    for i in range(4):
        PL.text(i-0.25,nt_mean_percs[i]+0.8,'%.1f' % nt_mean_percs[i])
    PL.xticks(range(4),nts)
    PL.ylim((0,26))
    PL.xlabel('PAM distal nucleotide\nadjacent to the cut site')
    PL.ylabel('I1 repeated left nucleotide\nas percent of total mutated reads')
    PL.show(block=False)
    saveFig('i1_rtp_nt_%s' % label)

def computeI1Repeats(data, label=''):
    merged_data = mergeWithIndelData(data)
    return (merged_data['I1_Rpt Left Reads - NonAmb'] + merged_data['Ambiguous Rpt Reads'])*100.0/merged_data['Total reads'], merged_data['Total reads'].median()

def computeFractionWithI1Repeats(data, label=''):
    merged_data = mergeWithIndelData(data)
    pie_label = 'Oligos with I1 Repeats'
    perc_with_11Rpt = len(merged_data.loc[(merged_data['I1_Rpt Left Reads - NonAmb'] + merged_data['Ambiguous Rpt Reads'])> 0.0])*100.0/len(merged_data)
    pie_data = {pie_label:perc_with_11Rpt}
    PL.figure()
    PL.pie([perc_with_11Rpt, 1-perc_with_11Rpt], labels=[pie_label,'Oligos without I1 Repeats'], autopct='%.1f', labeldistance=1.05, startangle=90.0, counterclock=False)
    PL.title(label)
    PL.show(block=False)
    return pie_data, [pie_label], merged_data['Total reads'].median()

def computePieDataNorm(data,label=''):
    pie_labels = ['I1_Rpt Left Reads','I1_Rpt Right Reads','I1_NonRpt Reads']
    data['I1 reads'] = data[pie_labels].sum(axis=1)
    pie_data = {x + ' Norm': (data[x]*100.0/data['I1 reads']).mean(axis=0) for x in pie_labels}
    return pie_data, [x + ' Norm' for x in pie_labels], data['Total reads'].median()

def computePieDataI1to3(data,label=''):
    pie_labels = []
    for ilen in [1,2,3]:
        pie_labels.extend(['I%d_Rpt Left Reads' % ilen,'I%d_Rpt Right Reads' % ilen,'I%d_NonRpt Reads'  % ilen])
    pie_data = {x: (data[x]*100.0/data['Total reads']).mean(axis=0) for x in pie_labels}
    return pie_data, pie_labels, data['Total reads'].median()

def computePieDataI1to3Norm(data,label=''):
    pie_labels = []
    for ilen in [1,2,3]:
        pie_labels.extend(['I%d_Rpt Left Reads' % ilen,'I%d_Rpt Right Reads' % ilen,'I%d_NonRpt Reads'  % ilen])
    data['I1 reads'] = data[pie_labels].sum(axis=1)
    pie_data = {x + ' Norm': (data[x]*100.0/data['I1 reads']).mean(axis=0) for x in pie_labels}
    return pie_data, [x + ' Norm' for x in pie_labels], data['Total reads'].median()

def plotMergedI1Repeats(all_result_outputs, label=''):
    merged_data = mergeSamples(all_result_outputs, ['I1_Rpt Left Reads - NonAmb','Total reads'], data_label='i1IndelData', merge_on=['Oligo Id','Repeat Nucleotide Left'])
    nt_mean_percs, nts = [], ['A','T','G','C']
    for nt in nts:
        nt_data  = merged_data.loc[merged_data['Repeat Nucleotide Left'] == nt]  
        nt_mean_percs.append((nt_data['I1_Rpt Left Reads - NonAmb Sum']*100.0/nt_data['Total reads Sum']).mean())
    PL.figure(figsize=(3,3))
    PL.bar(range(4),nt_mean_percs)
    for i in range(4):
        PL.text(i-0.25,nt_mean_percs[i]+0.8,'%.1f' % nt_mean_percs[i])
    PL.xticks(range(4),nts)
    PL.ylim((0,26))
    PL.xlabel('PAM distal nucleotide\nadjacent to the cut site')
    PL.ylabel('I1 repeated left nucleotide\nas percent of total mutated reads')
    PL.show(block=False)
    saveFig('i1_rtp_nt')

def plotMergedPieDataWithAmbig(all_result_outputs, label='', norm='I1 Total'):
    pie_labels = ['I1_Rpt Left Reads - NonAmb','Ambiguous Rpt Reads','I1_Rpt Right Reads - NonAmb','I1_NonRpt Reads']
    merged_data = mergeSamples(all_result_outputs, pie_labels + [norm], data_label='i1IndelData', merge_on=['Oligo Id'])
    labels = ['Repeated\nleft nucleotide', 'Ambiguous\n(Left = Right)', 'Repeated\nright nucleotide', 'Non-repeated\nnucleotide']
    pie_data = {x: (merged_data[x + ' Sum']*100.0/merged_data[norm + ' Sum']).mean(axis=0) for x in pie_labels}
    PL.figure(figsize=(3,3))
    PL.pie([pie_data[x] for x in pie_labels], labels=labels, autopct='%.1f', labeldistance=1.05, startangle=120.0, counterclock=False)
    PL.title('Single nucleotide insertions (I1)')
    PL.show(block=False)
    saveFig('ambig_pie')

def plotDominantPieDataWithAmbig(all_result_outputs, label=''):
    pie_labels = ['I1_Rpt Left Reads - NonAmb','Ambiguous Rpt Reads','I1_Rpt Right Reads - NonAmb','I1_NonRpt Reads']
    mci_merged_data = mergeSamples(all_result_outputs, [], data_label='i1IndelData')
    mci_merged_data['Equal MCI'] = (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 2']) & (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 3'])
    mci_common_i1 = mci_merged_data.loc[mci_merged_data['Equal MCI'] & (mci_merged_data['MCI Type'] == 'I1')]
    labels = ['Repeated\nleft nucleotide', 'Ambiguous\n(Left = Right)', 'Repeated\nright nucleotide', 'Non-repeated\nnucleotide']
    pie_data = []
    for label in pie_labels:
        is_rpt = lambda row: row['MCI Reads'] == row[label]
        pie_data.append(sum(mci_common_i1.apply(is_rpt,axis=1))*100.0/len(mci_common_i1))
    PL.figure(figsize=(3,3))
    PL.pie(pie_data, labels=labels, autopct='%.1f', labeldistance=1.05, startangle=180.0, counterclock=False)
    PL.title('Dominant single nucleotide insertions (I1)\n%d from %d gRNAs' % (len(mci_common_i1), len(mci_merged_data)))
    PL.show(block=False)
    saveFig('I1_dom_pie_3_rep')

def plotDominantBars(all_result_outputs, label=''):
    pie_labels = ['I1_Rpt Left Reads - NonAmb','Ambiguous Rpt Reads','I1_Rpt Right Reads - NonAmb','I1_NonRpt Reads']
    mci_merged_data = mergeSamples(all_result_outputs, [], data_label='i1IndelData')
    mci_merged_data['Equal MCI'] = (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 2']) & (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 3'])
    mci_merged_data['Is Dominant I1'] = (mci_merged_data['Equal MCI'] & (mci_merged_data['MCI Type'] == 'I1'))
    
    oligo_data =  pd.read_csv(getHighDataDir() + '/ST_June_2017/data/self_target_oligos_details_with_pam_details.csv',sep='\t')
    remove_under = lambda x: x.replace('_','')
    oligo_data['Oligo Id'] = oligo_data['ID'].apply(remove_under)
    merged_mci_data = pd.merge(mci_merged_data, oligo_data[['Oligo Id','Guide']], how='inner',on='Oligo Id')

    nt_perc_i1, cnt_labels = [], []
    nts = 'ATGC'
    for nt in nts:
        is_nt = lambda guide: (guide[-4] == nt)
        nt_data = merged_mci_data.loc[merged_mci_data['Guide'].apply(is_nt)]
        nt_perc_i1.append(sum(nt_data['Is Dominant I1'])*100.0/len(nt_data))
        cnt_labels.append('%d/%d' % (sum(nt_data['Is Dominant I1']),  len(nt_data)))
    
    PL.figure()
    PL.bar(range(4), nt_perc_i1, width=0.8)
    for i, cnt in enumerate(cnt_labels):
        PL.text(i-0.3,nt_perc_i1[i]+5.0,cnt)
    PL.xticks(range(4), [x for x in nts])
    PL.xlabel('Nucleotide on Left of cut-site')
    PL.ylabel('Percent gRNAs with single nucleotide insertion\nas most common indel in all 3 replicates')
    PL.show(block=False)
    saveFig('I1_bar_3_rep')
    
def plotVertHistI1(all_result_outputs, label=''):
    plotVerticalHistSummary(all_result_outputs,label=label,data_label='I1RptPercs',plot_label='i1_rpt_verthist', hist_width=5000, hist_bins=100)

def plotBarSummaryI1AmbIndels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_i1_ambig',  data_label='AmbigPieData', stacked=True)

def plotBarSummaryI1Indels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_i1',  data_label='PieData', stacked=True)

def plotBarSummaryNormI1Indels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_norm_i1',  data_label='PieDataNorm', stacked=True)

def plotBarSummaryI1to3Indels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_I1to3',  data_label='PieDataI1to3', stacked=True)

def plotBarSummaryNormI1to3Indels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_norm_I1to3',  data_label='PieDataI1to3Norm', stacked=True)

def plotBarSummaryI1RptFracs(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_rpt_fracs',  data_label='FracWithI1Rpt', stacked=True)

def runAnalysis():
	
    spec = {'results_specs': [{'results_dir':getHighDataDir() + '/i1/i1_summaries',
                              'dirname_to_result_fn': lambda x: '%s.txt' % x,
                              'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4]},
                            {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
                              'dirname_to_result_fn': lambda x: '%s.txt' % x,
                              'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4]}],
            'py_func_load': loadI1andMCIData,
            'py_funcs_per_result': [(mergeWithIndelData, 'i1IndelData')], 
            'py_funcs_all_results': [plotDominantBars,plotDominantPieDataWithAmbig,plotMergedPieDataWithAmbig, plotMergedI1Repeats],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname': 'Oligo Id',
            'partitions': ['Non-Targeting'],
            'samples': ['K562 New']
            }
    analyseResultsPerPartition( spec ) 
   
if __name__ == '__main__':
    setHighDataDir('..')
    runAnalysis()
    import pdb; pdb.set_trace()