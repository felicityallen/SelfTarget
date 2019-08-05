import io, os, sys, csv, random
import pandas as pd
import numpy as np
import pylab as PL

from scipy.stats import pearsonr, spearmanr
from sklearn.manifold import TSNE

from selftarget.util import getPlotDir, analyseResultsPerPartition, mergeSamples
from selftarget.plot import plotBarSummary, plotCorrelations, saveFig
from selftarget.data import parseSampleName, setHighDataDir, getHighDataDir
from selftarget.indel import tokFullIndel

OLD_COLORS = ['C0','C9','C3','C1','darkgreen', 'C4','C7']
COLORS = ['C7','C0','C9','C3','C1','darkgreen', 'C2','#DDDDDD']

MIN_READS = 20

DEL_LABELS = ['D1','D2','D>2, MH','D>2, No MH']
INS_LABELS = ['I1','I>1']
ALL_LABELS = ['I+D'] + INS_LABELS + DEL_LABELS

def loadData(results_file, guideset=set()):
    data = pd.read_csv(results_file, sep='\t')
    if len(guideset) > 0: data = data.loc[data['Oligo Id'].isin(guideset)] 
    data['MCI Ratio'] = data['MCI Reads']*100.0/data['Total reads']
    data['All Deletions'] = data[DEL_LABELS].sum(axis=1)
    data['All Insertions'] = data[INS_LABELS].sum(axis=1)
    return data

def passData(data, label=''):
    return data

def perOligoMCI(data, label=''):
    data['Altered Sequence'] = data['Altered Sequence'].fillna('')
    data['Microhomology Sequence'] = data['Microhomology Sequence'].fillna('')
    po_data =  data.groupby('Oligo Id')['MCI Type','Most Common Indel','Altered Sequence','Microhomology Sequence'].sum()    #Note: concatenates strings if Oligos have multiple most common indels, so when matching replicates, has to be the same for the other replicates (unlikely)
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
    merged_data =  mergeSamples(all_result_outputs, ['Total reads'] + ALL_LABELS, data_label='perOligoCounts')
    for col in ALL_LABELS:
        merged_data[col + ' Perc'] = merged_data[col + ' Sum']*100.0/merged_data['Total reads Sum']
    merged_data.to_csv('data_dump_indel_pie.txt',sep='\t',columns=['Oligo Id'] + [col + ' Perc' for col in ALL_LABELS])
    pie_vals = [merged_data[col + ' Perc'].mean() for col in ALL_LABELS]
    PL.figure(figsize=(4,4))
    PL.pie(pie_vals, labels=ALL_LABELS, autopct='%.1f', labeldistance=1.05, startangle=90.0, counterclock=False, colors=COLORS)
    PL.title('Average distribution\n of mutations\n per gRNA')
    PL.show(block=False)
    saveFig('pie_chart_cats')

def plotMCIPie(all_result_outputs, label=''):
    mci_merged_data = mergeSamples(all_result_outputs, ['MCI Type','Most Common Indel'], data_label='perOligoMCI')
    mci_common = mci_merged_data.loc[(mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 2']) & (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 3'])]
    pie_vals, pie_labels = [], []
    for mci_type in ALL_LABELS:
        pie_vals.append(len(mci_common.loc[mci_common['MCI Type'] == mci_type]))
        pie_labels.append(mci_type)
    pie_vals.append(len(mci_merged_data)-len(mci_common))
    pie_labels.append('Inconsistent\nbetween\nreplicates')

    PL.figure(figsize=(4,4))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.05, startangle=90.0, counterclock=False, colors=COLORS)
    PL.title('Most frequent\nmutation per gRNA')
    PL.show(block=False)
    saveFig('pie_chart_cats_dominant')

def plotD1(all_result_outputs, label=''):
    mci_merged_data = mergeSamples(all_result_outputs, [], data_label='perOligoMCI')
    mci_merged_data['Equal MCI'] = (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 2']) & (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 3'])
    mci_common = mci_merged_data.loc[mci_merged_data['Equal MCI']]
    pie_vals, pie_labels = [], []
    dmci_data = mci_common.loc[(mci_common['MCI Type'] == 'D1')] #Note: type check discards equally most common indels

    spans_cutsite = lambda indel: tokFullIndel(indel)[2]['L'] < -1 and tokFullIndel(indel)[2]['R'] > 0
    for nt in 'ATGC':
        is_mh = lambda alt_seq: len(alt_seq) >= 2 and alt_seq == (len(alt_seq)*nt)
        num_repeat_nt = len(dmci_data.loc[dmci_data['Altered Sequence'].apply(is_mh) & dmci_data['Most Common Indel'].apply(spans_cutsite)])
        pie_vals.append(num_repeat_nt*100.0/len(dmci_data))
        print(num_repeat_nt)
        pie_labels.append('Removal of %s\nfrom %s|%s' % (nt,nt,nt))
    is_non_repeat = lambda seq: len(seq) < 2 or seq != (seq[0]*len(seq))
    num_non_repeat  = len(dmci_data.loc[dmci_data['Altered Sequence'].apply(is_non_repeat) | ~dmci_data['Most Common Indel'].apply(spans_cutsite)])
    pie_vals.append(num_non_repeat*100.0/len(dmci_data))
    print(num_non_repeat)
    pie_labels.append('Removal from non-repeat')
    PL.figure(figsize=(4,4))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=OLD_COLORS)
    PL.title('Size 1 deletions that are\n"most common" for their gRNA in all 3 replicates\n(%d gRNAs from %d total)' % (len(dmci_data), len(mci_merged_data)))
    PL.show(block=False)
    saveFig('pie_chart_D1')
    

    oligo_data =  pd.read_csv(getHighDataDir() + '/ST_June_2017/data/self_target_oligos_details_with_pam_details.csv',sep='\t')
    remove_under = lambda x: x.replace('_','')
    oligo_data['Oligo Id'] = oligo_data['ID'].apply(remove_under)
    merged_mci_data = pd.merge(mci_merged_data, oligo_data[['Oligo Id','Guide']], how='inner',on='Oligo Id')
    print(len(merged_mci_data))

    nt_dbl_perc_d1, cnt_labels = [], []
    is_d1 = lambda indel: (indel.split('_')[0] == 'D1')
    non_dbl_nt = lambda row: row['Guide'][-4] != row['Guide'][-3]    
    nts = 'ATGC'
    for nt in nts:
        double_nt = lambda row: row['Guide'][-4:-2] == (nt+nt)
        dbl_data = merged_mci_data.loc[merged_mci_data.apply(double_nt,axis=1)]
        num_dbl_d1 = sum(dbl_data['Most Common Indel'].apply(is_d1) & dbl_data['Equal MCI'] & (dbl_data['Oligo Id']!='Oligo28137')) #Oligo28137: Corner case where a guide has CT|T and loses the C
        nt_dbl_perc_d1.append(num_dbl_d1*100.0/len(dbl_data))
        cnt_labels.append('%d/%d' % (num_dbl_d1,len(dbl_data)))
        print(len(dbl_data))
    non_dbl_data = merged_mci_data.loc[merged_mci_data.apply(non_dbl_nt,axis=1)]
    print(len(non_dbl_data))
    num_non_dbl_d1 = sum(non_dbl_data['Most Common Indel'].apply(is_d1) & non_dbl_data['Equal MCI'])
    nt_dbl_perc_d1.append(num_non_dbl_d1*100.0/len(non_dbl_data))
    cnt_labels.append('%d/%d' % (num_non_dbl_d1,len(non_dbl_data)))
    
    PL.figure()
    PL.bar(range(5), nt_dbl_perc_d1, width=0.8)
    for i, cnt in enumerate(cnt_labels):
        PL.text(i-0.3,nt_dbl_perc_d1[i]+5.0,cnt)
    PL.xticks(range(5), ['%s' % x*2 for x in nts] + ['Other'])
    PL.ylim((0,40))
    PL.xlabel('Nucleotides on either side of cut site')
    PL.ylabel('Percent gRNAs with single nucleotide deletion\nas most common indel in all 3 replicates')
    PL.show(block=False)
    saveFig('D1_bar_3_rep')


def plotD2(all_result_outputs, label=''):

    #Merge replicates
    mci_merged_data = mergeSamples(all_result_outputs, [], data_label='perOligoMCI')
    mci_merged_data['Equal MCI'] = (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 2']) & (mci_merged_data['Most Common Indel']==mci_merged_data['Most Common Indel 3'])

    oligo_data =  pd.read_csv(getHighDataDir() + '/ST_June_2017/data/self_target_oligos_details_with_pam_details.csv',sep='\t')
    remove_under = lambda x: x.replace('_','')
    oligo_data['Oligo Id'] = oligo_data['ID'].apply(remove_under)
    mci_merged_data_guides = pd.merge(mci_merged_data, oligo_data[['Oligo Id','Guide']], how='inner',on='Oligo Id')
    mci_common = mci_merged_data_guides.loc[mci_merged_data['Equal MCI']]
    dmci_data = mci_common.loc[(mci_common['MCI Type'] == 'D2')] #Note: type check discards equally most common indels

    pie_vals, pie_labels = [], []
    is_left_rpt = lambda row: row['Guide'][-5] == row['Guide'][-3] and tokFullIndel(row['Most Common Indel'])[2]['R'] >= 1 and tokFullIndel(row['Most Common Indel'])[2]['L'] <= -3
    is_right_rpt = lambda row: row['Guide'][-4] == row['Guide'][-2] and tokFullIndel(row['Most Common Indel'])[2]['R'] >= 2 and tokFullIndel(row['Most Common Indel'])[2]['L'] <= -2
    is_left_only_rpt = lambda row: is_left_rpt(row) and not is_right_rpt(row)
    is_right_only_rpt = lambda row: is_right_rpt(row) and not is_left_rpt(row)
    is_both_rpt = lambda row: is_right_rpt(row) and is_left_rpt(row)
    
    lrpt_data = dmci_data.loc[dmci_data.apply(is_left_only_rpt, axis=1)]
    pie_labels.append('Y|XY->Y'); pie_vals.append(len(lrpt_data))
    rrpt_data = dmci_data.loc[dmci_data.apply(is_right_only_rpt, axis=1)]
    pie_labels.append('XY|X->X'); pie_vals.append(len(rrpt_data))
    rpt_data = dmci_data.loc[dmci_data.apply(is_both_rpt, axis=1)]
    pie_labels.append('XY|XY->XY'); pie_vals.append(len(rpt_data))
    
    is_r0 = lambda row: tokFullIndel(row['Most Common Indel'])[2]['R'] == 0
    ro_data = dmci_data.loc[dmci_data.apply(is_r0, axis=1)]
    pie_labels.append('Z|XY->Z'); pie_vals.append(len(ro_data))
    is_l1 = lambda row: tokFullIndel(row['Most Common Indel'])[2]['L'] == -1
    l1_data = dmci_data.loc[dmci_data.apply(is_l1, axis=1)]
    pie_labels.append('XY|Z->Z'); pie_vals.append(len(l1_data))
    
    seen_ids = set(rpt_data['Oligo Id']).union(set(ro_data['Oligo Id'])).union(set(l1_data['Oligo Id'])).union(set(lrpt_data['Oligo Id'])).union(set(rrpt_data['Oligo Id']))
    is_unseen = lambda id: id not in seen_ids
    unseen_data = dmci_data.loc[dmci_data['Oligo Id'].apply(is_unseen)]
    print(unseen_data)
    assert(len(unseen_data) == 0)
    #pie_labels.append('Other')
    #pie_vals.append(len(unseen_data))

    #pie_labels = [x for x in dmci_data['Most Common Indel'].unique()]
    #pie_vals  = [len(dmci_data.loc[dmci_data['Most Common Indel']==indel]) for indel in pie_labels]
    PL.figure(figsize=(4,4))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=COLORS)
    PL.title('Size 2 deletions that are\n"most common" for their gRNA in all 3 replicates\n(%d gRNAs from %d total)' % (len(dmci_data), len(mci_merged_data)))
    PL.show(block=False)
    saveFig('pie_chart_D2_indel_cats')

    PL.figure(figsize=(12,8))

    #XY|XY->XY
    PL.subplot(2,3,1)
    pie_vals, pie_labels = [], []
    for mh_str in [y+x for x in 'ATGC' for y in 'ATGC']:
        pie_labels.append(mh_str)
        is_mh_str = lambda guide: guide[-5:-3] == mh_str
        pie_vals.append(len(rpt_data.loc[rpt_data['Guide'].apply(is_mh_str)]))
    for dnt, cnt in zip(pie_labels, pie_vals):
        print(dnt, cnt*100/sum(pie_vals))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=COLORS)
    PL.title('XY|XY->XY\n(%d gRNAs)' % len(rpt_data))
    PL.show(block=False)

    #__|
    PL.subplot(2,3,2)
    pie_vals, pie_labels = [], []
    for mh_str in [y+x for x in 'ATGC' for y in 'ATGC']:
        pie_labels.append(mh_str)
        is_mh_str = lambda guide: guide[-5:-3] == mh_str
        pie_vals.append(len(ro_data.loc[ro_data['Guide'].apply(is_mh_str)]))
    for dnt, cnt in zip(pie_labels, pie_vals):
        print(dnt, cnt*100/sum(pie_vals))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=COLORS)
    PL.title('XY| -> __|\n(%d gRNAs)'  % len(ro_data))
    PL.show(block=False)

    #|__
    PL.subplot(2,3,3)
    pie_vals, pie_labels = [], []
    for mh_str in [y+x for x in 'ATGC' for y in 'ATGC']:
        pie_labels.append(mh_str)
        is_mh_str = lambda guide: guide[-3:-1] == mh_str
        pie_vals.append(len(l1_data.loc[l1_data['Guide'].apply(is_mh_str)]))
    for dnt, cnt in zip(pie_labels, pie_vals):
        print(dnt, cnt*100/sum(pie_vals))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=COLORS)
    PL.title('|XY -> |__\n(%d gRNAs)'  % len(l1_data))
    PL.show(block=False)

    #XY|X->X
    PL.subplot(2,3,4)
    pie_vals, pie_labels = [], []
    for nt in 'ATGC':
        pie_labels.append('%sN|%s -> %s' % (nt,nt,nt))
        is_mh_str = lambda guide: guide[-5] == nt
        pie_vals.append(len(lrpt_data.loc[lrpt_data['Guide'].apply(is_mh_str)]))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=COLORS)
    PL.title('XY|X->X\n(%d gRNAs)'  % len(lrpt_data))
    PL.show(block=False)

    #X|YX->X
    PL.subplot(2,3,5)
    pie_vals, pie_labels = [], []
    for nt in 'ATGC':
        pie_labels.append('%s|N%s -> %s' % (nt,nt,nt))
        is_mh_str = lambda guide: guide[-4] == nt
        pie_vals.append(len(rrpt_data.loc[rrpt_data['Guide'].apply(is_mh_str)]))
    PL.pie(pie_vals, labels=pie_labels, autopct='%.1f', labeldistance=1.1, counterclock=False, colors=COLORS)
    PL.title('X|YX->X\n(%d gRNAs)'  % len(rrpt_data))
    PL.show(block=False)
    PL.subplots_adjust(left=0.05,right=0.95,top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)
    saveFig('D2_nts_per_cat')

    PL.figure(figsize=(12,8))

    #XY|XY->XY
    PL.subplot(2,3,1)
    bar_heights, bar_labels, d2_dnt_counts, dnt_counts = [], [], [], []
    for dnt in [y+x for x in 'ATGC' for y in 'ATGC']:
        has_dnt = lambda guide: guide[-5:-3] == dnt and guide[-3:-1] == dnt 
        dnt_data = mci_merged_data_guides.loc[mci_merged_data_guides['Guide'].apply(has_dnt)]
        dnt_counts.append(len(set(rpt_data['Oligo Id']).intersection(set(dnt_data['Oligo Id']))))
        d2_dnt_counts.append(len(dnt_data))
        bar_heights.append(dnt_counts[-1]*100.0/d2_dnt_counts[-1] if d2_dnt_counts[-1] > 0 else 0)
        bar_labels.append(dnt)
        print(dnt, dnt_counts[-1]*100.0/d2_dnt_counts[-1] if d2_dnt_counts[-1] > 0 else 0)
    PL.bar(range(len(bar_labels)), bar_heights, width=0.8)
    for i, (hgt, d2cnt, cnt) in enumerate(zip(bar_heights, d2_dnt_counts, dnt_counts)):
        PL.text(i-0.3,hgt+15,'%d/%d' % (cnt,d2cnt), rotation='vertical')
    PL.xticks(range(len(bar_labels)), bar_labels, rotation='vertical')
    PL.ylim((0,90))
    PL.xlabel('XY')
    PL.title('XY|XY->XY')
    PL.ylabel('Percent gRNAs with XY|XY->XY deletion\nas most common indel in all 3 replicates')
    PL.show(block=False)

    #__|
    PL.subplot(2,3,2)
    bar_heights, bar_labels, d2_dnt_counts, dnt_counts = [], [], [], []
    for dnt in [y+x for x in 'ATGC' for y in 'ATGC']:
        has_dnt = lambda guide: guide[-5:-3] == dnt 
        dnt_data = mci_merged_data_guides.loc[mci_merged_data_guides['Guide'].apply(has_dnt)]
        dnt_counts.append(len(set(ro_data['Oligo Id']).intersection(set(dnt_data['Oligo Id']))))
        d2_dnt_counts.append(len(dnt_data))
        bar_heights.append(dnt_counts[-1]*100.0/d2_dnt_counts[-1] if d2_dnt_counts[-1] > 0 else 0)
        bar_labels.append(dnt)
        print(dnt, dnt_counts[-1]*100.0/d2_dnt_counts[-1] if d2_dnt_counts[-1] > 0 else 0)
    PL.bar(range(len(bar_labels)), bar_heights, width=0.8)
    for i, (hgt, d2cnt, cnt) in enumerate(zip(bar_heights, d2_dnt_counts, dnt_counts)):
        PL.text(i-0.3,hgt+1.5,'%d/%d' % (cnt,d2cnt), rotation='vertical')
    PL.xticks(range(len(bar_labels)), bar_labels, rotation='vertical')
    PL.ylim((0,8))
    PL.xlabel('XY')
    PL.title('XY| -> __|')
    PL.ylabel('Percent gRNAs with XY| -> __| deletion\nas most common indel in all 3 replicates')
    PL.show(block=False)

    #|__
    PL.subplot(2,3,3)
    bar_heights, bar_labels, d2_dnt_counts, dnt_counts = [], [], [], []
    for dnt in [y+x for x in 'ATGC' for y in 'ATGC']:
        has_dnt = lambda guide: guide[-3:-1] == dnt 
        dnt_data = mci_merged_data_guides.loc[mci_merged_data_guides['Guide'].apply(has_dnt)]
        dnt_counts.append(len(set(l1_data['Oligo Id']).intersection(set(dnt_data['Oligo Id']))))
        d2_dnt_counts.append(len(dnt_data))
        bar_heights.append(dnt_counts[-1]*100.0/d2_dnt_counts[-1] if d2_dnt_counts[-1] > 0 else 0)
        bar_labels.append(dnt)
        print(dnt, dnt_counts[-1]*100.0/d2_dnt_counts[-1] if d2_dnt_counts[-1] > 0 else 0)
    PL.bar(range(len(bar_labels)), bar_heights, width=0.8)
    for i, (hgt, d2cnt, cnt) in enumerate(zip(bar_heights, d2_dnt_counts, dnt_counts)):
        PL.text(i-0.3,hgt+1.5,'%d/%d' % (cnt,d2cnt), rotation='vertical')
    PL.xticks(range(len(bar_labels)), bar_labels, rotation='vertical')
    PL.ylim((0,8))
    PL.xlabel('XY')
    PL.title('|XY -> |__')
    PL.ylabel('Percent gRNAs with |XY -> |__ deletion\nas most common indel in all 3 replicates')
    PL.show(block=False)

    #XY|X->X
    PL.subplot(2,3,4)
    bar_heights, bar_labels, d2_nt_counts, nt_counts = [], [], [], []
    for nt in 'ATGC':
        has_nt = lambda guide: guide[-3] == nt and guide[-5] == nt
        nt_data = mci_merged_data_guides.loc[mci_merged_data_guides['Guide'].apply(has_nt)]
        nt_counts.append(len(set(lrpt_data['Oligo Id']).intersection(set(nt_data['Oligo Id']))))
        d2_nt_counts.append(len(nt_data))
        bar_heights.append(nt_counts[-1]*100.0/d2_nt_counts[-1] if d2_nt_counts[-1] > 0 else 0)
        bar_labels.append(nt)
    PL.bar(range(len(bar_labels)), bar_heights, width=0.8)
    for i, (hgt, d2cnt, cnt) in enumerate(zip(bar_heights, d2_nt_counts, nt_counts)):
        PL.text(i-0.3,hgt+0.05,'%d/%d' % (cnt,d2cnt))
    PL.xticks(range(len(bar_labels)), bar_labels)
    PL.ylim((0,5))
    PL.xlabel('X')
    PL.title('XY|X->X')
    PL.ylabel('Percent gRNAs with XY|X->X deletion\nas most common indel in all 3 replicates')
    PL.show(block=False)

    #X|YX->X
    PL.subplot(2,3,5)
    bar_heights, bar_labels, d2_nt_counts, nt_counts = [], [], [], []
    for nt in 'ATGC':
        has_nt = lambda guide: guide[-4] == nt and guide[-2] == nt
        nt_data = mci_merged_data_guides.loc[mci_merged_data_guides['Guide'].apply(has_nt)]
        nt_counts.append(len(set(rrpt_data['Oligo Id']).intersection(set(nt_data['Oligo Id']))))
        d2_nt_counts.append(len(nt_data))
        bar_heights.append(nt_counts[-1]*100.0/d2_nt_counts[-1] if d2_nt_counts[-1] > 0 else 0)
        bar_labels.append(nt)
    PL.bar(range(len(bar_labels)), bar_heights, width=0.8)
    for i, (hgt, d2cnt, cnt) in enumerate(zip(bar_heights, d2_nt_counts, nt_counts)):
        PL.text(i-0.3,hgt+0.05,'%d/%d' % (cnt,d2cnt))
    PL.xticks(range(len(bar_labels)), bar_labels)
    PL.ylim((0,5))
    PL.xlabel('X')
    PL.title('X|YX->X')
    PL.ylabel('Percent gRNAs with X|YX->X deletion\nas most common indel in all 3 replicates')
    PL.show(block=False)

    PL.subplots_adjust(left=0.05,right=0.95,top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)
    saveFig('D2_nts_per_cat_bars')

def plotBarSummaryPieIndels(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_indel_bar_plots', stacked=True, combine_reps=True, colors=COLORS)

def plotBarSummaryPieIndelsNoCombineReps(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_indel_bar_plots_cl', stacked=True, combine_reps=False, colors=COLORS, cell_line_order=['CHO', 'E14TG2A', 'BOB','RPE1', 'HAP1','K562','eCAS9','TREX2'])

def plotBarSummaryPieIndelsNoCombineRepsTP(all_result_outputs, label=''):
    plotBarSummary(all_result_outputs, label=label, plot_label='pie_indel_bar_plots_tp', stacked=True, combine_reps=False, colors=COLORS)


def runAnalysis():
	
    spec = {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(perOligoCounts,'perOligoCounts'),(perOligoMCI,'perOligoMCI'),(computePercentages,'PercData')],
            'py_funcs_all_results': [plotSumPie,plotMCIPie],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname':'Oligo Id',
            'partitions': ['Real Guides'],
            'samples': ['K562 New']
            }
    #analyseResultsPerPartition( spec ) 

    spec = {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(perOligoMCI,'perOligoMCI')],
            'py_funcs_all_results': [plotD2, plotD1],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname':'Oligo Id',
            'partitions': ['Non-Targeting'],
            'samples': ['K562 New']
            }
    #analyseResultsPerPartition( spec ) 

    spec = {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(computePieData, 'PieData')],
            'py_funcs_all_results': [plotBarSummaryPieIndelsNoCombineReps],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname':'Oligo Id',
            'partitions': ['Real Guides'],
            'samples': ['DPI7']
            }
    analyseResultsPerPartition( spec ) 

    spec = {'results_dir':getHighDataDir() + '/indel_details/indel_pie_summaries_per_oligo',
            'dirname_to_result_fn': lambda x: '%s.txt' % x,
            'result_to_dirname_fn': lambda x: x.split('/')[-1][:-4],
            'py_func_load': loadData,
            'py_funcs_per_result': [(computePieData, 'PieData')],
            'py_funcs_all_results': [plotBarSummaryPieIndelsNoCombineRepsTP],
            'check_output_fn': lambda x: True,
            'reads_colname': 'Total reads',
            'min_reads': MIN_READS,
            'id_colname':'Oligo Id',
            'partitions': ['Non-Targeting'],
            'samples': ['K562_All_TP']
            }
    #analyseResultsPerPartition( spec ) 

if __name__ == '__main__':
    setHighDataDir('..')
    runAnalysis()
    import pdb; pdb.set_trace()