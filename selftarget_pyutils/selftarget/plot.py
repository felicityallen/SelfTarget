import io, os, csv, sys, re, random
import numpy as np
import pylab as PL
import pandas as pd

from scipy.stats import pearsonr, spearmanr

from selftarget.util import getPlotDir
from selftarget.data import shortDirLabel, getSimpleName, parseSampleName

def saveFig(fig_fileprefix):
    PL.rcParams['svg.fonttype'] = 'none'
    PL.savefig(getPlotDir() + '/' + fig_fileprefix +'.svg', bbox_inches='tight')  
    PL.savefig(getPlotDir()  + '/' + fig_fileprefix +'.png', bbox_inches='tight')   

def sanitizeLabel(label):
    for badchar in [' ','(',')','>']:
        label = label.replace(badchar,'_')
    return label

def avPieSummaries(pie_summaries):
    pie_labels = pie_summaries[0][1]
    pie_datas =[x[0] for x in pie_summaries]
    pie_data = {x: np.mean([y[x] for y in pie_datas]) for x in pie_labels}
    return pie_data, pie_labels, np.mean([x[2] for x in pie_summaries])

def plotBarSummary(all_result_outputs, label='', data_label='PieData', plot_label='bar_plots', stacked=False, combine_reps=False, colors=['C0','C1','C2','C3','C4','C5','C6','C7','C8']):
    summaries = [(x[0][data_label], x[1]) for x in all_result_outputs]
    if combine_reps: 
        combined_summaries = []
        cell_lines = set([parseSampleName(x[1])[0] for x in summaries])
        #cell_lines = ['K562','CHO','HAP1','RPE1','BOB', 'E14TG2A','CHO','TREX2','eCAS9']
        cell_lines = ['CHO', 'E14TG2A', 'BOB','RPE1', 'HAP1','K562','eCAS9','TREX2']
        for cell_line in cell_lines:
            cell_line_summaries = [x[0] for x in summaries if (parseSampleName(x[1])[0] == cell_line)]
            mapping = {'BOB':'Human iPSC','E14TG2A':'Mouse ESC'}
            combined_summaries.append((avPieSummaries(cell_line_summaries),(cell_line if cell_line not in mapping else mapping[cell_line])))
        summaries = combined_summaries

    PL.figure(figsize=(6,4))
    pie_labels = summaries[0][0][1]
    N, M = len(pie_labels), len(summaries)
    width = 0.8 if stacked else 0.8/N
    bottoms = np.array([0.0] * M)
    for i, pie_label in enumerate(pie_labels):
        bar_pos = [i*width*int(not stacked)+j for j in np.arange(M)]
        bar_heights = [x[0][0][pie_label] for x in summaries]
        print(pie_label,bar_heights)
        PL.bar(bar_pos,bar_heights,width,bottom=bottoms, label=pie_label, color=colors[i])
        if stacked:
            bottoms += np.array(bar_heights)
    PL.legend(loc='center right')
    #PL.title(label)
    PL.xticks([x + N/2*width*int(not stacked) for x in np.arange(M)], ['%s' % (getSimpleName(x[1]) if not combine_reps else x[1]) for x in summaries], rotation='vertical')
    PL.xlim((-1,M*1.6))
    PL.subplots_adjust(left=0.15,right=0.95,top=0.95, bottom=0.25)
    PL.ylabel('Percent Mutated Reads')
    PL.show(block=False) 
    saveFig(plot_label)

def plotBoxPlotSummary(all_result_outputs, label='', data_label='', y_label='', plot_label=''):
   
    data_values = [x[0][data_label][0].values for x in all_result_outputs]
    sample_names = [shortDirLabel(x[1]) + ' (%d)' % x[0][data_label][1] for x in all_result_outputs]

    PL.figure(figsize=(12,8))
    for i,dvs in enumerate(data_values):
        PL.boxplot([dvs], positions=[i], showfliers=True, sym='.', widths=0.8)
    PL.xticks(range(len(sample_names)), sample_names, rotation='vertical')
    PL.ylabel(y_label)   
    PL.title(label)
    PL.show(block=False)
    PL.savefig(getPlotDir() + '/%s_%s.png' % (plot_label, sanitizeLabel(label)), bbox_inches='tight') 

def plotVerticalHistSummary(all_result_outputs, label='', data_label='', y_label='', plot_label='', hist_width=1000, hist_bins=100, oligo_id_str='Oligo ID', val_str = 'Cut Rate', total_reads_str= 'Total Reads'):
   
    datas = [x[0][data_label][0] for x in all_result_outputs]
    sample_names = [shortDirLabel(x[1]) for x in all_result_outputs]
    
    merged_data = pd.merge(datas[0],datas[1],how='inner',on=oligo_id_str, suffixes=['', ' 2'])
    for i, data in enumerate(datas[2:]):
        merged_data = pd.merge(merged_data, data,how='inner',on=oligo_id_str, suffixes=['', ' %d' % (i+3)])
    suffix = lambda i: ' %d' % (i+1) if i > 0 else ''

    xpos = [x*hist_width for x in range(len(sample_names))]

    PL.figure(figsize=(12,8))
    for i,label1 in enumerate(sample_names):
        dvs = merged_data[val_str + suffix(i)]
        PL.hist(dvs, bins=hist_bins, bottom=i*hist_width, orientation='horizontal')
    PL.xticks(xpos, sample_names, rotation='vertical')
    PL.ylabel(y_label)   
    PL.title(label)
    PL.show(block=False)
    PL.savefig(getPlotDir() + '/%s_%s.png' % (plot_label, label.replace(' ','_')), bbox_inches='tight') 

def plotCorrelations(all_result_outputs, label='', data_label='', y_label='', plot_label='', plot_scatters=False, oligo_id_str='Oligo ID', val_str = 'Cut Rate', total_reads_str= 'Total Reads', scatter_samples={}, sdims=(0,0), scatter_fig=None, add_leg=True):
    
    datas = [x[0][data_label][0] for x in all_result_outputs]
    sample_names = [shortDirLabel(x[1]) for x in all_result_outputs]
    
    merged_data = pd.merge(datas[0],datas[1],how='inner',on=oligo_id_str, suffixes=['', ' 2'])
    for i, data in enumerate(datas[2:]):
        merged_data = pd.merge(merged_data, data,how='inner',on=oligo_id_str, suffixes=['', ' %d' % (i+3)])
    suffix = lambda i: ' %d' % (i+1) if i > 0 else ''

    N = len(sample_names)
    if plot_scatters: 
        if scatter_fig is None: PL.figure()
        else: PL.figure(scatter_fig.number)
        s_dims = (N,N) if len(scatter_samples) == 0 else sdims
    pcorrs, scorrs = np.zeros((N,N)), np.zeros((N,N))
    for i,label1 in enumerate(sample_names):
        for j,label2 in enumerate(sample_names):
            dvs1, dvs2, ids = merged_data[val_str + suffix(i)], merged_data[val_str + suffix(j)], merged_data[oligo_id_str]
            pcorrs[i,j] = pearsonr(dvs1, dvs2)[0]
            scorrs[i,j] = spearmanr(dvs1, dvs2)[0]
            if plot_scatters:
                if (label1, label2) in scatter_samples: idx = scatter_samples[(label1, label2)]
                elif len(scatter_samples) == 0: idx = i*N+j+1
                else: continue
                PL.subplot(s_dims[0],s_dims[1],idx)
                trs1, trs2 = merged_data[total_reads_str + suffix(i)], merged_data[total_reads_str + suffix(j)]
                for thr in [20,50,100,500,1000]:
                    thr_dvs = [(dv1, dv2,id) for (dv1, dv2, tr1, tr2,id) in zip(dvs1,dvs2,trs1,trs2,ids) if (tr1 >= thr and tr2 >= thr) ]
                    pcorrs[i,j] = pearsonr([x[0] for x in thr_dvs], [x[1] for x in thr_dvs])[0]
                    PL.plot([x[0] for x in thr_dvs],[x[1] for x in thr_dvs],'.', label='>%d Reads' % thr)
                    #if thr == 1000:
                        #for (x,y,id) in thr_dvs:
                            #if random.random() < 0.01:
                    #            PL.text(x,y,id)
                    

                PL.plot([0,100],[0,100],'k--')
                PL.xlabel('K562 Replicate A')
                PL.ylabel('K562 Replicate B')
                if add_leg: PL.legend()
                PL.title('%s (%.2f)' % (y_label, pcorrs[i,j]))
    if plot_scatters: 
        PL.subplots_adjust(left=0.05,right=0.95,top=0.9, bottom=0.1, hspace=0.4)
        PL.show(block=False)
        saveFig(plot_label)
        #PL.savefig(getPlotDir() + '/scatter_%s_%s.png' % (plot_label, sanitizeLabel(label)), bbox_inches='tight') 
    PL.figure()
    PL.subplot(1,2,1)
    PL.imshow(pcorrs, cmap='hot', vmin = 0.0, vmax = 1.0, interpolation='nearest')
    PL.xticks(range(N),sample_names, rotation='vertical')
    PL.yticks(range(N),sample_names)
    PL.title(y_label + ': Pearson')
    PL.colorbar()
    PL.subplot(1,2,2)
    PL.imshow(scorrs, cmap='hot', vmin = 0.0, vmax = 1.0, interpolation='nearest')
    PL.xticks(range(N),sample_names, rotation='vertical')
    PL.yticks(range(N),sample_names)
    PL.title(y_label + ': Spearman')
    PL.colorbar()
    PL.show(block=False)
    PL.savefig(getPlotDir() + '/%s_%s.png' % (plot_label, sanitizeLabel(label)), bbox_inches='tight') 