import io, os, csv, sys, re, random
import numpy as np
import pylab as PL
import pandas as pd

from scipy.stats import pearsonr, spearmanr

from selftarget.util import getPlotDir
from selftarget.data import shortDirLabel, getSimpleName, parseSampleName

FIG_TYPE = 'both'

def setFigType(val):
    global FIG_TYPE
    FIG_TYPE = val

def saveFig(fig_fileprefix, bbox=True):
    PL.rcParams['svg.fonttype'] = 'none'
    for ftype in ['svg','png']:
        if FIG_TYPE == 'both' or FIG_TYPE==ftype:
            if bbox:
                PL.savefig(getPlotDir() + '/' + fig_fileprefix +'.'+ftype, bbox_inches='tight')
            else:
                PL.savefig(getPlotDir() + '/' + fig_fileprefix +'.'+ftype) 

def sanitizeLabel(label):
    for badchar in [' ','(',')','>']:
        label = label.replace(badchar,'_')
    return label

def avPieSummaries(pie_summaries):
    pie_labels = pie_summaries[0][1]
    pie_datas =[x[0] for x in pie_summaries]
    pie_data = {x: np.mean([y[x] for y in pie_datas]) for x in pie_labels}
    return pie_data, pie_labels, np.mean([x[2] for x in pie_summaries])

def plotBarSummary(all_result_outputs, label='', data_label='PieData', plot_label='bar_plots', stacked=False, combine_reps=False, colors=['C0','C1','C2','C3','C4','C5','C6','C7','C8'], legcol=1, figsize=(6,4), cell_line_order = []):
    summaries = [(x[0][data_label], x[1]) for x in all_result_outputs]
    mapping = {'BOB':'Human iPSC','E14TG2A':'Mouse ESC'}
    if combine_reps: 
        combined_summaries = []
        for cell_line in cell_line_order:
            cell_line_summaries = [x[0] for x in summaries if (parseSampleName(x[1])[0] == cell_line)]
            combined_summaries.append((avPieSummaries(cell_line_summaries),(cell_line if cell_line not in mapping else mapping[cell_line])))
        summaries = combined_summaries

    PL.figure(figsize=figsize)
    pie_labels = summaries[0][0][1]
    N, M = len(pie_labels), len(summaries)
    width = 0.8 if stacked else 0.8/N
    bottoms = np.array([0.0] * M)
    for i, pie_label in enumerate(pie_labels):
        bar_heights = [x[0][0][pie_label] for x in summaries]
        cell_lines = [parseSampleName(x[1])[0] for x in summaries]
        if combine_reps or len(cell_line_order)==0:
            bar_pos = [i*width*int(not stacked)+j for j in np.arange(M)]
        else:
            bar_pos, prev_cl, xticks, xlabels, ncl = [-1.1*width], cell_lines[0], [], [], 0
            for cl in cell_lines:
                if cl != prev_cl: 
                    bar_pos.append(bar_pos[-1] + width*1.5)
                    xticks.append((bar_pos[-1]+bar_pos[-ncl])*0.5)
                    xlabels.append(mapping[prev_cl] if prev_cl in mapping else prev_cl)
                    ncl = 0
                else: bar_pos.append(bar_pos[-1] + width*1.1)
                prev_cl = cl
                ncl += 1
            xticks.append((bar_pos[-1]+bar_pos[-2]-width*0.4)*0.5)
            xlabels.append(mapping[prev_cl] if prev_cl in mapping else prev_cl)
            bar_pos = bar_pos[1:]
        print(pie_label,bar_heights)
        PL.bar(bar_pos,bar_heights,width,bottom=bottoms, label=pie_label, color=colors[i])
        if stacked:
            bottoms += np.array(bar_heights)
    PL.legend(loc='center right', ncol=legcol)
    #PL.title(label)
    if combine_reps:
        PL.xticks([x + N/2*width*int(not stacked) for x in np.arange(M)], [x[1] for x in summaries], rotation='vertical')
    elif len(cell_line_order)==0:
        PL.xticks([x + N/2*width*int(not stacked) for x in np.arange(M)], ['%s' % (getSimpleName(x[1],include_dpi=True) if not combine_reps else x[1]) for x in summaries], rotation='vertical')
    else:
        PL.xticks(xticks, xlabels, rotation='vertical')
    PL.xlim((-1,M*1.6))
    PL.subplots_adjust(left=0.15,right=0.95,top=0.95, bottom=0.25)
    PL.ylabel('Percent Mutated Reads')
    PL.show(block=False) 
    saveFig(plot_label)

def plotBoxPlotSummary(all_result_outputs, label='', data_label='', y_label='', plot_label='', cl_order=[]):
    data_values = [x[0][data_label][0].values for x in all_result_outputs]
    #sample_names = [getSimpleName(x[1]) + '\n(Median reads = %d)' % x[0][data_label][1] for x in all_result_outputs]
    sample_names = [getSimpleName(x[1]) for x in all_result_outputs]
    if len(cl_order)>0:
        cell_lines = [' '.join(x.split()[:-2]) for x in sample_names]
        print(cell_lines)
        reordered_data, reordered_sample_names = [],[]
        for cell_line in cl_order:
            for i, cline in enumerate(cell_lines):
                if cline == cell_line:
                      reordered_data.append(data_values[i])
                      reordered_sample_names.append(sample_names[i])
        sample_names = reordered_sample_names
        data_values = reordered_data

    PL.figure(figsize=(5,5))
    for i,dvs in enumerate(data_values):
        print(np.median(dvs))
        PL.boxplot([dvs], positions=[i], showfliers=True, sym='.', widths=0.8)
    PL.xticks(range(len(sample_names)), sample_names, rotation='vertical')
    PL.xlim((-0.5,len(sample_names)-0.5))
    PL.ylim((0,5))
    PL.ylabel(y_label)   
    PL.title(label)
    PL.subplots_adjust(bottom=0.3)
    PL.show(block=False)
    saveFig( '%s_%s' % (plot_label, sanitizeLabel(label))) 

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

                PL.plot([0,100],[0,100],'k--')
                PL.xlabel('K562 Replicate A')
                PL.ylabel('K562 Replicate B')
                if add_leg: PL.legend()
                PL.title('%s (%.2f)' % (y_label, pcorrs[i,j]))
    if plot_scatters: 
        PL.subplots_adjust(left=0.05,right=0.95,top=0.9, bottom=0.1, hspace=0.4)
        PL.show(block=False)
        saveFig(plot_label)

    if not plot_scatters:
        PL.figure(figsize=(10,10))
        #PL.subplot(1,2,1)
        PL.imshow(pcorrs, cmap='hot', vmin = 0.0, vmax = 1.0, interpolation='nearest')
        PL.xticks(range(N),sample_names, rotation='vertical')
        PL.yticks(range(N),sample_names)
        PL.title(y_label + ': Pearson')
        PL.colorbar()
        #PL.subplot(1,2,2)
        #PL.imshow(scorrs, cmap='hot', vmin = 0.0, vmax = 1.0, interpolation='nearest')
        #PL.xticks(range(N),sample_names, rotation='vertical')
        #PL.yticks(range(N),sample_names)
        #PL.title(y_label + ': Spearman')
        #PL.colorbar()
        PL.show(block=False)
        PL.savefig(getPlotDir() + '/%s_%s.png' % (plot_label, sanitizeLabel(label)), bbox_inches='tight') 