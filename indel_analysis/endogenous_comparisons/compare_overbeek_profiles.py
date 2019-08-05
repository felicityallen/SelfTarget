import io, os, csv, sys, re
import numpy as np
import pandas as pd
import pylab as PL
import Bio.Seq
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr

from selftarget.indel import tokFullIndel
from selftarget.profile import KL, symmetricKL,symmetricClassKLTopNIndels, classSymmetricKL,symmetricClassKLTop5Indels, fetchIndelSizeCounts, printCompareProfiles, readSummaryToProfile, fetchRepresentativeCleanReads, getProfileCounts
from selftarget.oligo import getFileForOligoIdx
from selftarget.view import plotProfiles
from selftarget.data import getNullTargetPamDetails, getHighDataDir, setHighDataDir
from selftarget.util import getPlotDir, setPlotDir
from selftarget.plot import saveFig, setFigType

def plotInFrame(overbeek_inframes, ours_inframes, oof_sel_overbeek_ids, pred_results_dir):

    PL.figure(figsize=(4.2,4.2))
    data = pd.read_csv(pred_results_dir + '/old_new_kl_predicted_summaries.txt', sep='\t').fillna(-1.0)
    label1, label2 = 'New 2x800x In Frame Perc', 'New 1600x In Frame Perc'
    xdata, ydata = data[label1], data[label2]
    PL.plot(xdata,ydata, '.', label='Synthetic between library (R=%.2f)' %  pearsonr(xdata,ydata)[0], color='C0',alpha=0.15)
    PL.plot(overbeek_inframes, ours_inframes, '^', label='Synthetic vs Endogenous (R=%.2f)' % pearsonr(overbeek_inframes, ours_inframes)[0], color='C1')
    for (x,y,id) in zip(overbeek_inframes, ours_inframes, oof_sel_overbeek_ids):
        if abs(x-y) > 25.0: PL.text(x,y,id)
    PL.plot([0,100],[0,100],'k--')
    PL.ylabel('Percent In-Frame Mutations')
    PL.xlabel('Percent In-Frame Mutations')
    PL.legend()
    PL.xticks([],[])
    PL.yticks([],[])
    PL.show(block=False)
    saveFig('in_frame_full_scatter')

def loadMappings():
    f = io.open(getHighDataDir() + '/overbeek_to_oligo_mapping.txt')
    reader = csv.reader(f, delimiter='\t')
    mappings = {}
    for toks in reader:
        overbeek_id = 'Overbeek' + toks[0].split()[-1]
        oligo_id = toks[1].split('_')[0]
        old = (toks[2] == 'Old')
        if overbeek_id not in mappings:
            mappings[overbeek_id] = []
        mappings[overbeek_id].append((oligo_id, old))
    f.close()
    return mappings

def computePercAbove30(profile):
    
    above_count, below_count = 0, 0
    i, o, s = fetchIndelSizeCounts(profile)
    for x in s['D']:
        if x > 30:
            above_count += s['D'][x]
        else:
            below_count += s['D'][x]
    count_i = 0
    for x in s['I']:
        count_i += s['I'][x]

    return above_count*100.0/(below_count + above_count + count_i)


def compareOverbeekProfiles(selected_overbeek_id=None, pred_results_dir = '../indel_prediction/model_testing'):

    new_dirs = ['ST_June_2017/data/K562_800x_LV7A_DPI7/mapped_reads/Oligos_71',
               'ST_June_2017/data/K562_800x_LV7A_DPI10/mapped_reads/Oligos_71',
               'ST_June_2017/data/K562_800x_LV7B_DPI7/mapped_reads/Oligos_71',
               'ST_June_2017/data/K562_800x_LV7B_DPI10/mapped_reads/Oligos_71',
               'ST_June_2017/data/K562_1600x_LV7B_DPI5/mapped_reads/Oligos_71',
               'ST_Feb_2018/data/CAS9_12NA_1600X_DPI7/mapped_reads/Oligos_71'
               ]

    #Old Samples
    old_dirs = ['ST_June_2017/data/K562_1600x_6OA_DPI5/mapped_reads/Oligos_71',
               'ST_June_2017/data/K562_1600x_6OA_DPI7/mapped_reads/Oligos_71',
               'ST_April_2017/data/K562_800x_6OA_DPI3_Old7/mapped_reads/Oligos_71',
               'ST_April_2017/data/K562_800x_6OA_DPI7_Old8/mapped_reads/Oligos_71',
               'ST_April_2017/data/K562_800x_6OA_DPI10_Old9/mapped_reads/Oligos_71',
               'ST_April_2017/data/K562_800x_6OB_DPI3_Old10/mapped_reads/Oligos_71',
               'ST_April_2017/data/K562_800x_6OB_DPI7_Old11/mapped_reads/Oligos_71',
               'ST_April_2017/data/K562_800x_6OB_DPI10_Old12/mapped_reads/Oligos_71'
               ]
    remove_long_indels = False
    remove_wt, wt_thresh = True, 3.0
    mappings = loadMappings()
        
    all_overbeek_profiles, all_new_profiles, all_old_profiles, all_our_profiles, sel_overbeek_ids,oldnew_overbeek_ids, old_ids, new_ids = [],[],[],[], [],[],[],[]

    overbeek_inframes, ours_inframes, oof_sel_overbeek_ids = [], [], []

    kls, kls_old, kls_new, log_reads, overbeek_ids, above30_percentages, log_reads_new, log_reads_old, min_log_reads = [],[],[],[],[],[],[],[], []
    for idx in range(1,97):

        overbeek_id = 'Overbeek%d' % idx
        if selected_overbeek_id is not None and selected_overbeek_id != overbeek_id:
            continue
        if overbeek_id not in mappings:
            continue
        
        overbeek_filename = getHighDataDir() +'/overbeek_fastq_files/' + overbeek_id + '_mappedindelsummary.txt'

        p1, p1_new, p1_old, o1, rep_reads1, rep_reads2 = {}, {}, {}, {}, {},{}
        nreads2, nreads1, nreads_old, nreads_new, nnull_old, nnull_new, nnull1, nnull2 = 0,0,0,0,0,0,0,0
        
        #Read the overbreek profile
        numread2, perc_accept2, num_null2  = readSummaryToProfile(overbeek_filename, o1, oligoid=overbeek_id, remove_long_indels=remove_long_indels, remove_wt = False)
        if selected_overbeek_id is not None: 
            fetchRepresentativeCleanReads(getHighDataDir() +'/overbeek_fastq_files/' + overbeek_id + '_mappedindelprofiles.txt', rep_reads2, oligoid=overbeek_id)
            pam_loc2, pam_dir2 = getNullTargetPamDetails(getHighDataDir() +'/overbeek_control_fastq_files/' + overbeek_id + '_exptargets.txt', oligoid=overbeek_id)
        nreads2 += numread2
        nnull2 += num_null2

        if numread2 == 0: continue
        
        p1_new_reps, p1_old_reps = [{},{}],[{},{}]
        rr_new_reps, rr_old_reps = [{},{}],[{},{}]
        #Read all the new and old profiles
        pam_loc1, pam_dir1 = None, None
        for oligo_id, is_old in mappings[overbeek_id]:

            #Read all reads for all our K562 profiles
            oligo_idx = eval(oligo_id[5:])
            _, oligo_fileprefix = getFileForOligoIdx(oligo_idx, ext='')
            oligo_filename = oligo_fileprefix + '_mappedindelsummary.txt'
            read_filename = oligo_fileprefix + '_mappedindelprofiles.txt'
            exptarget_filename = oligo_fileprefix + '_exptargets.txt'
            if is_old:
                oligo_dirs, p1_old_new, null_oligo_dir = old_dirs, p1_old, 'ST_April_2017/data/NULL_Old/mapped_reads/Oligos_71'
                p1_reps, rr_reps = p1_old_reps, rr_old_reps
            else:
                oligo_dirs, p1_old_new,  null_oligo_dir = new_dirs, p1_new,'ST_April_2017/data/NULL_New/mapped_reads/Oligos_71'
                p1_reps, rr_reps = p1_new_reps, rr_new_reps
                
            for oligo_dir in [getHighDataDir()+ '/' + x for x in oligo_dirs]:
                nr1, pa1, nn1  = readSummaryToProfile(oligo_dir + '/' + oligo_filename, p1_old_new, oligoid=oligo_id, remove_long_indels=remove_long_indels, remove_wt = remove_wt, wt_thresh=wt_thresh)
                numread1, perc_accept1, num_null1  = readSummaryToProfile(oligo_dir + '/' + oligo_filename, p1, oligoid=oligo_id, remove_long_indels=remove_long_indels, remove_wt = remove_wt, wt_thresh=wt_thresh)
                if 'DPI7' in oligo_dir:
                    rep_idx = 0 if '800x' in oligo_dir else 1
                    nr_rep, pa_rep, nn_rep  = readSummaryToProfile(oligo_dir + '/' + oligo_filename, p1_reps[rep_idx], oligoid=oligo_id, remove_long_indels=remove_long_indels, remove_wt = remove_wt, wt_thresh=wt_thresh)
                if selected_overbeek_id is not None: 
                    fetchRepresentativeCleanReads(oligo_dir + '/' + read_filename, rep_reads1, oligoid=oligo_id)
                    if 'DPI7' in oligo_dir:
                        fetchRepresentativeCleanReads(oligo_dir + '/' + read_filename, rr_reps[rep_idx], oligoid=oligo_id)
                    if pam_loc1 is None:
                        pam_loc1, pam_dir1 = getNullTargetPamDetails(getHighDataDir()+ '/' + null_oligo_dir + '/' + exptarget_filename, oligoid=oligo_id )
                if is_old: nreads_old += numread1; nnull_old += num_null1
                else: nreads_new += numread1; nnull_new += num_null1
                nreads1 += numread1; nnull1 += num_null1

        kls.append(symmetricKL(p1, o1, True))
        kls_old.append(symmetricKL(p1_old, o1, True))
        kls_new.append(symmetricKL(p1_new, o1, True))
            
        log_reads.append(np.log10(nreads1 - nnull1 + 0.5))
        log_reads_old.append(np.log10(nreads_old - nnull_old + 0.5))
        log_reads_new.append(np.log10(nreads_new - nnull_new + 0.5))
        min_log_reads.append(min(log_reads_old[-1], log_reads_new[-1]))
        above30_percentages.append(computePercAbove30(o1))
        overbeek_ids.append(overbeek_id)

        if log_reads[-1] > 2.0:
            all_overbeek_profiles.append(o1)
            all_our_profiles.append(p1)
            sel_overbeek_ids.append(overbeek_id[8:])
            if above30_percentages[-1] < 50.0:
                oif, oof, _ = fetchIndelSizeCounts(o1)
                pif, pof, _ = fetchIndelSizeCounts(p1)
                overbeek_inframes.append(oif*100.0/(oif+oof))
                ours_inframes.append(pif*100.0/(pif+pof))
                oof_sel_overbeek_ids.append(overbeek_id)

        if min_log_reads[-1] > 2.0:
            all_new_profiles.append(p1_new)
            all_old_profiles.append(p1_old)
            oldnew_overbeek_ids.append(overbeek_id)
            old_ids.append([id for id,is_old in mappings[overbeek_id] if is_old][0])
            new_ids.append([id for id,is_old in mappings[overbeek_id] if not is_old][0])

        try:
            print(overbeek_id, [x for (x,y) in mappings[overbeek_id]], kls[-1], nreads2, nreads1)
        except KeyError:
            print('Could not find', overbeek_id)
            print(mappings)

        if selected_overbeek_id is not None:
            title = '%s (KL=%.1f)' % (overbeek_id, kls[-1])
            labels =  ['Conventional scaffold Rep A', 'Conventional scaffold  Rep B','Improved scaffold Rep A', 'Improved scaffold  Rep B', 'Endogenous Profile']
            plotProfiles([p1_old_reps[0],p1_old_reps[1],p1_new_reps[0],p1_new_reps[0], o1], [rr_old_reps[0],rr_old_reps[1],rr_new_reps[0],rr_new_reps[1],rep_reads2], [pam_loc1,pam_loc1,pam_loc1,pam_loc1,pam_loc2],[x=='REVERSE' for x in [pam_dir1, pam_dir1, pam_dir1,pam_dir1,pam_dir2]],labels, title=title)

    if selected_overbeek_id is None:

        plotInFrame(overbeek_inframes, ours_inframes, oof_sel_overbeek_ids, pred_results_dir)
          
        i = 1
        PL.figure(figsize=(5.5,5))
        for thr_l, thr_h in [(0.0,10.0), (10.0,20.0), (20.0,50.0), (50.0,90.0), (90.0,100.0)]:
            ydata = [kl for (kl,a30,id,reads) in zip(kls, above30_percentages, overbeek_ids, log_reads ) if a30 > thr_l and a30 <= thr_h]
            xdata = [reads for (kl,a30,id,reads) in zip(kls, above30_percentages, overbeek_ids, log_reads ) if a30 > thr_l and a30 <= thr_h]
            sel_ids = [id for (kl,a30,id,reads) in zip(kls, above30_percentages, overbeek_ids, log_reads ) if a30 > thr_l and a30 <= thr_h]
            PL.plot(xdata,ydata, 'o', label= '%d-%d%% Deletions > 30' % (thr_l, thr_h))
            for x,y,id in zip(xdata, ydata, sel_ids):
                if y > 3 and x > 2: 
                    PL.text(x,y,id)
        PL.legend()
        PL.plot([0,6],[0.77,0.77],'--', color='grey')
        PL.text(0.1,0.5,'Median between our replicates', color='grey')
        PL.ylabel('Symmetric KL Divergence',fontsize=12)
        PL.xlabel('Log10 Mutated Reads',fontsize=12)
        PL.xlim((0,5.5))
        PL.ylim((0,8))
        PL.show(block=False)
        saveFig('scatter_KL')
        i += 1

        print('Median=',np.median(kls), 'Mean KL=',np.mean(kls)) 
        print(len(overbeek_ids))

        #Compute pairwise KL between overbeek and ours
        N = len(sel_overbeek_ids)
        kl_mat = np.zeros((N,N))
        for i, o1 in enumerate(all_overbeek_profiles):
            for j,p1 in enumerate(all_our_profiles):
                kl_mat[i,j] = symmetricKL(o1, p1)
        PL.figure(figsize=(8,6))
        PL.imshow(kl_mat, cmap='hot_r',vmin = 0.0, vmax = 3.0, interpolation='nearest')
        PL.xticks(range(N),sel_overbeek_ids,rotation='vertical', fontsize=6)
        PL.yticks(range(N),sel_overbeek_ids,rotation='horizontal', fontsize=6)
        PL.xlabel('Synthetic Measurement', fontsize=12)
        PL.ylabel('Endogenous Measurement', fontsize=12)
        PL.title('KL', fontsize=12)
        PL.colorbar()
        PL.show(block=False) 
        saveFig('heatmap_KL')
            

if __name__ == '__main__':

    setHighDataDir('.')
    setPlotDir('plots')
    setFigType('png')
    compareOverbeekProfiles(selected_overbeek_id='Overbeek16', pred_results_dir='.')
    compareOverbeekProfiles(selected_overbeek_id=None, pred_results_dir='.')
