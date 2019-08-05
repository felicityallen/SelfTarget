import io, os, sys, csv, random
import pandas as pd
import numpy as np
import pylab as PL

from scipy.stats import pearsonr, spearmanr
from sklearn.manifold import TSNE

from selftarget.plot import saveFig
from selftarget.data import getHighDataDir, setHighDataDir

def renameCol(colname):
    replace_terms = [('Old v New','Between Library'),('Old v Predicted','Predicted Vs Add. Lib'),('New v Predicted','Predicted vs Main Lib'),('Comb v Predicted','Predicted vs Combined Libs'),('Class KL',''), ('KL','')]
    replace_terms = [('Old', 'Conventional Scaffold'),('New', 'Improved Scaffold'),('Combined', 'Combined Replicates'),('In Frame Perc', 'Percent In Frame Mutations'),('KL', '')]
    replace_terms = [('New 1600x', 'Replicate B'),('New 2x800x', 'Replicate A')]
    for term1, term2 in replace_terms:
         colname = colname.replace(term1, term2)
    return colname.strip()

def plotKLBoxes(data):
    cols = [x for x in data.columns if 'KL' in x and 'Class KL' not in x and 'Old' not in x and 'Conventional' not in x and 'Combined' not in x]
    cols.reverse()
    cols_label, max_kl = 'KL', 9
    PL.figure(figsize=(4,5))

    pt = data.loc[(data['Combined v Predicted KL'] > 0.75) & (data['Combined v Predicted KL'] < 0.8) & (data['Old v New KL'] > 0.75) & (data['Old v New KL'] < 0.8)]
    print(pt['Old Oligo Id'])

    PL.boxplot([data[col] for col in cols], positions=range(len(cols)),patch_artist=True,boxprops=dict(facecolor='C2'),medianprops=dict(linewidth=2.5, color='C1'),showfliers=False)
    PL.xticks(range(len(cols)),[renameCol(x) for x in cols], rotation='vertical')
    for i,col in enumerate(cols):
        PL.text(i-0.15, np.median(data[col])+0.02, '%.2f' % np.median(data[col]))
    PL.ylabel(cols_label)
    PL.subplots_adjust(left=0.1,right=0.95,top=0.95, bottom=0.5)
    PL.show(block=False)
    saveFig('kl_compare_old_new_predicted_%s' % cols_label.replace(' ',''))

def plotInFrameCorr(data):
    
    shi_data = pd.read_csv(getHighDataDir() + '/shi_deepseq_frame_shifts.txt',sep='\t')

    label1, label2 = 'New In Frame Perc', 'Predicted In Frame Per'
    PL.figure(figsize=(4,4))

    xdata, ydata = data[label1], data[label2]
    PL.plot(xdata,ydata, '.',alpha=0.15)
    PL.plot(shi_data['Measured Frame Shift'], shi_data['Predicted Frame Shift'], '^', color='orange')
    for x,y,id in zip(shi_data['Measured Frame Shift'], shi_data['Predicted Frame Shift'],shi_data['ID']):
        if x-y > 10:
            PL.text(x,y,id.split('/')[1][:-21])
    PL.plot([0,100],[0,100],'k--')
    PL.title('R=%.3f' % (pearsonr(xdata, ydata)[0]))
    PL.xlabel('percent in frame mutations (measured)')
    PL.ylabel('percent in frame mutations (predicted)')
    PL.ylim((0,80))
    PL.xlim((0,80))
    PL.show(block=False)
    saveFig('in_frame_corr_%s_%s' % (label1.replace(' ','_'),label2.replace(' ','_')))

def runAnalysis():
    data = pd.read_csv(getHighDataDir() + '/old_new_kl_predicted_summaries.txt', sep='\t').fillna(-1.0)
    plotKLBoxes(data)
    plotInFrameCorr(data)

if __name__ == '__main__':
    setHighDataDir('.')
    runAnalysis()
    import pdb; pdb.set_trace()
