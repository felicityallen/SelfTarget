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
    replace_terms = [('Old', 'Replicate A'),('New', 'Replicate B'),('Combined', 'Combined Replicates'),('In Frame Perc', 'Percent In Frame Mutations'),('KL', '')]
    for term1, term2 in replace_terms:
         colname = colname.replace(term1, term2)
    return colname.strip()

def plotKLBoxes(data):
    cols = [x for x in data.columns if 'KL' in x and 'Class KL' not in x]
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
    
    for label1, label2 in [('Combined in Frame Perc', 'Predicted In Frame Per')]:
        PL.figure(figsize=(4,4))
        xdata, ydata = data[label1], data[label2]
        PL.plot(xdata,ydata, '.')
        PL.plot([0,100],[0,100],'k--')
        PL.title('R=%.3f' % (pearsonr(xdata, ydata)[0]))
        PL.xlabel(renameCol(label1))
        PL.ylabel(renameCol(label2))
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
