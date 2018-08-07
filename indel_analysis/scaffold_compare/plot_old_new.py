import io, os, sys, csv, random
import pandas as pd
import numpy as np
import pylab as PL

from scipy.stats import pearsonr, spearmanr
from sklearn.manifold import TSNE

from selftarget.data import setHighDataDir, getHighDataDir
from selftarget.plot import saveFig

def renameCol(colname):
    replace_terms = [('Alt2',''),('Alt',''),('New v New', 'Within Library'),('Old v Old','Within Additional Library'),('Old v New','Between Library'),('Class KL',''), ('KL','')]
    for term1, term2 in replace_terms:
         colname = colname.replace(term1, term2)
    return colname.strip()

def runAnalysis():
	
    data = pd.read_csv(getHighDataDir() + '/old_new_kl_summaries.txt', sep='\t').fillna(-1.0)
    kl_cols = [x for x in data.columns if 'KL' in x and 'Class KL' not in x and 'Old v Old' not in x]
    max_kl = 9
    PL.figure(figsize=(2.5,4))
    bps= []
    box_types = [('C2','Within Library'),('C0','Between Library')]
    for i,(clr,box_type) in enumerate(box_types):
        col_box_data = [data[col] for col in kl_cols if renameCol(col) == box_type]
        pos = [2*x + i + 1 for x in range(len(col_box_data))]
        print('KL', box_type, np.median(col_box_data, axis=1))
        bps.append(PL.boxplot(col_box_data, positions=pos,patch_artist=True,boxprops=dict(facecolor=clr),showfliers=False))
    PL.xticks([1.5,3.5,5.5],['Same\ngRNA','Other\ngRNA','Other\ngRNA\n(Rpt)'])
    PL.plot([2.5, 2.5],[0, max_kl],'-', color='silver')
    PL.plot([4.5, 4.5],[0, max_kl],'-', color='silver')
    PL.xlim((0.5,6.5))
    PL.ylim((0,max_kl))
    PL.ylabel('KL')
    PL.subplots_adjust(left=0.1,right=0.95,top=0.95, bottom=0.25)
    PL.legend([bp["boxes"][0] for bp in bps],[x[1] for x in box_types], loc='upper left')
    PL.show(block=False)
    saveFig('kl_compare_old_new_KL')


if __name__ == '__main__':
    setHighDataDir('.')
    runAnalysis()
    import pdb; pdb.set_trace()