import io, os, sys, csv
import pandas as pd
import numpy as np
import pylab as PL

from selftarget.oligo import partitionGuides
from selftarget.util import getPlotDir
from selftarget.data import getSampleSelectors, sortSampleNames

def compilePieSummariesPerOligo(filename, out_file, is_first=True):
    data = pd.read_csv(filename, sep='\t')
    data['MH len'] = data['Microhomology Sequence'].apply(lambda x: len(x) if type(x) == str else 0 )
    is_dbl = lambda indel: ('I' in indel) and ('D' in indel)
    data['Is Dbl'] = data['Most Common Indel'].apply(is_dbl)
    data['I1'] = ((data['Type']=='I') & (data['Size']==1) & (~data['Is Dbl']))*data['MCI Reads']
    data['I>1'] = ((data['Type']=='I') & (data['Size']>1) & (~data['Is Dbl']))*data['MCI Reads']
    data['D1'] = ((data['Type']=='D') & (data['Size']==1) & (~data['Is Dbl']))*data['MCI Reads']
    data['D2'] = ((data['Type']=='D') & (data['Size']==2) & (~data['Is Dbl']))*data['MCI Reads']
    data['D>2, MH'] = ((data['Type']=='D') & (data['MH len']>=2) & (data['Size']>2) & (~data['Is Dbl']))*data['MCI Reads']
    data['D>2, No MH'] = ((data['Type']=='D') & (data['MH len']<2) & (data['Size']>2) & (~data['Is Dbl']))*data['MCI Reads']
    data['I+D'] = data['Is Dbl']*data['MCI Reads']

    if len(data) == 0: return

    cols = ['I1','I>1','D1','D2','D>2, MH','D>2, No MH','I+D']
    po_data = data.groupby('Oligo Id')[cols].sum()
    read_data = data.groupby('Oligo Id')[['Total reads']].mean()
    mci_data = data[data.groupby('Oligo Id')['MCI Reads'].transform(max)==data['MCI Reads']][['Oligo Id','MCI Reads', 'Most Common Indel','Altered Sequence','Microhomology Sequence'] + cols]
    mci_data['MCI Type'] = mci_data[cols].idxmax(axis=1)
    mci_data = mci_data.set_index('Oligo Id')
    read_data = read_data.merge(mci_data[['MCI Reads','Most Common Indel','MCI Type','Altered Sequence','Microhomology Sequence']],left_index=True, right_index=True, how='inner')
    po_data = po_data.merge(read_data,left_index=True, right_index=True,how='inner')
    po_data.to_csv(out_file, sep='\t', mode='a', header=is_first)

def compileAllPieSummariesPerOligo(results_subdir, out_file):
    filenames = os.listdir(results_subdir)
    fout = io.open(out_file, 'w'); fout.close() #Clear file
    for i,filename in enumerate(filenames):
        compilePieSummariesPerOligo(results_subdir + '/' + filename, out_file, i==0)

if __name__ == '__main__':

    if len(sys.argv) != 3:
        raise Exception('Usage: compile_pie_summaries_per_oligo.py highdir results_subdir')
    else:
        highdir = sys.argv[1]
        results_subdir = sys.argv[2]
        subdir = results_subdir.split('/')[-1]
        if not os.path.isdir(results_subdir): raise Exception('Not a directory - %s' % results_subdir)

        out_dir = highdir + '/indel_pie_summaries_per_oligo'
        if not os.path.isdir(out_dir): os.mkdir(out_dir)
        
        out_file = out_dir + '/' + subdir + '.txt'
        compileAllPieSummariesPerOligo(results_subdir, out_file)
