import io, os, sys, csv, shutil

from selftarget.util import setRunLocal, setPythonCmd, setPearExe, setIndelMapExe, runPerSubdir, runSubdir
from selftarget.util import setIndelGenI1Exe, setIndelMhExe, getIndelGenI1Exe, getIndelMhExe
from selftarget.data import setHighDataDir, getHighDataDir

from run_all_pear import runAllPear
from run_all_partition import runAllPartition
from run_all_map import runAllMap
from run_all_split_null_mappings import runAllSplitNullMappings
from run_all_mapped_split import runAllMappedSplit
from run_all_compile_nulls import runAllCompileNulls
from run_all_indelmap import runAllIndelMap

sys.path.append('..')
from combine_results_files import combineAllFiles

sys.path.append('../microhomology')
from run_all_collect_mh_frequencies_by_len import runAllCollectMHFrequenciesByLen
sys.path.append('../microhomology_mismatch')
from fetch_mh_mismatch_frequencies import fetchMhMismatchFrequencies

def printStatus(status):
    print('\n### ',status,' ###\n ')
    

#----------------------------------------------------------------------
# Copy all example data to results directory since script runs in place
#----------------------------------------------------------------------
shutil.copytree('/data/indel_processing_example', '/results/indel_processing_example')

setRunLocal(True)
if not os.path.isdir('/results/indel_processing_example'): os.mkdir('/results/indel_processing_example')
setHighDataDir('/results/indel_processing_example/')
setPythonCmd('python')
setPearExe('/usr/local/bin/pear')
setIndelMapExe('/usr/local/bin/indelmap')
setIndelGenI1Exe('/usr/local/bin/indelgen_i1')
setIndelMhExe('/usr/local/bin/indelmh')

#----------------------------------------------------------------
# Processing of raw reads to produce descriptions of indels
#----------------------------------------------------------------

#Note:  This provides a demonstration on a cut-down data set of just 4 oligos.
#       In practice, the dataset was much too be large to be run in one script like this. 
#       Instead individual steps were performed by calling each script in turn from the
#       command line to set off parallel jobs on a compute cluster.

printStatus('Combine paired-end reads')
runAllPear()                #Combines paired-end reads
printStatus('Divide reads')
runAllPartition(nump=1)     #Divides reads into smaller files (for parallel processing, in this case leave as one)
printStatus('Map reads to oligos')
runAllMap()                 #Maps reads to gRNA-targets (oligos)
printStatus('Re-organise mappings for plasmid library')
runAllSplitNullMappings()   #Copies gRNA-target (oligo) mappings for plasmid library (NULL) into per-oligo files
printStatus('Re-organise reads for all samples')
runAllMappedSplit()         #Copies mapped reads into per-oligo files
printStatus('Compile expanded template set from plasmid library')
runAllCompileNulls()         #Processes mappings of plasmid library to form expanded set of oligo templates (accounting for synthesis mutations)
printStatus('Re-map reads to expanded template set')
runAllIndelMap()            #Re-map all reads against expanded set of oligo templates, and format into summary files
printStatus('Indel summaries complete.')

#----------------------------------------------------------------
# Further processing of indels to compute summary information
#----------------------------------------------------------------

## I1
printStatus('I1')
cmd = getIndelGenI1Exe() + ' ' + getHighDataDir() +  'ST_June_2017/data/exp_target_pam_new.fasta'
print(cmd); os.system(cmd)
runPerSubdir('../i1/compile_i1.py', 'out_i1', None, extra_args=(getHighDataDir() + ' '))
combineAllFiles(getHighDataDir() + '/i1_summaries',True)

## Indel Category, Size and Most Common Indel Details
printStatus('Indel Details')
mis_dir = getHighDataDir() + '/more_indel_summaries'
runPerSubdir('../indel_details/compile_indel_details.py', 'out_details', None, extra_args=(getHighDataDir() + ' '))
idx, dirnames = 0, [mis_dir + '/' + x for x in os.listdir(mis_dir)]
idx = runSubdir(idx, dirnames, 'All', '../indel_details/compile_pie_summaries_per_oligo.py', 'out_pie_oligo', None, extra_args=(getHighDataDir() + ' '))
combineAllFiles(mis_dir,True)

## Microhomology
printStatus('Microhomology')
cmd = getIndelMhExe() + ' ' + getHighDataDir() +  'ST_June_2017/data/exp_target_pam_new.fasta ' + getHighDataDir() + 'ST_June_2017/data/exp_target_new_mh_indels.txt'
print(cmd); os.system(cmd)
runPerSubdir('../microhomology/fetch_mh_indel_frequencies.py', 'out_mh_indel', None, extra_args=(getHighDataDir() + ' ') )
runAllCollectMHFrequenciesByLen(input_dir=getHighDataDir() + '/mh_indel_frequencies', highdir=getHighDataDir(), scriptloc='../microhomology')

## Mismatch Microhomology:
printStatus('Mismatch-Microhomology')
print(getHighDataDir())
fetchMhMismatchFrequencies(getHighDataDir() + '/ST_June_2017/data/K562_800x_LV7A_DPI7', outdir= getHighDataDir() + '/mh_mismatch_indel_frequencies')
