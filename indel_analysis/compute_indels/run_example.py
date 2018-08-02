import io, os, sys, csv

from selftarget.util import setRunLocal, setPythonCmd, setPearExe, setIndelMapExe
from selftarget.data import setHighDataDir

from run_all_pear import runAllPear
from run_all_partition import runAllPartition
from run_all_map import runAllMap
from run_all_split_null_mappings import runAllSplitNullMappings
from run_all_mapped_split import runAllMappedSplit
from run_all_compile_nulls import runAllCompileNulls
from run_all_indelmap import runAllIndelMap

def printStatus(status):
    print('\n### ',status,' ###\n ')

setRunLocal(True)
setHighDataDir('example/')
setPythonCmd('python')
setPearExe('~/pear-0.9.10-bin-64/pear-0.9.10-bin-64')
setIndelMapExe('~/indelmap/bin/indelmap')

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