import io, os, sys, csv, shutil

from selftarget.util import setRunLocal, setPythonCmd, setIndelMapExe
from selftarget.data import setHighDataDir, getHighDataDir

from create_overbeek_fasta import extractReads
from create_overbeek_templates import createOverbeekTemplates
from compute_overbeek_indel_profiles import computeOverbeekIndelProfiles

def printStatus(status):
    print('\n### ',status,' ###\n ')

shutil.copytree('/data/endogenous_processing_example','/results/endogenous_processing_example')

#----------------------------------------------------------------------
# Configure environment
#----------------------------------------------------------------------
setRunLocal(True)
setHighDataDir('/results/endogenous_processing_example/')
setPythonCmd('python')
setIndelMapExe('/usr/local/bin/indelmap')

#----------------------------------------------------------------
# Processing of raw Van-Overbeek et al reads to produce descriptions of indels
#----------------------------------------------------------------
#Note:  This provides a demonstration on just 1 oligo, going from raw overbeek reads to indel descriptions.
#       Sam files are assumed to be already collected (for further details of this part see 
#       collect_overbeek_sams.py in same dir)

printStatus('Create fastq files from Van Overbeek sam files')
sam_dir, fastq_dir = getHighDataDir() + '/overbeek_sam_files', getHighDataDir() + '/overbeek_fastq_files'
if not os.path.isdir(fastq_dir): os.makedirs(fastq_dir)
extractReads(sam_dir + '/Overbeek_6.sam', fastq_dir + '/Overbeek6.fastq','chrX:66765045-66765067', 'Overbeek6')
sam_dir, fastq_dir = getHighDataDir() + '/overbeek_control_sam_files', getHighDataDir() + '/overbeek_control_fastq_files'
if not os.path.isdir(fastq_dir): os.makedirs(fastq_dir)
extractReads(sam_dir + '/Overbeek_6.sam', fastq_dir + '/Overbeek6.fastq','chrX:66765045-66765067', 'Overbeek6')

printStatus('Compute mutational profile from Van Overbeek data')
createOverbeekTemplates(selected_id='Overbeek6')
computeOverbeekIndelProfiles(highdir = getHighDataDir(), selected_id='Overbeek6')
