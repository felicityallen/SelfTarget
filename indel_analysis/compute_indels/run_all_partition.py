import io, os, sys
from selftarget.data import getAllDataDirs, getShortDir
from selftarget.util import getLogDir, runCmdCheckIdx 

if len(sys.argv) == 3:
    start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
else:
    start_idx, stop_idx = -1, -1

all_dir, out_dir = getAllDataDirs(), getLogDir()

nump = 50

idx = 0
for dirname in all_dir:

    print getShortDir(dirname), idx

    filenames = [x for x in os.listdir(dirname) if x.split('_')[-1] == 'pear.assembled.fastq']

    for filename in filenames:

        cmd = '~/run_python.sh partition_pear.py %s/%s %d' % (dirname,filename,nump)        
        idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_part')
        
 
