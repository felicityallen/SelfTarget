import io, os, sys
from selftarget.data import getAllDataDirs, getShortDir
from selftarget.util import getLogDir, runCmdCheckIdx, getPythonCmd

def runAllPartition(start_idx=0, stop_idx=1000000, nump=50):
    all_dir, out_dir = getAllDataDirs(), getLogDir()

    idx = 0
    for dirname in all_dir:

        print(getShortDir(dirname), idx)
        filenames = [x for x in os.listdir(dirname) if x.split('_')[-1] == 'pear.assembled.fastq']

        for filename in filenames:

            cmd = getPythonCmd() + ' partition_pear.py %s/%s %d' % (dirname,filename,nump)        
            idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_part')
        
if __name__ == '__main__':
    if len(sys.argv) == 3:
        start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
    else:
        start_idx, stop_idx = -1, -1
    runAllPartition(start_idx, stop_idx)
