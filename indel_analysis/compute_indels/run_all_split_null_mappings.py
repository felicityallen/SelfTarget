import io, sys, os
from selftarget.data import getAllNullDirs, getShortDir
from selftarget.util import getLogDir, runCmdCheckIdx

if __name__ == '__main__':

    start_idx, stop_idx, map_dir_ext = -1, -1, ''
    if len(sys.argv) >= 3:
        start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
    if len(sys.argv) >= 4: map_dir_ext = sys.argv[3]

    null_dirs, out_dir = getAllNullDirs(), getLogDir()
    
    idx = 0
    for dirname in null_dirs:

        print getShortDir(dirname), idx
       
        cmd = '~/run_python.sh split_null_mappings.py %s %s' % (dirname, map_dir_ext)
        idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_null_split')
                
                
                
