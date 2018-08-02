import io, sys, os
from selftarget.data import getAllNullDirs, getShortDir, getPamLocFile, getExpOligoFile
from selftarget.util import getLogDir, runCmdCheckIdx, getPythonCmd

def runAllCompileNulls(start_idx=0, stop_idx=100000, map_ext=''):

    alt_pam = ('alt_pam' in map_ext)  
    map_dir, out_label = '/mapped_reads%s/' % map_ext, 'out_null%s' % map_ext

    null_dirs, out_dir = getAllNullDirs(), getLogDir()

    idx = 0
    for dirname in null_dirs:

        print(getShortDir(dirname), idx)
        if not os.path.isdir(dirname + map_dir): continue
    
        subdirs = os.listdir(dirname + map_dir)
        for subdir in subdirs:
            if not os.path.isdir(dirname + map_dir + subdir):
                continue

            exp_oligo_file = getExpOligoFile('/'.join(subdir.split('/')[:-2]))
            cmd = getPythonCmd() + ' compile_mapped_null_profiles.py %s %d %s %s' % (dirname + map_dir + subdir, alt_pam, getPamLocFile(), exp_oligo_file)
            idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, out_label)
            
    
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: run_all_compile_nulls.py <start_idx> <stop_idx> <(opt)map_dir_ext>') 

    start_idx, stop_idx, map_ext = -1, -1, '' 
    if len(sys.argv) > 2:
        start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
    if len(sys.argv) >= 3: map_ext = sys.argv[3]   
    runAllCompileNulls(start_idx, stop_idx, map_ext)
