import io, os, sys
import itertools
from selftarget.data import getAllDataDirs, getShortDir, getSubdirs, isOldLib
from selftarget.util import getLogDir, runCmdCheckIdx, runSubdir

if __name__ == '__main__':
    
    all_dir, out_dir = getAllDataDirs(), getLogDir()
    idx = 0
    for old_lib in [False]: #[True,False]:
    
        lib_dirs = [x for x in all_dir if (isOldLib(x) == old_lib) and os.path.isdir(x + '/mapped_reads') and 'DPI7' in x and 'K562_1600x_LV7B_DPI7' not in x and '2A_TREX' not in x and 'K562_800x_7A_DPI7_may' not in x]
    
        for dirname1, dirname2 in itertools.combinations(lib_dirs,2):

            subdirs_1 = getSubdirs(dirname1, withpath=False)
            subdirs_2 = getSubdirs(dirname2, withpath=False)
            common_subdirs = set(subdirs_1).intersection(set(subdirs_2))
            
            label = '%s\t%s' % (getShortDir(dirname1), getShortDir(dirname2))
            extra_args = '%s %s ' % (dirname1, dirname2)
            idx = runSubdir(idx, common_subdirs, label, 'compare_pairwise.py', 'out_compare_pairwise', __file__, extra_args=extra_args)
            
