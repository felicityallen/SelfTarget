from selftarget.util import runSubdir
from selftarget.data import getAllDataDirs

if __name__ == '__main__':
    
    idx = 0
    idx = runSubdir(idx, getAllDataDirs(),'All', 'fetch_mh_mismatch_frequencies.py', 'out_get_mm_freq', __file__)
