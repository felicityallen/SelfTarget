from selftarget.util import runPerSubdir

if __name__ == '__main__':
    
    runPerSubdir('fetch_mh_indel_frequencies.py', 'out_mh_indel', __file__)
