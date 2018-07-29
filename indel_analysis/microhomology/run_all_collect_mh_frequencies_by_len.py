import io, sys, os
from selftarget.util import runSubdir

if __name__ == '__main__':
    
    idx = 0
    input_dir = 'mh_indel_frequencies'
    dirnames = [input_dir + '/' + x for x in os.listdir(input_dir)]
    for dirname in dirnames:
        for mh_len in range(2, 16):
            subdirs = [dirname + '/' + x for x in os.listdir(dirname)]
            idx = runSubdir(idx, subdirs, '%s MH Len=%d' % (dirname,mh_len), 'collect_mh_frequencies_by_len.py', 'out_mh_by_len', __file__, extra_args='%d ' % mh_len)
