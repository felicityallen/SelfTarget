import io, sys, os
from selftarget.util import runSubdir
    
def runAllCollectMHFrequenciesByLen(input_dir = 'mh_indel_frequencies', highdir='.', scriptloc='.'):

    idx = 0
    dirnames = [input_dir + '/' + x for x in os.listdir(input_dir)]
    for dirname in dirnames:
        for mh_len in range(2, 16):
            subdirs = [dirname + '/' + x for x in os.listdir(dirname)]
            idx = runSubdir(idx, subdirs, '%s MH Len=%d' % (dirname,mh_len), scriptloc + '/collect_mh_frequencies_by_len.py', 'out_mh_by_len', __file__, extra_args='%d %s ' % (mh_len, highdir))

if __name__ == '__main__':
    runAllCollectMHFrequenciesByLen()