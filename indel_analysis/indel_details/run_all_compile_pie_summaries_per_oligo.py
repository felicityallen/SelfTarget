import io, sys, os
from selftarget.util import runSubdir

if __name__ == '__main__':
    
    idx = 0
    input_dir = 'more_indel_summaries'
    dirnames = [input_dir + '/' + x for x in os.listdir(input_dir) if os.path.isdir(input_dir + '/' + x)]
    idx = runSubdir(idx, dirnames, 'All', 'compile_pie_summaries_per_oligo.py', 'out_pie_oligo', __file__, extra_args='. ')


