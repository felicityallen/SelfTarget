from selftarget.util import runPerSubdir

if __name__ == '__main__':
    
    runPerSubdir('compile_indel_details.py', 'out_mci', __file__, extra_args='. ')
