from selftarget.util import runPerSubdir

if __name__ == '__main__':
    
    runPerSubdir('compile_i1.py', 'out_i1', __file__, extra_args='. ')