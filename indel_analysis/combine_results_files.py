import io, os, sys, shutil

def combineFiles(dirname, has_hdr_line):
    if not os.path.isdir(dirname): return
    print(dirname)
    filenames = [dirname + '/' + x for x in os.listdir(dirname)]
    fout = io.open(dirname + '.txt', 'w')
    for i, filename in enumerate(filenames):
        f = io.open(filename)
        if has_hdr_line and i != 0:
            f.readline()
        for line in f: fout.write(line)
        f.close()
    fout.close()
   
def combineAllFiles(results_dir, has_hdr_line):
    dirnames = [x for x in os.listdir(results_dir) if '.txt' not in x]
    for dirname in dirnames:
        combineFiles(results_dir + '/' + dirname, has_hdr_line)
        shutil.rmtree(results_dir + '/' + dirname)
        
if __name__ == '__main__':
            
    if len(sys.argv) != 3:
        print('Usage: combine_result_files.py <results_dir> <has hdr line>')
    else:   
        results_dir = sys.argv[1]
        has_hdr_line = eval(sys.argv[2])
        combineAllFiles(results_dir, has_hdr_line)


