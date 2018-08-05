import io, os, sys

def combineSubdirs(dirname, has_hdr_line):
    if not os.path.isdir(dirname): return
    print(dirname)
    subdirs = [dirname + '/' + x for x in os.listdir(dirname) if os.path.isdir(dirname + '/' + x)]
    if len(subdirs) == 0: return
    all_files = set(os.listdir(subdirs[0]))
    for subdir in subdirs:
        all_files = all_files.union(set(os.listdir(subdir)))
    for filename in all_files:
        fout = io.open(dirname + '/' + filename, 'w')
        i = 0
        for subdir in subdirs:
            if not os.path.isfile(subdir + '/' + filename): continue
            f = io.open(subdir + '/' + filename)
            if has_hdr_line and i != 0:
                f.readline()
            for line in f: fout.write(line)
            f.close()
            i += 1
        fout.close()

if len(sys.argv) != 3:
    print('Usage: combine_results_subdirs.py <results_dir> <has hdr line>')
else:   
    results_dir = sys.argv[1]
    has_hdr_line = eval(sys.argv[2])

    dirnames = os.listdir(results_dir)
    for dirname in dirnames:
        combineSubdirs(results_dir + '/' + dirname, has_hdr_line)

