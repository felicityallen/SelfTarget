import io, os, sys, csv
from selftarget.data import getAllDataDirs, getExpOligoFile, getShortDir, getNullDir, getDirLabel, getSubdirs
from selftarget.util import getLogDir, runCmdCheckIdx 

overbeek_only = False
start_idx, stop_idx = -1, -1
queue = 'normal'
if len(sys.argv) < 3 or len(sys.argv) > 5:
    print 'Usage: run_all_indelmap.py start_idx stop_idx (opt)overbeek_only (opt)queue'
else:
    start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
    if len(sys.argv) >= 4: overbeek_only = eval(sys.argv[3])
    if len(sys.argv) >= 5: queue = sys.argv[4]
    
all_dir, out_dir = getAllDataDirs(), getLogDir()
if overbeek_only: print 'Computing for Overbeek guides only'

file_per_part = 1
max_files_per_dir = 20 
num_parts = (max_files_per_dir+file_per_part-1)/file_per_part

i, idx = 0, 0
for dirname in all_dir:

    i += 1

    print getShortDir(dirname), idx
    if not os.path.isdir(dirname + '/mapped_reads'): continue

    for subdir in getSubdirs(dirname, withpath=False):
        if overbeek_only and subdir != 'Oligos_71': continue
        cmd = './run_correct_indel.sh %s %s %s %d' % (dirname, getNullDir(dirname), subdir, file_per_part) 
        idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_correct_indelmap_%s' % getDirLabel(dirname), numj=num_parts, queue=queue)

