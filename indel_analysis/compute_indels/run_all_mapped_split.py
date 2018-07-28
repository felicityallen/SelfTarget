import io, os, sys
from selftarget.data import getAllDataDirs, getExpOligoFile, getShortDir
from selftarget.util import getLogDir, runCmdCheckIdx 

if len(sys.argv) < 3:
    print('Usage: run_all_mapped_split.py <start_idx> <stop_idx> <(opt)part([] runs all together, all runs all separately)> <(opt)map_dir_ext>') 

start_idx, stop_idx, part, map_dir_ext = -1, -1, '[]', '' 
if len(sys.argv) > 2:
    start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
if len(sys.argv) >= 4: part = sys.argv[3]
if len(sys.argv) >= 5: map_dir_ext = sys.argv[4]

if part == 'all': parts = ["[%d]" % x for x in range(50)]
else: parts = [part]

map_dir, out_label = '/mapped_reads%s' % map_dir_ext, 'out_split%s' % map_dir_ext

all_dir, out_dir = getAllDataDirs(), getLogDir()

for part in parts:

    idx = 0
    for dirname in all_dir:

        print getShortDir(dirname), idx
                
        extra_cmd = ''
        if part == '[]' or 0 in eval(part): extra_cmd = 'rm -rf %s/%s' % (dirname, map_dir)
        cmd = '~/run_python.sh split_mapped_reads_by_id.py %s %s %s/%s "%s"' % (dirname, getExpOligoFile(dirname), dirname, map_dir, part)
        idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, out_label, extra_cmd = extra_cmd, queue='normal')

