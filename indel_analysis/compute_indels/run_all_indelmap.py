import io, os, sys, csv
from selftarget.data import getAllDataDirs, getExpOligoFile, getShortDir, getNullDir, getDirLabel, getSubdirs
from selftarget.util import getLogDir, runCmdCheckIdx, getPythonCmd, getIndelMapExe

def runAllIndelMap(start_idx = 0, stop_idx = 10000000, overbeek_only = False, queue='normal', map_dir='/mapped_reads/', max_cut_dist=4, num_parts=1, order_by_incomplete=False):
    
    all_dir, out_dir = getAllDataDirs(), getLogDir()
    if overbeek_only: print('Computing for Overbeek guides only')

    completed_lookup = {}
    if order_by_incomplete:
        f = io.open('../quality_checks/status.log')
        completed_lookup = {toks[0]: min([eval(x) for x in toks[1:]]) != 0 for toks in csv.reader(f, delimiter='\t')}
        f.close()
    completed = [x for x in all_dir if getDirLabel(x) in completed_lookup and completed_lookup[getDirLabel(x)]]
    not_completed = [x for x in all_dir if x not in completed]

    max_files_per_dir = 20 
    file_per_part = int(max_files_per_dir/num_parts + 0.99)

    i, idx = 0, 0
    for dirname in not_completed + completed:

        if len(not_completed) == i: print('-------------------------------------------------')
        i += 1

        print(getShortDir(dirname), idx)
        if not os.path.isdir(dirname + map_dir): continue

        for subdir in getSubdirs(dirname, withpath=False):
            if overbeek_only and subdir != 'Oligos_71': continue
            args = (dirname, getNullDir(dirname), subdir, file_per_part, map_dir, max_cut_dist, getIndelMapExe(), getPythonCmd()) 
            cmd = getPythonCmd() + ' indelmap_subdir.py %s %s %s %d - 0 %s %d %s %s' % args
            idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_indelmap_%s' % getDirLabel(dirname), numj=num_parts, queue=queue)

if __name__ == '__main__':
    overbeek_only = False
    start_idx, stop_idx = -1, -1
    queue = 'normal'
    map_dir = '/mapped_reads/'
    max_cut_dist = 4
    fpp = 1
    if len(sys.argv) < 3 or len(sys.argv) > 9:
        print('Usage: run_all_indelmap.py start_idx stop_idx (opt)overbeek_only (opt)queue (opt)map_dir (opt)fpp (opt)max_cut_dist (opt)order_by_incomplete')
    else:
        start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
        if len(sys.argv) >= 4: overbeek_only = eval(sys.argv[3])
        if len(sys.argv) >= 5: queue = sys.argv[4]
        if len(sys.argv) >= 6: map_dir = sys.argv[5]
        if len(sys.argv) >= 7: fpp = eval(sys.argv[6])
        if len(sys.argv) >= 8: max_cut_dist = eval(sys.argv[7])
        if len(sys.argv) >= 9: order_by_incomplete = eval(sys.argv[8])

        file_per_part = fpp
        max_files_per_dir = 20 
        num_parts = int((max_files_per_dir+file_per_part-1)/file_per_part)

        runAllIndelMap(start_idx, stop_idx, overbeek_only, queue, map_dir, max_cut_dist, num_parts, order_by_incomplete)
    