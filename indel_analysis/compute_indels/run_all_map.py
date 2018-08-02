import io, os, sys
from selftarget.data import getAllDataDirs, getExpOligoFile, getShortDir
from selftarget.util import getLogDir, runCmdCheckIdx, getRunLocal, getIndelMapExe

def runAllMap(start_idx=0, stop_idx=100000000, recompute=True, unassembled_only=False, max_cut_dist=4, map_dir='mapping_files'):

    all_dir, out_dir = getAllDataDirs(), getLogDir()

    idx = 0
    for dirname in all_dir:

        exp_file = getExpOligoFile(dirname)
        print(getShortDir(dirname), idx)

        check_str = '_pear.unassembled_pear.assembled._' if unassembled_only else '_pear.assembled._'
        filenames = [x for x in os.listdir(dirname) if check_str in x]
        for filename in filenames:

                cmd_args = (dirname, filename, exp_file, dirname,map_dir,filename[:-6], max_cut_dist)
                cmd = getIndelMapExe() + ' %s/%s %s %s/%s/%s_mappings.txt 1 %d' % cmd_args
                extra_cmd = ''
                if not os.path.isdir(dirname + '/' + map_dir):
                    if getRunLocal(): os.mkdir(dirname + '/' + map_dir)
                    else: extra_cmd = 'mkdir %s' % (dirname + '/' + map_dir)
                if not recompute and os.path.isfile(dirname + '/' + map_dir + '/' + filename[:-6] + '_mappings.txt'): continue
                idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_map', extra_cmd = extra_cmd)
                

if __name__ == '__main__':

    start_idx, stop_idx, recompute, unassembled_only  = -1,-1,True, False
    map_dir, max_cut_dist = 'mapping_files', 4

    if len(sys.argv) >= 3:
        start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
    if len(sys.argv) >= 4:
        recompute = eval(sys.argv[3])
    if len(sys.argv) >= 5:
        unassembled_only = eval(sys.argv[4])
    if len(sys.argv) >= 6:
        map_dir = sys.argv[5]
    if len(sys.argv) >= 7:
        max_cut_dist = eval(sys.argv[6])
    else:
        print('Usage: run_all_map.py <start_idx> <stop_idx> <recompute_existing> <unassembled_only> <mapping_dir> <max_cut_dist>')
    runAllMap(start_idx, stop_idx, recompute, unassembled_only, max_cut_dist, map_dir)