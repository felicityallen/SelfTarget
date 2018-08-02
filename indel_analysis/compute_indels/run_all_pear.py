import io, os, sys
from selftarget.data import getAllDataDirs, getShortDir
from selftarget.util import runCmdCheckIdx, getLogDir, getPearExe

def runAllPear(start_idx=0, stop_idx=100000):
    all_dir, out_dir = getAllDataDirs(), getLogDir()

    idx = 0
    for dirname in all_dir:
        print(getShortDir(dirname), idx)
        
        r1_fasta_files = [x for x in os.listdir(dirname) if x[-9:] == '_R1.fastq' or x[-13:] == '_R1_001.fastq']
        for r1_file in r1_fasta_files:
            file_prefix = r1_file[:r1_file.index('_R1')]
            file_suffix = r1_file[len(file_prefix):].replace('1','2',1)
            if not os.path.isfile(dirname + '/' + file_prefix + file_suffix):
                print('Could not find matching R2 file:', dirname, file_prefix, file_suffix)
                continue
            cargs = (dirname, r1_file, dirname, file_prefix, file_suffix, dirname, file_prefix) 
            cmd = getPearExe() +' -f %s/%s -r %s/%s%s -o %s/%s_pear -n 20 -p 0.01' % cargs
            idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_pear')


if __name__ == '__main__':

    if len(sys.argv) > 1:
        start_idx = eval(sys.argv[1])
    else:
        start_idx = -1
    if len(sys.argv) > 2:
        stop_idx = eval(sys.argv[2])
    else:
        stop_idx = -1

    runAllPair(start_idx, stop_idx)
