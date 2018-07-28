import io, os, sys, csv
from selftarget.data import getAllDataDirs, getExpOligoFile, getShortDir, getNullDir, getDirLabel, getSubdirs
from selftarget.util import getLogDir, runCmdCheckIdx 
from selftarget.oligo import getOligoIdsFromMappedFastaFile, splitOligoFastaFile, getFileForOligoIdx

start_idx, stop_idx = -1, -1
queue = 'normal'
source = 'mapped_counts'
indelmap_only = False
if len(sys.argv) != 3 and len(sys.argv) != 4 and len(sys.argv) != 5 and len(sys.argv) != 6:
    print 'Usage: run_all_indelmap.py start_idx stop_idx (opt)queue (opt)indelmap_only (opt)source'
else:
    start_idx, stop_idx = eval(sys.argv[1]),eval(sys.argv[2])
    if len(sys.argv) >= 4: queue = sys.argv[3]
    if len(sys.argv) >= 5: indelmap_only = eval(sys.argv[4])
    if len(sys.argv) >= 6: source = sys.argv[5]
    
all_dir, out_dir = getAllDataDirs(), getLogDir()
idx = 0

if source == 'date_modified':

    f = io.open('../quality_checks/log_files/out_modified.log')
    for line in f:

        if 'mapped_reads' not in line: continue
        #e.g. /lustre/scratch117/cellgen/team227/fa9/self-target/indel_analysis/ST_June_2017/data/K562_800x_LV7B_DPI3/mapped_reads/Oligos_70/Oligos_70850-70899_mappedindelsummary.txt
        full_filename = line.split()[0]
        toks = full_filename.split('/')
        fasta_file = toks[-1][:-23] + '.fasta'
        filepath = '/'.join(toks[:-1])
        dirname = '/'.join(toks[:-3])
        nulldir = getNullDir(dirname)
        subdir = toks[-2]

        print idx, dirname, fasta_file
        cmd = '~/run_python.sh indelmap_subdir.py %s %s %s -1 1 %s %d' % (dirname, nulldir, subdir, fasta_file) 
        idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_incomplete_indelmap_%s' % getDirLabel(dirname), queue=queue)

    f.close()

elif source == 'mapped_counts':

    for dirname in all_dir:
        dirlabel = getDirLabel(dirname)
        nulldir = getNullDir(dirname)
        repeat_indelmaps = set()
        repeat_reformat = set()
        f = io.open('../quality_checks/mapped_read_summaries/%s.txt' % dirlabel)
        for row in csv.DictReader(f, delimiter='\t'):
            if eval(row['Mapping Files']) != eval(row['Split Fasta File']):
                print row
                print 'PROBLEM IN MAPPED SPLIT - RERUN!', dirname
                break
            elif eval(row['Split Fasta File']) != eval(row['Mapped Split']) or eval(row['Mapped Split Assigned']) != eval(row['Mapped Split']):
                oligo_idx = eval(row['ID'][5:])
                filepath, filename = getFileForOligoIdx(oligo_idx)
                repeat_indelmaps.add((filename, filepath))
                print row
            elif eval(row['Mapped Split Assigned']) != eval(row['Summary']):
                oligo_idx = eval(row['ID'][5:])
                filepath, filename = getFileForOligoIdx(oligo_idx)
                print row
                repeat_reformat.add((filename, filepath))
        f.close()

        print 'INDELMAP', idx, dirname
        for filename, subdir in repeat_indelmaps:
            cmd = '~/run_python.sh indelmap_subdir.py %s %s %s -1 1 %s 0' % (dirname, nulldir, subdir, filename) 
            idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_incomplete_indelmapcounts_%s' % getDirLabel(dirname), queue=queue)
        if not indelmap_only:
            print 'REFORMAT', idx, dirname
            for filename, subdir in repeat_reformat:
                if (filename, filepath) in repeat_indelmaps: continue
                cmd = '~/run_python.sh indelmap_subdir.py %s %s %s -1 1 %s 1' % (dirname, nulldir, subdir, filename) 
                idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, 'out_incomplete_indelmapcounts_%s' % getDirLabel(dirname), queue=queue)
            
