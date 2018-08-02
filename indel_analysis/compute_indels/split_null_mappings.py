import io, sys, os, csv
from Bio import SeqIO

from selftarget.oligo import loadPamLookup, loadOligosByBarcode, getFileForOligoIdx

def closeFiles(fhandles):
    for id in fhandles:
        fhandles[id].close()

def writeBatchToFile(read_by_file, output_dir):
    for (filedir, filename) in read_by_file:
        if not os.path.isdir(output_dir + '/' + filedir):
            os.mkdir(output_dir + '/' + filedir)
        mapfilename = filename[:-6] + '_mappings.txt'
        if mapfilename in os.listdir(output_dir + '/' + filedir):
            fout = io.open(output_dir + '/' + filedir + '/' + mapfilename, 'a')
        else:
            fout = io.open(output_dir + '/' + filedir + '/' + mapfilename, 'w')
        for line in read_by_file[(filedir, filename)]:
            fout.write(u'%s\n' % line)
        fout.close()

if __name__ == '__main__':

    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print('split_null_mappings.py <high_dir> <(opt)map_dir_ext>')
    else:
        highdir = sys.argv[1]
        map_dir_ext = sys.argv[2] if len(sys.argv) >= 3 else ''

        mapping_dir = '/mapping_files%s/' % map_dir_ext
        map_dir = '/mapped_reads%s/' % map_dir_ext 
        if not os.path.isdir(highdir + '/' + map_dir):
            os.mkdir(highdir + '/' + map_dir)

        mapping_files = os.listdir(highdir + mapping_dir)
        total, assigned = 0,0
        for mapfile in mapping_files:
            
            fastq_file = mapfile[:-13] + '.fastq'

            batch = 0
            read_by_file = {}
            batch_size = 10000
            
            f = io.open(highdir + mapping_dir + mapfile)
            rdr = csv.reader(f, delimiter='\t')
            for toks in rdr:
                if '@@@' in toks[0]: continue
                oligo_id = toks[1].split()[0]
                batch += 1
                total += 1
                if oligo_id == 'None':
                    continue
                if ',' in oligo_id:
                    continue # Ignore ambiguously mapped reads
                oligo_idx = eval(oligo_id[5:].split('_')[0])

                filepath, filename = getFileForOligoIdx( oligo_idx )
                if (filepath, filename) not in read_by_file:
                    read_by_file[(filepath,filename)] = []
                read_by_file[(filepath,filename)].append('\t'.join(toks))


                assigned += 1

                if batch >= batch_size:
                    writeBatchToFile(read_by_file, highdir + map_dir)
                    batch = 0
                    read_by_file = {}

            if batch > 0 :
                writeBatchToFile(read_by_file, highdir + map_dir)

        print('Total records:', total)
        print('Total assigned:', assigned)


