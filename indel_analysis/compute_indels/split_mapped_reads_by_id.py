import io
from Bio import SeqIO
import sys
import csv
import os
import os.path
import random

from selftarget.oligo import loadOligosByBarcode, getFileForOligoIdx

def loadMappings( filename ):
    lookup = {}
    f = io.open(filename)
    rdr = csv.reader(f, delimiter='\t')
    for toks in rdr:
        if toks[0][:3] == '@@@': continue
        lookup[toks[0]] = toks[1].split()[0]
    return lookup

def createDirectories(lookup, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for (bc1,bc2) in lookup:
        (oligo_id, oligo_idx) = lookup[(bc1,bc2)]
        filedir, filename = getFileForOligoIdx( oligo_idx )
        if filedir not in os.listdir(output_dir):
            os.mkdir(output_dir + '/' + filedir)

def closeFiles(fhandles):
    for id in fhandles:
        fhandles[id].close()

def initCounts(lookup):
    counts = {}
    for (bc1, bc2) in lookup:
        (oligo_id, oligo_idx) = lookup[(bc1,bc2)]
        counts[oligo_id] = 0
    return counts

def writeBatchToFile(read_by_file, output_dir):

    for (filedir, filename) in read_by_file:
        if filename in os.listdir(output_dir + '/' + filedir):
            fout = io.open(output_dir + '/' + filedir + '/' + filename, 'a')
        else:
            fout = io.open(output_dir + '/' + filedir + '/' + filename, 'w')
        for (read_id, read_seq, oligo_id) in read_by_file[(filedir, filename)]:
            fout.write(u'>%s.%s\n%s\n' % (oligo_id.split('_')[0],read_id, read_seq))
        fout.close()

if __name__ == '__main__':

    if len(sys.argv) != 4 and len(sys.argv) != 5 and len(sys.argv) != 6:
        print('split_mapped_reads_by_id.py <high_dir> <exp_oligo_file> <output_dir> <part> <(opt)selected_mapped_fasta>')
    else:
        highdir = sys.argv[1]
        exp_oligo_file = sys.argv[2]
        output_dir = sys.argv[3]
        part, selected_mapped_fasta = [], ''
        if len(sys.argv) >= 5: part = eval(sys.argv[4])
        if len(sys.argv) >= 6: selected_mapped_fasta = sys.argv[5]

        oligo_lookup = loadOligosByBarcode(exp_oligo_file)
        createDirectories(oligo_lookup, output_dir)
        counts = initCounts(oligo_lookup)

        if 'alt_pam' in output_dir: mapping_dir = '/mapping_files_alt_pam/'
        elif 'widencut' in output_dir: mapping_dir = '/mapping_files_widencut/'
        else: mapping_dir = '/mapping_files/'        

        mapping_files = os.listdir(highdir + mapping_dir)
        total, assigned = 0,0
        if len(part) > 0: mapping_files = [x for x in mapping_files if any([('_%d_' % pt) in x for pt in part])]
        for mapfile in mapping_files:
            
            fastq_file = mapfile[:-13] + '.fastq'
            lookup = loadMappings(highdir + mapping_dir + mapfile)
            batch = 0
            read_by_file = {}
            batch_size = 20000
            for record in SeqIO.parse(highdir + '/' + fastq_file,'fastq'):
                batch += 1
                total += 1
                if str(record.description) not in lookup:
                    print('Could not find', str(record.description),'in mapping file')
                    continue
                oligo_id = lookup[str(record.description)]
                if oligo_id == 'None':
                    continue
                if ',' in oligo_id:
                    continue # Ignore ambiguously mapped reads
                oligo_idx = eval(oligo_id[5:].split('_')[0])

                filepath, filename = getFileForOligoIdx( oligo_idx )
                if selected_mapped_fasta != '' and selected_mapped_fasta != filename: continue

                if (filepath, filename) not in read_by_file:
                    read_by_file[(filepath,filename)] = []
                read_by_file[(filepath,filename)].append((record.id, str(record.seq),oligo_id))
                assigned += 1
                counts[oligo_id] += 1

                if batch >= batch_size:
                    writeBatchToFile(read_by_file, sys.argv[3])
                    batch = 0
                    read_by_file = {}

            if batch > 0 :
                writeBatchToFile(read_by_file, sys.argv[3])

        print('Total records:', total)
        print('Total assigned:', assigned)
        for thresh in [0,10,100,500,1000]:
            print('Num counts > %d:' % thresh, sum([counts[x]>thresh for x in counts]))


