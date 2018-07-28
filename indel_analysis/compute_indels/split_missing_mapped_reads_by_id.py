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

def writeBatchToFile(read_by_file, output_dir):
    for (filedir, filename) in read_by_file:
        outfile = output_dir + '/' + filedir + '/' + filename
        
        if filename in os.listdir(output_dir + '/' + filedir):
            fout = io.open(outfile, 'a')
            read_ids = getReadIdsFromFile(outfile)
        else:
            fout = io.open(outfile, 'w')
            read_ids = set([])
        for (read_id, read_seq, oligo_id) in read_by_file[(filedir, filename)]:
            if read_id not in read_ids:
                import pdb; pdb.set_trace()
                fout.write(u'>%s.%s\n%s\n' % (oligo_id.split('_')[0],read_id, read_seq))
        fout.close()

def getReadIdsFromFile(outfile):
    f = io.open(outfile)
    read_ids = []
    for line in f:
        if line[0] == '>':
            read_id = line[1:-1].split('.')[-1]
            read_ids.append(read_id)
    f.close()
    return set(read_ids)
        
if __name__ == '__main__':

    if len(sys.argv) != 4 and len(sys.argv) != 5 and len(sys.argv) != 6:
        print 'split_missing_mapped_reads_by_id.py <high_dir> <exp_oligo_file> <output_dir> <part> <unassembled_only>'
    else:
        highdir = sys.argv[1]
        exp_oligo_file = sys.argv[2]
        output_dir = sys.argv[3]
        part = []
        if len(sys.argv) >= 5: part = eval(sys.argv[4])
        if len(sys.argv) >= 6: unassembled_only = eval(sys.argv[5])

        mapping_files = os.listdir(highdir + '/mapping_files')
        if len(part) > 0: mapping_files = [x for x in mapping_files if any([('_%d_' % pt) in x for pt in part])]
        if unassembled_only: mapping_files = [x for x in mapping_files if 'unassembled' in x]
        read_by_file = {}
        for mapfile in mapping_files:
            
            fastq_file = mapfile[:-13] + '.fastq'
            lookup = loadMappings(highdir + '/mapping_files/' + mapfile)

            for record in SeqIO.parse(highdir + '/' + fastq_file,'fastq'):
                if str(record.description) not in lookup:
                    print 'Could not find', str(record.description),'in mapping file'
                    continue
                oligo_id = lookup[str(record.description)]
                if oligo_id == 'None':
                    continue
                if ',' in oligo_id:
                    print 'Ambiguously mapped read:', str(record.description)
                    continue # Ignore ambiguously mapped reads
                oligo_idx = eval(oligo_id[5:].split('_')[0])

                filepath, filename = getFileForOligoIdx( oligo_idx )
                if (filepath, filename) not in read_by_file:
                    read_by_file[(filepath,filename)] = []
                read_by_file[(filepath,filename)].append((record.id, str(record.seq),oligo_id))

            writeBatchToFile(read_by_file, output_dir)

