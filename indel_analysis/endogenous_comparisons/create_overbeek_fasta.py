import io
import sys
from Bio import SeqIO
import os
import numpy as np
import re
import csv
import Bio.Seq

def extractReads( sam_file, output_fastq_file, location, oligo_id ):
	chr = location.split(':')[0]
	mid_idx = 0.5*sum([eval(x) for x in location.split(':')[1].split('-')])
	fout = io.open(output_fastq_file,'w')
	count = 0
	f = io.open(sam_file)
	for line in f:
		if line[0] == '@':
			continue
		toks = line.split('\t')
		if len(toks) < 11:
			continue
		if toks[2] != chr:
			continue #Ignore reads mapped to the wrong position
		read_loc = eval(toks[3])
		if abs(read_loc - mid_idx) > 2000:
			continue
		full_read = str(toks[9])
		qual_str = str(toks[10])
		if full_read == '':
			continue
		if sum([x == 'N' for x in full_read]) > 0.5*len(full_read):
			continue
		fout.write(u'@%s:%d\n%s\n+\n%s\n' % (oligo_id, count, full_read, qual_str))
		count += 1
	f.close()
	print(oligo_id, count)
	
	
def loadLocationSpacerLookup():
	f = io.open('overbeek_2016_guides_s1.txt')
	reader = csv.DictReader(f, delimiter='\t')
	lookup = {'Overbeek%d' % eval(row['Spacer ']): (row['Genomic location of spacer (hg19)'], row['Spacer sequence']) for row in reader}
	f.close()
	return lookup
	
if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('create_overbeek_fasta.py <samdir> <output_fastq_dir>')
	else:
		samdir = sys.argv[1]
		output_fastq_dir = sys.argv[2]
		lookup = loadLocationSpacerLookup()
		
		if not os.path.isdir(output_fastq_dir):
			os.mkdir(output_fastq_dir)
			
		f = io.open('overbeek_self_targets.csv')
		reader = csv.reader(f, delimiter='\t')
		for toks in reader:
			idx = eval(toks[-1].split()[-1])
			id, samid = 'Overbeek%d' % idx, 'Overbeek_%d' % idx,
			loc, spacer_seq = lookup[id]
			target_seq = toks[1][1:-1]
			
			sam_file = samdir + '/%s.sam' % samid
			output_fastq_file = output_fastq_dir + '/%s.fastq' % id
			extractReads( sam_file, output_fastq_file, loc, id )

