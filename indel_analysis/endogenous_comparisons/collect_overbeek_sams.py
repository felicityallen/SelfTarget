import io
import csv
import os

def splitLoc(loc):
    toks = loc.split(':')
    chrom = toks[0]
    loc_toks = toks[1].split('-')
    return chrom, eval(loc_toks[0]), eval(loc_toks[1])

def findNearbyLoc(loc, spacers):
    chr, start, fin = splitLoc(loc)
    for spacer_loc in spacers:
        schr, sstart, sfin = splitLoc(spacer_loc)
        if chr == schr and abs(sstart - start) <= 1 and abs(sfin - fin) <= 1:
            return spacer_loc
    return loc

#Load the Guide details
f = io.open('overbeek_2016_guides_s1.txt')
reader = csv.DictReader(f, delimiter='\t')
spacers = {}
for row in reader:
    guide_id = 'Overbeek_%d' % eval(row['Spacer '])
    guide_loc = row['Genomic location of spacer (hg19)']
    seq = row['Spacer sequence']
    spacers[guide_loc] = [guide_id, seq]
f.close()

#Load the SRA lookup
f = io.open('SraRunTable.txt')
reader = csv.DictReader(f, delimiter='\t')
srrs = {}
for row in reader:
    sample_id = row['Library_Name_s']
    sample_toks = sample_id.split('_')
    loc = sample_toks[1]
    if loc not in spacers:
        loc_before = loc
        loc = findNearbyLoc(loc, spacers)
    if loc not in spacers:
        continue
    offset = 0
    if sample_toks[-1][0] == 'R':
        offset -= 1
    cell_type = sample_toks[offset-2]
    timepoint = sample_toks[offset-1]
    srr = row['Run_s']
    if loc not in srrs:
        srrs[loc] = {}
    if cell_type not in srrs[loc]:
        srrs[loc][cell_type] = {}
    if timepoint not in srrs[loc][cell_type]:
        srrs[loc][cell_type][timepoint] = []
    srrs[loc][cell_type][timepoint].append(srr)
f.close()

#Pick which SRR to use for each guide (use lentviral d11 if available, else RNP 48hr)
srr_to_use, srr_to_use_ctrl = {}, {}
for spacer in spacers:
    if 'd11' in srrs[spacer]['K562']:
        srr_to_use[spacer] = srrs[spacer]['K562']['d11']
    elif '48hr' in srrs[spacer]['K562']:
        srr_to_use[spacer] = srrs[spacer]['K562']['48hr']
        print(spacers[spacer], '48hr')
    else:
        import pdb; pdb.set_trace()
    if 'WT' in srrs[spacer]['K562']:
        srr_to_use_ctrl[spacer] = srrs[spacer]['K562']['WT']
    else:
        raise Exception('No sample for spacer', spacers[spacer][0])
print(len(srr_to_use), len(srr_to_use_ctrl))

#Collect the sam_files
if 'overbeek_sam_files' not in os.listdir('.'):
    os.mkdir('overbeek_sam_files')
for spacer in spacers:
    for srr in srr_to_use[spacer]:
        if spacers[spacer][0] + '_' + srr  + '.sam' not in os.listdir('overbeek_sam_files'):
            cmd = 'sratoolkit.2.8.0-win64\\sratoolkit.2.8.0-win64\\bin\\sam-dump %s --output-file overbeek_sam_files\\%s_%s.sam' % (srr,spacers[spacer][0],srr)
            print(cmd); os.system(cmd)
for spacer in spacers:
    if spacers[spacer][0] + '.sam' not in os.listdir('overbeek_control_sam_files'):
        cmd = 'sratoolkit.2.8.0-win64\\sratoolkit.2.8.0-win64\\bin\\sam-dump %s --output-file overbeek_control_sam_files\\%s.sam' % (srr_to_use_ctrl[spacer],spacers[spacer][0])
        print(cmd); os.system(cmd)	

#Combine sam files for the same spacer
filenames = os.listdir('overbeek_sam_files')
for idx in range(1,97):
	id = 'Overbeek_%d' % idx
	fout = io.open('overbeek_sam_files/%s.sam' % id, 'w')
	for filename in filenames:
		if 'SRR' not in filename:
			continue
		if id + '_' not in filename:
			continue
		print(idx, filename)
		f = io.open('overbeek_sam_files/' + filename)
		for line in f:
			if line[0] == '@':
				continue
			fout.write(line)
	fout.close()
