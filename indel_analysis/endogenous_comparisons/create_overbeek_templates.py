import io, sys, os, re, csv
import numpy as np
import Bio.Seq
from Bio import SeqIO

from selftarget.data import getHighDataDir

def extractTemplateSequenceAndPamLoc( ctrl_sam_file, location, oligo_id, spacer_seq, primer ):
    chr = location.split(':')[0]
    mid_idx = 0.5*sum([eval(x) for x in location.split(':')[1].split('-')])
    perfect_reads = {}
    f = io.open(ctrl_sam_file)
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
        cigar_string = toks[5]
        if cigar_string[-1] != 'M':
            continue
        try:
            eval(cigar_string[:-1])		#Should be a number followed by M, (ie. all matches) else discard
        except:
            continue
        full_read = str(toks[9])
        if full_read not in perfect_reads:
            perfect_reads[full_read] = 0
        perfect_reads[full_read] += 1
    f.close()
    
    #Use the most common perfect sequence as the template
    counts = [(perfect_reads[seq],seq) for seq in perfect_reads]
    counts.sort(reverse=True)
    if len(counts) == 0:
        print('No perfect sequences for', oligo_id)
        return '',-1,''
        
    template_seq = counts[0][1]
    if counts[0][0] < 1000:
        print('Low count template sequence for', oligo_id, counts[0][0])

    #Find the Pam location and direction
    idx = template_seq.find(spacer_seq)
    if idx < 0:
        print('Could not find spacer sequence within template for', oligo_id)
        print(template_seq)
        print(spacer_seq)
        import pdb; pdb.set_trace()
    if spacer_seq[-2:] != 'GG' and spacer_seq[:2] == 'CC':
        pam_dir = 'REVERSE'
        pam_loc = idx + 3
    elif spacer_seq[-2:] == 'GG' and spacer_seq[:2] != 'CC':
        pam_dir = 'FORWARD'
        pam_loc = idx + 20
    elif spacer_seq[:-3] in primer:
        pam_dir = 'FORWARD'
        pam_loc = idx + 20
    elif Bio.seq.reverse_complement(spacer_seq[3:]) in primer:
        pam_dir = 'REVERSE'
        pam_loc = idx + 3
    else:
        raise Exception('Ambigous PAM location, please check more carefully...', oligo_id)
       
    return template_seq, pam_loc, pam_dir
    
    
def loadLocationSpacerLookup():
    f = io.open(getHighDataDir() + '/overbeek_2016_guides_s1.txt')
    reader = csv.DictReader(f, delimiter='\t')
    lookup = {'Overbeek%d' % eval(row['Spacer ']): (row['Genomic location of spacer (hg19)'], row['Spacer sequence'], row['sgRNA primer']) for row in reader}
    f.close()
    return lookup
    

def createOverbeekTemplates(selected_id = None):

    ctrl_samdir = getHighDataDir() + '/overbeek_control_sam_files'
    output_template_dir = getHighDataDir() + '/overbeek_template_files'
    if not os.path.isdir(output_template_dir):
        os.mkdir(output_template_dir)
        
    lookup = loadLocationSpacerLookup()
    
    f = io.open(getHighDataDir() + '/overbeek_self_targets.csv')
    reader = csv.reader(f, delimiter='\t')
    
    for toks in reader:
        idx = eval(toks[-1].split()[-1])
        id, samid = 'Overbeek%d' % idx, 'Overbeek_%d' % idx
        if selected_id is not None and selected_id != id:
            continue
        loc, spacer_seq, primer = lookup[id]
        fout = io.open(output_template_dir + '/%s_template.fasta' % id,'w')
        ctrl_sam_file = ctrl_samdir + '/%s.sam' % samid
        template_seq, pam_loc, pam_dir = extractTemplateSequenceAndPamLoc(ctrl_sam_file, loc, id, spacer_seq, primer )
        fout.write(u'>%s_%s %d %s\n%s\n' % (id, spacer_seq, pam_loc, pam_dir, template_seq))
        fout.close()

if __name__ == '__main__':
    createOverbeekTemplates()