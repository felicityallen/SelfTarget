import io
import sys
from Bio import SeqIO
import os
import numpy as np
import re
import csv
import Bio.Seq

from selftarget.util import getIndelMapExe

sys.path.append('../compute_indels')
from compile_mapped_null_profiles import compileMappedNull, convertToExpFile
from reformat_indel_profile import reformatIndelProfile

def loadFastqReads( filename, id, ftype='fastq'):
    lookup = {}
    for record in SeqIO.parse(filename,ftype):
        readid = str(record.id)
        if '.' not in readid:
            readid = id + '.' + readid
        lookup[readid] = str(record.seq)
    return lookup

def fetchOrigPamAndTemplate(template_file):
    f = io.open(template_file)
    toks = f.readline().split()
    id, pam_loc, pam_dir = toks[0][1:], eval(toks[1]), toks[2]
    seq = f.readline()[:-1]
    return {id.split('_')[0]: (pam_loc, pam_dir)}, seq

def filterMappings(mappings_file,  output_mappings_file):
    f = io.open(mappings_file)
    fout = io.open(output_mappings_file, 'w')
    for toks in csv.reader(f,delimiter='\t'):
        if toks[0][:3] == '@@@' or toks[1] == '': continue
        fout.write(u'\t'.join(toks) + '\n')
    fout.close()
    f.close()
    
def numMismatch(seq1,seq2):
    return sum([x != y for (x,y) in zip(seq1,seq2)])
    
def trimRead( read_seq, template_seq ):
    start_idx = read_seq.find(template_seq[:20])
    stop_idx = read_seq.find(template_seq[-20:])
    #Search with mismatches if needed
    if start_idx < 0:
        start_idx = read_seq.find(template_seq[:10])
        if start_idx < 0 or numMismatch(template_seq[:20], read_seq[start_idx:start_idx+20]) > 2:
            start_idx = read_seq.find(template_seq[10:20])-10
            if start_idx < 0 or numMismatch(template_seq[:20], read_seq[start_idx:start_idx+20]) > 2:
                start_idx = -1
    if stop_idx < 0:
        stop_idx = read_seq.find(template_seq[-10:])-10
        if stop_idx < 0 or numMismatch(template_seq[-20:], read_seq[stop_idx:stop_idx+20]) > 2:
            stop_idx = read_seq.find(template_seq[-20:-10])
            if stop_idx < 0 or numMismatch(template_seq[-20:], read_seq[stop_idx:stop_idx+20]) > 2:
                stop_idx = -1
    if start_idx >= 0 and stop_idx >= 0:
        return read_seq[start_idx:stop_idx+20]
    return ''
    
def trimReadsToTemplate(fastq_file, output_fasta, template_seq, id):
    fout = io.open(output_fasta, 'w')
    count, total = 0, 0
    for record in SeqIO.parse(fastq_file,'fastq'):
        trimmed_seq = trimRead( str(record.seq), template_seq )
        if trimmed_seq != '':
            fout.write(u'>%s\n%s\n' % (id + '.' + str(record.id),trimmed_seq))
            count += 1
        total += 1
    fout.close()

def computeOverbeekIndelProfiles(highdir='.', selected_id = None):

    nulldir = highdir + '/overbeek_control_fastq_files'
    testdir = highdir + '/overbeek_fastq_files'

    for idx in range(1,97):
        id = 'Overbeek%d' % idx
        if selected_id is not None and id != selected_id:
            continue
        
        fastq_file = testdir + '/%s.fastq' % id
        null_fastq_file = nulldir + '/%s.fastq' % id
        template_file = highdir + '/overbeek_template_files/%s_template.fasta' % id
        null_mappings_file = nulldir + '/%s_mappings.txt' % id
        mapped_file = testdir + '/%s_mappedindels.txt' % id
    
        #Compute the Null Profile and resulting expected templates
        cmd = getIndelMapExe() + ' %s %s %s 0' % (null_fastq_file, template_file, null_mappings_file[:-4] + '_unfilt.txt')
        print(cmd); os.system(cmd)
        filterMappings(null_mappings_file[:-4] + '_unfilt.txt', null_mappings_file)

        null_reads = loadFastqReads(null_fastq_file, id)
        pam_lookup, template_seq = fetchOrigPamAndTemplate(template_file)
        compileMappedNull(nulldir + '/%s' % id, null_reads, pam_lookup, {})
        convertToExpFile(nulldir + '/%s_nullsummary.txt' % id, nulldir + '/%s_exptargets.txt' % id, discard_long=False)
        
        #Compute the Indel Profile, by mapping the test reads against the null templates
        trimReadsToTemplate(fastq_file, fastq_file[:-6] + '_trimmed.fasta', template_seq, id)
        cmd = getIndelMapExe() + ' %s %s %s 0' % (fastq_file[:-6] + '_trimmed.fasta', nulldir + '/%s_exptargets.txt' % id, mapped_file[:-4] + '_unfilt.txt')
        print(cmd); os.system(cmd)
        filterMappings(mapped_file[:-4] + '_unfilt.txt', mapped_file)
        
        reads = loadFastqReads(fastq_file[:-6] + '_trimmed.fasta', id, ftype='fasta')
        reformatIndelProfile(testdir + '/%s' % id, reads)
        
if __name__ == '__main__':

    selected_id = None
    if len(sys.argv) == 2:
        selected_id = sys.argv[1]
    computeOverbeekIndelProfiles(selected_id=selected_id)



