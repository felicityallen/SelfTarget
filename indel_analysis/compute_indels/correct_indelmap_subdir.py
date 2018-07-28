import io
import sys
import os
import csv

def loadMappingsFile(filename):
    f = io.open(filename)
    #e.g D00570:327:HGHC7ADXY:2:1110:11908:26847 1:N:0:TACTTCGG  Oligo82002_GGGGATAATCCAACATGAGG D1_L-4C7R5      M-21D2S-22M15M18D2S19M26D1S33   0.04    8       1
    lookup = {toks[0].split()[0]: toks for toks in csv.reader(f, delimiter='\t')}
    f.close()
    return lookup

if __name__ == '__main__':

    if len(sys.argv) != 6:
        print 'Usage: correct_indelmap_subdir.py <dirname> <nulldir> <subdir> <part> <files_per_part>'
    else:
        dirname = sys.argv[1]
        null_dir = sys.argv[2]
        subdir = sys.argv[3]
        part = eval(sys.argv[4])-1
        fpp = eval(sys.argv[5])
        
        fasta_files = [x for x in os.listdir(dirname + '/mapped_reads/' + subdir) if x[-6:]=='.fasta']
        fasta_files.sort()
        for (idx, fasta_file) in enumerate(fasta_files):
            if (idx >= part*fpp and idx < (part+1)*fpp):

                prefix = dirname + '/mapped_reads/' + subdir + '/' + fasta_file[:-6]
                mindels_file = prefix + '_mappedindels.txt'
                mappings_file = prefix + '_mappings.txt'  
                    
                mappings = loadMappingsFile(mappings_file)
                f = io.open(mindels_file)
                fout = io.open(mindels_file + '__mod', 'w')
                for toks in csv.reader(f, delimiter='\t'):
                    if '@@@' in toks[0]: continue
                    #e.g Oligo82029.D00570:327:HGHC7ADXY:1:1104:18883:81216      Oligo82029:-:68.794     -       -       0.05    10      1
                    if toks[0].split('.')[0] == toks[1].split(':')[0]:
                        fout.write(u'\t'.join(toks) + '\n')
                    else:
                        oligo_id, read_id = toks[0].split('.')
                        mapping_toks = mappings[read_id]
                        fout.write(u'%s.%s\t%s:-:100.0\t%s\t%s\t%s\t%s\t1\n' % (oligo_id, read_id, oligo_id, mapping_toks[2],mapping_toks[3],mapping_toks[4],mapping_toks[5]))
                f.close()
                fout.close()
                
                cmd = 'mv %s %s' % (mindels_file + '__mod', mindels_file)
                print cmd
                os.system(cmd)
        
                cmd = 'python reformat_indel_profile.py %s' % prefix 
                print cmd
                os.system(cmd)
