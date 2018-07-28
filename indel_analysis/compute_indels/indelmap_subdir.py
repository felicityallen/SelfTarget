import io
import sys
import os

if __name__ == '__main__':

	if len(sys.argv) < 6 or len(sys.argv)>10:
		print 'Usage: indelmap_subdir.py <dirname> <nulldir> <subdir> <part> <files_per_part> <(opt)fastafile> <(opt)reformat_only> <map_dir> <maxcutdist>'
	else:
		dirname = sys.argv[1]
		null_dir = sys.argv[2]
		subdir = sys.argv[3]
		part = eval(sys.argv[4])-1
		fpp = eval(sys.argv[5])
                if len(sys.argv) >= 7:
                    selected_file = sys.argv[6]
                else: selected_file = ''
                if selected_file == '-': selected_file = ''		
                reformat_only = eval(sys.argv[7]) if len(sys.argv) >= 8 else False
                map_dir = sys.argv[8] if len(sys.argv) >= 9 else '/mapped_reads/'
                maxcutdist = eval(sys.argv[9]) if len(sys.argv) >= 10 else 4

                fasta_files = [x for x in os.listdir(dirname + map_dir + subdir) if x[-6:]=='.fasta']
		fasta_files.sort()
		for (idx, fasta_file) in enumerate(fasta_files):
			if (idx >= part*fpp and idx < (part+1)*fpp) or fasta_file == selected_file:

				prefix = subdir + '/' + fasta_file[:-6]
				nullexpfilename = map_dir + prefix + '_exptargets.txt'
				outfilename = map_dir + prefix + '_mappedindels.txt'  
                            
                                if not reformat_only:

                                    cmd = './run_indelmap.sh %s/%s %s/%s %s/%s %d' % (dirname, map_dir + prefix + fasta_file[-6:], null_dir, nullexpfilename, dirname,outfilename, maxcutdist) 
                                    print cmd
                                    os.system(cmd)
				
				cmd = 'python reformat_indel_profile.py %s' % dirname + '/' + map_dir + '/' + prefix 
				print cmd
				os.system(cmd)
