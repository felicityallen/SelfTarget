import io, sys, os

from selftarget.util import getIndelMapExe, getPythonCmd

if __name__ == '__main__':

    if len(sys.argv) < 6 or len(sys.argv)>11:
        print('Usage: indelmap_subdir.py <dirname> <nulldir> <subdir> <files_per_part> <(opt)fastafile> <(opt)reformat_only> <map_dir> <maxcutdist> <indelmap_exe> <py_cmd>')
    else:
        dirname = sys.argv[1]
        null_dir = sys.argv[2]
        subdir = sys.argv[3]
        fpp = eval(sys.argv[4])
        selected_file = sys.argv[5] if len(sys.argv) >= 6 else ''
        if selected_file == '-': selected_file = ''
        reformat_only = eval(sys.argv[6]) if len(sys.argv) >= 7 else False
        map_dir = sys.argv[7] if len(sys.argv) >= 8 else '/mapped_reads/'
        maxcutdist = eval(sys.argv[8]) if len(sys.argv) >= 9 else 4
        indelmap_exe = sys.argv[9] if len(sys.argv) >= 10 else getIndelMapExe()
        py_cmd = sys.argv[10] if len(sys.argv) >= 11 else getPythonCmd()

        part = eval(os.getenv("LSB_JOBINDEX"))-1  if "LSB_JOBINDEX" in os.environ else 0
       
        fasta_files = [x for x in os.listdir(dirname + map_dir + subdir) if x[-6:]=='.fasta']
        fasta_files.sort()
        for (idx, fasta_file) in enumerate(fasta_files):
            if (idx >= part*fpp and idx < (part+1)*fpp) or fasta_file == selected_file:

                prefix = subdir + '/' + fasta_file[:-6]
                nullexpfilename = map_dir + prefix + '_exptargets.txt'
                outfilename = map_dir + prefix + '_mappedindels.txt'  
                            
                if not reformat_only:
                    cmd = indelmap_exe + ' %s/%s %s/%s %s/%s 0 %d' % (dirname, map_dir + prefix + fasta_file[-6:], null_dir, nullexpfilename, dirname,outfilename, maxcutdist) 
                    print(cmd)
                    os.system(cmd)
                
                cmd = py_cmd + ' reformat_indel_profile.py %s' % dirname + '/' + map_dir + '/' + prefix 
                print(cmd)
                os.system(cmd)
