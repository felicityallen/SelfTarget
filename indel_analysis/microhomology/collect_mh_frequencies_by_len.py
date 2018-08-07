import io,sys,os,csv

from selftarget.oligo import loadAllOligoDetails
from selftarget.data import getHighDataDir, setHighDataDir

def collectMhOfLen(filename, mh_len, fout):

    det = loadAllOligoDetails(oligo_detail_dir=getHighDataDir() + '/ST_June_2017/data')
    oligo_details = {'Oligo' + x.split('_')[-1]:val for x,val in det.items()}

    indels_to_write = []
    max_reads, len_mh_max_reads, left_max_reads, right_max_reads, max_indel = 0, -1, -1, -1, ''
    f = io.open(filename)
    
    #Collect indels of the right length, write out with details of MH indel with max reads for that oligo
    for toks in csv.reader(f, delimiter='\t'):
        
        #Next Oligo (write out last)
        if toks[0][:3] == '@@@':
            if len(indels_to_write) > 0:
                oligo_line = u'%d\t%d\t%s\t%d\t%d\t%d\t%d' % (accpt_reads, accpt_nonnull_reads,  max_indel, max_reads, len_mh_max_reads, left_max_reads, right_max_reads)
                for indel_line in indels_to_write:
                    fout.write(u'%s\t%s\n' % (indel_line,oligo_line) )
            ctoks = toks[0][3:].split(':')
            oligo_id = ctoks[0]
            target = oligo_details[oligo_id]['Target']
            accpt_reads, accpt_nonnull_reads = eval(ctoks[1]),eval(ctoks[2])
            max_reads, len_mh_max_reads, left_max_reads, right_max_reads, max_indel = 0, -1, -1, -1, ''
            indels_to_write = []
            continue
        
        #MH details, collect MH's of correct length, and also track details of MH indel with max reads
        left, right, c_mh_len, indel, reads = eval(toks[0]), eval(toks[1]), eval(toks[2]), toks[3], eval(toks[-1])
        l_mh_seq, r_mh_seq = target[left:left+c_mh_len], target[right:right+c_mh_len]
        assert(l_mh_seq == r_mh_seq)
        gc_content = sum([x in ['G','C'] for x in l_mh_seq])*100.0/len(l_mh_seq)
        if reads > max_reads:
            max_reads, len_mh_max_reads, left_max_reads, right_max_reads, max_indel  = reads, c_mh_len, left, right, indel
        if c_mh_len != mh_len: continue
        indels_to_write.append(u'%s\t%s\t%d\t%d\t%d\t%.1f' % (oligo_id, indel, reads, left, right, gc_content))
    
    #Write last Oligo (if needed)
    if len(indels_to_write) > 0:
        oligo_line = u'%d\t%d\t%s\t%d\t%d\t%d\t%d' % (accpt_reads, accpt_nonnull_reads,  max_indel, max_reads, len_mh_max_reads, left_max_reads, right_max_reads)
        for indel_line in indels_to_write:
            fout.write(u'%s\t%s\n' % (indel_line,oligo_line) )
        
def collectMhFrequenciesOfLen(results_subdir, mh_len, outfile):
    fout = io.open(outfile, 'w')
    fout.write(u'Oligo ID\tIndel\tIndel Reads\tLeft Position\tRight Position\tGC Content\tTotal Reads\tNon-Null Reads\tMost Frequent MH Indel\tReads of most frequent MH\tLength of MH with most reads\tLeft Position of most frequent MH\tRight Position of most frequent MH\n')
    files = [results_subdir + '/' + x for x in os.listdir(results_subdir)]
    for filename in files:
        collectMhOfLen(filename, mh_len, fout)
    fout.close()

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print('Usage: collect_mh_frequencies_by_len.py results_mh_len highdir subdir')
    else:

        mh_len = eval(sys.argv[1])
        highdir = sys.argv[2]
        results_subdir = sys.argv[3]
        subdir = '/'.join(results_subdir.split('/')[-2:])
        setHighDataDir(highdir)
        
        if not os.path.isdir(results_subdir):
            raise Exception('No such directory:' + results_subdir)

        out_dir = highdir + '/mh_freqs_by_len' + '/' + subdir
        if not os.path.isdir(out_dir): os.makedirs(out_dir)
    
        collectMhFrequenciesOfLen(results_subdir, mh_len, out_dir + '/mh_indels_of_len_%d.txt' % mh_len)
        
