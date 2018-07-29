import io, os, sys, csv


def numMismatch(target1, target2):
    return sum([x != y for (x,y) in zip(target1, target2)])

def longestMatch(seq1, seq2):
    nm, max_nm, mh_end = 0, -1, -1
    for i, (x,y) in enumerate(zip(seq1, seq2)):
        if x != y: nm = 0
        else:
            nm += 1
            if nm > max_nm:
                max_nm = nm
                mh_end = i
    return max_nm, mh_end
    
def findMaxMicrohomology(target):
    mh_len, loc1, loc2 = 0, -1, -1
    for i in range(1,len(target)-1):
        max_nm, mh_end = longestMatch(target[:-i],target[i:])
        if max_nm > mh_len:
            mh_len = max_nm
            loc1, loc2 = mh_end-mh_len+1, i+mh_end-mh_len+1
            
    #Check!
    if target[loc1:loc1+mh_len] != target[loc2:loc2+mh_len]:
        print('Error:', target[loc1:loc1+mh_len], target[loc2:loc2+mh_len])
        import pdb; pdb.set_trace()
            
    return loc1, loc2, mh_len, target[loc1:loc1+mh_len]
    
def mergeMM(mh_mm1, mh_mm2):
    if mh_mm1 == mh_mm2: return []
    merged_mms = ['']
    for nt1,nt2 in zip(mh_mm1,mh_mm2):
        if nt1 != nt2: 
            merged_mms =  [x + nt1 for x in merged_mms] +  [x + nt2 for x in merged_mms]
        else:
            merged_mms = [x + nt1 for x in merged_mms]
    merged_mms.remove(mh_mm1)
    merged_mms.remove(mh_mm2)
    return merged_mms

def writeProperties(fout, row, best_match, fout_reads):

    loc1, loc2, mh_len, mh = findMaxMicrohomology(best_match['Target'])
    read1 = best_match['Target'][:loc1] + best_match['Target'][loc2:]

    read2 = row['Target'][:loc1+mh_len] + row['Target'][loc2+mh_len:]
    read3 = row['Target'][:loc1] + row['Target'][loc2:]
    mh_mm1 = row['Target'][loc1:loc1+mh_len]
    mh_mm2 = row['Target'][loc2:loc2+mh_len]
    merged_mms = mergeMM(mh_mm1, mh_mm2)

    fout_reads.write(u'>%s_%s_read1\n%s\n' % (row['ID'],best_match['ID'],read1))
    fout_reads.write(u'>%s_%s_read2\n%s\n' % (row['ID'],best_match['ID'],read2))
    fout_reads.write(u'>%s_%s_read3\n%s\n' % (row['ID'],best_match['ID'],read3))
    merged_reads = []
    for i, merged_mm in enumerate(merged_mms):
        merged_reads.append( row['Target'][:loc1] + merged_mm + row['Target'][loc2+mh_len:] )
        fout_reads.write(u'>%s_%s_read_merge%d\n%s\n' % (row['ID'],best_match['ID'],i,merged_reads[-1]))
    
    merged_str = '\t'.join(['%s\t%s' % (mm_str, read) for (mm_str, read) in zip(merged_mms, merged_reads)]) + '\t\t'*(2-len(merged_mms))
    str_args = (row['ID'],best_match['ID'],numMismatch(mh_mm1, mh_mm2), loc1, loc2, mh,mh_mm1,mh_mm2,read1, read2, read3, merged_str)
    fout.write(u'%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % str_args )

    
f = io.open('../ST_June_2017/data/self_target_oligos_details.csv')
comments_lookup = {}
fout = io.open('mh_mismatch_mappings.txt','w')
fout_reads = io.open('mh_mismatch_mappings_reads.fasta','w')
merged_str = '\t'.join(['Merged Mut %d MH\tMerged Mut %d Read' % (i+1,i+1) for i in range(2)])
fout.write(u'Oligo ID\tMapped Oligo Id\tNum Mismatches\tLeft Loc\tRight Loc\tOrig MH\tLeft Mut-MH\tRight Mut-MH\tExp Orig Read\tExp Left Read\tExp Right Read\t%s\n' % merged_str)
for row in csv.DictReader(f, delimiter='\t'):
    
    comment = row['Comments']
    
    #Collect Single MH Oligos (based on comments field)

    if comment[0] == '[' and 'PAM' not in comment:
        try:
            num_mh = len(eval(comment))
        except:
            num_mh = 1000
        if num_mh == 1:
            if comment not in comments_lookup:
               comments_lookup[comment] = []
            comments_lookup[comment].append(row)
        
    ctoks = comment.split()
    if ctoks[0] == 'modified':
        mh_str = ctoks[1]
        if mh_str not in comments_lookup: continue
        
        best_nmm, best_match = 100000, None
        for mh_row in comments_lookup[mh_str]:
            nmm = numMismatch(mh_row['Target'], row['Target']) 
            if nmm < best_nmm:
                best_nmm = nmm
                best_match = mh_row
                     
        writeProperties(fout, row, best_match, fout_reads)
        
f.close()
fout.close()
fout_reads.close()

 