import io, os, csv, sys, re
import numpy as np
import pandas as pd
from Bio import SeqIO
from selftarget.data import getExpOligoFile
import Bio

def loadPamLookup(filename, adj=20):
    lookup = {}
    f = io.open(filename)
    reader = csv.reader(f, delimiter='\t')
    for toks in reader:
            lookup[toks[0]] = (eval(toks[1])-adj, toks[2])
    f.close()
    return lookup
    
def loadOligosByBarcode( oligo_file ):
    lookup = {}
    for record in SeqIO.parse(oligo_file, 'fasta'):
        bc1, bc2 = str(record.seq[:10]),str(record.seq[-10:])
        oligo_id = record.id.split(' ')[0]
        oligo_idx = eval(record.id[5:].split('_')[0])
        lookup[(bc1,bc2[-3:])] = (record.id, oligo_idx)	#We want all BC1 and at least the last 3 of BC2 OR
        lookup[(bc1[:3],bc2)] = (record.id, oligo_idx)  #        all BC2 and at least the first 3 of BC1
    return lookup
    
#Removes the guide sequence from the ID
def getShortOligoId(full_oligo_id):
    return full_oligo_id.split('_')[0]
    
def getFileForOligoIdx( oligo_idx, ext='.fasta' ):
    first_oligo, last_oligo  = int(oligo_idx/50)*50, int((oligo_idx+50)/50)*50-1
    filedir = 'Oligos_' + str( int(oligo_idx / 1000) )
    filename = 'Oligos_' + str(first_oligo) + '-' + str(last_oligo) + ext
    return filedir, filename
    
def getSummaryFileSuffix(oligo_id):
    oligo_idx = getOligoIdxFromId(oligo_id)
    subdir, sumfilename = getFileForOligoIdx(oligo_idx, ext='_mappedindelsummary.txt')
    return subdir + '/' + sumfilename

def getOligoIdsFromFile(filename):
    oligo_ids = []
    f = io.open(filename)
    for line in f:
        if line[:3] == '@@@':
            oligo_ids.append(line[3:-1])
    f.close()
    return oligo_ids

def splitOligoFastaFile(filename, oligo_ids, filepath=None):
    if filepath is not None: out_filename = filepath + '/' + filename.split('/')[-1][:-6]
    else: out_filename = filename[:-6]
    out_filename +=  ('_' + oligo_ids[0] + '_split.fasta') 
    fout = io.open(out_filename,'w')
    f = io.open(filename)
    correct_id = False
    for line in f:
        if line[0] == '>':
            line_oligo_id = line[1:].split('.')[0]
            correct_id = (line_oligo_id in oligo_ids)
        if correct_id: fout.write(line) 
    f.close()
    return out_filename            
    
def getOligoIdxFromId(oligo_id):
    return eval(oligo_id.replace('_','')[5:])

def getOligoIdsFromMappedFastaFile(filename, return_counts=False):
    oligo_counts = {}
    f = io.open(filename)
    for line in f:
        if line[0] == '>':
            oligo_id = line[1:].split('.')[0]
            if oligo_id not in oligo_counts: oligo_counts[oligo_id] = 0
            oligo_counts[oligo_id] += 1 
    f.close()
    if not return_counts:
        return [x for x in oligo_counts]
    else:
        return oligo_counts

def getFullFilename( id, mapped_reads_dir = 'mapped_reads'):
    idx = eval(id.split('_')[0][5:])/1000
    return '%s/Oligos_%d/%s_indelprofile.txt' % (mapped_reads_dir, idx, id)
    
def loadAllOligoDetails(oligo_detail_dir='../ST_June_2017/data'):
    details = {}
    f = io.open(oligo_detail_dir + '/self_target_oligos_details_with_pam_details.csv')
    details = {row['ID']: row for row in csv.DictReader(f, delimiter='\t')}
    f.close()
    return details

def loadExpOligoLookup(subdir, exp_oligo_file=None):
    if exp_oligo_file is None:
        exp_oligo_file = getExpOligoFile('/'.join(subdir.split('/')[:-2]))
    lookup = {}
    for record in SeqIO.parse(exp_oligo_file, "fasta"):
        oligo_id, pam_loc, pam_dir = str(record.description).split()
        oligo_id = oligo_id.split('_')[0]
        seq = str(record.seq)
        oligo_idx = eval(oligo_id[5:])
        filedir, filename = getFileForOligoIdx(oligo_idx, ext='')
        if filename not in lookup: lookup[filename] = []
        lookup[filename].append((oligo_id, pam_loc, pam_dir, seq))
    return lookup

def hasMHLenNOrLonger(target, pam_loc, pam_dir, N):
    cut_site =  pam_loc-3 if pam_dir == 'FORWARD' else pam_loc+3
    for i in range(cut_site, len(target)-N):
        if target[i:i+N] in target[:cut_site]:
            return True
    return False

def partitionGuides(lib='Both', oligo_detail_dir='../ST_June_2017/data'):
    guide_matches_target = lambda row: (row['Guide'][1:] in row['Target']) or \
                                       (row['Guide'][1:] in Bio.Seq.reverse_complement(row['Target']))

    f = io.open(oligo_detail_dir + '/self_target_oligos_details_with_pam_details.csv')
    rdr = csv.DictReader(f, delimiter='\t')
    partitions = {'All':set(),'Non-Targeting no MH>3':set(), 'Real Guides':set(), 'PAM':set(),'Mismatch':set(), 'Non-Targeting with MH>3':set(), 'Non-Targeting':set(),'With MH>3':set(),'No MH>3':set()}
    non_targeting, real_guides, mismatch = [],[],[]
    for row in rdr:
        oligo_id = ''.join(row['ID'].split('_'))
        if lib == 'Old' and not eval(row['Old']): continue
        if lib == 'New' and eval(row['Old']): continue
        partitions['All'].add(oligo_id)
        if 'PAM' in row['Comments']:
            partitions['PAM'].add(oligo_id)
        elif 'mismatch' in row['Comments']:
            partitions['Mismatch'].add(oligo_id)
        elif not guide_matches_target(row):
            partitions['Mismatch'].add(oligo_id)
        else:
            has_mh = hasMHLenNOrLonger(row['Target'], eval(row['PAM Location']),row['PAM Direction'], 4)
            if has_mh:
                partitions['With MH>3'].add(oligo_id)
            else:
                partitions['No MH>3'].add(oligo_id)

            if 'Guide' in row['Comments']:
                partitions['Real Guides'].add(oligo_id)
            else:
                partitions['Non-Targeting'].add(oligo_id)
                if has_mh:
                    partitions['Non-Targeting with MH>3'].add(oligo_id)
                else:
                    partitions['Non-Targeting no MH>3'].add(oligo_id)
    return partitions    
    
def getOldLookup():
    f = io.open('../ST_June_2017/data/self_target_oligos_details.csv')
    rdr = csv.DictReader(f, delimiter='\t')
    lookup = {''.join(row['ID'].split('_')): eval(row['Old']) for row in rdr}
    f.close()
    return lookup
    
def getNullTargetPamDetails(exptargets_filename, oligoid=None):
    f = io.open(exptargets_filename)
    found = False
    for line in f:
        if line[0] == '>':
            id_str, pam_loc, pam_dir = line.split()
            id,indel,perc = id_str.split(':')
            if oligoid is not None and id[1:] != oligoid:
                continue
            if indel == '-':
                found = True
                break
    if found:
        return eval(pam_loc), pam_dir
    else:
        return None, None

def loadOldNewMapping(mapping_file_folder='../ST_June_2017/data'):
    data = pd.read_csv(mapping_file_folder + '/' + 'oligo_mapping_old_to_new.txt', sep=' ')
    data = data.loc[~data['New'].isin(["matchlessA","matchlessC","matchlessT"])]
    return data