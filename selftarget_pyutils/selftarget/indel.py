import io, os, csv, sys, re
import numpy as np

def tokFullIndel(indel):
    indel_toks = indel.split('_')
    indel_type, indel_details = indel_toks[0], ''
    if len(indel_toks) > 1:
        indel_details =  indel_toks[1]
    cigar_toks = re.findall(r'([CLRDI]+)(-?\d+)', indel_details)
    details, muts = {'I':0,'D':0,'C':0}, []
    for (letter,val) in cigar_toks:
        details[letter] = eval(val)
    if len(indel_toks) > 2 or (indel_type == '-' and len(indel_toks) > 1):
        mut_toks = re.findall(r'([MNDSI]+)(-?\d+)(\[[ATGC]+\])?', indel_toks[-1])
        for (letter,val,nucl) in mut_toks:
            if nucl == '':
                nucl = '[]'
            muts.append((letter, eval(val), nucl[1:-1]))
    if indel_type[0] == '-':
        isize = 0
    else:
        isize = eval(indel_type[1:])
    return indel_type[0],isize,details, muts
    
def computeReadLength(indel, oligo_indel):
    read_length = 79
    for ind in [indel, oligo_indel]:
        itype, isize, details, muts = tokFullIndel(ind)    
        if itype == 'I': read_length += isize
        elif itype == 'D': read_length -= isize
        for mut in muts:
            if mut[0] == 'I': read_length += mut[1]
            elif mut[0] == 'D': read_length -= mut[1]
    return read_length
 
def indelOutofGuideSeedPAM(indel):
    if indel == '-':
        return True
    itype, isize, details, muts = tokFullIndel(indel)
    if itype != '-':
        left, right = details['L'], details['R']
        if (left > -15 and left < 6) or (right > -14 and right < 7):
            return False
    for (letter, pos, nucl) in muts:
        if letter == 'M' or letter == 'S':
            if pos > -14 and pos < 7 and not (letter == 'M' and pos == 0):
                return False
    return True	
    
def isDetectable(indel, target, pam_det):
    pam_loc, pam_dir = pam_det
    itype, isize, details, muts = tokFullIndel(indel)
    if pam_dir == 'REVERSE':
        space_right, space_left = pam_loc-3, len(target)-pam_loc-3
    else:
        space_right, space_left = len(target)-pam_loc-3, pam_loc-3
    return ((abs(details['L']) < space_left) and (details['R'] < space_right))