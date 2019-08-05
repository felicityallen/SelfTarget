import io, os, csv, sys, re

SELFTARGET_ANALYSIS = '/lustre/scratch117/cellgen/team227/fa9/self-target/indel_analysis/'

def setHighDataDir(dirname):
    global SELFTARGET_ANALYSIS
    SELFTARGET_ANALYSIS = dirname

def getHighDataDir():
    return SELFTARGET_ANALYSIS

def getAllDataDirs():
    null_dir, old_dir, new_dir = [], [], []
    for data_dir in ['ST_April_2017/data/','ST_June_2017/data/','ST_Feb_2018/data/']:
        data_dir = SELFTARGET_ANALYSIS + data_dir
        null_dir.extend([data_dir + x for x in os.listdir(data_dir) if ('NULL' in x)])
        old_dir.extend([data_dir + x for x in os.listdir(data_dir) if ('DPI' in x  and ('Old' in x or '6O' in x or '12O' in x))])
        new_dir.extend([data_dir + x for x in os.listdir(data_dir) if ('DPI' in x  and ('Old' not in x and '6O' not in x and '12O' not in x))])
        old_dir.sort(); new_dir.sort()
    all_dir = null_dir + old_dir + new_dir
    all_dir = [x for x in all_dir if os.path.isdir(x)]
    return all_dir
    
def isNullDir(dirname):
    return 'NULL' in dirname or 'WT' in dirname

def getAllNullDirs():
    null_dir = set([getNullDir(x) for x in getAllDataDirs()])
    return [x for x in null_dir]

def isOldLib(full_dirname):
    dirname = full_dirname.split('/')[-1]
    if ('O' in dirname and 'BOB' not in dirname and 'CHO' not in dirname) or 'Old' in dirname or 'old' in dirname:
        return True
    return False

def getWTDir(full_dirname):
    highdir, dirname = '/'.join(full_dirname.split('/')[:-3]), '/'.join(full_dirname.split('/')[-3:])
    if isOldLib(dirname):
        return highdir + '/ST_Feb_2018/data/WT_12OA_DPI7'
    else: return highdir + '/ST_Feb_2018/data/WT_12NA_DPI7'
        
def getExpOligoFile(dirname):
    exp_data_dir = SELFTARGET_ANALYSIS + 'ST_June_2017/data/'
    if 'BOTH' in dirname:
        exp_oligo = exp_data_dir + 'exp_target_pam_both.fasta'
    elif isOldLib(dirname):
        exp_oligo = exp_data_dir + 'exp_target_pam_old.fasta'
    else:
        exp_oligo = exp_data_dir + 'exp_target_pam_new.fasta'
    return exp_oligo
    
def getNullTargetPamDetails(exptargets_filename, oligoid=None):
    pam_loc, pam_dir = None, None
    f = io.open(exptargets_filename)
    for line in f: 
        if line[0] != '>': continue
        id_str, a_pam_loc, a_pam_dir = line.split()
        id,indel,perc = id_str.split(':')
        if oligoid is not None and id[1:] != oligoid:
            continue
        if indel == '-':
            pam_loc, pam_dir = eval(a_pam_loc), a_pam_dir
            break
    return pam_loc, pam_dir

def getNullDir(dirname):
    if isOldLib(dirname):
        null_dir = 'ST_April_2017/data/NULL_Old'
    else:
        null_dir = 'ST_April_2017/data/NULL_New'
    return SELFTARGET_ANALYSIS + null_dir
    
def getShortDir(dirname):
    return dirname[len(SELFTARGET_ANALYSIS):]
    
def getDirLabel(dirname):
    dirlabel = dirname.split('/')[-3] + '_' + dirname.split('/')[-1]
    #dirlabel = dirlabel if dirlabel != 'ST_Feb_2018_K562_800x_7A_DPI7_may' else 'ST_Feb_2018_K562_800x_7A_DPI16_may'
    return dirlabel

def shortDirLabel(dir_label):
 return ' '.join(dir_label.split('_')[3:])
    
def getPamLocFile():
    return SELFTARGET_ANALYSIS + 'ST_June_2017/data/oligos_for_customarray_Dec2016_pamlocations.txt'

def parseSampleName(a_dirname):
    dirname = a_dirname[:-8] + 'DPI16_may' if 'K562_800x_7A_DPI7_may' in a_dirname else a_dirname
    toks = dirname.split('_') if '_' in dirname else dirname.split('-')
    dpi = eval([x for x in toks if 'DPI' in x][0][3:]) if 'DPI' in dirname else 0
    date = toks[1] if toks[1] in ['Feb','June','April'] else toks[1]
    cov = ''
    for x in ['1600x','800x','500x','1600X']: 
        if x in dirname: cov = x
    cellline = ''
    for x in ['K562','RPE1','CAS9','eCAS9','WT','TREX2','2A_TREX2','HAP1','E14TG2A','CHO','BOB', 'SIM']: 
        if x in dirname: cellline = x
    if cellline == 'CAS9': cellline = 'K562'
    if cellline == '2A_TREX2': cellline = 'TREX2_2A'
    if cellline == 'WT': cellline = ' WT'
    #if cellline == 'K562' and str.upper(cov) == '1600X': cellline = 'K562_1600x'
    month = ''
    for i,x in enumerate(['may','dec']):
        if x in dirname: month = str(i)+x
    virus = ''
    for x in ['12NA','12NB','7A','7B','6OA','6OB','12OA','12OB']:
        if x in dirname: virus = x
    return cellline, dpi, virus, cov, date, month, a_dirname

def sortSampleNames(dirnames):
    cl_tp_labels = [parseSampleName(dirname) for dirname in dirnames]
    cl_tp_labels.sort()

    return [x[-1] for x in cl_tp_labels]
    
def getSubdirs(dirname, withpath=True, mapped_dir='mapped_reads'):
    subdirs = [x for x in os.listdir(dirname + ('/%s/' % mapped_dir)) if os.path.isdir(dirname + ('/%s/' % mapped_dir) + x)]
    if withpath: return [dirname + ('/%s/' % mapped_dir) + x for x in subdirs]
    else: return subdirs
    
def getIndelSummaryFiles(subdir, withpath=True):
    if withpath: return [subdir + '/' + x for x in os.listdir(subdir) if 'mappedindelsummary' in x]
    else: return [x for x in os.listdir(subdir) if 'mappedindelsummary' in x]
    
def getDirNameFromSubdir(subdir):
    return '/'.join(subdir.split('/')[:-2])
    
def createResultDirectory(outdir, subdir, with_subdir = False):
    #e.g. to make outdir/ST_April_2017_K562_etc/Oligos_35
    if not os.path.isdir(outdir): os.mkdir(outdir)
    dirname = getDirNameFromSubdir(subdir)
    full_outdir = outdir + '/' + getDirLabel(dirname)
    if not os.path.isdir(full_outdir): os.mkdir(full_outdir)
    if with_subdir:
        full_outdir += ('/' + subdir.split('/')[-1])
        if not os.path.isdir(full_outdir): os.mkdir(full_outdir)
    return full_outdir

def getSimpleName(dirname, include_dpi=False):
    cellline, dpi, virus, cov, date, month, a_dirname = parseSampleName(dirname)
    if cellline == 'BOB': cellline = 'Human iPSC'
    if 'E14' in cellline: cellline = 'Mouse ESC'
    rep = 'D' if 'K562_1600x' in dirname else 'C' if '_CAS9' in dirname else ('A' if 'A' in virus else 'B')
    out = cellline + ' Rep ' + rep
    if include_dpi: out += ' DPI %d' % dpi 
    return out

def getSampleSelectors(include_wt = False):

    one_sel = lambda x: 'ST_Feb_2018_CAS9_12NA_1600X_DPI7' in x
    test_sel = lambda x: parseSampleName(x)[0] == 'K562' and not isOldLib(x) and 'DPI7' in x and ('_may' not in x and '1600x' not in x)
    test_exp = lambda x: 'CHO' in x or 'BOB' in x
    good_sel = lambda x: not (('RPE' in x and '7A' in x) or ('RPE' in x and 'dec' not in x) or ('DPI16' in x) or ('K562_800x_7A_DPI7_may' in x) or ('DPI7' in x and 'K562' in x and '1600' in x and '7B' in x) or ('CAS9' in x and '1600' in x and 'DPI3' in x) or ('WT' in x and 'DPI20' in x) or ('2A_TREX' in x and '12NB' in x and 'DPI3' in x))
    all_sel = lambda x: ('NULL' not in x or include_wt) and ('WT' not in x or include_wt) and good_sel(x)
    old_sel = lambda x: isOldLib(x) and all_sel(x)
    new_sel = lambda x: not isOldLib(x) and all_sel(x)
    tp_sel = lambda x: not isOldLib(x) and parseSampleName(x)[0] == 'K562' and all_sel(x)
    cl_sel = lambda x: not isOldLib(x) and parseSampleName(x)[1] == 7 and all_sel(x) and '2A_TREX' not in x
    old_k562 = lambda x: isOldLib(x) and parseSampleName(x)[0] == 'K562' and all_sel(x)
    k562_dpi7_new = lambda x: cl_sel(x) and parseSampleName(x)[0] == 'K562'
    k562_dpi7_old = lambda x: old_sel(x) and parseSampleName(x)[0] == 'K562'
    trex_dpi7 = lambda x: parseSampleName(x)[0] == 'TREX2' and parseSampleName(x)[1] == 7 and all_sel(x)
    k562_trex_new = lambda x: k562_dpi7_new(x) or trex_dpi7(x)

    return {'DPI7': cl_sel,
            'K562 New': k562_dpi7_new,
            'Old Scaffold': old_sel, 
            'K562_All_TP':tp_sel, 
            'New Scaffold': new_sel,
            'K562 TREX': k562_trex_new,
			'All': all_sel,
            'Concat New': lambda x: 'data_' in x and 'N' in x,
            'Concat Old': lambda x: 'data_' in x and 'O' in x}



