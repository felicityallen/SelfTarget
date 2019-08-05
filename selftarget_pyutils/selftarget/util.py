import io, os, csv, sys, re
from Bio import SeqIO
from selftarget.data import getHighDataDir, getAllDataDirs, isNullDir, getShortDir, getSubdirs, getSampleSelectors, sortSampleNames, parseSampleName
from selftarget.oligo import partitionGuides
import pandas as pd

RUN_LOCAL = False
PYTHON_CMD = '~/run_python.sh'
LOG_DIR = 'log_files'
PLOT_DIR = 'plots'
PICKLE_DIR = 'pickle'
PEAR_EXE = '~/pear-0.9.10-bin-64/pear-0.9.10-bin-64'
INDELMAP_EXE = '/lustre/scratch117/cellgen/team227/fa9/indelmap/bin/indelmap'
INDELMH_EXE = '/lustre/scratch117/cellgen/team227/fa9/indelmap/bin/indelmh'
INDELGENI1_EXE = '/lustre/scratch117/cellgen/team227/fa9/indelmap/bin/indelgen_i1'
INDELGEN_EXE = '/lustre/scratch117/cellgen/team227/fa9/indelmap/bin/indelgen'

def setIndelGenExe(exe):
    global INDELGEN_EXE
    INDELGEN_EXE = exe

def setIndelMapExe(exe):
    global INDELMAP_EXE
    INDELMAP_EXE = exe

def setIndelGenI1Exe(exe):
    global INDELGENI1_EXE
    INDELGENI1_EXE = exe

def setIndelMhExe(exe):
    global INDELMH_EXE
    INDELMH_EXE = exe

def setRunLocal(val):
    global RUN_LOCAL
    RUN_LOCAL = val

def setPythonCmd(cmd):
    global PYTHON_CMD
    PYTHON_CMD = cmd

def setPlotDir(plot_dir):
    global PLOT_DIR
    PLOT_DIR = plot_dir

def setPearExe(pear_exe):
    global PEAR_EXE
    PEAR_EXE = pear_exe

def getPythonCmd():
    return PYTHON_CMD

def getRunLocal():
    return RUN_LOCAL

def getLogDir():
    out_dir = LOG_DIR
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    return out_dir    
    
def getIndelMapExe():
    return INDELMAP_EXE

def getIndelMhExe():
    return INDELMH_EXE

def getIndelGenExe():
    return INDELGEN_EXE

def getIndelGenI1Exe():
    return INDELGENI1_EXE

def getPlotDir():
    out_dir = PLOT_DIR
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    return out_dir 
    
def getPearExe():
    return PEAR_EXE

def getPickleDir():
    out_dir = PICKLE_DIR
    if out_dir not in os.listdir('.'):
        os.mkdir(out_dir)
    return out_dir 

def getCommonKeys(dicts):
    if len(dicts) == 0: return set()
    common_keys = set(dicts[0].keys())
    for a_dict in dicts[1:]:
        common_keys = common_keys.intersection(set(a_dict.keys()))
    return common_keys 
    
def mergeSamples(all_result_outputs, cols_to_sum, merge_on='Oligo Id',data_label='Data'):
    merge_cols = [merge_on] if type(merge_on) != list else merge_on
    if len(cols_to_sum) > 0: datas = [x[0][data_label][merge_cols + cols_to_sum] for x in all_result_outputs]
    else: datas = [x[0][data_label] for x in all_result_outputs]
    if len(datas) > 1:
        merged_data = pd.merge(datas[0],datas[1],how='outer',on=merge_cols, suffixes=['', ' 2']).fillna(0)
    else:
        merged_data = datas[0]
    for i, data in enumerate(datas[2:]):
        merged_data = pd.merge(merged_data, data,how='outer',on=merge_on, suffixes=['', ' %d' % (i+3)]).fillna(0)
    suffix = lambda i: ' %d' % (i+1) if i > 0 else ''
    for col in cols_to_sum:
        merged_data[col + ' Sum'] = merged_data[[col + suffix(i) for i in range(len(datas))]].sum(axis=1)
    return merged_data

def runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, out_prefix, extra_cmd = '', queue='normal', numj=0):
    if idx >= start_idx and idx <= stop_idx:
        if extra_cmd != '':
            print(extra_cmd)
            os.system(extra_cmd)
        if RUN_LOCAL:
            if numj > 1:
                import pdb; pdb.set_trace()
                print('Num Job>1, this probably wont work!!')
                for j in range(numj):
                    cmd_j = "%s %d > %s/%s.txt" % (cmd, j, out_dir, out_prefix)
                    print(cmd_j); os.system(cmd_j)
            else:
                #cmd = "%s > %s/%s.txt" % (cmd, out_dir, out_prefix)
                print(cmd); os.system(cmd)
        else:
            jcmd = (' -J%s[1-%d]' % (out_prefix, numj)) if numj > 0 else ''
            cmd = "bsub -R'select[mem>8000] rusage[mem=8000]' -M8000 -G teamparts -n 1%s -o %s/%s.log -q %s " % (jcmd, out_dir, out_prefix, queue) + cmd
            print(cmd); os.system(cmd)
    return idx + 1
    
def loadFastaReadsById( filename ):
    lookup = {}
    for record in SeqIO.parse(filename,'fasta'):
        lookup[str(record.id)] = str(record.seq)
    return lookup

    
def loadFileToDict(filename, id_col='Oligo ID'):
    f = io.open(filename)
    rdr = csv.DictReader(f, delimiter='\t')
    if id_col is None: id_col = rdr.fieldnames[0]
    lookup = {row[id_col]: row for row in rdr}
    f.close()
    return lookup
  
def startup():
    if not RUN_LOCAL:
        start_idx = eval(sys.argv[1]) if len(sys.argv) > 1 else -1
        stop_idx = eval(sys.argv[2]) if len(sys.argv) > 2 else -1
        queue = sys.argv[3] if len(sys.argv) > 3 else 'normal'
        grouped = eval(sys.argv[4]) if len(sys.argv) > 4 else 0 
        if len(sys.argv) == 1:
            print('Usage: <python_script.py> start_idx stop_idx (opt)queue (opt)grouped(1 or 0) (opt)other_args\n')
    else:
        start_idx, stop_idx, queue, grouped = 0, 100000, 'normal', 0
    return start_idx, stop_idx, queue, grouped
  
def runSubdir(idx, subdirs, label,  python_script, out_label, caller, extra_args=''):
    #Runs a python script on a list of subdirs, either one job per dirname (GROUPED), or
    #one job per subdir (INDIVIDUAL)
    INDIVIDUAL, GROUPED, RECURSING = 0,1,2
    start_idx, stop_idx, queue, grouped = startup()
    out_dir = getLogDir()
 
    if grouped != RECURSING: print(idx, label)
    
    first_idx = idx
    for subdir in subdirs:
        if grouped == RECURSING:
            if idx >= start_idx and idx <= stop_idx:
                cmd = PYTHON_CMD + ' %s %s%s' % (python_script, extra_args, subdir) 
                print(cmd)
                os.system(cmd)
            idx += 1
        elif grouped == INDIVIDUAL:
            cmd = PYTHON_CMD + ' %s %s%s' % (python_script, extra_args, subdir) 
            idx = runCmdCheckIdx(cmd, idx, start_idx, stop_idx, out_dir, out_label,queue=queue)
        else: idx += 1
    if grouped == GROUPED:
        if first_idx >= start_idx and idx <= stop_idx+1:
            cmd = PYTHON_CMD + ' %s %d %d %s %d' % (caller,first_idx,idx-1,queue,RECURSING) 
            runCmdCheckIdx(cmd, 0, 0, 0, out_dir, out_label,queue=queue)
    return idx
  
def runPerSubdir(python_script, out_label, caller, extra_args='', include_null=False):
    idx = 0
    for dirname in [x for x in getAllDataDirs() if include_null or not isNullDir(x)]:
        if not os.path.isdir(dirname + '/mapped_reads'):
            print(getShortDir(dirname)), 'No mapped_reads directory'
        else:
            subdirs = getSubdirs(dirname)
            idx = runSubdir(idx, subdirs, getShortDir(dirname), python_script, out_label, caller, extra_args=extra_args)

def saveToPickle(object, filename):
    f = io.open(filename, 'wb')
    import pickle; pickle.dump(object, f)
    f.close()

def loadFromPickle(filename):
    f = io.open(filename, 'rb')
    import pickle; object = pickle.load(f); f.close()
    return object

def defaultLoadData(results_file, guideset=set(), oligo_id_str='Oligo Id'):

    data = pd.read_csv(results_file, sep='\t')
    if len(guideset) > 0: data = data.loc[data[oligo_id_str].isin(guideset)] 
    return data

def getCommonGuideset(results_list, part_guideset, data_function=defaultLoadData, id_colname='Oligo Id', reads_colname='Total reads', min_reads=0):
    common_guides = set([x for x in part_guideset])
    results_to_skip = []
    for result in results_list:
        if len(result) == 1:    data = data_function(result[0])
        else:   data = data_function(result)
        data = data.loc[data[reads_colname] > min_reads]
       
        guides = set([x for x in data[id_colname]])
        if len(guides) < 0.2*len(common_guides) and result != results_list[0]:
            results_to_skip.append(result)
            print('Skipping %s due to lack of oligo measurements' % result)
            continue
        common_guides = common_guides.intersection(guides)
    return common_guides, results_to_skip

# spec = { py_func_load: Python function to load data from each result, 
#                py_funcs_per_result: List of python functions to call on each result, 
#                py_funcs_all_results: List of python functions to call on all results combined, 
#                results_dir: Directory in which to find results, 
#                dirname_to_result_fn: function to go from dirname to result (file or directory as needed by load function),
#                result_to_dirname_fn: function invert above,
#                reads_colname: column name in data resulting from py_func_load on which to apply min_read constraint,
#                min_reads: minimum number of reads for each oligo in all samples,
#                check_output_fn: function to check if per result output is ok to proceed to combined }
def analyseResultsPerPartition( spec ):
    partitions = partitionGuides(oligo_detail_dir=getHighDataDir()+ '/ST_June_2017/data')
    samples_selectors = getSampleSelectors(include_wt=('include_wt' in spec and spec['include_wt']))
    
    #Backwards compatibility for single results_spec: move to listed spec
    if 'results_specs' not in spec: 
        spec['results_specs'] = [{x:spec[x] for x in ['results_dir','result_to_dirname_fn','dirname_to_result_fn']}]
    
    #Find common samples
    all_dir_sets = [[sp['result_to_dirname_fn'](x) for x in os.listdir(sp['results_dir'])] for sp in spec['results_specs']]
    all_dirs = set(all_dir_sets[0])
    for dirs in all_dir_sets[1:]:
        all_dirs = all_dirs.intersection(set(dirs))
    all_dirs = sortSampleNames([x for x in all_dirs])
    
    #Process various partitions of guides and samples according to provided functions
    for part_desc in spec['partitions']:
        print(part_desc)
        for sel_desc in spec['samples']:
            print(sel_desc)
            selector = samples_selectors[sel_desc]

            #Collect results to be processed
            all_result_outputs = []
            get_result = lambda x, sp: sp['results_dir'] + '/' + sp['dirname_to_result_fn'](x)

            results_list = [[get_result(x, sp) for sp in spec['results_specs']] for x in all_dirs if selector(x)]
            
            #Find which guides are common to all samples
            if 'use_common_only' not in spec or spec['use_common_only']:
                common_guides, results_to_skip = getCommonGuideset(results_list, partitions[part_desc], data_function=spec['py_func_load'],reads_colname=spec['reads_colname'],min_reads=spec['min_reads'],id_colname=spec['id_colname'])
            else: common_guides, results_to_skip = set(), []

            #Run all analysis functions on common guides
            for result in [x for x in results_list if x not in results_to_skip]:
                dirname = spec['results_specs'][0]['result_to_dirname_fn'](result[0])
                
                #Load Data
                if len(spec['results_specs']) == 1:
                    #If only one data source, flatten (for backwards compatibility)
                    data = spec['py_func_load'](result[0], common_guides)
                else:
                    data = spec['py_func_load'](result, common_guides)
                print(dirname)
                
                #Run analysis per results
                per_result_outputs = {x[1]:x[0](data, label=dirname) for x in spec['py_funcs_per_result']}
                if spec['check_output_fn'](per_result_outputs):
                    all_result_outputs.append((per_result_outputs, dirname))
            
            #Run analysis across all results        
            plot_label = part_desc + ' ' + sel_desc + ' (%d Oligos)' % len(common_guides)
            for cfunc in spec['py_funcs_all_results']:
                print(cfunc.__name__)
                cfunc(all_result_outputs, label=plot_label)
              

