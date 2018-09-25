import io, os, csv, sys

from predictor.model import computePredictedProfile, readTheta, setFeaturesDir, setReadsDir
from predictor.features import calculateFeaturesForGenIndelFile, readFeaturesData
from predictor.predict import predictMutationsBulk, predictMutationsSingle

if __name__ == '__main__':
      
    if len(sys.argv) == 3: #Batch mode
    
        batch_file = sys.argv[1]
        output_prefix = sys.argv[2]
    
        if not os.path.isfile(batch_file):
            raise Exception('Count not find batch file ' + batch_file)
            
        predictMutationsBulk(batch_file, output_prefix)
    
    elif len(sys.argv) == 4:    #Single mode
        target_seq = sys.argv[1]
        if sum([x not in 'ATGC' for x in target_seq]) > 0:
            raise Exception('Invalid target sequence, expecting string containing only A,T,G,C:\n%s' % target_seq)
        try:
            pam_idx = eval(sys.argv[2]) 
        except:
            raise('Could not parse PAM index, expected an integer %s' % pam_idx)
        output_prefix = sys.argv[3]
            
        predictMutationsSingle(target_seq, pam_idx, output_prefix)
    
    else:
        err_str = 'FORECasT: Invalid inputs. Usage:\n\nSingle gRNA: python FORECasT.py <guide_sequence> <PAM index (0 based)> <output_prefix>'
        err_str += '\n\nBatch gRNA: python FORECasT.py <batch_filename> <output_prefix>\n'
        raise Exception(err_str)
    
