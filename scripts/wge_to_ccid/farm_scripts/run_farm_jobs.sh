#!/usr/bin/env bash
source activate forecast
export PYTHONPATH=/nfs/users/nfs_c/cellgeni-su/anton/SelfTarget/:$PYTHONPATH
export CRISPR_ANALYSER=/lustre/scratch117/cellgen/cellgeni/forecast/CRISPR-Analyser/bin/crispr_analyser
bsub -J "myarray[33142-58884]" -o "output-%J.txt" -e "error-%J.txt" -M 3000 -R "select[mem>=3000] rusage[mem=3000]" ./wge_to_ccid.py map -f /lustre/scratch117/cellgen/cellgeni/forecast/mouse/input.\$LSB_JOBINDEX
