#!/usr/bin/env sh
set -e

dirname=$1
nulldir=$2
subdir=$3
fpp=$4

python correct_indelmap_subdir.py $dirname $nulldir $subdir $LSB_JOBINDEX $fpp

