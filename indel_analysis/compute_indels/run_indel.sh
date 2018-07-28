#!/usr/bin/env sh
set -e

dirname=$1
nulldir=$2
subdir=$3
fpp=$4

python indelmap_subdir.py $dirname $nulldir $subdir $LSB_JOBINDEX $fpp $5 $6 $7 $8 $9 $10 $11 $12

