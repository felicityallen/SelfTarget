#!/usr/bin/env sh
set -e

pearfile=$1
oligofilewithpam=$2
outputfile=$3
maxcutdist=$4

/lustre/scratch117/cellgen/team227/fa9/indelmap/bin/indelmap $pearfile $oligofilewithpam $outputfile 0 $maxcutdist
