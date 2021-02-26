#!/bin/bash
#
# This scripts plot a raw 3D fp-32 dataset using nimage application.
# Usage:


if [ $# -lt 4 ]; then
    echo "Missing parameters!"
    echo "Usage:"
    echo "\t #> nplot.sh <filepath> <nx> <ny> <nz> [perc] [title]" 
    exit
fi

file=$1
nx=$2
ny=$3
nz=$4
title=$1
perc=96
if [ $# -gt 4 ]; then
  perc=$5
fi
if [ $# -gt 5 ]; then
 title=$6
fi


./data/cwp/bin/nimage n1=$nz n2=$ny perc=$perc title=$title < $file &

