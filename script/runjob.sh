#!/bin/bash
export MALLOC_CHECK_=2
export OMP_NUM_THREADS=`getconf _NPROCESSORS_ONLN`
mapby="none"
binfile=$1
inputpath=$2
if [ "$#" -eq 3 ]; then
    mapby=$3
fi
mpirun --bind-to $mapby $binfile $inputpath
