#!/bin/bash

bin=$1
json=$2

export OMP_NUM_THREADS=72
srun $bin $json