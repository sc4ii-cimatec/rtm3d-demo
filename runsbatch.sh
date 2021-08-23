#!/bin/bash
partition=CPUlongB
account="cenpes-lde"
jobname=$1
json=$2
cnodes=1
tnodes=1
cpus=72
outfile="data/output/slog_$jobname.log"
bin="bin/RTM3D.bin"


pwd; hostname; date

echo "Running 3DRTM $SLURM_CPUS_ON_NODE CPU cores"


sbatch -p ${partition} -A $account -J $jobname --nodes $cnodes \
  --ntasks-per-node=$tnodes -c $cpus --export='ALL' \
  -o $outfile ./runsrun.sh $bin $json
