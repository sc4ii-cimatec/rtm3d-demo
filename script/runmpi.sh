#!/bin/bash
export MALLOC_CHECK_=2
rootdir=`pwd`
advisor=advixe-cl

if [ ! -d "$rootdir/script" ]; then
    echo ">> Missing 'script/' dir! Must run perf test from root folder. "
    exit 1
fi
if [ "$#" -lt 2 ]; then
    echo "Usage: "
    echo "> ./runmpi.sh input.json build_flag"
    echo ""
    exit 1
fi

if [[ -z "${SPARTION}" ]]; then
  SPARTION="CPUlongB"
fi
if [[ -z "${JOBNAME}" ]]; then
  JOBNAME="3D-aRTM"
fi

# input file
inputJson=$1
nprocs=$2
build=$3
reportdir=$4
binfile=bin/RTM3D_MPI.bin

#build host
if [ "$build" = "y" -o "$build" = "Y" ]; then
    rm $binfile > /dev/null 2>&1
    if [ "$build" = "Y" ]; then
        # clean build
	    $rootdir/script/buildhost.sh "on" "off" "y"
    else
        # just build
        $rootdir/script/buildhost.sh "on" "off" "n"
    fi
	if [ ! -f "$binfile" ]; then
	    exit 1
	fi   
fi

#mname=`cat $inputJson | grep mname | cut -d: -f2 | cut -d\" -f2`
OUTPUTFILE="$reportdir/job-output-NP$nprocs.log"
# run
export OMP_NUM_THREADS=72
echo "> sbatch -p ${SPARTION} -J $JOBNAME --nodes $nprocs --ntasks-per-node=1\
 --export='ALL' -o $OUTPUTFILE $rootdir/script/runjob.sh $binfile $inputJson"
sbatch -p ${SPARTION} -J $JOBNAME --nodes $nprocs --ntasks-per-node=1 \
--export='ALL' -o $OUTPUTFILE $rootdir/script/runjob.sh $binfile $inputJson

