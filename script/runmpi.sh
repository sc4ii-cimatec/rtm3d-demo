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
if [[ -z "${ACCOUNT}" ]]; then
  ACCOUNT="cenpes-lde"
fi

# input file
inputJson=$1
nprocs=$2
ntask=$3
build=$4
reportdir=$5
acc=$6

binfile=bin/RTM3D_MPI.bin
if [ "$acc" == "gpu" ];then
    binfile=bin/RTM3D_MPI_GPU.bin
fi

#build host
if [ "$build" = "y" -o "$build" = "Y" ]; then
    rm $binfile > /dev/null 2>&1
    if [ "$build" = "Y" ]; then
        # clean build
	    $rootdir/script/buildhost.sh "on" "$acc" "y"
    else
        # just build
        $rootdir/script/buildhost.sh "on" "$acc" "n"
    fi
	if [ ! -f "$binfile" ]; then
	    exit 1
	fi   
fi

#mname=`cat $inputJson | grep mname | cut -d: -f2 | cut -d\" -f2`
OUTPUTFILE="$reportdir/job-output-NP$nprocs.log"
# run
export OMP_NUM_THREADS=72
echo "> sbatch -p ${SPARTION} -A $ACCOUNT -J $JOBNAME --nodes $nprocs --ntasks-per-node=$ntask\
 --export='ALL' -o $OUTPUTFILE $rootdir/script/runjob.sh $binfile $inputJson"
sbatch -p ${SPARTION} -J $JOBNAME -A $ACCOUNT --nodes $nprocs --ntasks-per-node=$ntask \
--export='ALL' -o $OUTPUTFILE $rootdir/script/runjob.sh $binfile $inputJson

