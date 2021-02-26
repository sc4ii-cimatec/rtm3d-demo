#!/bin/bash
export MALLOC_CHECK_=2

PWRPID=0
rootdir=`pwd`
startPowerMonitor(){
    if [ "$#" -lt 1 ]; then
        echo ">> !! Missing power report file! Abort !! <<"
        exit 1
    fi

    reportfile=$1
    # queries only for gpu 0
    nvidia-smi -i 0 --query-gpu=index,timestamp,power.draw,clocks.sm --format=csv -l 1 > $reportfile &
    local pid=$!
    PWRPID=$pid
}

if [ ! -d "$rootdir/script" ]; then
    echo ">> Missing 'script/' dir! Must run perf test from root folder. "
    exit 1
fi
if [ "$#" -lt 2 ]; then
    echo "Usage: "
    echo "> ./rungpu.sh input.json build_flag"
    echo ""
    exit 1
fi

# input file
inputJson=$1
build=$2
reportdir=$3
binfile=bin/RTM3D_GPU.bin
#build host
if [ "$build" = "y" -o "$build" = "Y" ]; then
    rm $binfile > /dev/null 2>&1
    if [ "$build" = "Y" ]; then
        # clean build
	    $rootdir/script/buildhost.sh "off" "gpu" "y"
    else
        # just build
        $rootdir/script/buildhost.sh "off" "gpu" "n"
    fi
	if [ ! -f "$binfile" ]; then
	    exit 1
	fi   
fi

powercsv=$reportdir/gpupower.csv
startPowerMonitor $powercsv
echo "> ./$binfile $inputJson"
./$binfile $inputJson
kill -9 $PWRPID > /dev/null 2>&1
