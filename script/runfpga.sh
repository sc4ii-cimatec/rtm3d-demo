#!/bin/bash
export MALLOC_CHECK_=2

PWRPID=0
rootdir=`pwd`
sdaflow="hw"
startPowerMonitor(){
    if [ "$#" -lt 1 ]; then
        echo ">> !! Missing power report file! Abort !! <<"
        exit 1
    fi

    reportfile=$1
    # queries only for gpu 0
    $rootdir/script/fpgapwr.sh > $reportfile &
    local pid=$!
    PWRPID=$pid
}

source $XILINX_XRT/setup.sh > /dev/null 2>&1
# sdaflow="sw_emu"


if [ ! -d "$rootdir/script" ]; then
    echo ">> Missing 'script/' dir! Must run perf test from root folder. "
    exit 1
fi

if [ "$#" -lt 2 ]; then
    echo "Usage: "
    echo "> ./runfpga.sh input.json build_flag [kernel.xclbin]"
    echo ""
    exit 1
fi

# input file
inputJson=$1
build=$2
reportdir=$3
#kernel="kernel/$sdaflow/rtmforward_maxY128_maxZ512_b16_nPEZ4_nPEX2_nFSM2.xclbin"
kernel="kernel/$sdaflow/rtmforward_maxY32_maxZ64_b16_nPEZ2_nPEX2_nFSM4.xclbin"
if [ "$#" -eq 4 ]; then
    kernel=$4
fi

binfile=bin/RTM3D_FPGA.bin
#build host
if [ "$build" = "y" -o "$build" = "Y" ]; then
    rm $binfile > /dev/null 2>&1
    if [ "$build" = "Y" ]; then
        # clean build
	    $rootdir/script/buildhost.sh "off" "fpga" "y"
    else
        # just build
        $rootdir/script/buildhost.sh "off" "fpga" "n"
    fi
	if [ ! -f "$binfile" ]; then
	    exit 1
	fi   
fi

# run
XPLATFORM="xilinx_u280_xdma_201920_3"
if [ "$sdaflow" != "hw" ]; then
    export XCL_EMULATION_MODE=$sdaflow
fi

powercsv=$reportdir/fpgapower.csv
startPowerMonitor $powercsv
echo "> ./$binfile $inputJson $kernel"
./$binfile $inputJson $kernel
kill -9 $PWRPID > /dev/null 2>&1
# reset fpga device
# echo "> Reset FPGA board:"
# xbutil reset