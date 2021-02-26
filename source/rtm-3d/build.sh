#!/bin/sh

help(){
    echo "> Usage: build.sh <build parameters> "
    echo "> build.sh --target=<host/kernel> --mpi=<on/off> --acc=<fpga/gpu/off>" 
    echo ">          --define='-D<MACRO_NAME[=VALUE]'"
    exit 1
}
# TARGET=sw_emu 
# DEVICE=xilinx_u280-es1_xdma_201910_1 
# KNAME=vadd
# SWTYPE=cpu
if [ $# -eq 0 ]; then
    echo "> Missing parameters!"
    help
fi

target="host"
mpi="off"
acc="none"
buildcmd=""
buildmacros=""
clean="n"
while [ $# -gt 0 ]; do
  case "$1" in
    --clean=*)
      clean="${1#*=}"
      ;;
    --target=*)
      host="${1#*=}"
      ;;
    --mpi=*)
      mpi="${1#*=}"
      ;;
    --acc=*)
      acc="${1#*=}"
      ;;
    --define=*)
      buildmacros="$buildmacros${1#*=}"
      ;;
    *)
      help
  esac
  shift
done

if [ $clean = "y" ]; then
    make -f HOSTMakefile.mk clean
fi

if [ $target = "host" ]; then
    if [ $mpi = "on" ]; then
        mpi="ON"
    elif [ $mpi = "off" ]; then
        mpi="OFF"
    else   
        echo "> Invalid MPI: $mpi"
        help
        exit 1
    fi
    if [ $acc = "fpga" ]; then
        acc="FPGA"
    elif [ $acc = "gpu" ]; then
        acc="GPU"
    elif [ $acc = "off" ]; then
        acc="NONE"
    else
        echo "> Invalid ACC: $acc."
        help
        exit 1
    fi
    if [ "$buildmacros" = "" ]; then
        RTM3D_BUILD_MACROS=""
    else
        #RTM3D_BUILD_MACROS="DEFS='$buildmacros'"
        RTM3D_BUILD_MACROS="'$buildmacros'"
    fi
    export RTM3D_BUILD_MACROS="$buildmacros"
    buildcmd="make -f HOSTMakefile.mk MPI=$mpi ACC=$acc $macros all"
# elif [ $1 = "host-only" ]; then
#     SWTYPE=$1
#     make $SWTYPE
# elif [ $1 = "fpga-only" ]; then
#     SWTYPE=$1
#     TARGET=$2
#     KNAME=$3
#     make clean TARGET=$2 DEVICE=$DEVICE KNAME=$KNAME
#     make $SWTYPE TARGET=$2 DEVICE=$DEVICE KNAME=$KNAME
else
    echo "> Invalid target parameter: $target"
    exit 1
fi

nThreads=40
mkdir -p "build"
LOGFILE="build/build.log"
echo "> $buildcmd"
`$buildcmd > $LOGFILE 2>&1`
output=`cat $LOGFILE | grep "BUILD SUCCESSFUL"`
#echo $output
if [ "$output" ]; then
    echo "> Build Successful!"
    #rm $LOGFILE
else
    echo "> Build Failed! Check '$LOGFILE' for details."
fi


