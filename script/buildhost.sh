#!/bin/bash
rootdir=`pwd`
srcroot="$rootdir/source/rtm-3d"
outdir="bin/"
acc="off"
buildclean="y"
mpibuild=$1
acc=$2

if [ "$#" -gt 2 ]; then
buildclean=$3 # clean build
fi

echo "> Building Host Program..."
#check source root
if [ ! -f "$srcroot/build.sh" ]; then
    echo "> Invalid SOURCE ROOT folder ($srcroot). "
    exit 1
fi
cd $srcroot
rm *.bin > /dev/null 2>&1
rm $rootdir/build.log > /dev/null 2>&1
./build.sh --clean=$buildclean --target=host --mpi=$mpibuild --acc=$acc --define="$buildmacros"
binfile=`ls *.bin`
if [ -f "$binfile" ]; then
    mkdir -p $rootdir/bin
    mv $binfile $rootdir/bin
else

    cp build/build.log $rootdir
fi
cd $rootdir 
exit 0