#!/bin/bash
USER=xilinx.fpga 
PASSKEY="W39@Ld&"
ROOTDIR=`pwd`
DATURL="ftp://200.9.65.22/rtm3d-demo-data.tar.xz"
DATFILE="rtm3d-demo-data.tar.xz"
DATDIR="$ROOTDIR/data"
mkdir -p $DATDIR
cd $DATDIR
echo "> Initializing DEMO environment..."
echo "> Downloading demonstration data package..."
echo "wget $DATURL --user=$USER --password=$PASSKEY"
wget $DATURL --user=$USER --password=$PASSKEY

echo "> Extracting demonstration data..."
tar -xvf $DATFILE 

echo "> Build CWP library for nplot command. Please, answer 'y' when prompted."
cd $DATDIR/cwp 
./buildcwp.sh -y

cd $ROOTDIR
echo "> All set for demonstration."