#!/bin/bash

ROOTDIR=`pwd`
daturl="https://ucc44bb8f1bdd85e039bee9420cb.dl.dropboxusercontent.com/cd/0/get/BJuRvOuXE8wEiALjizy8JMr81XD_iSnS8qH2lvepOhMnFxCvikPvZubaRxZY263xKiSiMG03UoDJHG5g6JsIFRfpMDaDEVmun6nrdqzQ9OdpLP1Iq2XtgIAAD15WeHdFz8A/file?_download_id=52453586887120844318846471753034249484744742225488098457927462395&_notify_domain=www.dropbox.com&dl=1"
datfile="rtm3d-demo-data.tar.xz"
datdir="$ROOTDIR/data"
mkdir -p $datdir
cd $datdir
echo "> Initializing DEMO environment..."
echo "> Downloading demonstration data package..."
wget --no-check-certificate $daturl -O $datfile

echo "> Extracting demonstration data..."
tar -xvf $datfile 

echo "> Build CWP library for nplot command. Please, answer 'y' when prompted."
cd $datdir/cwp 
./buildcwp.sh -y

cd $ROOTDIR
echo "> All set for demonstration."