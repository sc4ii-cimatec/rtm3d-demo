#!/bin/bash
echo "index,timestamp,power.draw"
index=0
while true
do
    tstamp=`date +"%Y/%m/%d %H:%M:%S.%3N"`
    pwr=`$XILINX_XRT/bin/xbutil dump | grep \"power\": | cut -d: -f2 | tr -d \" | tr -d \ `
    echo "$index, $tstamp, $pwr"
    sleep 0.5
done