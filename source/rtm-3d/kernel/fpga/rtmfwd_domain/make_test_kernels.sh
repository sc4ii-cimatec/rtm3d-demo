#!/bin/bash
mlist=("3LAYERS_64x74x84" "MONSERRAT_122x101x194" "SEG-EAGE-214x214x210" "OVERTHRUST_402x242x432")
setModelParam(){
    mname=$1
    NX=64;NY=74;NZ=84;NT=512;BLEN=16;
    if [ "$mname" = "3LAYERS_64x74x84" ]; then
        NX=64;NY=74;NZ=84;NT=512;BLEN=12;
    elif [ "$mname" = "MONSERRAT_122x101x194" ]; then
        NX=122;NY=101;NZ=194;NT=2304;BLEN=15;
	elif [ "$mname" = "SEG-EAGE-214x214x210" ]; then
		NX=214;NY=214;NZ=210;NT=2368;BLEN=15;
	elif [ "$mname" = "OVERTHRUST_402x242x432" ]; then
		NX=402;NY=242;NZ=432;NT=6192;BLEN=15;
	fi
}

copy_bin(){
	BIN_NAME="rtmforward_maxY${maxY}_maxZ${maxZ}_b${blen}_nPEZ${nPEZ}_nPEX${nPEX}_nFSM${nFSM}"
	BUILD_DIR="./build/${sda_flow}/${BIN_NAME}"
	BIN_DIR="./xclbin/${sda_flow}/"
	mkdir -p ${BIN_DIR}
	cp "$BUILD_DIR/${BIN_NAME}.xclbin" "${BIN_DIR}"
}

kbuild(){
	echo ">> ******************************************************"
	echo ">>"
	echo ">> make ${sda_flow} RTM_numFSMs=$nFSM RTM_nPEZ=$nPEZ RTM_nPEX=$nPEX RTM_maxZ=$maxZ RTM_maxY=$maxY RTM_MaxB=$blen"
	echo ">>"
	make ${sda_flow} RTM_numFSMs=$nFSM RTM_nPEZ=$nPEZ RTM_nPEX=$nPEX RTM_maxZ=$maxZ RTM_maxY=$maxY RTM_MaxB=$blen
	echo ">> ******************************************************"
}
    
sda_flow="sw_emu"
if [ "$#" -eq 1 ]; then
sda_flow=$1
fi
fsmlist=("2" "4")
	
export XILINX_XRT="/opt/xilinx/xrt"
export XILINX_VITIS="/tools/Xilinx/Vitis/2020.1"
export XILINX_VIVADO="/tools/Xilinx/Vivado/2020.1"
export PATH="$XILINX_VIVADO/bin:XILINX_VITIS/bin:$PATH"
export PYTHONPATH="/opt/xilinx/xrt/python"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu:/opt/xilinx/xrt/lib:/usr/local/cuda-10.2/lib64"
export LIBRARY_PATH="/opt/xilinx/xrt/lib:/usr/lib/x86_64-linux-gnu:/usr/local/cuda-10.2/lib64$LD_LIBRARY_PATH"
source /opt/xilinx/xrt/setup.sh

rm build/*.vitis.log > /dev/null 2>&1
mkdir -p build

dimlistY=("32" "64" "128" "256")
dimlistZ=("512" "256")
for maxDimZ in "${dimlistZ[@]}"
do
	for maxDimY in "${dimlistY[@]}"
	do
		for nFSM in "${fsmlist[@]}"
		do
			# build 256x256x16
			maxY=$maxDimY; maxZ=$maxDimZ;
			nPEX=2;nPEZ=4;
			blen=16;
			echo ">> ******************************************************"
			echo ">> Build Kernel: ($maxDimY x $maxDimZ) nFSM=$nFSM "
			kbuild >> build/kernel_${maxDimY}x${maxDimZ}_${nFSM}.vitis.log 2>&1
			copy_bin
			echo ">> ******************************************************"
		done
	done
done
echo "> All Done"
