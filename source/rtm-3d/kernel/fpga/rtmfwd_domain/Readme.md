# RTM3D Kernels Compilation

FPGA kernels must be compiled before running. Three 
compilation alternatives are available: sw_emu, hw_emu, 
and hw. Output .xclbin files and compilation
reports are generated into "build/out_SDA_FLOW/". 

## Software emulation
By default, a small dataset:
```
	make sw_emu
```
For other dataset (as large as the Pluto model):
```
	make sw_emu RTM_x=400 RTM_y=212 RTM_z=256 RTM_time=1280
```
## Hardware emulation
By default, a small dataset:
```
	make hw_emu
```
For other dataset (as large as the Pluto model):
```
	make hw_emu RTM_x=400 RTM_y=212 RTM_z=256 RTM_time=1280
```
## Hardware
Build bitstream
```
	make hw
```
Run on hardware (Alveo U280)
``` 
or
	make run RTM_x=400 RTM_y=212 RTM_z=256 RTM_time=1280 RTM_verify=0
```
