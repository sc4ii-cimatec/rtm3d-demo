##############################################
# FPGA Host Flags and Kernel compilation
##############################################
include ./opencl.mk

ifndef XILINX_VITIS
XILINX_VITIS:=undefined
endif

ifndef XILINX_XRT
XILINX_XRT:=undefined
endif
########################################################################
# Check Xilinx Environment Variables
check-xilinx: check-xilinx-xrt check-xilinx-vitis check-xilinx-vivado

# Checks for XILINX_XRT
check-xilinx-xrt:
ifeq ($(XILINX_XRT), undefined)
	$(error XILINX_XRT not defined. Please define XILINX_XRT variable poiting to \
	a valid XRT installation )
else
	$(ECHO) "> XILINX_XRT installation found at: $(XILINX_XRT)"
endif

#Checks for XILINX_VITIS
check-xilinx-vitis:
ifeq ($(XILINX_VITIS), undefined)
	$(error XILINX_VITIS not defined. Please define XILINX_VITIS variable poiting to \
	a valid VITIS installation)
else
	$(ECHO) "> XILINX_VITIS installation found at: $(XILINX_VITIS)"
endif

#Checks for XILINX_VIVADO
check-xilinx-vivado:
ifeq ($(XILINX_VIVADO), undefined)
	$(error XILINX_VIVADO not defined. Please define XILINX_VIVADO variable poiting to \
	a valid VIVADO installation)
else
	$(ECHO) "> XILINX_VIVADO installation found at: $(XILINX_VIVADO)"
endif
########################################################################

.PHONY: all
all:

.PHONY: fpga-all
fpga-all: ACC_FPGA_KBUILD

.PHONY: fpga-clean
fpga-clean:
	$(RM) $(ACC_OBJS) *.xo

# HOST FLAGS
GCC_VERSION=6.2.0
GCC_PATH=${XILINX_VIVADO}/tps/lnx64
XCC = $(GCC_PATH)/gcc-$(GCC_VERSION)/bin/g++
XOPENCL_LIB_PATH:=${XILINX_XRT}/lib
X_LDIR := -L"${XILINX_XRT}/lib" 
X_IDIR := -I"kernel/fpga/include/" -I$(XILINX_XRT)/include $(opencl_CXXFLAGS)
X_LFLAGS = $(opencl_LDFLAGS) -lOpenCL -lxilinxopencl -lstdc++ -lrt -pthread -Wl,--rpath=${XOPENCL_LIB_PATH}

ACC_FPGA_FLAGS:
	$(ECHO) "> Setting HOST XCL Environment..."
	$(eval HOST_LDIR += $(X_LDIR))
	$(eval HOST_IDIR += $(X_IDIR))
	$(eval HOST_LDFLAGS += $(X_LFLAGS))
#	$(eval HOST_CPP := $(XCC))

ACC_FPGA_KBUILD:
	$(ECHO) "> Building fpga kernel..."
