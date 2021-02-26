#
# Copyright 2019-2020 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
############################## Help Section ##############################
.PHONY: help

help::
	$(ECHO) "Makefile Usage:"
	$(ECHO) "  make all TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to generate the design for specified Target and Shell."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make clean "
	$(ECHO) "      Command to remove the generated non-hardware files."
	$(ECHO) ""
	$(ECHO) "  make cleanall"
	$(ECHO) "      Command to remove all the generated files."
	$(ECHO) ""
	$(ECHO) "  make sd_card TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to prepare sd_card files."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make run TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to run application in emulation."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make build TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to build xclbin application."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  make host DEVICE=<FPGA platform> HOST_ARCH=<aarch32/aarch64/x86>"
	$(ECHO) "      Command to build host application."
	$(ECHO) "      By default, HOST_ARCH=x86. HOST_ARCH is required for SoC shells"
	$(ECHO) ""
	$(ECHO) "  NOTE: For SoC shells, ENV variable SYSROOT needs to be set."
	$(ECHO) ""

############################## Setting up Project Variables ##############################
MK_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
XF_PROJ_ROOT ?= $(shell pwd)
CUR_DIR := $(patsubst %/,%,$(dir $(MK_PATH)))
XFLIB_DIR = $(XF_PROJ_ROOT)

TARGET ?= sw_emu
HOST_ARCH := x86
SYSROOT := ${SYSROOT}
DEVICE ?= xilinx_u280_xdma_201920_3

ifneq ($(findstring u280, $(DEVICE)), u280)
$(error [ERROR]: This project is not supported for $(DEVICE).)
endif
include ./params.mk
include ./utils.mk

BIN_NAME:=${RTM_kernel}_maxY${RTM_maxY}_maxZ${RTM_maxZ}_b${RTM_MaxB}_nPEZ${RTM_nPEZ}_nPEX${RTM_nPEX}_nFSM${RTM_numFSMs}
BUILD_DIR := $(XF_PROJ_ROOT)/build
XDEVICE := $(call device2xsa, $(DEVICE))
TEMP_DIR := $(BUILD_DIR)/tmp/_x_$(TARGET)
XL_BUILD_DIR := $(BUILD_DIR)/$(TARGET)/$(BIN_NAME)
TEMP_REPORT_DIR := $(XL_BUILD_DIR)/reports/
BUILD_REPORT_DIR := $(TEMP_REPORT_DIR)/
EMCONFIG_DIR := $(XL_BUILD_DIR)

# Setting tools
VPP := v++

############################## Setting up Host Variables ##############################
#Include Required Host Source Files
HOST_SRCS += $(CUR_DIR)/main.cpp
HOST_SRCS += $(XFLIB_DIR)/include/xcl2/xcl2.cpp

CXXFLAGS += -I$(XFLIB_DIR)/include/
CXXFLAGS += -I$(XFLIB_DIR)/include/
CXXFLAGS += -I$(XFLIB_DIR)/include/
CXXFLAGS += -I$(XFLIB_DIR)/include/xcl2
CXXFLAGS += -I$(XFLIB_DIR)/../xf_blas/L1/tests/sw/include/



# Host compiler global settings
CXXFLAGS += -I$(XILINX_XRT)/include -I$(XILINX_VIVADO)/include -std=c++14 -O3 -Wall -Wno-unknown-pragmas -Wno-unused-label
LDFLAGS += -L$(XILINX_XRT)/lib -lOpenCL -lpthread -lrt -Wno-unused-label -Wno-narrowing -DVERBOSE
CXXFLAGS += -fmessage-length=0 -g
CXXFLAGS +=-I$(CUR_DIR)/src/ 


EXE_NAME := host.exe
EXE_FILE := $(XL_BUILD_DIR)/$(EXE_NAME)
HOST_ARGS := 

ifneq ($(HOST_ARCH), x86)
	LDFLAGS += --sysroot=$(SYSROOT)
endif

############################## Setting up Kernel Variables ##############################
# Kernel compiler global settings
VPP_FLAGS += -t $(TARGET) --platform $(XPLATFORM) --save-temps
LDCLFLAGS += --optimize 2 --jobs 8
VPP_FLAGS += -I$(XFLIB_DIR)/include/
VPP_FLAGS += -I$(XFLIB_DIR)/include/
VPP_FLAGS += -I$(XFLIB_DIR)/../xf_blas/include/

rtmforward_VPP_FLAGS += --kernel_frequency 250

# Kernel linker flags
LDCLFLAGS_rtmforward += --kernel_frequency 250 --config $(XF_PROJ_ROOT)/conn_u280.ini --config $(XF_PROJ_ROOT)/opts.ini
############################## Declaring Binary Containers ##############################
BINARY_CONTAINERS += $(XL_BUILD_DIR)/$(BIN_NAME).xclbin
BINARY_CONTAINER_rtmforward_OBJS += $(TEMP_DIR)/$(BIN_NAME).xo

############################## Setting Targets ##############################
CP = cp -rf

.PHONY: all clean cleanall docs emconfig
all: check_vpp check_platform | $(EXE_FILE) $(BINARY_CONTAINERS) emconfig

.PHONY: host
host: $(EXE_FILE) | check_xrt

.PHONY: xclbin
xclbin: check_vpp | $(BINARY_CONTAINERS)

.PHONY: build
build: xclbin

############################## Setting Rules for Binary Containers (Building Kernels) ##############################
$(TEMP_DIR)/$(BIN_NAME).xo: $(XFLIB_DIR)/src/rtmforward_domain_HBC.cpp
	$(ECHO) "Compiling Kernel: rtmforward"
	mkdir -p $(TEMP_DIR)
	cd $(BUILD_DIR); $(VPP) $(rtmforward_VPP_FLAGS) $(VPP_FLAGS) --temp_dir $(TEMP_DIR) --report_dir $(TEMP_REPORT_DIR) -c -k ${RTM_kernel} -I'$(<D)' -o'$@' '$<'

$(XL_BUILD_DIR)/$(BIN_NAME).xclbin: $(BINARY_CONTAINER_rtmforward_OBJS)
	mkdir -p $(XL_BUILD_DIR)
	cd $(BUILD_DIR); $(VPP) $(VPP_FLAGS) --temp_dir $(XL_BUILD_DIR) --report_dir $(BUILD_REPORT_DIR) -l $(LDCLFLAGS) $(LDCLFLAGS_rtmforward) -o'$@' $(+)

############################## Setting Rules for Host (Building Host Executable) ##############################
$(EXE_FILE): $(HOST_SRCS) | check_xrt
	mkdir -p $(XL_BUILD_DIR)
	cd $(BUILD_DIR); $(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(XPLATFORM) --od $(EMCONFIG_DIR)

############################## Setting Essential Checks and Running Rules ##############################
run: all
	XCL_EMULATION_MODE=$(TARGET) make data_gen
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	$(CP) $(EMCONFIG_DIR)/emconfig.json .
	XCL_EMULATION_MODE=$(TARGET) $(EXE_FILE) $(HOST_ARGS)
else
	$(EXE_FILE) $(HOST_ARGS)
endif

############################## Cleaning Rules ##############################
cleanh:
	-$(RMDIR) $(EXE_FILE) vitis_* TempConfig system_estimate.xtxt *.rpt .run/
	-$(RMDIR) src/*.ll _xocc_* .Xil dltmp* xmltmp* *.log *.jou *.wcfg *.wdb

cleank:
	-$(RMDIR) $(XL_BUILD_DIR)/*.xclbin _vimage *xclbin.run_summary qemu-memory-_* emulation/ _vimage/ pl* start_simulation.sh *.xclbin
	-$(RMDIR) _x_temp.*/_x.* _x_temp.*/.Xil _x_temp.*/profile_summary.* 
	-$(RMDIR) _x_temp.*/dltmp* _x_temp.*/kernel_info.dat _x_temp.*/*.log 
	-$(RMDIR) _x_temp.* 

cleanall: cleanh cleank
	-$(RMDIR) $(XL_BUILD_DIR)  build_dir.* emconfig.json *.html $(TEMP_DIR) $(BUILD_DIR)/reports *.csv *.run_summary $(CUR_DIR)/*.raw
	-$(RMDIR) $(XFLIB_DIR)/common/data/*.xe2xd* $(XFLIB_DIR)/common/data/*.orig*

sw_emu:
	+make TARGET=sw_emu build
hw_emu:
	+make TARGET=hw_emu build
hw:
	+make TARGET=hw build

clean: cleanh
