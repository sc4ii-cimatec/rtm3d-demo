##################################################
# Host Build Process
##################################################
include ./utils.mk
include ./opencl.mk
include ./GPUMakefile.mk
include ./FPGAMakefile.mk

help::
	$(ECHO) "> Makefile Usage:"
	$(ECHO) ">  make all MPI<ON/OFF> ACC=<FPGA/GPU/NONE> "
	$(ECHO) ">      Builds the HOST program."
	$(ECHO) "> "
	$(ECHO) "> "

# Build custom definitions
ifndef RTM3D_BUILD_MACROS
RTM3D_BUILD_MACROS:=
endif

# Host Env
BASE_DIR=$(shell pwd)
SRCEXT 		:= cpp
OBJEXT 		:= o
DEPEXT      := d
HOST_BUILD_DIR:=${BASE_DIR}/build
HOST_SRC_DIR=${BASE_DIR}/src
HOST_SRCS := $(shell find $(HOST_SRC_DIR) -type f -name *.$(SRCEXT))
HOST_OBJS  := $(patsubst $(HOST_SRC_DIR)/%,$(HOST_BUILD_DIR)/%,$(HOST_SRCS:.$(SRCEXT)=.$(OBJEXT)))
HOST_CPP := g++
HOST_IDIR = -I${BASE_DIR}/include -I${BASE_DIR}/kernel/gpu/ -I${BASE_DIR}/kernel/fpga/

# Host compiler global settings
HOST_CXXFLAGS += -fmessage-length=0 -ggdb -rdynamic -std=c++11 -O3 -fopenmp
HOST_LDFLAGS += -lrt -lstdc++ -lm -fopenmp
HOST_LDIR := 
HOST_DFLAGS += -DOMP $(RTM3D_BUILD_MACROS)
HOST_BIN_NAME :=

# MPI-related rules
MPI_RULE:=
ifeq ($(MPI),ON)
MPI_RULE=MPI
endif

# Accelerator-related rules
ACC_RULE:=ACC_NONE
ifeq ($(ACC),FPGA)
ACC_RULE:=ACC_FPGA
endif
ifeq ($(ACC),GPU)
ACC_RULE:=ACC_GPU
endif

# Build all
.PHONY: all
all: BUILD_DIR $(MPI_RULE) $(ACC_RULE) HOST_BNAME HOST_BUILD ACC_GPU_CLEAN

# Clean all
.PHONY: clean
clean:
	@$(ECHO) "> Cleaning old object and binary files... "
	@$(RM) *.bin
	@$(RM) -rf $(HOST_BUILD_DIR)/

#Copy Resources from Resources Directory to Target Directory
BUILD_DIR:
	@mkdir -p $(HOST_BUILD_DIR)

#Pull in dependency info for *existing* .o files
-include $(HOST_OBJS:.$(OBJEXT)=.$(DEPEXT))

# Compile
$(HOST_BUILD_DIR)/%.$(OBJEXT): $(HOST_SRC_DIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	@echo "Host DFLAGS: $(HOST_DFLAGS)"
	@$(HOST_CPP) $(HOST_CXXFLAGS) $(HOST_DFLAGS) $(HOST_IDIR) -c -o $@ $<
	@$(HOST_CPP) $(HOST_CXXFLAGS) $(HOST_DFLAGS) $(HOST_IDIR) -MM $(HOST_SRC_DIR)/$*.$(SRCEXT) > $(HOST_BUILD_DIR)/$*.$(DEPEXT)
	@cp -f $(HOST_BUILD_DIR)/$*.$(DEPEXT) $(HOST_BUILD_DIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(HOST_BUILD_DIR)/$*.$(OBJEXT):|' < $(HOST_BUILD_DIR)/$*.$(DEPEXT).tmp > $(HOST_BUILD_DIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(HOST_BUILD_DIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(HOST_BUILD_DIR)/$*.$(DEPEXT)
	@rm -f $(HOST_BUILD_DIR)/$*.$(DEPEXT).tmp

#Link
HOST_TARGET: ${HOST_OBJS}
	@echo $(HOST_CPP) -o ${HOST_BIN_NAME} $^ $(HOST_LDIR) $(HOST_LDFLAGS) $(ACC_OBJS)
	@$(HOST_CPP) -o ${HOST_BIN_NAME} $^ $(HOST_LDIR) $(HOST_LDFLAGS) $(ACC_OBJS)

# Build Host software
HOST_BUILD: HOST_TARGET
	@$(ECHO) "> HOST_BIN_NAME: ${HOST_BIN_NAME}"
	@$(ECHO) "> BUILD SUCCESSFUL"

# Sets host's binary file name
HOST_BNAME:
	$(eval HOST_BIN_NAME=RTM3D$(HOST_BIN_NAME).bin)

# Set MPI specific flags
MPI:
	$(eval HOST_DFLAGS += "-DRTM_MPI")
	$(eval HOST_CPP := mpic++)
	$(eval HOST_BIN_NAME :=$(HOST_BIN_NAME)_MPI)

# Set FPGA-related flags
ACC_FPGA: check-xilinx ACC_FPGA_FLAGS
	$(eval HOST_DFLAGS += "-DRTM_ACC_FPGA")
	$(eval HOST_BIN_NAME :=$(HOST_BIN_NAME)_FPGA)
	

# Set GPU-related flags and build .cu files
ACC_GPU: check-cuda ACC_GPU_FLAGS ACC_GPU_KBUILD

# No accelerator used
ACC_NONE:
	$(eval HOST_DFLAGS += "-DRTM_ACC_NONE")
