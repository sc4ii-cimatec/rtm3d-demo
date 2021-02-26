##########################################################
# As GPU Kernel codes are linked into the HOST program,
# this file builds all .cu cuda kernels and generates
# their corresponding object files for later linkage.
# GPU Kernel Source files (.cu) must be stored at 
# 'kernel/gpu/' folder.
##########################################################
ACC_SRCS :=
ACC_OBJS :=
ifndef GPU_ARCH
	GPU_ARCH:=sm_70
endif
GPUCC:=nvcc
GPUCFLAGS:=  -O3 -w --ftz=false --prec-div=true --fmad=false -arch=$(GPU_ARCH)
GPUDFLAGS:= -DRTM_ACC_GPU
GPULDFLAGS := -lcudart -lcuda
GPULDIR := -L$(CUDA_ROOT)/lib64 -L/usr/lib/x86_64-linux-gnu/
GPUIDIR:=-I$(CUDA_ROOT)/include -I"./include" -I"./kernel/gpu/"

ifndef CUDA_ROOT
	CUDA_ROOT:=/usr/local/cuda
endif
.PHONY: all
all:

.PHONY: check-cuda
check-cuda: 
	$(ECHO) "> Checking CUDA environment variables... "

.PHONY: gpu-all
gpu-all: ACC_GPU_KBUILD

.PHONY: gpu-clean
gpu-clean:
	$(RM) $(ACC_OBJS) *.o

ACC_GPU_KBUILD:
	$(eval ACC_SRCS += ./kernel/gpu/rtmgpu_finitediff.cu)
	$(eval ACC_OBJS += rtmgpu_finitediff.o)
	$(ECHO) "> BUILDING GPU KERNEL: $(GPUCC) -c $(GPUIDIR) $(GPUDFLAGS) $(GPUCFLAGS) $(ACC_SRCS) $(GPULDIR) $(GPULDFLAGS)"
	$(GPUCC) -c $(GPUIDIR) $(GPUDFLAGS) $(GPUCFLAGS) $(ACC_SRCS) $(GPULDIR) $(GPULDFLAGS)

ACC_GPU_FLAGS:
	$(eval HOST_DFLAGS += "-DRTM_ACC_GPU")
	$(eval HOST_LDIR += -L"$(CUDA_ROOT)/lib64" $(GPULDIR))
	$(eval HOST_LDFLAGS += -lcudart -lcuda)
	$(eval HOST_IDIR += -I$(CUDA_ROOT)/include)
	$(eval HOST_BIN_NAME:=$(HOST_BIN_NAME)_GPU)

ACC_GPU_CLEAN: gpu-clean