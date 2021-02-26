#+-------------------------------------------------------------------------------
# The following parameters are assigned with default values. These parameters can
# be overridden through the make command line
#+-------------------------------------------------------------------------------
PROFILE := no
DEBUG := no

ifndef CUDA_ROOT
	CUDA_ROOT:=/usr/local/cuda
endif


# Cleaning stuff
RM = @rm -f
RMDIR = @rm -rf

ECHO:= @echo

.PHONY: all
all: