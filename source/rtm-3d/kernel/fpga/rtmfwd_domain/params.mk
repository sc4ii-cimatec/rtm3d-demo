PY_SCRIPT=${XFLIB_DIR}/L1/tests/sw/python/operation.py

RTM_kernel = rtmforward
RTM_dataType      = float

ifndef RTM_x
RTM_x  = 64
endif
ifndef RTM_y
RTM_y  = 74
endif
ifndef RTM_z
RTM_z  = 84
endif
ifndef RTM_maxZ
RTM_maxZ  = 280
endif
ifndef RTM_maxY
RTM_maxY  = 180
endif
ifndef RTM_MaxB
RTM_MaxB  = 16
endif
RTM_NXB = $(RTM_MaxB)
RTM_NYB = $(RTM_MaxB)
RTM_NZB = $(RTM_MaxB)

RTM_order  = 8

ifndef RTM_numFSMs
RTM_numFSMs = 2
endif
ifndef RTM_nPEZ
RTM_nPEZ = 3
endif
ifndef RTM_nPEX
RTM_nPEX = 4
endif
ifndef RTM_time
RTM_time = 4
endif

RTM_verify=1
RTM_device = 0

MACROS += -D RTM_dataType=$(RTM_dataType) \
		  -D RTM_numFSMs=$(RTM_numFSMs) \
		  -D RTM_maxY=$(RTM_maxY) \
		  -D RTM_maxZ=$(RTM_maxZ) \
		  -D RTM_order=$(RTM_order) \
		  -D RTM_MaxB=$(RTM_MaxB) \
		  -D RTM_nPEZ=$(RTM_nPEZ) \
		  -D RTM_nPEX=$(RTM_nPEX) \
		  -D RTM_NXB=$(RTM_NXB) \
		  -D RTM_NYB=$(RTM_NYB) \
		  -D RTM_NZB=$(RTM_NZB) \
		  -D RTM_parEntries=$(RTM_parEntries)

CXXFLAGS += ${MACROS}
rtmforward_VPP_FLAGS += ${MACROS} 

DATA_DIR = ./$(BUILD_DIR)/dataset_z${RTM_z}_y${RTM_y}_x${RTM_x}_t${RTM_time}/
HOST_ARGS = $(BINARY_CONTAINERS) $(RTM_z) $(RTM_y) $(RTM_x) $(RTM_time) ${DATA_DIR} ${RTM_verify} ${RTM_device}

data_gen:
	mkdir -p ${DATA_DIR} 
	python3 ${PY_SCRIPT} --func testForward --path ${DATA_DIR} --z ${RTM_z} --y ${RTM_y} --x ${RTM_x} --time ${RTM_time} --nxb ${RTM_NXB} --nyb ${RTM_NYB} --nzb ${RTM_NZB} --order ${RTM_order} --verify ${RTM_verify}

run_hw: data_gen
	$(EXE_FILE) $(HOST_ARGS)
