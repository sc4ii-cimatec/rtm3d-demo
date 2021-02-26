#ifndef _RTM_GPU_
#define _RTM_GPU_
#include <stdio.h>
#include <stdlib.h>
#include<unistd.h>
#include <execinfo.h>
#include <RTMBase.hpp>
#include <RTMDeviceGrid.hpp>

#ifdef RTM_ACC_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#endif

using namespace std;

#ifdef RTM_ACC_GPU
#define CUDACHECK(cmd) do {  \
  cudaError_t e = cmd;                              \
  if( e != cudaSuccess ) {                          \
    printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); \
    printf("!> CUDA ERROR:\n!> %s:%d '%s (%d)'\n",            \
        __FILE__,__LINE__,cudaGetErrorString(e),e);   \
    printf("!> Stack trace: \n"); \
    void *array[10]; \
    size_t size;      \
    size = backtrace(array, 10); \
    char ** symbols = backtrace_symbols (array, size); \
    for (int i=0; i<size; i++){ printf("!>\t+ %s \n", symbols[i]);}\
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); \
    exit(EXIT_FAILURE);                             \
  }                                                 \
} while(0)

#define CUDANGRIDS(dim, bsize) (dim%bsize!=0)?((dim/bsize)+1):(dim/bsize)

#endif


#define DEFAULT_BLOCK_SIZE 32
#define SMALL_BLOCK_SIZE 10

/**
 *     CUDA Device Properties struct
 *     struct cudaDeviceProp {
 *       char name[256];
 *       size_t totalGlobalMem;
 *       size_t sharedMemPerBlock;
 *       int regsPerBlock;
 *       int warpSize;
 *       size_t memPitch;
 *       int maxThreadsPerBlock;
 *       int maxThreadsDim[3];
 *       int maxGridSize[3];
 *       size_t totalConstMem;
 *       int major;
 *       int minor;
 *       int clockRate;
 *       size_t textureAlignment;
 *       int deviceOverlap;
 *       int multiProcessorCount;
 *       int kernelExecTimeoutEnabled;
 *       int integrated;
 *       int canMapHostMemory;
 *       int computeMode;
 *   }
 * 
 * */
template<typename HostPtr_type, typename DevPtr_type=HostPtr_type>
class GPUDeviceGrid : public DeviceGrid<HostPtr_type,DevPtr_type>
{
public:
    GPUDeviceGrid(/* args */):DeviceGrid<HostPtr_type,DevPtr_type>{}{}

    virtual void createDeviceBuffer(uint32_t bufferFlags=0);
    virtual void removeDeviceBuffer();
    virtual void moveFromDevice();
    virtual void moveToDevice();
    virtual void devMemSet(HostPtr_type val);
};

#endif