#include <iostream>
#include <vector>
#include <Misc.hpp>
#include <RTM.hpp>
#include <RTMAcc.hpp>

#ifdef RTM_ACC_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#include <rtmgpu.hpp>
#endif

template<>
void GPUDeviceGrid<RTMData_t, RTMDevPtr_t>::createDeviceBuffer(uint32_t bufferFlags)
{
#ifdef RTM_ACC_GPU
    CUDACHECK(cudaMalloc(&DEV_PTR, BUFFER_LENGTH));
    RTMGPUPlatform::subMemory(BUFFER_LENGTH);
#endif
}
template<>
void GPUDeviceGrid<RTMData_t, RTMDevPtr_t>::removeDeviceBuffer()
{
#ifdef RTM_ACC_GPU
    if(DEV_PTR!=nullptr){
        //CUDACHECK(cudaFree(DEV_PTR));
        cudaFree(DEV_PTR);
        RTMGPUPlatform::addMemory(BUFFER_LENGTH);
        DEV_PTR=nullptr;
    }
#endif
}
template<>
void GPUDeviceGrid<RTMData_t, RTMDevPtr_t>::moveFromDevice()
{
#ifdef RTM_ACC_GPU
    CUDACHECK(cudaMemcpy(HOST_PTR, DEV_PTR, BUFFER_LENGTH, cudaMemcpyDeviceToHost));
#endif
}
template<>
void GPUDeviceGrid<RTMData_t, RTMDevPtr_t>::moveToDevice()
{
#ifdef RTM_ACC_GPU
    // printf(">> HOST_PTR=%p; DEV_PTR=%p; len=%ld MB\n", 
    // HOST_PTR, DEV_PTR, BUFFER_LENGTH/1000000);fflush(stdout);
    CUDACHECK(cudaMemcpy(DEV_PTR, HOST_PTR, BUFFER_LENGTH, cudaMemcpyHostToDevice));
#endif
}

template<>
void GPUDeviceGrid<RTMData_t, RTMDevPtr_t>::devMemSet(RTMData_t val){
#ifdef RTM_ACC_GPU
    int ival = (int) val;
    //printf(">> %s: HOST_PTR=%p; DEV_PTR=%p; len=%ld MB\n", __func__, HOST_PTR, DEV_PTR, length_in_bytes/1000000);
    if(DEV_PTR!=NULL)
        CUDACHECK(cudaMemset(DEV_PTR,ival,BUFFER_LENGTH));
#endif
}
