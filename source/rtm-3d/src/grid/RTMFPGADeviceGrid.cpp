#include <iostream>
#include <vector>
#include <Misc.hpp>
#include <RTM.hpp>
#include <RTMFPGA.hpp>

#ifdef RTM_ACC_FPGA
#include "xcl2/xcl2.hpp"
#endif

template<>
void FPGADeviceGrid<RTMData_t, RTMDevPtr_t>::createDeviceBuffer(uint32_t bufferFlags)
{
#ifdef RTM_ACC_FPGA
#endif
}
template<>
void FPGADeviceGrid<RTMData_t, RTMDevPtr_t>::removeDeviceBuffer()
{
#ifdef RTM_ACC_FPGA

#endif
}
template<>
void FPGADeviceGrid<RTMData_t, RTMDevPtr_t>::moveFromDevice()
{
#ifdef RTM_ACC_FPGA

#endif
}
template<>
void FPGADeviceGrid<RTMData_t, RTMDevPtr_t>::moveToDevice()
{
#ifdef RTM_ACC_FPGA

#endif
}

template<>
void FPGADeviceGrid<RTMData_t, RTMDevPtr_t>::devMemSet(RTMData_t val)
{
#ifdef RTM_ACC_FPGA

#endif
}
