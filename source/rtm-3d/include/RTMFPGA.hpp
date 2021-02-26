/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef RTM_FPGA_HPP
#define RTM_FPGA_HPP
#include <RTMDeviceGrid.hpp>
#include <chrono>
#include <vector>

using namespace std;

#ifdef  RTM_ACC_FPGA
#include "xcl2/xcl2.hpp"
#include "xrtm/fpga.hpp"
#include "xrtm/types.hpp"
#include "xrtm/binFiles.hpp"
#include "xrtm/utils.hpp"
#define CLBuffer_t cl::Buffer
#else
#define CLBuffer_t void *
#endif

#ifndef RTM_FPGA_nFSM
#define RTM_FPGA_nFSM 2
#endif
#ifndef RTM_FPGA_nPEX
#define RTM_FPGA_nPEX 2
#endif
#ifndef RTM_FPGA_nPEZ
#define RTM_FPGA_nPEZ 4
#endif
#ifndef RTM_FPGA_stOrder
#define RTM_FPGA_stOrder 8
#endif

#define RTM_FPGA_V2DT2_BASE 24
#define RTM_FPGA_P0_BASE    0
#define RTM_FPGA_P1_BASE    6
#define RTM_FPGA_PP0_BASE   12
#define RTM_FPGA_PP1_BASE   18

#ifdef  RTM_ACC_FPGA
#define RTMFPGADevice FPGA
#else
#define RTMFPGADevice void
#endif

template<typename HostPtr_type, typename DevPtr_type=CLBuffer_t>
class FPGADeviceGrid : public DeviceGrid<HostPtr_type, DevPtr_type>
{
public:
    FPGADeviceGrid(/* args */):DeviceGrid<HostPtr_type,DevPtr_type>{}{}

    virtual void createDeviceBuffer(uint32_t bufferFlags=0);
    virtual void removeDeviceBuffer();
    virtual void moveFromDevice();
    virtual void moveToDevice();
    virtual void devMemSet(HostPtr_type val);
};


#endif