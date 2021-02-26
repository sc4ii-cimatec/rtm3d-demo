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
#ifndef XF_RTM_FPGA_HPP
#define XF_RTM_FPGA_HPP

#include "utils.hpp"
#include "binFiles.hpp"
#include "types.hpp"
#include "xcl2/xcl2.hpp"
#include <chrono>
#include <vector>

using namespace std;

/**
 * @file fpga.hpp
 * @brief class FPGA, ForwardKernel and BackwardKernel are defined here.
 */



/**
 * @brief class FPGA is used to manage FPGA info
 */
class FPGA {
   public:
    FPGA(string p_xclbinFile, size_t p_deviceId = 0) {
        m_deviceId = p_deviceId;
        m_devices = xcl::get_xil_devices();

        m_Context = cl::Context({m_devices[m_deviceId]});
        m_CommandQueue = cl::CommandQueue(m_Context, m_devices[m_deviceId], CL_QUEUE_PROFILING_ENABLE);
        m_bins = xcl::import_binary_file(p_xclbinFile);
        m_Program = cl::Program(m_Context, {m_devices[m_deviceId]}, m_bins);
    }

    cl::Kernel createKernel(const string& p_name) const { return cl::Kernel(m_Program, p_name.c_str()); }
    const cl::Context& getContext() const { return m_Context; }
    const cl::CommandQueue& getCommandQueue() const { return m_CommandQueue; }

    const cl::Device operator[](const size_t p_Idx) const { return m_devices[p_Idx]; }

    template <typename T>
    vector<cl::Buffer> createDeviceBuffer(cl_mem_flags p_flags,
                                          const vector<HostBuffer_t<T> >& p_buffer,
                                          size_t p_base_channel) const {
        int p_hbm_pc = p_buffer.size();
        vector<cl::Buffer> l_buffer(p_hbm_pc);
        for (int i = 0; i < p_hbm_pc; i++) {
            cl_int err;
            size_t l_bufferSize = sizeof(T) * p_buffer[i].size();
            cl_mem_ext_ptr_t ext;
            ext.flags = XCL_MEM_TOPOLOGY | (i + p_base_channel);
            ext.obj = (void*)p_buffer[i].data();
            ext.param = nullptr;
            l_buffer[i] =
                cl::Buffer(m_Context, p_flags | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR, l_bufferSize, &ext, &err);
            if (err != CL_SUCCESS) {
                printf("Failed to allocate device buffer!\n");
                throw std::bad_alloc();
            }
        }
        return l_buffer;
    }

    template <typename T>
    cl::Buffer createDeviceBuffer(cl_mem_flags p_flags, const HostBuffer_t<T>& p_buffer) const {
        size_t l_bufferSize = sizeof(T) * p_buffer.size();
        cl_int err;
        cl::Buffer buffer =
            cl::Buffer(m_Context, p_flags | CL_MEM_USE_HOST_PTR, l_bufferSize, (void*)p_buffer.data(), &err);
        if (err != CL_SUCCESS) {
            printf("Failed to allocate device buffer!\n");
            throw std::bad_alloc();
        }
        return buffer;
    }

   private:
    size_t m_deviceId;
    vector<cl::Device> m_devices;
    cl::Context m_Context;
    cl::CommandQueue m_CommandQueue;
    cl::Program::Binaries m_bins;
    cl::Program m_Program;
};

/**
 * @brief class ForwardKernel is used to manage and run forward kernel on FPGA
 */

template <typename t_DataType, size_t t_Order, size_t t_PEZ, size_t t_PEX>
class ForwardKernel {
    typedef WideData<t_DataType, t_PEX> t_WideTypeX;
    typedef WideData<t_WideTypeX, t_PEZ> t_WideType;

    static const size_t t_SizePerChannel = 256 * 1024 * 1024;
    static const size_t t_ChannelSize = t_SizePerChannel / sizeof(t_WideType);

   public:
    ForwardKernel(const FPGA* fpga,
                  size_t p_z,
                  size_t p_y,
                  size_t p_x,
                  size_t p_zb,
                  size_t p_yb,
                  size_t p_xb,
                  size_t p_time,
                  size_t p_v2dt2_base = 0,
                  size_t p_p0_base = 0,
                  size_t p_p1_base = 0,
                  size_t p_pp0_base = 0,
                  size_t p_pp1_base = 0) {
        static const string t_KernelName = "rtmforward";
        m_fpga = fpga;
        m_Kernel = m_fpga->createKernel(t_KernelName);

        m_z = p_z;
        m_y = p_y;
        m_x = p_x;
        m_cube = m_x * m_y * m_z;
        m_zb = p_zb;
        m_yb = p_yb;
        m_xb = p_xb;
        m_time = p_time;

        m_src.resize(p_time);// = new HostBuffer_t<t_DataType>(p_time);
        m_v2dt2.resize(m_cube / t_PEX / t_PEZ);

        m_v2dt2_base = p_v2dt2_base;
        m_pp0_base = p_pp0_base;
        m_pp1_base = p_pp1_base;
        m_p0_base = p_p0_base;
        m_p1_base = p_p1_base;

        l_cubeSize = m_cube / t_PEX / t_PEZ;
        p_numChannels = (t_ChannelSize - 1 + l_cubeSize) / t_ChannelSize;
    }

    void createDeviceBuffers(size_t upbSize=1)
    {
        l_pp0.resize(p_numChannels);
        l_pp1.resize(p_numChannels);
        l_p0.resize(p_numChannels);
        l_p1.resize(p_numChannels);
        l_v2dt2.resize(p_numChannels);
        m_upb = new HostBuffer_t<t_DataType>(upbSize);

        for (int i = 0; i < p_numChannels; i++) {
            size_t l_vectorSize = l_cubeSize >= t_ChannelSize ? t_ChannelSize : l_cubeSize;
            l_cubeSize -= t_ChannelSize;

            l_p0[i].resize(l_vectorSize);
            l_p1[i].resize(l_vectorSize);
            l_pp0[i].resize(l_vectorSize);
            l_pp1[i].resize(l_vectorSize);
            l_v2dt2[i].resize(l_vectorSize);

            copy(m_v2dt2.begin() + i * t_ChannelSize, m_v2dt2.begin() + i * t_ChannelSize + l_vectorSize,
                 l_v2dt2[i].begin());
        }

        d_coefx = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_coefx);
        d_coefy = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_coefy);
        d_coefz = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_coefz);
        d_srcs = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_src);
        d_taperx = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_taperx);
        d_tapery = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_tapery);
        d_taperz = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_ONLY, m_taperz);
        d_upb = m_fpga->createDeviceBuffer<t_DataType>(CL_MEM_READ_WRITE, *m_upb);
        d_v2dt2 = m_fpga->createDeviceBuffer<t_WideType>(CL_MEM_READ_ONLY, l_v2dt2, m_v2dt2_base);
        d_p0 = m_fpga->createDeviceBuffer<t_WideType>(CL_MEM_READ_WRITE, l_p0, m_p0_base);
        d_p1 = m_fpga->createDeviceBuffer<t_WideType>(CL_MEM_READ_WRITE, l_p1, m_p1_base);
        d_pp0 = m_fpga->createDeviceBuffer<t_WideType>(CL_MEM_READ_WRITE, l_pp0, m_pp0_base);
        d_pp1 = m_fpga->createDeviceBuffer<t_WideType>(CL_MEM_READ_WRITE, l_pp1, m_pp1_base);

        inBufVec.push_back(d_coefx);
        inBufVec.push_back(d_coefy);
        inBufVec.push_back(d_coefz);
        inBufVec.push_back(d_taperx);
        inBufVec.push_back(d_tapery);
        inBufVec.push_back(d_taperz);
        inBufVec.push_back(d_srcs);
        for (int i = 0; i < p_numChannels; i++) {
            inBufVec.push_back(d_p0[i]);
            inBufVec.push_back(d_pp0[i]);
            inBufVec.push_back(d_v2dt2[i]);
            /**** ALU ***/
            inBufVec.push_back(d_p1[i]);
            inBufVec.push_back(d_pp1[i]);

            outBufVec.push_back(d_p0[i]);
            outBufVec.push_back(d_pp0[i]);
            outBufVec.push_back(d_p1[i]);
            outBufVec.push_back(d_pp1[i]);
            /**** ALU ***/
        }
        outBufVec.push_back(d_upb);

    }
    /**
     * @brief run kernel, memory selection, HBC mode
     *
     * @param p_sel, memory port selection
     *
     * @param p_shot, shot id
     * @param p_sy, shot y coordinate
     * @param p_sx, shot x coordinate
     *
     * @param p_p, seismic snapshot
     * @param p_upb, upper boundary data
     */
    double run(const bool p_sel,
               size_t p_shot,
               size_t p_sy,
               size_t p_sx,
               bool saveSnaps,
               HostBuffer_t<t_WideType>& p_pp,
               HostBuffer_t<t_WideType>& p_p) 
    {
        cl::Context m_Context = m_fpga->getContext();
        cl::CommandQueue m_CommandQueue = m_fpga->getCommandQueue();
        int fArg = 0;
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(m_z));
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(m_y));
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(m_x));
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(m_time));
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(m_zb));
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(p_sy));
        m_Kernel.setArg(fArg++, static_cast<uint32_t>(p_sx));
        m_Kernel.setArg(fArg++, d_srcs);
        m_Kernel.setArg(fArg++, d_coefz);
        m_Kernel.setArg(fArg++, d_coefy);
        m_Kernel.setArg(fArg++, d_coefx);
        m_Kernel.setArg(fArg++, d_taperz);
        m_Kernel.setArg(fArg++, d_tapery);
        m_Kernel.setArg(fArg++, d_taperx);
        m_Kernel.setArg(fArg++, d_v2dt2[0]);
        m_Kernel.setArg(fArg++, d_p0[0]);
        m_Kernel.setArg(fArg++, d_p1[0]);
        m_Kernel.setArg(fArg++, d_p0[0]);
        m_Kernel.setArg(fArg++, d_p1[0]);
        m_Kernel.setArg(fArg++, d_pp0[0]);
        m_Kernel.setArg(fArg++, d_pp1[0]);
        m_Kernel.setArg(fArg++, d_pp0[0]);
        m_Kernel.setArg(fArg++, d_pp1[0]);
        m_Kernel.setArg(fArg++, d_upb);
        m_CommandQueue.finish();

        m_CommandQueue.enqueueMigrateMemObjects(inBufVec, 0 /* 0 means from host*/);
        m_CommandQueue.finish();
        auto start = chrono::high_resolution_clock::now();
        m_CommandQueue.enqueueTask(m_Kernel);
        m_CommandQueue.finish();
        auto finish = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = finish - start;
        m_CommandQueue.enqueueMigrateMemObjects(outBufVec, CL_MIGRATE_MEM_OBJECT_HOST);
        m_CommandQueue.finish();
        if(saveSnaps){
            for (int i = 0; i < p_numChannels; i++) {
                if (p_sel) {
                    p_p.insert(p_p.end(), l_p0[i].begin(), l_p0[i].end());
                    p_pp.insert(p_pp.end(), l_pp0[i].begin(), l_pp0[i].end());
                } else {
                    p_p.insert(p_p.end(), l_p1[i].begin(), l_p1[i].end());
                    p_pp.insert(p_pp.end(), l_pp1[i].begin(), l_pp1[i].end());
                }
            }
        }
        return elapsed.count();
    }

    double runModeling(const bool p_sel,
               size_t p_shot,
               size_t p_sy,
               size_t p_sx) 
    {
        HostBuffer_t<t_WideType>p_pp;
        HostBuffer_t<t_WideType>p_p;
        return run(p_sel, p_shot, p_sy, p_sx, false, p_pp, p_p);
    }

    HostBuffer_t<t_DataType> & getUPBBuffer()
    {
        return *m_upb;
    }

    void setTaper(HostBuffer_t<t_DataType> &tp){
        m_taperx.assign(tp.begin(), tp.end());
        m_tapery.assign(tp.begin(), tp.end());
        m_taperz.assign(tp.begin(), tp.end());
    }
    void setCoefs(HostBuffer_t<t_DataType> &cx, 
        HostBuffer_t<t_DataType> &cy, 
        HostBuffer_t<t_DataType> &cz)
    {   
        m_coefx.assign(cx.begin(), cx.end());
        m_coefy.assign(cy.begin(), cy.end());
        m_coefz.assign(cz.begin(), cz.end());
    }
    void setSourceFunction(HostBuffer_t<t_DataType> & src){
        m_src.assign(src.begin(), src.end());
    }
    void setV2Dt2(HostBuffer_t<t_DataType> & p_v2dt2){
        converter<t_PEX, t_PEZ>(m_x, m_y, m_z, p_v2dt2, m_v2dt2.data());
    }

   private:
    const FPGA* m_fpga;
    cl::Kernel m_Kernel;

    size_t m_x, m_y, m_z, m_xb, m_yb, m_zb, m_time, m_shots, m_cube;

    size_t m_v2dt2_base;
    size_t m_pp0_base;
    size_t m_pp1_base;
    size_t m_p0_base;
    size_t m_p1_base;


    HostBuffer_t<t_WideType> m_v2dt2;
    HostBuffer_t<t_DataType> m_taperx, m_tapery, m_taperz;
    HostBuffer_t<t_DataType> m_coefx, m_coefy, m_coefz;
    HostBuffer_t<t_DataType> m_src;
    HostBuffer_t<t_DataType> *m_upb ;

    size_t l_cubeSize;
    size_t p_numChannels;

    
    vector<HostBuffer_t<t_WideType> >l_pp0;
    vector<HostBuffer_t<t_WideType> >l_pp1;
    vector<HostBuffer_t<t_WideType> >l_p0;
    vector<HostBuffer_t<t_WideType> >l_p1;
    vector<HostBuffer_t<t_WideType> >l_v2dt2;

    vector<cl::Memory> inBufVec, outBufVec;

    cl::Buffer d_coefx;
    cl::Buffer d_coefy;
    cl::Buffer d_coefz;
    cl::Buffer d_srcs;
    cl::Buffer d_upb;
    cl::Buffer d_taperx;
    cl::Buffer d_tapery;
    cl::Buffer d_taperz;
    vector<cl::Buffer> d_v2dt2;
    vector<cl::Buffer> d_p0;
    vector<cl::Buffer> d_p1;
    vector<cl::Buffer> d_pp0;
    vector<cl::Buffer> d_pp1;

};

#endif
