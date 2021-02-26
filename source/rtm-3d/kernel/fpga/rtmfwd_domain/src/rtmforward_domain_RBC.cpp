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

/**
 * @file rtmforward.cpp
 * @brief It defines the forward kernel function
 */

#include "rtmforward_domain_rbc.hpp"

template <typename t_Domain>
void forward(bool p_sel,
             const t_Domain& p_domain,
             RTM_TYPE s[RTM_numFSMs],
             const unsigned int p_t,
             const RTM_dataType* p_src,
             const RTM_type* p_v2dt2,
             RTM_type* p_pi0,
             RTM_type* p_pi1,
             RTM_type* p_po0,
             RTM_type* p_po1,
             RTM_type* p_ppi0,
             RTM_type* p_ppi1,
             RTM_type* p_ppo0,
             RTM_type* p_ppo1) {
    const int l_entries = s[0].getCube();

    hls::stream<RTM_type> l_vt_in, l_vt_out;
    hls::stream<RTM_type> l_pp_in, l_pp_out;

    hls::stream<RTM_type> l_p[RTM_numFSMs + 1];
#pragma HLS ARRAY_PARTITION variable = l_p complete dim = 1

#pragma HLS DATAFLOW

    p_domain.template memSelStream(p_sel, p_pi1, p_pi0, l_p[0]);
    p_domain.template memSelStream(p_sel, p_ppi1, p_ppi0, l_pp_in);

    p_domain.template mem2stream(p_v2dt2, l_vt_in);

#if RTM_numFSMs == 1
    s[0].forward(p_src[0], l_vt_in, l_vt_out, l_pp_in, l_p[0], l_pp_out, l_p[1]);
#else

    hls::stream<RTM_type> l_vt[RTM_numFSMs - 1];
#pragma HLS ARRAY_PARTITION variable = l_vt complete dim = 1
#pragma HLS stream depth = RTM_TYPE::t_FifoDepth variable = l_vt
#pragma HLS RESOURCE variable = l_vt core = fifo_uram

    hls::stream<RTM_type> l_pp[RTM_numFSMs - 1];
#pragma HLS ARRAY_PARTITION variable = l_pp complete dim = 1
#pragma HLS stream depth = RTM_TYPE::t_FifoDepth variable = l_pp
#pragma HLS RESOURCE variable = l_pp core = fifo_uram

    s[0].forward(p_src[0], l_vt_in, l_vt[0], l_pp_in, l_p[0], l_pp[0], l_p[1]);
    for (int i = 1; i < RTM_numFSMs - 1; i++) {
#pragma HLS UNROLL
        s[i].forward(p_src[i], l_vt[i - 1], l_vt[i], l_pp[i - 1], l_p[i], l_pp[i], l_p[i + 1]);
    }
    s[RTM_numFSMs - 1].forward(p_src[RTM_numFSMs - 1], l_vt[RTM_numFSMs - 2], l_vt_out, l_pp[RTM_numFSMs - 2],
                               l_p[RTM_numFSMs - 1], l_pp_out, l_p[RTM_numFSMs]);
#endif
    dataConsumer(l_entries, l_vt_out);

    p_domain.template streamSelMem(p_sel, l_p[RTM_numFSMs], p_po0, p_po1);
    p_domain.template streamSelMem(p_sel, l_pp_out, p_ppo0, p_ppo1);
}

extern "C" void rtmforward(const unsigned int p_z,
                           const unsigned int p_y,
                           const unsigned int p_x,
                           const unsigned int p_t,
                           const unsigned int p_srcz,
                           const unsigned int p_srcy,
                           const unsigned int p_srcx,
                           const RTM_dataType* p_src,
                           const RTM_dataType* p_coefz,
                           const RTM_dataType* p_coefy,
                           const RTM_dataType* p_coefx,
                           const RTM_type* p_v2dt2,
                           RTM_type* p_pi0,
                           RTM_type* p_pi1,
                           RTM_type* p_po0,
                           RTM_type* p_po1,
                           RTM_type* p_ppi0,
                           RTM_type* p_ppi1,
                           RTM_type* p_ppo0,
                           RTM_type* p_ppo1) {
    SCALAR(p_x)
    SCALAR(p_y)
    SCALAR(p_z)
    SCALAR(p_t)
    SCALAR(p_srcx)
    SCALAR(p_srcy)
    SCALAR(p_srcz)
    SCALAR(return )

    POINTER(p_coefx, gmemParam)
    POINTER(p_coefy, gmemParam)
    POINTER(p_coefz, gmemParam)
    POINTER(p_src, gmemParam)

    POINTER(p_pi0, gmem_pi0)
    POINTER(p_pi1, gmem_pi1)
    POINTER(p_po0, gmem_po0)
    POINTER(p_po1, gmem_po1)

    POINTER(p_ppi0, gmem_ppi0)
    POINTER(p_ppi1, gmem_ppi1)
    POINTER(p_ppo0, gmem_ppo0)
    POINTER(p_ppo1, gmem_ppo1)
    POINTER(p_v2dt2, gmem_v2dt2)

    RTM_TYPE l_s[RTM_numFSMs];
#pragma HLS ARRAY_PARTITION variable = l_s complete dim = 1
    RTM_dataType l_src[RTM_numFSMs];
#pragma HLS ARRAY_PARTITION variable = l_src complete dim = 1

    DOMAIN_TYPE l_domain(p_x, p_y, p_z);

    for (int i = 0; i < RTM_numFSMs; i++) {
        l_s[i].setCoef(p_coefz, p_coefy, p_coefx);
        l_s[i].setSrc(p_srcz, p_srcy, p_srcx);
    }

    for (int t = 0; t < p_t / RTM_numFSMs; t++) {
        bool cont = l_domain.reset();
        while (cont) {
            for (int i = 0; i < RTM_numFSMs; i++) {
                l_s[i].setDomain(l_domain);
                l_s[i].setDim(l_domain.m_z, l_domain.m_extDim, l_domain.m_x);
                l_src[i] = p_src[t * RTM_numFSMs + i];
            }
            forward(t & 0x01, l_domain, l_s, t, l_src, p_v2dt2, p_pi0, p_pi1, p_po0, p_po1, p_ppi0, p_ppi1, p_ppo0,
                    p_ppo1);
            cont = l_domain.next();
        }
    }
}
