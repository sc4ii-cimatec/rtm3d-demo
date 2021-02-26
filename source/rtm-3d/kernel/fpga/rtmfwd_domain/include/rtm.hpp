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

#ifndef XF_RTM_RTM3D_HPP
#define XF_RTM_RTM3D_HPP

/**
 * @file rtm.hpp
 * @brief class RTM3D derived from Stencil3D is defined here, it provides L1 primary
 * streaming modules for 3D-RTM kernels
 */

#ifndef __SYNTHESIS__
#include <cassert>
#endif
#include <cstring>
#include "hls_stream.h"
#include "stencil.hpp"
#include "dataMover.hpp"

namespace xf {
namespace rtm {

/**
 * @brief RTM3D class defines the basic operations for 2D RTM3D
 *
 * @tparam t_DataType the basic wavefield datatype
 * @tparam t_Order is the spatial discretization order
 * @tparam t_MaxDim is the maximum height this kernel can process
 * @tparam t_PE is the number of processing elements
 */

template <typename t_DataType, int t_Order, int t_MaxDimZ = 128, int t_MaxDimY = 128, int t_MaxB = 40, int t_PE = 1>
class RTM3D : public Stencil3D<t_DataType, t_Order, t_MaxDimZ, t_MaxDimY, t_PE> {
   public:
    RTM3D(int p_NZB = t_MaxB, int p_NYB = t_MaxB, int p_NXB = t_MaxB) {
#ifndef __SYNTHESIS__
        assert(t_MaxB >= p_NXB);
        assert(t_MaxB >= p_NYB);
        assert(t_MaxB >= p_NZB);
        assert(p_NXB % t_PE == 0);
        assert(p_NYB % t_PE == 0);
        assert(p_NZB % t_PE == 0);
#endif
        m_NXB = p_NXB;
        m_NYB = p_NYB;
        m_NZB = p_NZB;
    }
    typedef Stencil3D<t_DataType, t_Order, t_MaxDimZ, t_MaxDimY, t_PE> t_StencilType;
    using t_StencilType::t_NumData;
    using t_StencilType::t_FifoDepth;
    using t_StencilType::t_HalfOrder;
    typedef typename t_StencilType::t_PairType t_PairType;
    typedef typename t_StencilType::t_WideType t_WideType;
    typedef typename t_WideType::t_TypeInt t_InType;
    typedef typename t_PairType::t_TypeInt t_PairInType;

    typedef blas::WideType<t_InType, t_HalfOrder / t_PE> t_UpbType;
    typedef typename t_UpbType::t_TypeInt t_UpbInType;

    void setBoundaryDim(int p_NZB, int p_NYB, int p_NXB) {
#ifndef __SYNTHESIS__
        assert(this->m_x >= 2 * p_NXB);
        assert(this->m_y >= 2 * p_NYB);
        assert(this->m_z >= 2 * p_NZB);
        assert(t_MaxB >= p_NXB);
        assert(t_MaxB >= p_NYB);
        assert(t_MaxB >= p_NZB);
        assert(p_NZB % t_PE == 0);
#endif
        m_NXB = p_NXB;
        m_NYB = p_NYB;
        m_NZB = p_NZB;
    }

    inline int getZB() {
#pragma HLS INLINE
        return m_NZB;
    }
    inline int getYB() {
#pragma HLS INLINE
        return m_NYB;
    }
    inline int getXB() {
#pragma HLS INLINE
        return m_NXB;
    }

    void setTaper(const t_DataType* p_taperz, const t_DataType* p_tapery, const t_DataType* p_taperx) {
        for (int i = 0; i < m_NXB; i++) {
#pragma HLS PIPELINE
            m_taperx[i] = p_taperx[i];
        }

        for (int i = 0; i < m_NYB; i++) {
#pragma HLS PIPELINE
            m_tapery[i] = p_tapery[i];
        }

        for (int i = 0; i < m_NZB / t_PE; i++) {
            t_WideType l_w;
            for (int pe = 0; pe < t_PE; pe++) {
#pragma HLS PIPELINE
                l_w.unshift(p_taperz[i * t_PE + pe]);
            }
            m_taperz[i] = l_w;
        }
    }
    void setSrc(int p_srcz, int p_srcy, int p_srcx) {
        m_srcz = p_srcz;
        m_srcy = p_srcy;
        m_srcx = p_srcx;
    }

   private:
    /**
     * @brief taper function applys to the absorbing boundary
     *
     * @param p_in is the stream of input wavefield
     * @param p_out is the stream of output wavefield
     */
    void taper(hls::stream<t_PairInType>& p_in, hls::stream<t_PairInType>& p_out) {
        for (int i = 0, j = 0, k = 0, t = 0; t < this->m_cube; t++) {
#pragma HLS PIPELINE
            t_PairType l_in = p_in.read();
            t_PairType l_out;
            t_DataType tx, ty;
            t_WideType tz;

            if (k < m_NZB / t_PE) {
                tz = m_taperz[k];
                if (i < m_NXB)
                    tx = m_taperx[i];
                else if (i >= this->m_x - m_NXB)
                    tx = m_taperx[this->m_x - i - 1];
                else
                    tx = 1;

                if (j < m_NYB)
                    ty = m_tapery[j];
                else if (j >= this->m_y - m_NYB)
                    ty = m_tapery[this->m_y - j - 1];
                else
                    ty = 1;
            } else {
                tz = 1;
                ty = 1;
                tx = 1;
            }

            for (int n = 0; n < t_NumData; n++) {
#pragma HLS UNROLL
                t_WideType l_w = l_in[n];
                t_WideType l_wo;
                for (int pe = 0; pe < t_PE; pe++) {
#pragma HLS UNROLL
                    l_wo[pe] = l_w[pe] * tz[pe] * tx * ty;
                }
                l_out[n] = l_wo;
            }
            p_out.write(l_out);

            if (k == this->m_zPE - 1 && j == this->m_y - 1) {
                k = 0;
                j = 0;
                i++;
            } else if (k == this->m_zPE - 1) {
                k = 0;
                j++;
            } else
                k++;
        }
    }

    void extractUPB(hls::stream<t_PairInType>& p_pin,
                    hls::stream<t_PairInType>& p_pout,
                    hls::stream<t_UpbInType>& p_upb) {
        t_UpbType l_upb;
        for (int i = 0, j = 0, k = 0, t = 0; t < this->m_cube; t++) {
#pragma HLS PIPELINE
            t_PairType l_val = p_pin.read();
            p_pout.write(l_val);
            l_upb.unshift(l_val[1]);

            if (k == m_NZB / t_PE - 1) {
                p_upb.write(l_upb);
            }

            if (k == this->m_zPE - 1 && j == this->m_y - 1) {
                k = 0;
                j = 0;
                i++;
            } else if (k == this->m_zPE - 1) {
                k = 0;
                j++;
            } else
                k++;
        }
    }

    void addSrc(const t_DataType p_src, hls::stream<t_PairInType>& p_p0, hls::stream<t_PairInType>& p_p1) {
        for (int i = 0, j = 0, k = 0, t = 0; t < this->m_cube; t++) {
#pragma HLS PIPELINE
            t_PairType l_val = p_p0.read();
            t_PairType l_oVal;
            l_oVal[0] = l_val[0];
            t_WideType l_wo = l_val[1];
            if (i == m_srcx && j == m_srcy && k == m_srcz / t_PE) {
                l_wo[0] = ((t_WideType)l_val[1])[0] + p_src;
            }

            l_oVal[1] = l_wo;
            p_p1.write(l_oVal);

            if (k == this->m_zPE - 1 && j == this->m_y - 1) {
                k = 0;
                j = 0;
                i++;
            } else if (k == this->m_zPE - 1) {
                k = 0;
                j++;
            } else
                k++;
        }
    }
    void copy(hls::stream<t_PairInType>& p_p0, hls::stream<t_PairInType>& p_p1) {
        for (int t = 0; t < this->m_cube; t++) {
#pragma HLS PIPELINE
            t_PairType l_val = p_p0.read();
            p_p1.write(l_val);
        }
    }

   public:
    /**
     * @brief forward defines the forward streaming module
     *
     * @param p_src is the source wavefield at given time stamp
     * @param p_v2dt2 is the pow(v * dt, 2)
     * @param p_vt is a copy of p_v2dt2
     * @param p_p0 is the stream of input wavefield p(t-1) and p(t)
     * @param p_p1 is the stream of output wavefield p(t) and p(t+1)
     */

    void forward(const t_DataType p_src,
                 hls::stream<t_InType>& p_v2dt2,
                 hls::stream<t_InType>& p_vt,
                 hls::stream<t_PairInType>& p_p0,
                 hls::stream<t_PairInType>& p_p1) {
#pragma HLS DATAFLOW
        hls::stream<t_PairInType> l_p0, l_p1;
        copy(p_p0, l_p0);
        this->propagate(p_v2dt2, p_vt, l_p0, l_p1);
        addSrc(p_src, l_p1, p_p1);
    }
    /**
     * @brief forward defines the forward streaming module
     *
     * @param p_src is the source wavefield at given time stamp
     * @param p_upb is the upper boundary to be saved
     * @param p_v2dt2 is the pow(v * dt, 2)
     * @param p_vt is a copy of p_v2dt2
     * @param p_p0 is the stream of input wavefield p(t-1) and p(t)
     * @param p_p1 is the stream of output wavefield p(t) and p(t+1)
     */

    void forward(const t_DataType p_src,
                 hls::stream<t_UpbInType>& p_upb,
                 hls::stream<t_InType>& p_v2dt2,
                 hls::stream<t_InType>& p_vt,
                 hls::stream<t_PairInType>& p_p0,
                 hls::stream<t_PairInType>& p_p1) 
    {
#pragma HLS DATAFLOW
        hls::stream<t_PairInType> l_p0, l_p1, l_upb;
        taper(p_p0, l_p0);
        this->propagate(p_v2dt2, p_vt, l_p0, l_p1);
        addSrc(p_src, l_p1, l_upb);
        extractUPB(l_upb, p_p1, p_upb);
    }

    template <unsigned int t_NumStream, typename T>
    static void saveUpb(unsigned int p_n, int p_t, hls::stream<T> p_s[t_NumStream], T* p_mem) {
        int l_totalSize = p_n * t_NumStream;
        int l_index[t_NumStream];
#pragma HLS ARRAY_PARTITION variable = l_index complete dim = 1
        for (int i = 0; i < t_NumStream; i++)
#pragma HLS UNROLL
            l_index[i] = 0;

        while (l_totalSize > 0) {
#pragma HLS PIPELINE
            T l_val;
            for (int i = 0; i < t_NumStream; i++) {
#pragma HLS UNROLL
                if (p_s[i].read_nb(l_val)) {
                    l_totalSize--;
                    p_mem[p_t * p_n * t_NumStream + i * p_n + l_index[i]] = l_val;
                    l_index[i]++;
                    break;
                }
            }
        }
    }

   protected:
    int m_NXB, m_NYB, m_NZB;
    int m_srcx, m_srcy, m_srcz;
    int m_recz, m_recy, m_recx;

    t_DataType m_taperx[t_MaxB];
    t_DataType m_tapery[t_MaxB];
    t_InType m_taperz[t_MaxB / t_PE];
};
}
}
#endif
