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
 * @file vecMoverB1.hpp
 * @brief common data movers for vectors used in BLAS L1 routines.
 *
 * This file is part of Vitis BLAS Library.
 */

#ifndef XF_BLAS_VECMOVERB1_HPP
#define XF_BLAS_VECMOVERB1_HPP

#include "hls_stream.h"
#include "ap_int.h"
#include "ap_shift_reg.h"

namespace xf {

namespace blas {

template <typename t_DataType, typename t_DesDataType = t_DataType>
void mem2stream(unsigned int p_n, const t_DataType* p_in, hls::stream<t_DesDataType>& p_out) {
    for (unsigned int i = 0; i < p_n; ++i) {
#pragma HLS PIPELINE
        t_DesDataType l_val = p_in[i];
        p_out.write(l_val);
    }
} // end mem2stream

template <typename t_DataType, typename t_DesDataType = t_DataType>
void stream2mem(unsigned int p_n, hls::stream<t_DataType>& p_in, t_DesDataType* p_out) {
    for (unsigned int i = 0; i < p_n; ++i) {
#pragma HLS PIPELINE
        t_DesDataType l_val = p_in.read();
        p_out[i] = l_val;
    }
} // end stream2mem


} // namespace blas

} // namespace xf
#endif
