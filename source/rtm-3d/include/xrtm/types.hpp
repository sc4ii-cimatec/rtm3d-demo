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
#ifndef XF_RTM_TYPES_HPP
#define XF_RTM_TYPES_HPP

#include "Misc.hpp"

/**
 * @file types.hpp
 * @brief class WideData provides wide datatype for host code
 * @tparam t_DataType the basic datatype
 * @tparam t_DataSize is the number of data in the object
 */

template <typename t_DataType, int t_DataSize>
class WideData {
   private:
    t_DataType m_data[t_DataSize];

   public:
    WideData() {
        for (int i = 0; i < t_DataSize; i++) m_data[i] = 0;
    }

    WideData(const WideData& wd) {
        for (int i = 0; i < t_DataSize; i++) m_data[i] = wd[i];
    }

    WideData(const t_DataType l_value) {
        for (int i = 0; i < t_DataSize; i++) m_data[i] = l_value;
    }

    t_DataType& operator[](int index) { return m_data[index]; }

    const t_DataType& operator[](int index) const { return m_data[index]; }
};

template <size_t t_PEX, size_t t_PEZ, typename t_DataType>
void converter(
    size_t l_x, size_t l_y, size_t l_z, 
    const WideData<WideData<t_DataType, t_PEX>, t_PEZ>* mem, 
    t_DataType * vec) {
    size_t index = 0;
    for (size_t i = 0; i < l_x / t_PEX; i++) {
        for (size_t j = 0; j < l_y; j++) {
            for (size_t k = 0; k < l_z / t_PEZ; k++) {
                for (size_t pez = 0; pez < t_PEZ; pez++) {
                    for (size_t pex = 0; pex < t_PEX; pex++) {
                        size_t offset0 = (i * t_PEX + pex) * l_y * l_z + j * l_z + k * t_PEZ + pez;
                        vec[offset0] = mem[index][pez][pex];
                    }
                }
                index++;
            }
        }
    }
}

template <size_t t_PEX, size_t t_PEZ, typename t_DataType>
void converter(
    size_t l_x, size_t l_y, size_t l_z, 
    const WideData<WideData<t_DataType, t_PEX>, t_PEZ>* mem, 
    HostBuffer_t<t_DataType> & vec) {
    size_t index = 0;
    for (size_t i = 0; i < l_x / t_PEX; i++) {
        for (size_t j = 0; j < l_y; j++) {
            for (size_t k = 0; k < l_z / t_PEZ; k++) {
                for (size_t pez = 0; pez < t_PEZ; pez++) {
                    for (size_t pex = 0; pex < t_PEX; pex++) {
                        size_t offset0 = (i * t_PEX + pex) * l_y * l_z + j * l_z + k * t_PEZ + pez;
                        vec.at(offset0) = mem[index][pez][pex];
                    }
                }
                index++;
            }
        }
    }
}

template <size_t t_PEX, size_t t_PEZ, typename t_DataType>
void converter(
    size_t l_x, size_t l_y, size_t l_z, 
    t_DataType * vec, 
    WideData<WideData<t_DataType, t_PEX>, t_PEZ>* mem) {
    size_t index = 0;
    for (size_t i = 0; i < l_x / t_PEX; i++) {
        for (size_t j = 0; j < l_y; j++) {
            for (size_t k = 0; k < l_z / t_PEZ; k++) {
                for (size_t pez = 0; pez < t_PEZ; pez++) {
                    for (size_t pex = 0; pex < t_PEX; pex++) {
                        size_t offset0 = (i * t_PEX + pex) * l_y * l_z + j * l_z + k * t_PEZ + pez;
                        mem[index][pez][pex] = vec[offset0];
                    }
                }
                index++;
            }
        }
    }
}
template <size_t t_PEX, size_t t_PEZ, typename t_DataType>
void converter(
    size_t l_x, size_t l_y, size_t l_z, 
    const WideData<WideData<t_DataType, t_PEX>, t_PEZ>* mem, 
    vector<t_DataType>& vec) {
    size_t index = 0;
    for (size_t i = 0; i < l_x / t_PEX; i++) {
        for (size_t j = 0; j < l_y; j++) {
            for (size_t k = 0; k < l_z / t_PEZ; k++) {
                for (size_t pez = 0; pez < t_PEZ; pez++) {
                    for (size_t pex = 0; pex < t_PEX; pex++) {
                        size_t offset0 = (i * t_PEX + pex) * l_y * l_z + j * l_z + k * t_PEZ + pez;
                        vec[offset0] = mem[index][pez][pex];
                    }
                }
                index++;
            }
        }
    }
}

template <size_t t_PEX, size_t t_PEZ, typename t_DataType>
void converter(
    size_t l_x, size_t l_y, size_t l_z, 
    const vector<t_DataType>& vec, 
    WideData<WideData<t_DataType, t_PEX>, t_PEZ>* mem) {
    size_t index = 0;
    for (size_t i = 0; i < l_x / t_PEX; i++) {
        for (size_t j = 0; j < l_y; j++) {
            for (size_t k = 0; k < l_z / t_PEZ; k++) {
                for (size_t pez = 0; pez < t_PEZ; pez++) {
                    for (size_t pex = 0; pex < t_PEX; pex++) 
                    {
                        size_t offset0 = (i * t_PEX + pex) * l_y * l_z + j * l_z + k * t_PEZ + pez;
                        mem[index][pez][pex] = vec[offset0];
                    }
                }
                index++;
            }
        }
    }
}

template <size_t t_PEX, size_t t_PEZ, typename t_DataType>
void converter(
    size_t l_x, size_t l_y, size_t l_z, 
    const HostBuffer_t<t_DataType> & vec, 
    WideData<WideData<t_DataType, t_PEX>, t_PEZ>* mem) {
    size_t index = 0;
    for (size_t i = 0; i < l_x / t_PEX; i++) {
        for (size_t j = 0; j < l_y; j++) {
            for (size_t k = 0; k < l_z / t_PEZ; k++) {
                for (size_t pez = 0; pez < t_PEZ; pez++) {
                    for (size_t pex = 0; pex < t_PEX; pex++) 
                    {
                        size_t offset0 = (i * t_PEX + pex) * l_y * l_z + j * l_z + k * t_PEZ + pez;
                        mem[index][pez][pex] = vec[offset0];
                    }
                }
                index++;
            }
        }
    }
}


template <size_t t_PEX, size_t t_PEZ, size_t t_Order, typename t_DataType>
void converter_upb(size_t p_x, size_t p_y, size_t p_time, HostBuffer_t<t_DataType>& pin, t_DataType* out) 
{
    size_t index = 0;
    for (size_t i = 0; i < p_time; i++){
        for (size_t k = 0; k < p_x / t_PEX; k++){
            for (size_t j = 0; j < p_y; j++){
                for (size_t po = 0; po < t_Order / 2; po++){
                    for (size_t pe = 0; pe < t_PEX; pe++){
                        size_t offset0 = (i * p_x * p_y * t_Order / 2) + ((k * t_PEX + pe) * p_y * t_Order / 2) +
                            (j * t_Order / 2 + po);
                        out[offset0] = pin[index++];
                    }
                }
            }
        }
    }
}

template <size_t t_PEX, size_t t_PEZ, size_t t_Order, typename t_DataType>
void converter_upb(size_t p_x, size_t p_y, size_t p_time, HostBuffer_t<t_DataType>& pin, HostBuffer_t<t_DataType>& pout) 
{
    size_t index = 0;
    for (size_t i = 0; i < p_time; i++){
        for (size_t k = 0; k < p_x / t_PEX; k++){
            for (size_t j = 0; j < p_y; j++){
                for (size_t po = 0; po < t_Order / 2; po++){
                    for (size_t pe = 0; pe < t_PEX; pe++){
                        size_t offset0 = (i * p_x * p_y * t_Order / 2) + ((k * t_PEX + pe) * p_y * t_Order / 2) +
                            (j * t_Order / 2 + po);
                        pout.at(offset0) = pin.at(index++);
                    }
                }
            }
        }
    }
}

#endif
