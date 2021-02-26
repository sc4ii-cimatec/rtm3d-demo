/**
 * @file RTM.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-05-07
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef RTM_H
#define RTM_H

#include <cstdlib>
#include <algorithm>  // all_of, find, for_each
#include <cassert>    // assert
#include <cstddef>    // nullptr_t, ptrdiff_t, size_t
#include <functional> // hash, less
#include <iosfwd>     // istream, ostream
#include <iterator>   // random_access_iterator_tag
#include <string>     // string, stoi, to_string
#include <utility>    // declval, forward, move, pair, swap
#include <iostream>
#include <algorithm>
#include <exception>
#include <vector>
#include <math.h>
#include <tuple>
#include <Misc.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMParam.hpp>
#include <RTMException.hpp>
#include <RTMGridPartition.hpp>
#include <RTMKernel.hpp>
#include <RTMController.hpp>
#include <RTMUtil.hpp>
#include <RTMSeismic.hpp>
#include <RTMProcess.hpp>
#include <RTMAcc.hpp>
#include <RTMPlatform.hpp>
#include <RTMCPUPlatform.hpp>
#include <RTMGPUPlatform.hpp>
#include <RTMFPGAPlatform.hpp>

#if defined(RTM_ACC_GPU) || defined (RTM_ACC_FPGA)
#define RTM_ACC
#endif


/**
 * @brief Function static inline RTMGRID_SWAP()
 * @param A P/PR grid
 * @param B PP/PPR grid
 */
static inline void RTMGRID_SWAP(RTMCube<RTMData_t, RTMDevPtr_t> ** A, RTMCube<RTMData_t, RTMDevPtr_t> ** B){
    RTMCube<RTMData_t, RTMDevPtr_t> * swap = *B;
	*B = *A;
	*A = swap;
}

#endif