/**
 * @file Misc.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-05-06
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#ifndef MISC_H
#define MISC_H

#include <cstdlib>
#include <chrono>
#include <string.h>
#include <algorithm>
#include <vector>
#include <chrono>

/**
 * @brief Functions declaration  
 * @see  Misc.cpp
 */
typedef std::chrono::high_resolution_clock 	Clock;
typedef std::chrono::milliseconds 			milliseconds;
typedef std::chrono::seconds	 			seconds;
typedef Clock::time_point 					timepoint;


template <typename T>
struct rtm_aligned_allocator {
    using value_type = T;
    T* allocate(std::size_t num) {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, 4096, num * sizeof(T))) throw std::bad_alloc();
        return reinterpret_cast<T*>(ptr);
    }
    void deallocate(T* p, std::size_t num) { free(p); }
};

template <typename T>
using HostBuffer_t = std::vector<T, rtm_aligned_allocator<T> >;

template<typename T>
HostBuffer_t<T> & SLICE(HostBuffer_t<T> &v, int m, int n);

/**
 * @brief Maximum value between two values
 * @param a  Float value
 * @param b  Float value
 * @return   Maximum float value between a and b
 */
float max(float a, float b);

/**
 * @brief Minimum value between two values
 * @param a Float value
 * @param b Float value
 * @return  Minimum float value between a and b
 */
float min(float a, float b);

/**
 * @brief Function used for calculating the Random border condition
 * @param v      Float value 
 * @param maxvel Float value used for calculating the lowest allowed velocity value
 * @param minvel Float value
 * @param dist   Float value
 * @param range  Float value
 * @return       Float velocity value
 */
float randbetween(float v, float maxvel, float minvel, int dist, int range);

/**
 * @brief Function used for calculating the Random border condition
 * @param v      Float value 
 * @param maxvel Float value used for calculating the lowest allowed velocity value
 * @param minvel Float value
 * @param dist   Float value
 * @param range  Float value
 * @return       Float velocity value
 */
float randparabolicval(float v, float maxvel, float minvel, int dist, int range);

/**
 * @brief Function used for calculating a random velocity value 
 * @param v   Float value
 * @param k   Float value
 * @param nkb Float value
 * @return float, velocity value random
 */
float randlinearval(float v, int k, int nkb);

/**
 * @brief Function used for calculating the laplacian-coefficient 
 * @param coef  Array of pointers
 * @param order Compilation parameter of RTM
 */
void makeo2 (float *coef,int order);

/**
 * @brief Selection of laplacian coefficient 
 * @param coefs  Array of pointers
 * @param _order Integer value can be 2, 4, 6 and 8
 */
void laplacian_coefs(float * coefs, int _order);

timepoint tic();
timepoint toc();

/**
 * @brief Calculation of time in seconds
 * @param t0 Timepoint value
 * @param t1 Timepoint value
 * @return float Returns the time in seconds
 */
float elapsed_s(timepoint t0, timepoint t1);

/**
 * @brief Calculation of time in milliseconds
 * @param t0 Timepoint value
 * @param t1 Timepoint value
 * @return float  Returns the time in milliseconds
 */
float elapsed_ms(timepoint t0, timepoint t1);


/**
 * @brief Swaps two pointers
 * @param A P/PR grid
 * @param B PP/PPR grid
 */
void SWAP_PTR(void ** A, void ** B);

void TO_UPPER(std::string &data);

void TO_LOWER(std::string &data);

#endif