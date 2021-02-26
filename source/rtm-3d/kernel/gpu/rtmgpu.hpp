#ifndef _RTMGPU_FINITEDIFF_HPP
#define _RTMGPU_FINITEDIFF_HPP

#ifdef RTM_ACC_GPU
#include <cuda.h>
#include <assert.h>

void RTMSTEP_WRAPPER(dim3 dimGrid, dim3 dimBlock,
					 uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen,
					 float *__restrict__ P, float *__restrict__ PP, float *__restrict__ V2DT2,
					 float *__restrict__ coefs, bool rtm2D);

void RTMSEISM_WRAPPER(dim3 dimGrid, dim3 dimBlock, 
	uint32_t rcvOffsetX, uint32_t rcvDistX, int rcvCountX, 
	uint32_t rcvOffsetY, uint32_t rcvDistY, int rcvCountY, uint32_t rcvDepthZ,
	uint32_t pStartX, uint32_t pEndX,
	uint32_t pStartY, uint32_t pEndY,
	uint32_t border_length, uint32_t nze, uint32_t nt, uint32_t it, bool modeling,
	float *__restrict__ SEISM, float *__restrict__ PPR);

void RTMSOURCE_WRAPPER(dim3 dimGrid, dim3 dimBlock,
	uint32_t sx, uint32_t sy, uint32_t sz, float eval, 
	uint32_t pStartX, uint32_t pEndX,
	uint32_t pStartY, uint32_t pEndY,
	uint32_t nze,float * __restrict__ PP);

void RTMTAPER_WRAPPER(dim3 dimGrid, dim3 dimBlock,
	uint32_t pStartX, uint32_t pEndX,uint32_t pStartY, uint32_t pEndY, 
	uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen, 
	float * __restrict__  taper, float * __restrict__  P, bool upperOnly=false);

void RTMWUPB_WRAPPER(dim3 dimGrid, dim3 dimBlock, uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen, uint32_t it, bool rw,
	float *__restrict__ PP,float *__restrict__ UPB);

void RTMIMG_WRAPPER(dim3 dimGrid, dim3 dimBlock, 
	uint32_t nxe, uint32_t nye, uint32_t nze,
	float * __restrict__ IMG, float * __restrict__ PS, float * __restrict__ PR);

void RTM_FREQIMG_WRAPPER(dim3 dimGrid, dim3 dimBlock, 
		uint64_t iw,uint64_t lw,
		uint64_t gStartX, uint64_t gEndX, 
		uint64_t gStartY, uint64_t gEndY, 
		uint64_t gStartZ, uint64_t gEndZ,
		uint64_t nxe, uint64_t nye, uint64_t nze, uint64_t blen,   
		float *__restrict__ wList,
		float *__restrict__ IMG,
		float *__restrict__ PSRe, float *__restrict__ PSIm,
		float *__restrict__ PRRe, float *__restrict__ PRIm);

void RTM_UPDATEFREQ_WRAPPER(dim3 dimGrid, dim3 dimBlock, 
		uint64_t it, uint64_t iw,uint64_t lw, uint64_t nt,
		uint64_t gStartX, uint64_t gEndX, uint64_t gStartY, uint64_t gEndY, 
		uint64_t gStartZ, uint64_t gEndZ,
		uint64_t nxe, uint64_t nye, uint64_t nze, uint64_t blen, 
		float *__restrict__ kernelRe, float *__restrict__ kernelIm,
		float *__restrict__ PS, float *__restrict__ PR,
		float *__restrict__ PSRe, float *__restrict__ PSIm,
		float *__restrict__ PRRe, float *__restrict__ PRIm);

void RTMSTEP_MULTIWAVE_WRAPPER(dim3 dimGrid, dim3 dimBlock,
	uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze,
	float *__restrict__ P0, float *__restrict__ PP0,
	float *__restrict__ P1, float *__restrict__ PP1, 
	float *__restrict__ V2DT2,
	float *__restrict__ coefs, bool rtm2D);

#endif // RTM_ACC_GPU
#endif