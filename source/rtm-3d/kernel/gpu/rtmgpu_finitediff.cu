#ifdef RTM_ACC_GPU
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>
#include <rtmgpu.hpp>

/**	<<<(NXE,NYE,NZE),(10,10,10)>>> */
__global__ void kgpu_updateFreqContributions(uint64_t it, uint64_t iw,uint64_t lw, uint64_t nt,
	uint64_t gStartX, uint64_t gEndX, uint64_t gStartY, uint64_t gEndY, uint64_t gStartZ, uint64_t gEndZ,
	uint64_t nxe, uint64_t nye, uint64_t nze, uint64_t blen, 
	float *__restrict__ kernelRe, float *__restrict__ kernelIm,
	float *__restrict__ PS, float *__restrict__ PR,
	float *__restrict__ PSRe, float *__restrict__ PSIm,
	float *__restrict__ PRRe, float *__restrict__ PRIm)
{
	uint64_t ix = (uint64_t)((blockIdx.x * blockDim.x) + threadIdx.x); // Global X index
	uint64_t iy = (uint64_t)((blockIdx.y * blockDim.y) + threadIdx.y); // Global Y index
	uint64_t iz = (uint64_t)((blockIdx.z * blockDim.z) + threadIdx.z); // Global Z index
	uint64_t gnx = gEndX-gStartX;
	uint64_t gny = gEndY-gStartY;
	uint64_t gnz = gEndZ-gStartZ;

	uint64_t px = gStartX + blen + ix;
	uint64_t py = gStartY + blen + iy;
	uint64_t pz = gStartZ + blen + iz;

	if (ix < gnx)
	{
		if (iy < gny )
		{
			if (iz < gnz)
			{
				// printf("ix=%lu iy=%lu iz=%lu iw=%lu lw=%lu \n", ix, iy, iz, iw, lw);
				uint64_t koffset = iw*nt + it;
				uint64_t poffset = px*nye*nze + py*nze + pz;
				uint64_t foffset = (lw*gnx*gny*gnz) + (ix*gny*gnz) + (iy*gnz) + iz;
				PSRe[foffset] += (PS[poffset]*kernelRe[koffset]);
				PRRe[foffset] += (PR[poffset]*kernelRe[koffset]);
				PSIm[foffset] += (PS[poffset]*kernelIm[koffset]);
				PRIm[foffset] += (PR[poffset]*kernelIm[koffset]);
			}
		}
	}
}
__global__ void kgpu_freqDomainImageCondition(uint64_t iw,uint64_t lw,
	uint64_t gStartX, uint64_t gEndX, 
	uint64_t gStartY, uint64_t gEndY, 
	uint64_t gStartZ, uint64_t gEndZ,
	uint64_t nxe, uint64_t nye, uint64_t nze, uint64_t blen,  
	float *__restrict__ w2List,
	float *__restrict__ IMG,
	float *__restrict__ PSRe, float *__restrict__ PSIm,
	float *__restrict__ PRRe, float *__restrict__ PRIm)
{
	uint64_t ix = (uint64_t)((blockIdx.x * blockDim.x) + threadIdx.x); // Global X index
	uint64_t iy = (uint64_t)((blockIdx.y * blockDim.y) + threadIdx.y); // Global Y index
	uint64_t iz = (uint64_t)((blockIdx.z * blockDim.z) + threadIdx.z); // Global Z index
	uint64_t gnx = gEndX-gStartX;
	uint64_t gny = gEndY-gStartY;
	uint64_t gnz = gEndZ-gStartZ;

	uint64_t px = gStartX + blen + ix;
	uint64_t py = gStartY + blen + iy;
	uint64_t pz = gStartZ + blen + iz;
	if (ix < gnx)
	{
		if (iy < gny )
		{
			if (iz < gnz)
			{
				// printf("ix=%lu iy=%lu iz=%lu iw=%lu lw=%lu \n", ix, iy, iz, iw, lw);
				uint64_t poffset = px*nye*nze + py*nze + pz;
				uint64_t foffset = (lw*gnx*gny*gnz) + (ix*gny*gnz) + (iy*gnz) + iz;
				IMG[poffset] += (w2List[iw])*(PSRe[foffset]*PRRe[foffset] - PSIm[foffset]*PRIm[foffset]);
			}
		}
	}
}

/**	<<<(NYE,NZE,1),(32,32,1)>>> */
__global__ void kgpu_rtmStep2D(
	uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen,
	float *__restrict__ P, float *__restrict__ PP, float *__restrict__ V2DT2,
	float *__restrict__ coefs, bool rtm2D)
{
	uint32_t half_order = st_order / 2;
	uint32_t ix = blen; // Global X index
	uint32_t iy = blockIdx.x * blockDim.x + threadIdx.x; // Global Y index
	uint32_t iz = blockIdx.y * blockDim.y + threadIdx.y; // Global Z index
	uint32_t kx, ky, kz;
	uint32_t cx_base = 0;
	uint32_t cy_base = (st_order) + 1;
	uint32_t cz_base = (2 * (st_order + 1));

	if (ix >= nxe || iy >= nye || iz >= nze){
		return;
	}
	if (ix >= half_order && ix < (nxe - half_order))
	{
		if (iy >= half_order && iy < (nye - half_order))
		{
			if (iz >= half_order && iz < (nze - half_order))
			{
				uint32_t offset = ix * nye * nze + iy * nze + iz;
				float pp = PP[offset];
				float p = P[offset];
				float v2dt2 = V2DT2[offset];
				float acmy = 0.0, acmz = 0.0, lapl = 0.0;
				#pragma unroll (8)
				for (uint32_t io = 0; io <= st_order; io++)
				{
					ky = (ix * nye * nze) + (iy - half_order + io) * nze + iz;
					kz = ix * nye * nze + iy * nze + (iz - half_order + io);
					acmy += P[ky] * coefs[cy_base + io];
					acmz += P[kz] * coefs[cz_base + io];
				}
				lapl = acmz + acmy;
				PP[offset] = (2 * p) - pp + (v2dt2 * lapl);
			}
		}
	}
}


/**	<<<(NXE,NYE,NZE),(32,32,1)>>> */
__global__ void kgpu_rtmStep(
	uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze,
	float *__restrict__ P, float *__restrict__ PP, float *__restrict__ V2DT2,
	float *__restrict__ coefs, bool rtm2D)
{
	uint32_t half_order = st_order / 2;
	uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index
	uint32_t iz = blockIdx.z * blockDim.z + threadIdx.z; // Global Z index
	uint32_t kx, ky, kz;
	uint32_t cx_base = 0;
	uint32_t cy_base = (st_order) + 1;
	uint32_t cz_base = (2 * (st_order + 1));

	if (ix >= nxe || iy >= nye || iz >= nze){
		return;
	}
	if (ix >= half_order && ix < (nxe - half_order))
	{
		if (iy >= half_order && iy < (nye - half_order))
		{
			if (iz >= half_order && iz < (nze - half_order))
			{
				uint32_t offset = ix * nye * nze + iy * nze + iz;
				float pp = PP[offset];
				float p = P[offset];
				float v2dt2 = V2DT2[offset];
				float acmx = 0.0, acmy = 0.0, acmz = 0.0, lapl = 0.0;
				#pragma unroll (8)
				for (uint32_t io = 0; io <= st_order; io++)
				{
					kx = ((ix - half_order + io) * nye * nze) + iy * nze + iz;
					ky = (ix * nye * nze) + (iy - half_order + io) * nze + iz;
					kz = ix * nye * nze + iy * nze + (iz - half_order + io);
					acmx += P[kx] * coefs[cx_base + io];
					acmy += P[ky] * coefs[cy_base + io];
					acmz += P[kz] * coefs[cz_base + io];
				}
				if (rtm2D)
				{
					acmx = 0.0;
				}
				lapl = acmz + acmy + acmx;
				PP[offset] = (2 * p) - pp + (v2dt2 * lapl);
			}
		}
	}
}

/**	<<<(NXE,NYE,NZE),(32,32,1)>>> */
__global__ void kgpu_rtmStepMultiWave(
	uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze,
	float *__restrict__ P0, float *__restrict__ PP0,
	float *__restrict__ P1, float *__restrict__ PP1,
	float *__restrict__ V2DT2,
	float *__restrict__ coefs, bool rtm2D)
{
	uint32_t half_order = st_order / 2;
	uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index
	uint32_t iz = blockIdx.z * blockDim.z + threadIdx.z; // Global Z index
	uint32_t kx, ky, kz, io;
	uint32_t cx_base = 0;
	uint32_t cy_base = 0;
	uint32_t cz_base = 0;
	uint32_t offset = 0;
	float acmx_0 = 0.0, acmy_0 = 0.0, acmz_0 = 0.0, lapl_0 = 0.0;
	float acmx_1 = 0.0, acmy_1 = 0.0, acmz_1 = 0.0, lapl_1 = 0.0;

	if (ix >= nxe || iy >= nye || iz >= nze){
		return;
	}
	if (ix >= half_order && ix < (nxe - half_order))
	{
		if (iy >= half_order && iy < (nye - half_order))
		{
			if (iz >= half_order && iz < (nze - half_order))
			{
				cx_base = 0;
				cy_base = (st_order) + 1;
				cz_base = (2 * (st_order + 1));
				offset = ix * nye * nze + iy * nze + iz;
				acmx_0 = 0.0, acmy_0 = 0.0, acmz_0 = 0.0, lapl_0 = 0.0;
				acmx_1 = 0.0, acmy_1 = 0.0, acmz_1 = 0.0, lapl_1 = 0.0;
				#pragma unroll (8)
				for (io = 0; io <= st_order; io++)
				{
					kx = ((ix - half_order + io) * nye * nze) + iy * nze + iz;
					ky = (ix * nye * nze) + (iy - half_order + io) * nze + iz;
					kz = ix * nye * nze + iy * nze + (iz - half_order + io);
					acmx_0 += P0[kx] * coefs[cx_base + io];
					acmy_0 += P0[ky] * coefs[cy_base + io];
					acmz_0 += P0[kz] * coefs[cz_base + io];

					acmx_1 += P1[kx] * coefs[cx_base + io];
					acmy_1 += P1[ky] * coefs[cy_base + io];
					acmz_1 += P1[kz] * coefs[cz_base + io];
				}
				if (rtm2D)
				{
					acmx_0 = 0.0;
					acmx_1 = 0.0;
				}
				lapl_0 = acmz_0 + acmy_0 + acmx_0;
				lapl_1 = acmz_1 + acmy_1 + acmx_1;
				PP0[offset] = (2 * P0[offset]) - PP0[offset] + (V2DT2[offset] * lapl_0);
				PP1[offset] = (2 * P1[offset]) - PP1[offset] + (V2DT2[offset] * lapl_1);
			}
		}
	}
}



/**	<<<(NXE,NYE,NZE),(32,32,1)>>> */
__global__ void kgpu_taperAllBorders(uint32_t pStartX, uint32_t pEndX, uint32_t pStartY, uint32_t pEndY,
									 uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen,
									 float *__restrict__ taper, float *__restrict__ P)
{
	// local grid coordinates
	uint32_t lx = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t ly = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index
	uint32_t lz = blockIdx.z * blockDim.z + threadIdx.z; // Global Z index

	// global grid coordinates
	uint32_t ix = lx + pStartX;
	uint32_t iy = ly + pStartY;
	uint32_t iz = lz;

	// local grid dimensions
	uint32_t gxe = pEndX - pStartX;
	uint32_t gye = pEndY - pStartY;
	uint32_t gze = nze;

	// check if thread is within grid limits
	if (lx >= gxe || ly >= gye || lz >= gze)
		return;

	//printf("l(%d,%d,%d); i(%d,%d,%d), blen=%d n(%d,%d,%d)\n",lx, ly, lz, ix, iy, iz, blen, nxe, nye, nze);
	uint32_t poffset = (lx * gye * gze) + (ly * gze) + lz;
	float val0 = P[poffset];
	/********************************************************/
	if (iy < blen)
	{ //left
		val0 *= taper[blen - 1 - iy];
	}
	else if (iy >= (nye - blen))
	{ // right
		val0 *= taper[iy - (nye - blen)];
	}
	if (ix < blen)
	{ // front
		val0 *= taper[blen - 1 - ix];
	}
	else if (ix >= (nxe - blen))
	{ // back
		val0 *= taper[ix - (nxe - blen)];
	}
	if (iz < blen)
	{ // top
		val0 *= taper[blen - 1 - iz];
	}
	else if (iz > (nze - blen))
	{ // bottom
		val0 *= taper[iz - ((nze - blen))];
	}
	/********************************************************/
	P[poffset] = val0;
}

/**	<<<(NXE,NYE,BLEN),(10,10,10)>>> */
__global__ void kgpu_taperUpperBorders(
	uint32_t pStartX, uint32_t pEndX, uint32_t pStartY, uint32_t pEndY,
	uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen,
	float *__restrict__ taper, float *__restrict__ P)
{
	// local grid coordinates
	uint32_t lx = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t ly = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index
	uint32_t lz = blockIdx.z * blockDim.z + threadIdx.z; // Global Z index

	// global grid coordinates
	uint32_t ix = lx + pStartX;
	uint32_t iy = ly + pStartY;
	uint32_t iz = lz;

	// local grid dimensions
	uint32_t gxe = pEndX - pStartX;
	uint32_t gye = pEndY - pStartY;
	uint32_t gze = nze;

	if (lx >= gxe || ly >= gye || lz >= gze)
		return;

	uint32_t poffset = (lx * gye * gze) + (ly * gze) + lz;
	float val0 = P[poffset];
	/********************************************************/
	if (iy < blen)
	{ //left
		val0 *= taper[blen - 1 - iy];
	}
	else if (iy >= (nye - blen))
	{ // right
		val0 *= taper[iy - (nye - blen)];
	}
	if (ix < blen)
	{ // front
		val0 *= taper[blen - 1 - ix];
	}
	else if (ix >= (nxe - blen))
	{ // back
		val0 *= taper[ix - (nxe - blen)];
	}
	if (iz < blen)
	{ // top
		val0 *= taper[blen - 1 - iz];
	}
	/********************************************************/

	P[poffset] = val0;
}

__global__ void kgpu_applySrceEnergy(uint32_t sx, uint32_t sy, uint32_t sz, float eval,
									 uint32_t pStartX, uint32_t pEndX,
									 uint32_t pStartY, uint32_t pEndY,
									 uint32_t nze, float *__restrict__ PP)
{
	uint32_t nxe = pEndX - pStartX;
	uint32_t nye = pEndY - pStartY;
	// check if it is inside the process limits
	if (sx >= pStartX && sx < pEndX)
	{
		if (sy >= pStartY && sy < pEndY)
		{
			uint32_t poffset = (sx - pStartX) * nye * nze + (sy - pStartY) * nze + sz;
			//printf(">> kSrce: P(%d,%d,%d) <= %.10f + %.10f (P=%p)\n", sx, sy,sz, PP[poffset], eval, PP);
			PP[poffset] += eval;
		}
	}
}

/** <<< (NX,NY,1), (32,32) >>> */
__global__ void kgpu_rwRcvrEnergy(uint32_t rcvOffsetX, uint32_t rcvDistX, uint32_t rcvCountX,
								  uint32_t rcvOffsetY, uint32_t rcvDistY, uint32_t rcvCountY, uint32_t rcvDepthZ,
								  uint32_t pStartX, uint32_t pEndX,
								  uint32_t pStartY, uint32_t pEndY,
								  uint32_t border_length, uint32_t nze, uint32_t nt, uint32_t it, bool modeling,
								  float *__restrict__ SEISM, float *__restrict__ PPR)
{
	// coordinates in the receiver grid matrix
	uint32_t rgx = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t rgy = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index

	if (rgx >= rcvCountX || rgy >= rcvCountY)
		return;

	// cordinates in the pressure grid global grid
	uint32_t rz = rcvDepthZ;
	uint32_t ry = rgy * rcvDistY + rcvOffsetY + border_length;
	uint32_t rx = rgx * rcvDistX + rcvOffsetX + border_length;
	uint32_t nxe = pEndX - pStartX;
	uint32_t nye = pEndY - pStartY;

	// check if it is inside the process limits
	if (rx >= pStartX && rx < pEndX)
	{
		if (ry >= pStartY && ry < pEndY)
		{
			uint32_t poffset = (rx - pStartX) * nye * nze + (ry - pStartY) * nze + rz;
			uint32_t soffset = (rgx * rcvCountY * nt) + (rgy * nt) + it;
			if (modeling)
			{
				SEISM[soffset] = PPR[poffset];
			}
			else
			{
				PPR[poffset] += SEISM[soffset];
			}
		}
	}
}

/**	<<<(NXE,NYE,1),(32,32,1)>>> */
__global__ void kgpu_rwUPB(uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze,
						   uint32_t blen, uint32_t it, bool rw,
						   float *__restrict__ PP, float *__restrict__ UPB)
{
	uint32_t h_order = st_order / 2;
	uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index

	if (ix < nxe)
	{
		if (iy < nye)
		{
			for (uint32_t iz = blen - h_order; iz < blen; iz++)
			{
				uint32_t uoffset = it * nxe * nye * h_order + ix * nye * h_order + iy * h_order + (iz - (blen - h_order));
				uint32_t poffset = ix * nye * nze + iy * nze + iz;
				if (rw)
				{
					UPB[uoffset] = PP[poffset];
				}
				else
				{
					PP[poffset] = UPB[uoffset];
				}
			}
		}
	}
}

/**	<<<(NXE,NYE,NZE),(32,32,1)>>> */
__global__ void kernel_imgcondition(uint32_t nxe, uint32_t nye, uint32_t nze,
									float *__restrict__ IMG, float *__restrict__ PS, float *__restrict__ PR)
{
	uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x; // Global X index
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y; // Global Y index
	uint32_t iz = blockIdx.z * blockDim.z + threadIdx.z; // Global Z index
	uint32_t offset = ix * nye * nze + iy * nze + iz;

	if (ix < nxe)
	{
		if (iy < nye)
		{
			if (iz < nze)
			{
				IMG[offset] += PS[offset] * PR[offset];
			}
		}
	}
}

void RTMSTEP_WRAPPER(dim3 dimGrid, dim3 dimBlock,
					 uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen,
					 float *__restrict__ P, float *__restrict__ PP, float *__restrict__ V2DT2,
					 float *__restrict__ coefs, bool rtm2D)
{
	if(rtm2D){
		kgpu_rtmStep2D<<<dimGrid, dimBlock>>>(st_order, nxe, nye, nze, blen, P, PP, V2DT2, coefs, rtm2D);
	}else{
		kgpu_rtmStep<<<dimGrid, dimBlock>>>(st_order, nxe, nye, nze, P, PP, V2DT2, coefs, rtm2D);
	}
}

void RTMSTEP_MULTIWAVE_WRAPPER(dim3 dimGrid, dim3 dimBlock,
	uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze,
	float *__restrict__ P0, float *__restrict__ PP0,
	float *__restrict__ P1, float *__restrict__ PP1, 
	float *__restrict__ V2DT2,
	float *__restrict__ coefs, bool rtm2D)
{
	kgpu_rtmStepMultiWave<<<dimGrid, dimBlock>>>(st_order, nxe, nye, nze, P0, PP0, 
		P1, PP1, V2DT2, coefs, rtm2D);
}

void RTMSEISM_WRAPPER(dim3 dimGrid, dim3 dimBlock,
					  uint32_t rcvOffsetX, uint32_t rcvDistX, int rcvCountX,
					  uint32_t rcvOffsetY, uint32_t rcvDistY, int rcvCountY, uint32_t rcvDepthZ,
					  uint32_t pStartX, uint32_t pEndX,
					  uint32_t pStartY, uint32_t pEndY,
					  uint32_t border_length, uint32_t nze, uint32_t nt, uint32_t it, bool modeling,
					  float *__restrict__ SEISM, float *__restrict__ PPR)
{

	kgpu_rwRcvrEnergy<<<dimGrid, dimBlock>>>(rcvOffsetX, rcvDistX, rcvCountX,
											 rcvOffsetY, rcvDistY, rcvCountY, rcvDepthZ,
											 pStartX, pEndX,
											 pStartY, pEndY,
											 border_length, nze, nt, it, modeling,
											 SEISM, PPR);
}

void RTMSOURCE_WRAPPER(dim3 dimGrid, dim3 dimBlock,
					   uint32_t sx, uint32_t sy, uint32_t sz, float eval,
					   uint32_t pStartX, uint32_t pEndX,
					   uint32_t pStartY, uint32_t pEndY,
					   uint32_t nze, float *__restrict__ PP)
{
	kgpu_applySrceEnergy<<<dimGrid, dimBlock>>>(sx, sy, sz, eval,
												pStartX, pEndX,
												pStartY, pEndY,
												nze, PP);
}

void RTMTAPER_WRAPPER(dim3 dimGrid, dim3 dimBlock,
					  uint32_t pStartX, uint32_t pEndX, uint32_t pStartY, uint32_t pEndY,
					  uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen,
					  float *__restrict__ taper, float *__restrict__ P, bool upperOnly)
{
	if (upperOnly)
	{
		kgpu_taperUpperBorders<<<dimGrid, dimBlock>>>(
			pStartX, pEndX, pStartY, pEndY,
			nxe, nye, nze, blen,
			taper, P);
	}
	else
	{
		kgpu_taperAllBorders<<<dimGrid, dimBlock>>>(pStartX, pEndX, pStartY, pEndY,
													nxe, nye, nze, blen,
													taper, P);
	}
}

void RTMWUPB_WRAPPER(dim3 dimGrid, dim3 dimBlock, uint32_t st_order, uint32_t nxe, uint32_t nye, uint32_t nze, uint32_t blen, uint32_t it, bool rw,
					 float *__restrict__ PP, float *__restrict__ UPB)
{

	kgpu_rwUPB<<<dimGrid, dimBlock>>>(st_order, nxe, nye,
									  nze, blen, it, rw,
									  PP, UPB);
}

void RTMIMG_WRAPPER(dim3 dimGrid, dim3 dimBlock,
					uint32_t nxe, uint32_t nye, uint32_t nze,
					float *__restrict__ IMG, float *__restrict__ PS, float *__restrict__ PR)
{
	kernel_imgcondition<<<dimGrid, dimBlock>>>(nxe, nye, nze,
											   IMG, PS, PR);
}


void RTM_FREQIMG_WRAPPER(dim3 dimGrid, dim3 dimBlock, 
		uint64_t iw,uint64_t lw,
		uint64_t gStartX, uint64_t gEndX, 
		uint64_t gStartY, uint64_t gEndY, 
		uint64_t gStartZ, uint64_t gEndZ,
		uint64_t nxe, uint64_t nye, uint64_t nze, uint64_t blen,   
		float *__restrict__ wList,
		float *__restrict__ IMG,
		float *__restrict__ PSRe, float *__restrict__ PSIm,
		float *__restrict__ PRRe, float *__restrict__ PRIm)
{

	kgpu_freqDomainImageCondition<<<dimGrid, dimBlock>>>(iw,lw,
		gStartX, gEndX, 
		gStartY, gEndY, 
		gStartZ, gEndZ,
		nxe, nye, nze, blen,   
		wList,IMG,PSRe, PSIm,PRRe, PRIm);
}

void RTM_UPDATEFREQ_WRAPPER(dim3 dimGrid, dim3 dimBlock, 
		uint64_t it, uint64_t iw,uint64_t lw, uint64_t nt,
		uint64_t gStartX, uint64_t gEndX, uint64_t gStartY, uint64_t gEndY, 
		uint64_t gStartZ, uint64_t gEndZ,
		uint64_t nxe, uint64_t nye, uint64_t nze, uint64_t blen,
		float *__restrict__ kernelRe, float *__restrict__ kernelIm,
		float *__restrict__ PS, float *__restrict__ PR,
		float *__restrict__ PSRe, float *__restrict__ PSIm,
		float *__restrict__ PRRe, float *__restrict__ PRIm)
{
	kgpu_updateFreqContributions<<<dimGrid, dimBlock>>>(it,iw,lw,nt,
		gStartX, gEndX, gStartY, gEndY,gStartZ, gEndZ, 
		nxe, nye, nze, blen, 
		kernelRe, kernelIm,PS, PR,PSRe, PSIm,PRRe, PRIm);
}

#endif // RTM_ACC_GPU