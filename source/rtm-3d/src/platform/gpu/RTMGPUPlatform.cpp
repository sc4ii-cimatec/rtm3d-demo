#include <assert.h>
#include <iostream>
#include <fstream>
#include <Misc.hpp>
#include <RTM.hpp>
#include <RTMGrid.hpp>
#include <RTMGPU.hpp>
#include <RTMGPUPlatform.hpp>

#ifdef RTM_ACC_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#include <rtmgpu.hpp>
#endif

void RTMGPUPlatform::destroyRTMPlatform()
{
#ifdef RTM_ACC_GPU
    CUDACHECK(cudaDeviceReset());
#endif
}
void RTMGPUPlatform::initRTMPlatform()
{
#ifdef RTM_ACC_GPU
    RTM_PRINT("Initializing CUDA GPU Platform...", rtmParam->verbose);
    CUDACHECK(cudaGetDeviceCount(&nDevices));
    deviceID=pLimits->lRank%nDevices;
    // each process uses a different GPU when
    // multiple devices are available
    CUDACHECK(cudaSetDevice(deviceID));
    CUDACHECK(cudaDeviceReset());
    CUDACHECK(cudaGetDeviceProperties(&deviceProperties,deviceID));
    // printf(">>> [P%d] GPU Properties: \n",  pLimits->pRank);
    // printf(">>> [P%d] nDevices= %d \n",  pLimits->pRank, nDevices);
    // printf(">>> [P%d] devId   = %d (lRank=%d) \n",  pLimits->pRank, deviceID, pLimits->pRank);
    //printf(">>> \t\t + name                    : %s \n", deviceProperties.name);
    //printf(">>> \t\t + totalGlobalMem          : %d GB \n",deviceProperties.totalGlobalMem/1000000000);
    //printf(">>> \t\t + sharedMemPerBlock       : %d \n", deviceProperties.sharedMemPerBlock);
    //printf(">>> \t\t + regsPerBlock            : %d \n", deviceProperties.regsPerBlock);
    //printf(">>> \t\t + warpSize                : %d \n", deviceProperties.warpSize);
    //printf(">>> \t\t + memPitch                : %d \n", deviceProperties.memPitch);
    //printf(">>> \t\t + maxThreadsPerBlock      : %d \n", deviceProperties.maxThreadsPerBlock);
    //printf(">>> \t\t + maxThreadsDim           : (%d,%d,%d) \n", deviceProperties.maxThreadsDim[0], deviceProperties.maxThreadsDim[1], deviceProperties.maxThreadsDim[2]);
    //printf(">>> \t\t + maxGridSize             : (%d,%d,%d) \n", deviceProperties.maxGridSize[0], deviceProperties.maxGridSize[1], deviceProperties.maxGridSize[2]);
    //printf(">>> \t\t + totalConstMem           : %d \n", deviceProperties.totalConstMem);
    //printf(">>> \t\t + major                   : %d \n", deviceProperties.major);
    //printf(">>> \t\t + minor                   : %d \n", deviceProperties.minor);
    //printf(">>> \t\t + clockRate               : %d \n", deviceProperties.clockRate);
    //printf(">>> \t\t + deviceOverlap           : %d \n", deviceProperties.deviceOverlap);
    //printf(">>> \t\t + multiProcessorCount     : %d \n", deviceProperties.multiProcessorCount);
    //printf(">>> \t\t + kernelExecTimeoutEnabled: %d \n", deviceProperties.kernelExecTimeoutEnabled);
    //printf(">>> \t\t + integrated              : %d \n", deviceProperties.integrated);
    //printf(">>> \t\t + canMapHostMemory        : %d \n", deviceProperties.canMapHostMemory);
    //printf(">>> \t\t + computeMode             : %d \n", deviceProperties.computeMode);
    
#endif
}

void RTMGPUPlatform::rtmStepMultipleWave(RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
                         RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid, 
                         RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
                         const RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid)
{
    grtmStepMultiWave(P0Grid->getDevPtr(), PP0Grid->getDevPtr(),P1Grid->getDevPtr(), PP1Grid->getDevPtr(),  
    stencil->getDevPtr(), v2dt2Grid.getDevPtr());
}

void RTMGPUPlatform::rtmUpdateFreqContributions(int it, int iw, int lw, 
                              RTMCube<RTMData_t, RTMDevPtr_t> *PSGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid,
                              RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSImGrid,
                              RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRImGrid,
                              RTMPlane<RTMData_t,RTMDevPtr_t> * kernelRe, RTMPlane<RTMData_t,RTMDevPtr_t> * kernelIm)
{
    grtmUpdateFreqContribution((uint64_t) it, (uint64_t) iw, (uint64_t)lw,
        PSReGrid->getStartX(), PSReGrid->getEndX(),
        PSReGrid->getStartY(), PSReGrid->getEndY(),
        PSReGrid->getStartZ(), PSReGrid->getEndZ(),
		kernelRe->getDevPtr(), kernelIm->getDevPtr(),
		PSGrid->getDevPtr(), PRGrid->getDevPtr(),
		PSReGrid->getDevPtr(), PSImGrid->getDevPtr(),
		PRReGrid->getDevPtr(), PRImGrid->getDevPtr());
}

void RTMGPUPlatform::rtmFreqDomainImageCondition(int iw, int lw, RTMVector<RTMData_t,RTMDevPtr_t> * w2List, 
                                    RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                    RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSImGrid,
                                    RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRImGrid)
{
    grtmFreqImgCondition((uint64_t)iw, (uint64_t)lw,
        PSReGrid->getStartX(), PSReGrid->getEndX(),
        PSReGrid->getStartY(), PSReGrid->getEndY(),
        PSReGrid->getStartZ(), PSReGrid->getEndZ(),
        w2List->getDevPtr(),imgGrid->getDevPtr(),
		PSReGrid->getDevPtr(), PSImGrid->getDevPtr(),
		PRReGrid->getDevPtr(), PRImGrid->getDevPtr());
}

void RTMGPUPlatform::rtmApplySource(RTMCube<RTMData_t, RTMDevPtr_t> * PGrid,
                                     RTMSeismicSource<RTMData_t,RTMDevPtr_t> *srcGrid, 
                                     uint32_t it)
{
    /* source position does not consider the extended grid */
    int sxe = srcGrid->getX() + rtmParam->blen;
    int sye = srcGrid->getY() + rtmParam->blen;
    int sze = srcGrid->getZ() + rtmParam->blen;
    RTMData_t srcVal = (*srcGrid)[it];
    grtmSource(sxe, sye, sze, srcVal, PGrid->getDevPtr());
}
void RTMGPUPlatform::rtmRestoreReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> * PPGrid, RTMReceiverGrid<RTMData_t,RTMDevPtr_t> *rcvGrid, uint32_t it)
{
    grtmSeism(PPGrid->getDevPtr(), rcvGrid->getDevPtr(), rcvGrid->getNZ(), it, false);
}

void RTMGPUPlatform::rtmSaveReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> * PPGrid, RTMReceiverGrid<RTMData_t,RTMDevPtr_t> *rcvGrid, uint32_t it)
{
    grtmSeism(PPGrid->getDevPtr(), rcvGrid->getDevPtr(), rcvGrid->getNZ(), it, true);
}

void RTMGPUPlatform::rtmSaveUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, 
                                        RTMGridCollection<RTMData_t,RTMDevPtr_t> *upbGrid,uint32_t it)
{
    grtmUPB(it, PGrid->getDevPtr(), upbGrid->getDevPtr(), true);
}
void RTMGPUPlatform::rtmRestoreUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, 
                                            RTMGridCollection<RTMData_t,RTMDevPtr_t> *upbGrid,uint32_t it)
{
    grtmUPB(it, PGrid->getDevPtr(), upbGrid->getDevPtr(), false);
}

void RTMGPUPlatform::rtmStep(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid, 
                            RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
                             const RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid)
{
    grtmStep(PGrid->getDevPtr(), PPGrid->getDevPtr(),  stencil->getDevPtr(), v2dt2Grid.getDevPtr());
}

/**
 * Cross-correlation image condition
 * OBS: IMG, PS and PR grids must have the same dimensions
 * */
void RTMGPUPlatform::rtmImageCondition(RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                       RTMCube<RTMData_t, RTMDevPtr_t> *PSGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid)
{
    int ix, iy, iz;
    int nx = imgGrid->getNX();
    int ny = imgGrid->getNY();
    int nz = imgGrid->getNZ();

    assert((imgGrid->getNX() == PSGrid->getNX()) && (imgGrid->getNX() == PRGrid->getNX()));
    assert((imgGrid->getNY() == PSGrid->getNY()) && (imgGrid->getNY() == PRGrid->getNY()));
    assert((imgGrid->getNZ() == PSGrid->getNZ()) && (imgGrid->getNZ() == PRGrid->getNZ()));

    grtmImgCondition(imgGrid->getDevPtr(), PSGrid->getDevPtr(), PRGrid->getDevPtr());
}

void RTMGPUPlatform::rtmTaperAllBorders(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid, RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper)
{
    grtmTaperBorders(rtmGrid->getDevPtr(), rtmTaper->getDevPtr(), false);
}

void RTMGPUPlatform::rtmTaperUpperBorders(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid, RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper)
{
    grtmTaperBorders(rtmGrid->getDevPtr(), rtmTaper->getDevPtr(), true);
}


void RTMGPUPlatform::grtmStep(RTMDevPtr_t * devP, RTMDevPtr_t * devPP,  RTMDevPtr_t * devCoefs, RTMDevPtr_t * devV2DT2)
{
#ifdef RTM_ACC_GPU
    uint32_t gridx;
    uint32_t gridy;
    uint32_t gridz;
    bool rtm2D = rtmParam->nx==1;

    if (rtm2D){
        gridx = CUDANGRIDS(plen_y, SMALL_BLOCK_SIZE);
        gridy = CUDANGRIDS(plen_z, SMALL_BLOCK_SIZE);
        gridz = 1;
    }else{
        gridx = CUDANGRIDS(plen_x, SMALL_BLOCK_SIZE);
        gridy = CUDANGRIDS(plen_y, SMALL_BLOCK_SIZE);
        gridz = CUDANGRIDS(plen_z, SMALL_BLOCK_SIZE);
    }
	dim3 dimGrid(gridx, gridy, gridz);	
	dim3 dimBlock(SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE);

    // printf(">> PLEN_X=%d; GridNXE = %d \n", plen_x, gridx);
    // printf(">> PLEN_Y=%d; GridNYE = %d \n", plen_y, gridy);
    // printf(">> PLEN_Z=%d; GridNZE = %d \n", plen_z, gridz);
	
    RTMSTEP_WRAPPER(dimGrid, dimBlock, st_order, plen_x, plen_y, plen_z, rtmParam->blen,
    devP, devPP, devV2DT2, devCoefs, rtm2D);
	CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}

void RTMGPUPlatform::grtmStepMultiWave(
RTMDevPtr_t * devP0, RTMDevPtr_t * devPP0, 
RTMDevPtr_t * devP1, RTMDevPtr_t * devPP1,  
RTMDevPtr_t * devCoefs, RTMDevPtr_t * devV2DT2)
{
#ifdef RTM_ACC_GPU
    uint32_t gridx = CUDANGRIDS(plen_x, SMALL_BLOCK_SIZE);
    uint32_t gridy = CUDANGRIDS(plen_y, SMALL_BLOCK_SIZE);
    uint32_t gridz = CUDANGRIDS(plen_z, SMALL_BLOCK_SIZE);
    bool rtm2D = rtmParam->nx==1;
	dim3 dimGrid(gridx, gridy, gridz);	
	dim3 dimBlock(SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE);

    // printf(">> PLEN_X=%d; GridNXE = %d \n", plen_x, gridx);
    // printf(">> PLEN_Y=%d; GridNYE = %d \n", plen_y, gridy);
    // printf(">> PLEN_Z=%d; GridNZE = %d \n", plen_z, gridz);
	
    RTMSTEP_MULTIWAVE_WRAPPER(dimGrid, dimBlock, st_order, plen_x, plen_y, plen_z, 
    devP0, devPP0, devP1, devPP1, devV2DT2, devCoefs, rtm2D);
	CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}



void RTMGPUPlatform::grtmSeism(RTMDevPtr_t * devPPR, RTMDevPtr_t * devSeism, int _nt, int _it, bool modeling){
#ifdef RTM_ACC_GPU
    uint32_t rcvOffsetX = rtmParam->receiver_start_x; 
    uint32_t rcvDistX = rtmParam->receiver_distance_x; 
    uint32_t rcvCountX = rtmParam->receiver_count_x; 
	uint32_t rcvOffsetY = rtmParam->receiver_start_y; 
    uint32_t rcvDistY = rtmParam->receiver_distance_y; 
    uint32_t rcvCountY = rtmParam->receiver_count_y; 
    uint32_t rcvDepthZ = rtmParam->receiver_depth_z + rtmParam->blen;
	uint32_t pStartX = pLimits->processArea.xStart; 
    uint32_t pEndX = pLimits->processArea.xEnd;
	uint32_t pStartY = pLimits->processArea.yStart; 
    uint32_t pEndY = pLimits->processArea.yEnd;
    uint32_t blen = rtmParam->blen;
	uint32_t nze = rtmParam->nz + 2*rtmParam->blen; 
    uint32_t nt = _nt; 
    uint32_t it = _it; 
    uint32_t gNXE = CUDANGRIDS(rtmParam->receiver_count_x, DEFAULT_BLOCK_SIZE);
    uint32_t gNYE = CUDANGRIDS(rtmParam->receiver_count_y, DEFAULT_BLOCK_SIZE);
	dim3 dimGrid(gNXE, gNYE, 1);	
	dim3 dimBlock(DEFAULT_BLOCK_SIZE, DEFAULT_BLOCK_SIZE, 1);

    RTMSEISM_WRAPPER(dimGrid, dimBlock, 
	rcvOffsetX, rcvDistX, rcvCountX, 
    rcvOffsetY, rcvDistY, rcvCountY, rcvDepthZ,
	pStartX,    pEndX,    pStartY,   pEndY,
	blen, nze, nt, it, 
    modeling,
    devSeism, devPPR);
	CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}

void RTMGPUPlatform::grtmSource(uint32_t sx, uint32_t sy, uint32_t sz, RTMData_t eval, RTMDevPtr_t * devPP)
{
#ifdef RTM_ACC_GPU
    dim3 dimGrid(1, 1, 1);	
	dim3 dimBlock(1, 1, 1);
    uint32_t pStartX = pLimits->processArea.xStart; 
    uint32_t pEndX = pLimits->processArea.xEnd;
	uint32_t pStartY = pLimits->processArea.yStart; 
    uint32_t pEndY = pLimits->processArea.yEnd;
	uint32_t nze = rtmParam->nz + 2*rtmParam->blen; 
    
    if (sx >= pStartX && sx < pEndX){
		if (sy >= pStartY && sy < pEndY){
            RTMSOURCE_WRAPPER(dimGrid, dimBlock,
            sx, sy, sz, eval, 
            pStartX, pEndX,
            pStartY, pEndY,
            nze,devPP);
            CUDACHECK( cudaPeekAtLastError() );
	        CUDACHECK( cudaDeviceSynchronize() );
        }
    }
#endif
}

void RTMGPUPlatform::grtmTaperBorders( RTMDevPtr_t *  P, RTMDevPtr_t * TAPER, bool upperBorderOnly)
{
#ifdef RTM_ACC_GPU
    uint32_t gridx = CUDANGRIDS(plen_x, SMALL_BLOCK_SIZE);
    uint32_t gridy = CUDANGRIDS(plen_y, SMALL_BLOCK_SIZE);
    uint32_t gridz;
    if(upperBorderOnly){
        gridz = CUDANGRIDS(rtmParam->blen, SMALL_BLOCK_SIZE);
    }else{
        gridz = CUDANGRIDS(plen_z, SMALL_BLOCK_SIZE);
    }
    dim3 dimGrid(gridx, gridy, gridz);	
	dim3 dimBlock(SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE);
    
    uint32_t pStartX = pLimits->processArea.xStart; 
    uint32_t pEndX = pLimits->processArea.xEnd;
	uint32_t pStartY = pLimits->processArea.yStart; 
    uint32_t pEndY = pLimits->processArea.yEnd;
    uint32_t nxe = rtmParam->nx + 2*rtmParam->blen;
    uint32_t nye = rtmParam->ny + 2*rtmParam->blen;
	uint32_t nze = rtmParam->nz + 2*rtmParam->blen;
    uint32_t blen = rtmParam->blen;
    
    RTMTAPER_WRAPPER(dimGrid, dimBlock,
                    pStartX, pEndX,pStartY, pEndY, 
                    nxe, nye, nze, blen, 
                    TAPER, P, upperBorderOnly);
    CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}

void RTMGPUPlatform::grtmSwapPtr(RTMDevPtr_t ** devA, RTMDevPtr_t ** devB)
{
#ifdef RTM_ACC_GPU
    // RTMDevPtr_t * DEV_SWAP_PTR = *devA;
    // *devA = *devB;
    // *devB = DEV_SWAP_PTR;
#endif
}

void RTMGPUPlatform::grtmImgCondition(RTMDevPtr_t * IMG, RTMDevPtr_t * PS, RTMDevPtr_t * PR)
{
#ifdef RTM_ACC_GPU
    uint32_t gridx = CUDANGRIDS(plen_x, SMALL_BLOCK_SIZE);
    uint32_t gridy = CUDANGRIDS(plen_y, SMALL_BLOCK_SIZE);
    uint32_t gridz = CUDANGRIDS(plen_z, SMALL_BLOCK_SIZE);
	dim3 dimGrid(gridx, gridy, gridz);	
	dim3 dimBlock(SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE);
    RTMIMG_WRAPPER(dimGrid, dimBlock, 
	plen_x, plen_y, plen_z,
	IMG, PS, PR);
    CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}

 void RTMGPUPlatform::grtmUPB(uint32_t it, RTMDevPtr_t * devPP, RTMDevPtr_t * devUPB, bool rw){
#ifdef RTM_ACC_GPU
	uint32_t nze = rtmParam->nz + 2*rtmParam->blen;
    uint32_t gridx = CUDANGRIDS(plen_x, DEFAULT_BLOCK_SIZE);
    uint32_t gridy = CUDANGRIDS(plen_y, DEFAULT_BLOCK_SIZE);
	dim3 dimGrid(gridx, gridy, 1);	
	dim3 dimBlock(DEFAULT_BLOCK_SIZE, DEFAULT_BLOCK_SIZE, 1);

    RTMWUPB_WRAPPER(dimGrid, dimBlock, rtmParam->stencil_order, plen_x, plen_y, nze, 
    rtmParam->blen, it, rw,	devPP,devUPB);
    CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
 }

void RTMGPUPlatform::grtmFreqImgCondition(uint64_t iw,  uint64_t lw,
        uint64_t iStartX, uint64_t iEndX, uint64_t iStartY, uint64_t iEndY, 
        uint64_t iStartZ, uint64_t iEndZ,
        RTMDevPtr_t * w2List,RTMDevPtr_t * IMG,
		RTMDevPtr_t * PSRe, RTMDevPtr_t * PSIm,
		RTMDevPtr_t * PRRe, RTMDevPtr_t * PRIm){
#ifdef RTM_ACC_GPU 
    uint64_t nxe = rtmParam->nx + 2*rtmParam->blen;
    uint64_t nye = rtmParam->ny + 2*rtmParam->blen;
	uint64_t nze = rtmParam->nz + 2*rtmParam->blen;
    uint64_t blen = rtmParam->blen;

    if(iStartX >= rtmParam->nx || iStartY >= rtmParam->ny || iStartZ>=rtmParam->nz){
        return; // invalid area
    }
    // this is necessary to avoid border regions
    uint32_t iNX = iEndX - iStartX;
    uint32_t iNY = iEndY - iStartY;
    uint32_t iNZ = iEndZ - iStartZ;
    uint32_t gridx = CUDANGRIDS(iNX, SMALL_BLOCK_SIZE);
    uint32_t gridy = CUDANGRIDS(iNY, SMALL_BLOCK_SIZE);
    uint32_t gridz = CUDANGRIDS(iNZ, SMALL_BLOCK_SIZE);
	
    // printf("P[%d]: gridx = %d gridy=%d gridz=%d inY=%d\n", pLimits->pRank, gridx, gridy, gridz, iNY);
    dim3 dimGrid(gridx, gridy, gridz);
	dim3 dimBlock(SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE);
    
    RTM_FREQIMG_WRAPPER(dimGrid, dimBlock, iw,lw,
		iStartX, iEndX, 
		iStartY, iEndY, 
		iStartZ, iEndZ,
		nxe, nye, nze, blen,  
		w2List,
		IMG,
		PSRe, PSIm,
		PRRe, PRIm);
    
    CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}

void RTMGPUPlatform::grtmUpdateFreqContribution(uint64_t it, uint64_t iw, uint64_t lw,
        uint64_t iStartX, uint64_t iEndX, uint64_t iStartY, uint64_t iEndY, 
        uint64_t iStartZ, uint64_t iEndZ,
		RTMDevPtr_t * kernelRe, RTMDevPtr_t * kernelIm,
		RTMDevPtr_t * PS, RTMDevPtr_t * PR,
		RTMDevPtr_t * PSRe, RTMDevPtr_t * PSIm,
		RTMDevPtr_t * PRRe, RTMDevPtr_t * PRIm){
#ifdef RTM_ACC_GPU
    uint64_t nxe = rtmParam->nx + 2*rtmParam->blen;
    uint64_t nye = rtmParam->ny + 2*rtmParam->blen;
	uint64_t nze = rtmParam->nz + 2*rtmParam->blen;
    uint64_t blen = rtmParam->blen;
    if(iStartX >= rtmParam->nx || iStartY >= rtmParam->ny || iStartZ>=rtmParam->nz){
        return; // invalid area
    }
    // this is necessary to avoid border regions
    uint32_t iNX = iEndX - iStartX;
    uint32_t iNY = iEndY - iStartY;
    uint32_t iNZ = iEndZ - iStartZ;
    uint32_t gridx = CUDANGRIDS(iNX, SMALL_BLOCK_SIZE);
    uint32_t gridy = CUDANGRIDS(iNY, SMALL_BLOCK_SIZE);
    uint32_t gridz = CUDANGRIDS(iNZ, SMALL_BLOCK_SIZE);
	dim3 dimGrid(gridx, gridy, gridz);
	dim3 dimBlock(SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE, SMALL_BLOCK_SIZE);


    // printf(" pNX = %d pNY = %d pNZ=%d \n pStartX = %d pEndX = %d \n pStartY = %d pEndY = %d \n",
	//  iNX, iNY, iNZ, iStartX, iEndX, iStartY, iEndY);
    RTM_UPDATEFREQ_WRAPPER(dimGrid, dimBlock, 
		it, iw,lw, rtmParam->nt,
		iStartX, iEndX, 
		iStartY, iEndY, 
		iStartZ, iEndZ,
		nxe, nye, nze, blen,
		kernelRe, kernelIm,
		PS, PR,
		PSRe, PSIm,
		PRRe, PRIm);
    CUDACHECK( cudaPeekAtLastError() );
	CUDACHECK( cudaDeviceSynchronize() );
#endif
}