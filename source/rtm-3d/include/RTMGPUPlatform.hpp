#ifndef _RTM_GPU_PLATFORM_H
#define _RTM_GPU_PLATFORM_H

#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMParam.hpp>
#include <RTMGPU.hpp>
#include <RTMGridPartition.hpp>
#include <RTMPlatform.hpp>

#ifdef RTM_ACC_GPU
#include <cuda.h>
#include <cuda_runtime.h>
#include <rtmgpu.hpp>
#endif

using namespace std;
using namespace rtmparam;


/**
 * @brief Class RTMGPUPlatform
 * @see RTMKernel.cpp
 */
class RTMGPUPlatform : public RTMPlatform<RTMData_t, RTMDevPtr_t>
{
protected:
    int nt;
    int ntstep;
    int blen;
    int st_order;
    int plen_x; // process limits length in x
    int plen_y; // process limits length in y
    int plen_z; // process limits length in z
    int nDevices;
    int deviceID;

#ifdef RTM_ACC_GPU
    cudaDeviceProp deviceProperties;
#endif

public: // constructors
    RTMGPUPlatform(RTMParam &_rtmParam, RTMProcessLimits &_pLimits) : 
    RTMPlatform{_rtmParam, _pLimits} 
    {
        nt = _rtmParam.nt;
        ntstep = _rtmParam.ntstep;
        blen = _rtmParam.blen;
        st_order = _rtmParam.stencil_order;
        plen_x = _pLimits.processArea.xEnd-_pLimits.processArea.xStart;
        plen_y = _pLimits.processArea.yEnd-_pLimits.processArea.yStart;
        plen_z = _rtmParam.nz + 2*blen;

    };

public:
    /**************************************************************************************/
    /** This class implements the following interface from base class RTMPlatform for GPU */
    virtual void initRTMPlatform();
    virtual void destroyRTMPlatform();
    virtual void rtmTaperAllBorders(RTMCube<RTMData_t, RTMDevPtr_t> *grid,
                                    RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper);
    virtual void rtmTaperUpperBorders(RTMCube<RTMData_t, RTMDevPtr_t> *grid, RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper);
    virtual void rtmStep(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid, RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
                         const RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid);

    virtual void rtmStepMultipleWave(RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
                         RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid, 
                         RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
                         const RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid);

    virtual void rtmImageCondition(RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                   RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid);
    virtual void rtmApplySource(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMSeismicSource<RTMData_t,RTMDevPtr_t> *srcGrid, uint32_t it);
    virtual void rtmRestoreReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMReceiverGrid<RTMData_t,RTMDevPtr_t> *rcvGrid, uint32_t it);
    virtual void rtmSaveReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMReceiverGrid<RTMData_t,RTMDevPtr_t> *rcvGrid, uint32_t it);
    virtual void rtmSaveUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *upbGrid, uint32_t it);
    virtual void rtmRestoreUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *upbGrid, uint32_t it);
    virtual void rtmUpdateFreqContributions(int it, int iw, int lw, 
                              RTMCube<RTMData_t, RTMDevPtr_t> *PSGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid,
                              RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSImGrid,
                              RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRImGrid,
                              RTMPlane<RTMData_t,RTMDevPtr_t> * kernelRe, RTMPlane<RTMData_t,RTMDevPtr_t> * kernelIm);

     virtual void rtmFreqDomainImageCondition(int iw, int lw, RTMVector<RTMData_t,RTMDevPtr_t> * wList,
                              RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                              RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PSImGrid,
                              RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t,RTMDevPtr_t> *PRImGrid);
    
    /**************************************************************************************/

    static size_t getAvailableMemory(){
#ifdef RTM_ACC_GPU
        size_t AVAILABLE_MEMORY;
        size_t TOTAL_MEMORY;
        cudaMemGetInfo(&AVAILABLE_MEMORY, &TOTAL_MEMORY);
        return AVAILABLE_MEMORY;
#else  
        return 0;
#endif
    }

    static void subMemory(size_t _mem){
#ifdef RTM_ACC_GPU
        // size_t AVAILABLE_MEMORY;
        // size_t TOTAL_MEMORY;
        // cudaMemGetInfo(&AVAILABLE_MEMORY, &TOTAL_MEMORY);
        // if (AVAILABLE_MEMORY > _mem){
        //     AVAILABLE_MEMORY-=_mem;
        // }else{
        //     char msg[256];
        //     sprintf(msg, "[-GPU Error: requested mem=%lu(MB) available=%lu (MB)]",
        //     _mem/1000000, AVAILABLE_MEMORY/1000000);
        //     string s(msg);
        //     RTMException ex(s);
        //     throw ex;
        // }
#endif
    }
    static void addMemory(size_t _mem){
#ifdef RTM_ACC_GPU
        // size_t AVAILABLE_MEMORY;
        // size_t TOTAL_MEMORY;
        // cudaMemGetInfo(&AVAILABLE_MEMORY, &TOTAL_MEMORY);
        // if ((AVAILABLE_MEMORY+_mem) <= TOTAL_MEMORY){
        //     AVAILABLE_MEMORY+=_mem;
        // }else{
        //     AVAILABLE_MEMORY=TOTAL_MEMORY;
        // }
#endif
    }
private:
    /** kernel wrappers **/
    void grtmStep(RTMDevPtr_t * devP, RTMDevPtr_t * devPP, RTMDevPtr_t * coefs, RTMDevPtr_t * devV2DT2);
    void grtmStepMultiWave(RTMDevPtr_t * devP0, RTMDevPtr_t * devPP0, 
                           RTMDevPtr_t * devP1, RTMDevPtr_t * devPP1, RTMDevPtr_t * coefs, RTMDevPtr_t * devV2DT2);
    void grtmSeism(RTMDevPtr_t * devPPR, RTMDevPtr_t * devSeism, int _nt, int _it, bool modeling);
    void grtmSource(uint32_t sx, uint32_t sy, uint32_t sz, RTMData_t eval, RTMDevPtr_t * devPP);
    void grtmTaperBorders( RTMDevPtr_t *  P, RTMDevPtr_t * TAPER, bool upperBorderOnly=false);
    void grtmSwapPtr(RTMDevPtr_t ** devA, RTMDevPtr_t ** devB);
    void grtmImgCondition(RTMDevPtr_t * IMG, RTMDevPtr_t * PS, RTMDevPtr_t * PR);
    void grtmUPB(uint32_t it, RTMDevPtr_t * devPP, RTMDevPtr_t * devUPB, bool rw);
    void grtmFreqImgCondition(uint64_t iw,  uint64_t lw,
        uint64_t iStartX, uint64_t iEndX, uint64_t iStartY, uint64_t iEndY, 
        uint64_t iStartZ, uint64_t iEndZ,
        RTMDevPtr_t * w2List,RTMDevPtr_t * IMG,
		RTMDevPtr_t * PSRe, RTMDevPtr_t * PSIm,
		RTMDevPtr_t * PRRe, RTMDevPtr_t * PRIm);
    void grtmUpdateFreqContribution(uint64_t it, uint64_t iw, uint64_t lw,
        uint64_t iStartX, uint64_t iEndX, uint64_t iStartY, uint64_t iEndY, 
        uint64_t iStartZ, uint64_t iEndZ,
		RTMDevPtr_t * kernelRe, RTMDevPtr_t * kernelIm,
		RTMDevPtr_t * PS, RTMDevPtr_t * PR,
		RTMDevPtr_t * PSRe, RTMDevPtr_t * PSIm,
		RTMDevPtr_t * PRRe, RTMDevPtr_t * PRIm);
};

#endif