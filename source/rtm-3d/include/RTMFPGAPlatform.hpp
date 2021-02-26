#ifndef _RTM_FPGA_PLATFORM_H
#define _RTM_FPGA_PLATFORM_H

#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMParam.hpp>
#include <RTMFPGA.hpp>
#include <RTMGridPartition.hpp>
#include <RTMPlatform.hpp>
#include <RTMCPUPlatform.hpp>

#ifdef  RTM_ACC_FPGA
#include "xcl2/xcl2.hpp"
#endif

using namespace std;
using namespace rtmparam;


/**
 * @brief Class RTMFPGAPlatform
 * @see RTMKernel.cpp
 */
class RTMFPGAPlatform : public RTMCPUPlatform
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

    RTMFPGADevice * fpgaDevice = nullptr;
#ifdef  RTM_ACC_FPGA
    typedef ForwardKernel<RTMData_t, RTM_FPGA_stOrder, RTM_FPGA_nPEZ, RTM_FPGA_nPEX> RTMForwardKernel;
    typedef WideData<RTMData_t, RTM_FPGA_nPEX> RTM_wideTypeX;
    typedef WideData<RTM_wideTypeX, RTM_FPGA_nPEZ> RTM_wideType;
#else
    typedef int RTMForwardKernel;
    typedef int RTM_wideTypeX;
    typedef int RTM_wideType;
#endif
    RTMForwardKernel * fwdKernel;

public: // constructors
    RTMFPGAPlatform(RTMParam &_rtmParam, RTMProcessLimits &_pLimits) : 
    RTMCPUPlatform{_rtmParam, _pLimits} 
    {
        nt = _rtmParam.nt;
        ntstep = _rtmParam.ntstep;
        blen = _rtmParam.blen;
        st_order = _rtmParam.stencil_order;
        plen_x = _pLimits.processArea.xEnd-_pLimits.processArea.xStart;
        plen_y = _pLimits.processArea.yEnd-_pLimits.processArea.yStart;
        plen_z = _rtmParam.nz + 2*blen;
    };

    virtual void initRTMPlatform();
    virtual void destroyRTMPlatform();

public:
    /* fpga specific functions */
    void rtmSeismicModeling(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> * shotDescriptor,
        RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
        RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper,
        RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid);
    
    void rtmForwardPropagation(
        RTMShotDescriptor<RTMData_t, RTMDevPtr_t> * shotDescriptor,
        RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
        RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper,
        RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid,
        RTMCube<RTMData_t, RTMDevPtr_t> *snap0Grid,
        RTMCube<RTMData_t, RTMDevPtr_t> *snap1Grid, 
        RTMGridCollection<RTMData_t,RTMDevPtr_t> *upbGrid);

    void rtmForwardPropagation(
        RTMShotDescriptor<RTMData_t, RTMDevPtr_t> * shotDescriptor,
        RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
        RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper,
        const RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid,
        RTMCube<RTMData_t, RTMDevPtr_t> *snap0Grid,
        RTMCube<RTMData_t, RTMDevPtr_t> *snap1Grid);
private:
    
};

#endif