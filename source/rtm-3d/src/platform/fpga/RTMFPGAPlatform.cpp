
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <Misc.hpp>
#include <RTMGrid.hpp>
#include <RTM.hpp>
#include <RTMController.hpp>
#include <RTMFPGAPlatform.hpp>

#ifdef RTM_ACC_FPGA
using namespace cl;
#include "xcl2/xcl2.hpp"
#endif
using namespace std;

void RTMFPGAPlatform::initRTMPlatform()
{
#ifdef RTM_ACC_FPGA
    if (fpgaDevice!=nullptr){
        delete fpgaDevice;
    }
    RTM_PRINT("Initializing Xilinx FPGA Platform...", rtmParam->verbose);

    /* check model dimensions*/
    assert(ntstep % RTM_FPGA_nFSM == 0);
    assert(plen_z % RTM_FPGA_nPEZ == 0);
    assert(plen_x % RTM_FPGA_nPEX == 0);

    /* creates fpga device, context and command queue */
    fpgaDevice = new RTMFPGADevice(rtmParam->fpga_xclbin);

    fwdKernel = new RTMForwardKernel(fpgaDevice, plen_z, plen_y, plen_x, 
        blen, blen, blen,ntstep, 
        RTM_FPGA_V2DT2_BASE,
        RTM_FPGA_P0_BASE,
        RTM_FPGA_P1_BASE,
        RTM_FPGA_PP0_BASE,
        RTM_FPGA_PP1_BASE);

#endif
}

void RTMFPGAPlatform::destroyRTMPlatform()
{
#ifdef RTM_ACC_FPGA
    if (fpgaDevice!=nullptr){
        delete fpgaDevice;
        fpgaDevice=nullptr;
    }
#endif
}

void RTMFPGAPlatform::rtmSeismicModeling(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> * shotDesc,
        RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
        RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper,
        RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid)
{
#ifdef RTM_ACC_FPGA
    RTM_PRINT("Running FPGA Seismic Modeling Kernel......", rtmParam->verbose);
    size_t blen = rtmParam->blen;
    size_t sx = shotDesc->getSource()->getX();
    size_t sxe = sx+blen;
    size_t sy = shotDesc->getSource()->getY();
    size_t sye = sy+blen;
    size_t sz = shotDesc->getSource()->getZ();
    size_t sze = sz+blen;
    size_t nxe = plen_x;
    size_t nye = plen_y;
    size_t nz = rtmParam->nz;
    size_t nt = rtmParam->nt;
    size_t itstep = 0, lt = 0, it = 0, ix = 0, iy = 0, iz = 0;
    size_t ntstep = rtmParam->ntstep;
    size_t hf_order = rtmParam->stencil_order/2;

    RTMGridCollection<RTMData_t, RTMDevPtr_t> *upbGrid = 
    new RTMGridCollection<RTMData_t, RTMDevPtr_t>(ntstep, nxe, nye,hf_order);
    upbGrid->setBorderLength(blen);
    upbGrid->reset();

    /*load kernel data*/
    bool selF = (ntstep / RTM_FPGA_nFSM) % 2 == 0;
    fwdKernel->setTaper(rtmTaper->getGridBuffer());
    fwdKernel->setCoefs(stencil->getStencilCoefVector(RTMDim::Xdim).getGridBuffer(),
        stencil->getStencilCoefVector(RTMDim::Ydim).getGridBuffer(),
        stencil->getStencilCoefVector(RTMDim::Zdim).getGridBuffer());
    fwdKernel->setV2Dt2(v2dt2Grid.getGridBuffer());
    fwdKernel->createDeviceBuffers(upbGrid->size());
    
    /* Forward propagation */
    RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid = shotDesc->getReceiverGrid();
    for (itstep = 0; itstep < nt; itstep += ntstep)
    {
        // update srce vector
        fwdKernel->setSourceFunction(SLICE(shotDesc->getSource()->getGridBuffer(), itstep, itstep+ntstep-1));
        fwdKernel->runModeling(selF, 0, sye, sxe);
        
        // save upb tmp file
        converter_upb<RTM_FPGA_nPEX, RTM_FPGA_nPEZ, RTM_FPGA_stOrder, RTMData_t>
        (plen_x, plen_y,ntstep,fwdKernel->getUPBBuffer(), upbGrid->getGridBuffer());
        
        /** copy seismogram to receiver grid **/
        for (it=0; it<ntstep; it++){
            for (ix=0; ix<rtmParam->receiver_count_x; ix++){
                for(iy=0; iy<rtmParam->receiver_count_y; iy++){
                    size_t ux = ix*rtmParam->receiver_distance_x + rtmParam->receiver_start_x + blen;
                    size_t uy = iy*rtmParam->receiver_distance_y + rtmParam->receiver_start_y + blen;
                    rcvGrid->set(upbGrid->get(it, ux, uy, hf_order-1), ix, iy, it);
                }
            }
        }
        // save shot seismic traces
        string seismFile;
        RTM_SEISMOGRAM_NAME(seismFile, rtmParam->datdir, sx, sy, sz,
                            rtmParam->receiver_count_x, rtmParam->receiver_count_y,
                            itstep, rtmParam->ntstep);
        RTM_PRINT("", rtmParam->verbose);
        RTM_PRINT("Saving shot seismic traces to '" + seismFile + "' ...", rtmParam->verbose);
        rcvGrid->saveToFile(seismFile);
    }
    delete upbGrid;
#endif
}

/* HBC Forward Propagation */
void RTMFPGAPlatform::rtmForwardPropagation(
        RTMShotDescriptor<RTMData_t, RTMDevPtr_t> * shotDesc,
        RTMStencil<RTMData_t,RTMDevPtr_t> *stencil,
        RTMTaperFunction<RTMData_t,RTMDevPtr_t> *rtmTaper,
        RTMVelocityModel<RTMData_t,RTMDevPtr_t> &v2dt2Grid,
        RTMCube<RTMData_t, RTMDevPtr_t> *snap0Grid,
        RTMCube<RTMData_t, RTMDevPtr_t> *snap1Grid, 
        RTMGridCollection<RTMData_t,RTMDevPtr_t> *upbGrid)
{
#ifdef RTM_ACC_FPGA
    RTM_PRINT("Running FPGA Foward Kernel......", rtmParam->verbose);
    HostBuffer_t<RTM_wideType> l_snap0, l_snap1;
    int sx = shotDesc->getSource()->getX()+blen;
    int sy = shotDesc->getSource()->getY()+blen;
    int sz = shotDesc->getSource()->getZ()+blen;
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    size_t itstep = 0, lt = 0, it = 0, ix = 0, iy = 0, iz = 0;
    size_t ntstep = rtmParam->ntstep;

    /*load kernel data*/
    bool selF = (ntstep / RTM_FPGA_nFSM) % 2 == 0;
    fwdKernel->setTaper(rtmTaper->getGridBuffer());
    fwdKernel->setCoefs(stencil->getStencilCoefVector(RTMDim::Xdim).getGridBuffer(),
        stencil->getStencilCoefVector(RTMDim::Ydim).getGridBuffer(),
        stencil->getStencilCoefVector(RTMDim::Zdim).getGridBuffer());
    fwdKernel->setV2Dt2(v2dt2Grid.getGridBuffer());
    fwdKernel->createDeviceBuffers(upbGrid->size());
    
    /* Forward propagation */
    for (itstep = 0; itstep < nt; itstep += ntstep)
    {
        // update srce vector
        fwdKernel->setSourceFunction(SLICE(shotDesc->getSource()->getGridBuffer(), 
            itstep, itstep+ntstep-1));
        fwdKernel->run(selF, 0, sy, sx, ((itstep+ntstep)>=nt), l_snap0, l_snap1);
        
        // save upb tmp file
        converter_upb<RTM_FPGA_nPEX, RTM_FPGA_nPEZ, RTM_FPGA_stOrder, RTMData_t>
        (plen_x, plen_y,ntstep,fwdKernel->getUPBBuffer(), upbGrid->getGridBuffer());
        string upbFile;
        int unxe = upbGrid->getNX();
        int unye = upbGrid->getNY();
        RTM_HBCUPB_NAME(upbFile, rtmParam->tmpdir, sx-blen, sy-blen, sz-blen, unxe, unye,
                        rtmParam->stencil_order/2, ntstep, itstep, pLimits->pRank, pLimits->nProcesses);
        
        RTM_PRINT("Saving UPB File at '" + upbFile + "' ...", rtmParam->verbose);
        upbGrid->saveToFile(upbFile);
    }
    converter<RTM_FPGA_nPEX, RTM_FPGA_nPEZ, RTMData_t>
        (plen_x, plen_y, plen_z, l_snap1.data(), snap1Grid->getGridBuffer());
    converter<RTM_FPGA_nPEX, RTM_FPGA_nPEZ, RTMData_t>
        (plen_x, plen_y, plen_z, l_snap0.data(), snap0Grid->getGridBuffer());
#endif 
}


