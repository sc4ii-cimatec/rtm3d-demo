#include <assert.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <Misc.hpp>
#include <RTM.hpp>
#include <RTMGrid.hpp>
#include <RTMCPUPlatform.hpp>

void RTMCPUPlatform::destroyRTMPlatform()
{
    // empty
}

void RTMCPUPlatform::initRTMPlatform()
{
    RTM_PRINT("Initializing CPU Platform...", rtmParam->verbose);
}

void RTMCPUPlatform::rtmUpdateFreqContributions(int it, int iw,int lw, 
                              RTMCube<RTMData_t, RTMDevPtr_t> *PSGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid,
                              RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSImGrid,
                              RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRImGrid,
                              RTMPlane<RTMData_t, RTMDevPtr_t> * kernelRe, RTMPlane<RTMData_t, RTMDevPtr_t> * kernelIm)
{
    int ix, iy, iz;
    int blen = rtmParam->blen;
    int nz = rtmParam->nz;

    int startX = PSReGrid->getStartX();
    int endX = PSReGrid->getEndX();
    int startY = PSReGrid->getStartY();
    int endY = PSReGrid->getEndY();

    // printf("startX=%d endX=%d startY=%d endY=%d \n", startX, endX, startY, endY);

    RTMData_t kre = kernelRe->get(iw,it);
    RTMData_t kim = kernelIm->get(iw,it);
    #pragma omp parallel for private(ix, iy, iz) collapse(2)
    for (ix=startX; ix<endX; ix++){
        for (iy=startY; iy<endY; iy++){
            uint64_t pOffset = PSGrid->getOffset(ix+blen,iy+blen,blen);
            uint64_t fOffset = PSReGrid->getOffset(lw, ix-startX,iy-startY,0);
            for (iz=0; iz<nz; iz++){
                RTMData_t psVal  = PSGrid->getByOffset(pOffset);
                RTMData_t prVal  = PRGrid->getByOffset(pOffset);
                RTMData_t sReVal = PSReGrid->getByOffset(fOffset);
                RTMData_t sImVal = PSImGrid->getByOffset(fOffset);
                RTMData_t rReVal = PRReGrid->getByOffset(fOffset);
                RTMData_t rImVal = PRImGrid->getByOffset(fOffset);
                sReVal += (kre*psVal);
                sImVal += (kim*psVal);
                rReVal += (kre*prVal);
                rImVal += (kim*prVal);
                PSReGrid->setByOffset(fOffset, sReVal);
                PSImGrid->setByOffset(fOffset, sImVal);
                PRReGrid->setByOffset(fOffset, rReVal);
                PRImGrid->setByOffset(fOffset, rImVal);
                pOffset++;
                fOffset++;
            }
        }
    }
}

void RTMCPUPlatform::rtmFreqDomainImageCondition(int iw,int lw, RTMVector<RTMData_t, RTMDevPtr_t> * w2List, 
                                    RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                    RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSImGrid,
                                    RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRImGrid)
{
    int ix, iy, iz;
    int blen = rtmParam->blen;
    int nz = rtmParam->nz;

    int startX = PSReGrid->getStartX();
    int endX = PSReGrid->getEndX();
    int startY = PSReGrid->getStartY();
    int endY = PSReGrid->getEndY();

    RTMData_t w2Val = w2List->get(iw);

    #pragma omp parallel for private(ix, iy, iz) collapse(2)
    for (ix=startX; ix<endX; ix++){
        for (iy=startY; iy<endY; iy++){
            uint64_t pOffset = imgGrid->getOffset(ix+blen,iy+blen,blen);
            uint64_t fOffset = PSReGrid->getOffset(lw, ix-startX,iy-startY,0);
            for (iz=0; iz<nz; iz++){
                RTMData_t imgVal = imgGrid->getByOffset(pOffset);
                RTMData_t sReVal = PSReGrid->getByOffset(fOffset);
                RTMData_t sImVal = PSImGrid->getByOffset(fOffset);
                RTMData_t rReVal = PRReGrid->getByOffset(fOffset);
                RTMData_t rImVal = PRImGrid->getByOffset(fOffset);
                // Im_f = Im_f + (w(iw)*w(iw))*((GSR(:,:,iw)*GRR(:,:,iw))-(GSI(:,:,iw)*GRI(:,:,iw)))
		        imgVal += (w2Val)*(sReVal*rReVal - sImVal*rImVal);
                imgGrid->setByOffset(pOffset, imgVal);
                pOffset++;
                fOffset++;
            }
        }
    }
}

void RTMCPUPlatform::rtmApplySource(RTMCube<RTMData_t, RTMDevPtr_t> * PGrid,
                                     RTMSeismicSource<RTMData_t, RTMDevPtr_t> *srcGrid, 
                                     uint32_t it)
{
    /* source position does not consider the extended grid */
    int sxe = srcGrid->getX() + rtmParam->blen;
    int sye = srcGrid->getY() + rtmParam->blen;
    int sze = srcGrid->getZ() + rtmParam->blen;
    RTMData_t srcVal = (*srcGrid)[it];
    PGrid->applyEnergy(sxe, sye, sze, srcVal);
}
void RTMCPUPlatform::rtmRestoreReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> * PPGrid, RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid, uint32_t it)
{
    int ix, iy;
    int blen = rtmParam->blen;
    #pragma omp parallel for private(ix, iy) collapse (2)
    for (ix=0; ix<rcvGrid->getNX(); ix++){
        for (iy=0; iy<rcvGrid->getNY(); iy++){
            int rx = (rcvGrid->getOffsetX() + ix*rcvGrid->getDistanceX())+blen;
            int ry = (rcvGrid->getOffsetY() + iy*rcvGrid->getDistanceY())+blen;
            int rz = rcvGrid->getOffsetZ()+blen;
            if((rx >= PPGrid->getStartX() && rx<PPGrid->getEndX())&&
            (ry >= PPGrid->getStartY() && ry<PPGrid->getEndY())){
                RTMData_t rval = rcvGrid->get(ix, iy, it);
                PPGrid->set(rval, rx-PPGrid->getStartX(),ry-PPGrid->getStartY(), rz);
            }
        }
    }
}

void RTMCPUPlatform::rtmSaveReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> * PPGrid, RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid, uint32_t it)
{
    int ix, iy;
    int blen = rtmParam->blen;
    #pragma omp parallel for private(ix, iy) collapse (2)
    for (ix=0; ix<rcvGrid->getNX(); ix++){
        for (iy=0; iy<rcvGrid->getNY(); iy++){
            int rx = (rcvGrid->getOffsetX() + ix*rcvGrid->getDistanceX())+blen;
            int ry = (rcvGrid->getOffsetY() + iy*rcvGrid->getDistanceY())+blen;
            int rz = rcvGrid->getOffsetZ()+blen;
            if((rx >= PPGrid->getStartX() && rx<PPGrid->getEndX())&&
            (ry >= PPGrid->getStartY() && ry<PPGrid->getEndY())){
                RTMData_t rval = PPGrid->get(rx-PPGrid->getStartX(),ry-PPGrid->getStartY(), rz);
                rcvGrid->set(rval, ix, iy, it);
            }
        }
    }

}

void RTMCPUPlatform::rtmSaveUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, 
                                        RTMGridCollection<RTMData_t, RTMDevPtr_t> *upbGrid,uint32_t it)
{
    int ix, iy, iz;
    int blen = rtmParam->blen;
    int hf_order = rtmParam->stencil_order/2;
    #pragma omp parallel for private(ix, iy, iz)
    for (ix = 0; ix < PGrid->getNX(); ix++){
        for (iy = 0; iy < PGrid->getNY(); iy++)
        {
            for (iz = blen - hf_order; iz < blen; iz++)
            {
                RTMData_t val0 = PGrid->get(ix,iy,iz);
                upbGrid->set(val0, it, ix, iy, iz - (blen - hf_order));
            }
        }
    }
}
void RTMCPUPlatform::rtmRestoreUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, 
                                            RTMGridCollection<RTMData_t, RTMDevPtr_t> *upbGrid,uint32_t it)
{
    int ix, iy, iz;
    int blen = rtmParam->blen;
    int hf_order = rtmParam->stencil_order/2;
    #pragma omp parallel for private(ix, iy, iz) collapse(3)
    for (ix = 0; ix < PGrid->getNX(); ix++)
    {
        for (iy = 0; iy < PGrid->getNY(); iy++)
        {
            for (iz = blen - hf_order; iz < blen; iz++)
            {
                RTMData_t val0 = upbGrid->get(it,ix,iy,iz-(blen-hf_order));
                PGrid->set(upbGrid->get(it,ix,iy,iz-(blen-hf_order)), ix,iy,iz);
            }
        }
    }
}
/**
 * rtmStep2D
 * 
 * This propagation method considers NX dimension as equal to 1.
 * Therefore the velocities cube becomes a plane (NYxNZ). The DX derivative
 * is not considered in the laplacian computation
 * */
void RTMCPUPlatform::rtmStep2D(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid,
                               RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                               const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{

    int ny = PGrid->getNY();
    int nz = PGrid->getNZ();
    int half_order = rtmParam->stencil_order / 2;
    int st_order = rtmParam->stencil_order;

    int startY = half_order;
    int startZ = half_order;
    int endY = ny - half_order;
    int endZ = nz - half_order;

    int fx = rtmParam->blen; // ignore borders
    int iy = 0, iz = 0;
    #pragma omp parallel for private(iy, iz) collapse(2)
    for (iy = startY; iy < (endY); iy++)
    {
        for (iz = startZ; iz < (endZ); iz++)
        {
            uint64_t pOffset = PGrid->getOffset(fx, iy, iz);
            RTMData_t acmy = 0.0;
            RTMData_t acmz = 0.0;
            RTMData_t lapl = 0.0;
            int ky = 0, kz = 0;
            for (int io = 0; io <= st_order; io++)
            {
                ky = iy - half_order + io;
                kz = iz - half_order + io;
                acmy += PGrid->get(fx, ky, iz) * stencil->getStencilCoef(RTMDim::Ydim, io);
                acmz += PGrid->get(fx, iy, kz) * stencil->getStencilCoef(RTMDim::Zdim, io);
            }
            lapl = acmy + acmz;

            RTMData_t pp = PPGrid->getByOffset(pOffset);
            RTMData_t p = PGrid->getByOffset(pOffset);
            RTMData_t v2dt2 = v2dt2Grid.getByOffset(pOffset);
            RTMData_t npp = 2 * p - pp + (v2dt2)*lapl;
            PPGrid->setByOffset(pOffset, npp);
            if (npp != npp || p != p || pp != pp || lapl != lapl)
            {
                char msg[1024];
                sprintf(msg, "*******************************************************\n");
                sprintf(msg, "> RTM ERROR: Numerical Dispersion at %s \n", __func__);
                sprintf(msg, "%s> PointPosition: (%d,%d,%d)\n", msg, fx, iy, iz);
                sprintf(msg, "%s> Grid Size    : (%d,%d,%d)\n", msg, 2 * fx + 1, ny, nz);
                sprintf(msg, "%s> P            : %.15f \n", msg, p);
                sprintf(msg, "%s> PP           : %.15f \n", msg, pp);
                sprintf(msg, "%s> V2DT2        : %.15f \n", msg, v2dt2);
                sprintf(msg, "%s> LAPL         : %.15f \n", msg, lapl);
                sprintf(msg, "%s*******************************************************\n", msg);
                string s(msg);
                RTMException ex(s);
                throw ex;
            }
        }
    }
}
void RTMCPUPlatform::rtmStep2DMultiWave(
    RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
    RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid,
    RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
    const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    int ny = P0Grid->getNY();
    int nz = P0Grid->getNZ();
    int half_order = rtmParam->stencil_order / 2;
    int st_order = rtmParam->stencil_order;

    int startY = half_order;
    int startZ = half_order;
    int endY = ny - half_order;
    int endZ = nz - half_order;

    int fx = rtmParam->blen; // ignore borders
    int iy = 0, iz = 0;
    #pragma omp parallel for private(iy, iz) collapse(2)
    for (iy = startY; iy < (endY); iy++)
    {
        for (iz = startZ; iz < (endZ); iz++)
        {
            RTMData_t acmy_0 = 0.0, acmy_1=0.0;
            RTMData_t acmz_0 = 0.0, acmz_1=0.0;
            RTMData_t lapl_0 = 0.0, lapl_1=0.0;
            RTMData_t pp0, p0, npp0, pp1, p1, v2dt2, npp1;
            int ky = 0, kz = 0;
            uint64_t pOffset = P0Grid->getOffset(fx, iy, iz);
            for (int io = 0; io <= st_order; io++)
            {
                ky = iy - half_order + io;
                kz = iz - half_order + io;
                acmy_0 += P0Grid->get(fx, ky, iz) * stencil->getStencilCoef(RTMDim::Ydim, io);
                acmz_0 += P0Grid->get(fx, iy, kz) * stencil->getStencilCoef(RTMDim::Zdim, io);
                acmy_1 += P1Grid->get(fx, ky, iz) * stencil->getStencilCoef(RTMDim::Ydim, io);
                acmz_1 += P1Grid->get(fx, iy, kz) * stencil->getStencilCoef(RTMDim::Zdim, io);
            }
            lapl_0 = acmy_0 + acmz_0;
            lapl_1 = acmy_1 + acmz_1;

            v2dt2 = v2dt2Grid.getByOffset(pOffset);
            pp0 = PP0Grid->getByOffset(pOffset);
            p0 = P0Grid->getByOffset(pOffset);
            pp1 = PP1Grid->getByOffset(pOffset);
            p1 = P1Grid->getByOffset(pOffset);
            npp0 = 2 * p0 - pp0 + (v2dt2)*lapl_0;
            npp1 = 2 * p1 - pp1 + (v2dt2)*lapl_1;
            PP0Grid->setByOffset(pOffset, npp0);
            PP1Grid->setByOffset(pOffset, npp1);
            if ((npp0 != npp0 || p0 != p0 || pp0 != pp0 || lapl_0 != lapl_0)||
                (npp1 != npp1 || p1 != p1 || pp1 != pp1 || lapl_1 != lapl_1))
            {
                char msg[1024];
                sprintf(msg, "*******************************************************\n");
                sprintf(msg, "> RTM ERROR: Numerical Dispersion at %s \n", __func__);
                sprintf(msg, "%s> PointPosition: (%d,%d,%d)\n", msg, fx, iy, iz);
                sprintf(msg, "%s> Grid Size    : (%d,%d,%d)\n", msg, 2 * fx + 1, ny, nz);
                sprintf(msg, "%s*******************************************************\n", msg);
                string s(msg);
                RTMException ex(s);
                throw ex;
            }
        }
    }
}


void RTMCPUPlatform::rtmStep3D(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid,
                               RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                               const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    uint32_t nx = PGrid->getNX();
    uint32_t ny = PGrid->getNY();
    uint32_t nz = PGrid->getNZ();
    uint32_t half_order = rtmParam->stencil_order / 2;
    uint32_t st_order = rtmParam->stencil_order;
    uint32_t startX = half_order;
    uint32_t startY = half_order;
    uint32_t startZ = half_order;
    uint32_t endX = nx - half_order;
    uint32_t endY = ny - half_order;
    uint32_t endZ = nz - half_order;
    size_t ix = 0, iy = 0;
    #pragma omp parallel for private(ix, iy) collapse(2)
    for (ix = startX; ix < (endX); ix++)
    {
        for (iy = startY; iy < (endY); iy++)
        {
            size_t pOffset = PGrid->getOffset(ix, iy, startZ);
            for (size_t iz = startZ; iz < (endZ); iz++)
            {
                // int tid = omp_get_thread_num();
                // int nthreads = omp_get_num_threads(); 
                // printf("Thread = %d/%d\n", tid, nthreads);
                RTMData_t acmx = 0.0;
                RTMData_t acmy = 0.0;
                RTMData_t acmz = 0.0;
                RTMData_t lapl = 0.0;
                size_t kx = 0, ky = 0, kz = 0;
                for (size_t io = 0; io <= st_order; io++)
                {
                    kx = ix - half_order + io;
                    ky = iy - half_order + io;
                    kz = iz - half_order + io;
                    acmx += PGrid->get(kx, iy, iz) * stencil->getStencilCoef(RTMDim::Xdim, io);
                    acmy += PGrid->get(ix, ky, iz) * stencil->getStencilCoef(RTMDim::Ydim, io);
                    acmz += PGrid->get(ix, iy, kz) * stencil->getStencilCoef(RTMDim::Zdim, io);
                }
                lapl = acmx + acmy + acmz;

                RTMData_t pp = PPGrid->getByOffset(pOffset);
                RTMData_t p = PGrid->getByOffset(pOffset);
                RTMData_t v2dt2 = v2dt2Grid.getByOffset(pOffset);
                RTMData_t npp = 2 * p - pp + (v2dt2)*lapl;
                PPGrid->setByOffset(pOffset, npp);
                pOffset++;
                if (npp != npp || p != p || pp != pp || lapl != lapl)
                {
                    char msg[1024];
                    sprintf(msg, "\n*******************************************************\n");
                    sprintf(msg, "%s> RTM ERROR: Numerical Dispersion at %s \n", msg, __func__);
                    sprintf(msg, "%s> PointPosition: (%d,%d,%d)\n", msg, ix, iy, iz);
                    sprintf(msg, "%s> Grid Size    : (%d,%d,%d)\n", msg, nx, ny, nz);
                    sprintf(msg, "%s> P            : %.15f \n", msg, p);
                    sprintf(msg, "%s> PP           : %.15f \n", msg, pp);
                    sprintf(msg, "%s> V2DT2        : %.15f \n", msg, v2dt2);
                    sprintf(msg, "%s> LAPL         : %.15f \n", msg, lapl);
                    sprintf(msg, "%s> NPP          : %.15f \n", msg, npp);
                    sprintf(msg, "%s> acmx=%.3f; acmy=%.3f; acmz=%.3f; \n", msg, acmx, acmy, acmz);
                    sprintf(msg, "%s*******************************************************\n", msg);
                    string s(msg);
                    RTMException ex(s);
                    throw ex;
                }
            }
        }
    }
}

void RTMCPUPlatform::rtmStep3DMultiWave(
    RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
    RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid,
    RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
    const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    uint32_t nx = P0Grid->getNX();
    uint32_t ny = P0Grid->getNY();
    uint32_t nz = P0Grid->getNZ();
    uint32_t half_order = rtmParam->stencil_order / 2;
    uint32_t st_order = rtmParam->stencil_order;
    uint32_t startX = half_order;
    uint32_t startY = half_order;
    uint32_t startZ = half_order;
    uint32_t endX = nx - half_order;
    uint32_t endY = ny - half_order;
    uint32_t endZ = nz - half_order;

    uint32_t ix = 0, iy = 0, iz = 0;
    #pragma omp parallel for private(ix, iy, iz) collapse(2)
    for (ix = startX; ix < (endX); ix++)
    {
        for (iy = startY; iy < (endY); iy++)
        {
            uint64_t pOffset = P0Grid->getOffset(ix, iy, startZ);
            for (iz = startZ; iz < (endZ); iz++)
            {   
                RTMData_t acmx_0 = 0.0, acmx_1=0.0;
                RTMData_t acmy_0 = 0.0, acmy_1=0.0;
                RTMData_t acmz_0 = 0.0, acmz_1=0.0;
                RTMData_t lapl_0 = 0.0, lapl_1=0.0;
                RTMData_t pp0, p0, npp0, pp1, p1, v2dt2, npp1;
                size_t kx = 0, ky = 0, kz = 0;
                for (size_t io = 0; io <= st_order; io++)
                {
                    kx = ix - half_order + io;
                    ky = iy - half_order + io;
                    kz = iz - half_order + io;

                    acmx_0 += P0Grid->get(kx, iy, iz) * stencil->getStencilCoef(RTMDim::Xdim, io);
                    acmy_0 += P0Grid->get(ix, ky, iz) * stencil->getStencilCoef(RTMDim::Ydim, io);
                    acmz_0 += P0Grid->get(ix, iy, kz) * stencil->getStencilCoef(RTMDim::Zdim, io);
                    acmx_1 += P1Grid->get(kx, iy, iz) * stencil->getStencilCoef(RTMDim::Xdim, io);
                    acmy_1 += P1Grid->get(ix, ky, iz) * stencil->getStencilCoef(RTMDim::Ydim, io);
                    acmz_1 += P1Grid->get(ix, iy, kz) * stencil->getStencilCoef(RTMDim::Zdim, io);
                }
                lapl_0 = acmx_0 + acmy_0 + acmz_0;
                lapl_1 = acmx_1 + acmy_1 + acmz_1;

                v2dt2 = v2dt2Grid.getByOffset(pOffset);
                pp0 = PP0Grid->getByOffset(pOffset);
                p0 = P0Grid->getByOffset(pOffset);
                npp0 = 2 * p0 - pp0 + (v2dt2)*lapl_0;
                pp1 = PP1Grid->getByOffset(pOffset);
                p1 = P1Grid->getByOffset(pOffset);
                npp1 = 2 * p1 - pp1 + (v2dt2)*lapl_1;
                PP0Grid->setByOffset(pOffset, npp0);
                PP1Grid->setByOffset(pOffset, npp1);
                if ((npp0 != npp0 || p0 != p0 || pp0 != pp0 || lapl_0 != lapl_0)||
                (npp1 != npp1 || p1 != p1 || pp1 != pp1 || lapl_1 != lapl_1))
                {
                    char msg[1024];
                    sprintf(msg, "*******************************************************\n");
                    sprintf(msg, "> RTM ERROR: Numerical Dispersion at %s \n", __func__);
                    sprintf(msg, "%s> PointPosition: (%d,%d,%d)\n", msg, ix, iy, iz);
                    sprintf(msg, "%s> Grid Size    : (%d,%d,%d)\n", msg, nx, ny, nz);
                    sprintf(msg, "%s*******************************************************\n", msg);
                    string s(msg);
                    RTMException ex(s);
                    throw ex;
                }
                pOffset++;
            }
        }
    }
}

void RTMCPUPlatform::rtmStep(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid, RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                             const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    try{
        int nx = rtmParam->nx;
        if (nx == 1)
        { // 2D vel model
            rtmStep2D(PGrid, PPGrid, stencil, v2dt2Grid);
        }
        else if (nx > 1)
        { // 3D vel model
            rtmStep3D(PGrid, PPGrid, stencil, v2dt2Grid);
        }
    }catch(RTMException &e){
		throw e;
    }
}

void RTMCPUPlatform::rtmStepMultipleWave(RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
                         RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid, 
                         RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                         const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    // for now...
    int nx = rtmParam->nx;
    if (nx == 1)
    { // 2D vel model
        rtmStep2DMultiWave(P0Grid, PP0Grid,P1Grid, PP1Grid, stencil, v2dt2Grid);
    }
    else if (nx > 1)
    { // 3D vel model
        rtmStep3DMultiWave(P0Grid, PP0Grid,P1Grid, PP1Grid, stencil, v2dt2Grid);
    }
}

/**
 * Cross-correlation image condition
 * OBS: IMG, PS and PR grids must have the same dimensions
 * */
void RTMCPUPlatform::rtmImageCondition(RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                       RTMCube<RTMData_t, RTMDevPtr_t> *PSGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid)
{
    int ix, iy, iz;
    int nx = imgGrid->getNX();
    int ny = imgGrid->getNY();
    int nz = imgGrid->getNZ();

    assert((imgGrid->getNX() == PSGrid->getNX()) && (imgGrid->getNX() == PRGrid->getNX()));
    assert((imgGrid->getNY() == PSGrid->getNY()) && (imgGrid->getNY() == PRGrid->getNY()));
    assert((imgGrid->getNZ() == PSGrid->getNZ()) && (imgGrid->getNZ() == PRGrid->getNZ()));

//printf("IMG:\n startX:%d endX=%d\n startY:%d endY=%d\n startZ:%d endZ=%d \n",startX, endX, startY, endY, startZ, endZ);
#pragma omp parallel for private(ix, iy, iz) collapse(2)
    for (ix = 0; ix < nx; ix++)
    {
        for (iy = 0; iy < ny; iy++)
        {
            uint64_t pOffset = PSGrid->getOffset(ix, iy, rtmParam->blen);
            for (iz = rtmParam->blen; iz < nz-rtmParam->blen; iz++)
            {
                RTMData_t val0 = imgGrid->getByOffset(pOffset) + (PSGrid->getByOffset(pOffset) * PRGrid->getByOffset(pOffset));
                imgGrid->set(val0, ix, iy, iz);
                pOffset++;
            }
        }
    }
}

void RTMCPUPlatform::rtmTaperAllBorders(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper)
{
    if (rtmParam->nx == 1)
    { // 2D scenario
        taper2D(rtmGrid, rtmTaper);
    }
    else
    { // 3D
        taper3D(rtmGrid, rtmTaper);
    }
}

void RTMCPUPlatform::rtmTaperUpperBorders(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper)
{
    if (rtmParam->nx == 1)
    { // 2D scenario
        taperUpper2D(rtmGrid, rtmTaper);
    }
    else
    { // 3D
        taperUpper3D(rtmGrid, rtmTaper);
    }
}

void RTMCPUPlatform::taperUpper3D(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid,
                                  RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper)
{
    unsigned int offset0, offset1, ix, iy, iz, dist;
    unsigned int lx, ly, lz, nx, ny, nz, blen;
    RTMData_t val0;

    nx = rtmParam->nx + 2 * rtmParam->blen;
    ny = rtmParam->ny + 2 * rtmParam->blen;
    nz = rtmParam->nz + 2 * rtmParam->blen;
    blen = rtmParam->blen;

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;
// top cube
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (ix = blen; ix < (nx - blen); ix++)
    {
        for (iy = blen; iy < (ny - blen); iy++)
        {
            for (iz = 0; iz < blen; iz++)
            {
                // top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);
            }
        }
    }

// now the corners...
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (ix = 0; ix < blen; ix++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            for (iy = 0; iy < ny; iy++)
            {

                // front-top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                if (iy < blen)
                {
                    val0 *= (*rtmTaper)[blen - 1 - iy];
                }
                else if (iy >= (ny - blen))
                {
                    val0 *= (*rtmTaper)[iy - (ny - blen)];
                }
                val0 *= (*rtmTaper)[blen - 1 - ix];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // back-top
                val0 = rtmGrid->getWithGlobalCoordinate(nx - blen + ix, iy, iz);
                if (iy < blen)
                {
                    val0 *= (*rtmTaper)[blen - 1 - iy];
                }
                else if (iy >= (ny - blen))
                {
                    val0 *= (*rtmTaper)[iy - (ny - blen)];
                }
                val0 *= (*rtmTaper)[ix];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate((nx - blen) + ix, iy, iz, val0);
            }
        }
    }
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (iy = 0; iy < blen; iy++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            for (ix = blen; ix < nx - blen; ix++)
            {
                // left-top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - iy];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // right-top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, ny - blen + iy, iz);
                val0 *= (*rtmTaper)[iy];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);
            }
        }
    }
}

void RTMCPUPlatform::taperUpper2D(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid,
                                  RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper)
{
    unsigned int offset0, offset1, ix, iy, iz, dist;
    unsigned int lx, ly, lz, nx, ny, nz, blen;
    RTMData_t val0, MAXVAL, MINVAL;

    nx = rtmParam->nx + 2 * rtmParam->blen;
    ny = rtmParam->ny + 2 * rtmParam->blen;
    nz = rtmParam->nz + 2 * rtmParam->blen;
    blen = rtmParam->blen;

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;
    // top and bottom cube

    ix = blen;
#pragma omp parallel for private(iy, iz, val0) collapse(2)
    for (iy = blen; iy < (ny - blen); iy++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            // top
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
            val0 *= (*rtmTaper)[blen - 1 - iz];
            rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);
        }
    }

    ix = blen;
#pragma omp parallel for private(iy, iz, val0) collapse(2)
    for (iy = 0; iy < blen; iy++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            // left-top
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
            val0 *= (*rtmTaper)[blen - 1 - iy];
            val0 *= (*rtmTaper)[blen - 1 - iz];
            rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

            // right-top
            val0 = rtmGrid->getWithGlobalCoordinate(ix, ny - blen + iy, iz);
            val0 *= (*rtmTaper)[iy];
            val0 *= (*rtmTaper)[blen - 1 - iz];
            rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);
        }
    }
}

void RTMCPUPlatform::taper3D(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid,
                             RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper)
{
    unsigned int offset0, offset1, ix, iy, iz, dist;
    unsigned int lx, ly, lz, nx, ny, nz, blen;
    RTMData_t val0, val1, MAXVAL, MINVAL;

    nx = rtmParam->nx + 2 * rtmParam->blen;
    ny = rtmParam->ny + 2 * rtmParam->blen;
    nz = rtmParam->nz + 2 * rtmParam->blen;
    blen = rtmParam->blen;

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;
// top and bottom cube
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (ix = blen; ix < (nx - blen); ix++)
    {
        for (iy = blen; iy < (ny - blen); iy++)
        {
            for (iz = 0; iz < blen; iz++)
            {
                // top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // bottom
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, (nz - blen) + iz);
                val0 *= (*rtmTaper)[iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, (nz - blen) + iz, val0);
            }
        }
    }

// left and right cubes
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (ix = blen; ix < nx - blen; ix++)
    {
        for (iz = blen; iz < (nz - blen); iz++)
        {
            for (iy = 0; iy < blen; iy++)
            {
                // left
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - iy];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // right
                val0 = rtmGrid->getWithGlobalCoordinate(ix, (ny - blen) + iy, iz);
                val0 *= (*rtmTaper)[iy];
                rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);
            }
        }
    }

// front and back cubes
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (iy = blen; iy < ny - blen; iy++)
    {
        for (iz = blen; iz < (nz - blen); iz++)
        {
            for (ix = 0; ix < blen; ix++)
            {
                // front
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - ix];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // back
                val0 = rtmGrid->getWithGlobalCoordinate((nx - blen) + ix, iy, iz);
                val0 *= (*rtmTaper)[ix];
                rtmGrid->setWithGlobalCoordinate((nx - blen) + ix, iy, iz, val0);
            }
        }
    }
// now the corners...
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (ix = 0; ix < blen; ix++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            for (iy = 0; iy < ny; iy++)
            {

                // front-top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                if (iy < blen)
                {
                    val0 *= (*rtmTaper)[blen - 1 - iy];
                }
                else if (iy >= (ny - blen))
                {
                    val0 *= (*rtmTaper)[iy - (ny - blen)];
                }
                val0 *= (*rtmTaper)[blen - 1 - ix];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // back-top
                val0 = rtmGrid->getWithGlobalCoordinate(nx - blen + ix, iy, iz);
                if (iy < blen)
                {
                    val0 *= (*rtmTaper)[blen - 1 - iy];
                }
                else if (iy >= (ny - blen))
                {
                    val0 *= (*rtmTaper)[iy - (ny - blen)];
                }
                val0 *= (*rtmTaper)[ix];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate((nx - blen) + ix, iy, iz, val0);

                // front-bottom
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, nz - blen + iz);
                if (iy < blen)
                {
                    val0 *= (*rtmTaper)[blen - 1 - iy];
                }
                else if (iy >= (ny - blen))
                {
                    val0 *= (*rtmTaper)[iy - (ny - blen)];
                }
                val0 *= (*rtmTaper)[blen - 1 - ix];
                val0 *= (*rtmTaper)[iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, (nz - blen) + iz, val0);

                // back-bottom
                val0 = rtmGrid->getWithGlobalCoordinate(nx - blen + ix, iy, nz - blen + iz);
                if (iy < blen)
                {
                    val0 *= (*rtmTaper)[blen - 1 - iy];
                }
                else if (iy >= (ny - blen))
                {
                    val0 *= (*rtmTaper)[iy - (ny - blen)];
                }
                val0 *= (*rtmTaper)[ix];
                val0 *= (*rtmTaper)[iz];
                rtmGrid->setWithGlobalCoordinate((nx - blen) + ix, iy, (nz - blen) + iz, val0);
            }
        }
    }
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (iy = 0; iy < blen; iy++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            for (ix = blen; ix < nx - blen; ix++)
            {
                // left-top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - iy];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // right-top
                val0 = rtmGrid->getWithGlobalCoordinate(ix, ny - blen + iy, iz);
                val0 *= (*rtmTaper)[iy];
                val0 *= (*rtmTaper)[blen - 1 - iz];
                rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);

                // left-bottom
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, nz - blen + iz);
                val0 *= (*rtmTaper)[blen - 1 - iy];
                val0 *= (*rtmTaper)[iz];
                rtmGrid->setWithGlobalCoordinate(ix, iy, (nz - blen) + iz, val0);

                // right-bottom
                val0 = rtmGrid->getWithGlobalCoordinate(ix, ny - blen + iy, nz - blen + iz);
                val0 *= (*rtmTaper)[iy];
                val0 *= (*rtmTaper)[iz];
                rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, (nz - blen) + iz, val0);
            }
        }
    }
#pragma omp parallel for private(ix, iy, iz, val0) collapse(3)
    for (ix = 0; ix < blen; ix++)
    {
        for (iy = 0; iy < blen; iy++)
        {
            for (iz = blen; iz < nz - blen; iz++)
            {
                // front-left
                val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - ix];
                val0 *= (*rtmTaper)[blen - 1 - iy];
                rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

                // front-right
                val0 = rtmGrid->getWithGlobalCoordinate(blen, ny - blen + iy, iz);
                val0 *= (*rtmTaper)[blen - 1 - ix];
                val0 *= (*rtmTaper)[iy];
                rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);

                // back-left
                val0 = rtmGrid->getWithGlobalCoordinate(nx - blen + ix, iy, iz);
                val0 *= (*rtmTaper)[ix];
                val0 *= (*rtmTaper)[blen - 1 - iy];
                rtmGrid->setWithGlobalCoordinate((nx - blen) + ix, iy, iz, val0);

                // back-right
                val0 = rtmGrid->getWithGlobalCoordinate(nx - blen + ix, ny - blen + iy, iz);
                val0 *= (*rtmTaper)[ix];
                val0 *= (*rtmTaper)[iy];
                rtmGrid->setWithGlobalCoordinate((nx - blen) + ix, (ny - blen) + iy, iz, val0);
            }
        }
    }
}

void RTMCPUPlatform::taper2D(RTMCube<RTMData_t, RTMDevPtr_t> *rtmGrid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper)
{
    unsigned int offset0, offset1, ix, iy, iz, dist;
    unsigned int lx, ly, lz, nx, ny, nz, blen;
    RTMData_t val0, val1;

    nx = rtmParam->nx + 2 * rtmParam->blen;
    ny = rtmParam->ny + 2 * rtmParam->blen;
    nz = rtmParam->nz + 2 * rtmParam->blen;
    blen = rtmParam->blen;

    lx = nx - 2 * blen;
    ly = ny - 2 * blen;
    lz = nz - 2 * blen;
    // top and bottom cube

    ix = blen;
#pragma omp parallel for private(iy, iz, val0) collapse(2)
    for (iy = blen; iy < (ny - blen); iy++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            // top
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
            val0 *= (*rtmTaper)[blen - 1 - iz];
            rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

            // bottom
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, (nz - blen) + iz);
            val0 *= (*rtmTaper)[iz];
            rtmGrid->setWithGlobalCoordinate(ix, iy, (nz - blen) + iz, val0);
        }
    }

    // left and right cubes
    ix = blen;
#pragma omp parallel for private(iy, iz, val0) collapse(2)
    for (iz = blen; iz < (nz - blen); iz++)
    {
        for (iy = 0; iy < blen; iy++)
        {
            // left
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
            val0 *= (*rtmTaper)[blen - 1 - iy];
            rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

            // right
            val0 = rtmGrid->getWithGlobalCoordinate(ix, (ny - blen) + iy, iz);
            val0 *= (*rtmTaper)[iy];
            rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);
        }
    }

    ix = blen;
#pragma omp parallel for private(iy, iz, val0) collapse(2)
    for (iy = 0; iy < blen; iy++)
    {
        for (iz = 0; iz < blen; iz++)
        {
            // left-top
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, iz);
            val0 *= (*rtmTaper)[blen - 1 - iy];
            val0 *= (*rtmTaper)[blen - 1 - iz];
            rtmGrid->setWithGlobalCoordinate(ix, iy, iz, val0);

            // right-top
            val0 = rtmGrid->getWithGlobalCoordinate(ix, ny - blen + iy, iz);
            val0 *= (*rtmTaper)[iy];
            val0 *= (*rtmTaper)[blen - 1 - iz];
            rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, iz, val0);

            // left-bottom
            val0 = rtmGrid->getWithGlobalCoordinate(ix, iy, nz - blen + iz);
            val0 *= (*rtmTaper)[blen - 1 - iy];
            val0 *= (*rtmTaper)[iz];
            rtmGrid->setWithGlobalCoordinate(ix, iy, (nz - blen) + iz, val0);

            // right-bottom
            val0 = rtmGrid->getWithGlobalCoordinate(ix, ny - blen + iy, nz - blen + iz);
            val0 *= (*rtmTaper)[iy];
            val0 *= (*rtmTaper)[iz];
            rtmGrid->setWithGlobalCoordinate(ix, (ny - blen) + iy, (nz - blen) + iz, val0);
        }
    }
}
