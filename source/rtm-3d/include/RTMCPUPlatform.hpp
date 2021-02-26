#ifndef _RTM_CPU_PLATFORM_H
#define _RTM_CPU_PLATFORM_H

#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMParam.hpp>
#include <RTMPlatform.hpp>

using namespace std;
using namespace rtmparam;

/**
 * @brief Class RTMPlatform
 * @see RTMKernel.cpp
 */
class RTMCPUPlatform : public RTMPlatform<RTMData_t, RTMDevPtr_t>
{
public: // constructors
     RTMCPUPlatform(RTMParam &_rtmParam, RTMProcessLimits &_pLimits) : RTMPlatform<RTMData_t, RTMDevPtr_t>{_rtmParam, _pLimits} {};

public:
     /**************************************************************************************/
     /** This class implements the following interface from base class RTMPlatform for CPU */
     virtual void initRTMPlatform();
     virtual void destroyRTMPlatform();
     virtual void rtmTaperAllBorders(RTMCube<RTMData_t, RTMDevPtr_t> *grid,
                                     RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper);
     virtual void rtmTaperUpperBorders(RTMCube<RTMData_t, RTMDevPtr_t> *grid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper);
     virtual void rtmStep(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid, RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                          const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);
     virtual void rtmStepMultipleWave(RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
                                      RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid,
                                      RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                                      const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);
     virtual void rtmImageCondition(RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                    RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid);
     virtual void rtmApplySource(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMSeismicSource<RTMData_t, RTMDevPtr_t> *srcGrid, uint32_t it);
     virtual void rtmRestoreReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid, uint32_t it);
     virtual void rtmSaveReceiverData(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid, uint32_t it);
     virtual void rtmSaveUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *upbGrid, uint32_t it);
     virtual void rtmRestoreUpperBorder(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *upbGrid, uint32_t it);

     virtual void rtmUpdateFreqContributions(int it, int iw, int lw,
                                             RTMCube<RTMData_t, RTMDevPtr_t> *PSGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PRGrid,
                                             RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSImGrid,
                                             RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRImGrid,
                                             RTMPlane<RTMData_t, RTMDevPtr_t> *kernelRe, RTMPlane<RTMData_t, RTMDevPtr_t> *kernelIm);

     virtual void rtmFreqDomainImageCondition(int iw, int lw, RTMVector<RTMData_t, RTMDevPtr_t> *wList,
                                              RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid,
                                              RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PSImGrid,
                                              RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRReGrid, RTMGridCollection<RTMData_t, RTMDevPtr_t> *PRImGrid);

     /**************************************************************************************/

private:
     /**
         * @brief Function taper3D
         * @param grid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         */
     void taper3D(RTMCube<RTMData_t, RTMDevPtr_t> *grid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper);

     /**
         * @brief Function taper2D
         * @param grid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         */
     void taper2D(RTMCube<RTMData_t, RTMDevPtr_t> *grid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper);

     /**
         * @brief Function taperUpper3D
         * @param grid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         */
     void taperUpper3D(RTMCube<RTMData_t, RTMDevPtr_t> *grid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper);

     /**
         * @brief Function taperUpper2D
         * @param grid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         */
     void taperUpper2D(RTMCube<RTMData_t, RTMDevPtr_t> *grid, RTMTaperFunction<RTMData_t, RTMDevPtr_t> *rtmTaper);

     /**
         * @brief Function rtmStep3D
         * @param PGrid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         * @param PPGrid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         * @param v2dt2 reference to RTMVelocityModel<RTMData_t,RTMDevPtr_t>
         * @param stencil reference to RTMStencil<RTMData_t,RTMDevPtr_t>
         */
     void rtmStep3D(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid,
                    RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                    const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);

     /**
         * @brief Function rtmStep2D
         * @param PGrid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         * @param PPGrid reference to  RTMCube<RTMData_t, RTMDevPtr_t>
         * @param v2dt2 reference to RTMVelocityModel<RTMData_t,RTMDevPtr_t>
         * @param stencil reference to RTMStencil<RTMData_t,RTMDevPtr_t>
         */
     void rtmStep2D(RTMCube<RTMData_t, RTMDevPtr_t> *PGrid, RTMCube<RTMData_t, RTMDevPtr_t> *PPGrid,
                    RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                    const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);

     void rtmStep3DMultiWave(RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
                             RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid,
                             RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                             const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);

     void rtmStep2DMultiWave(RTMCube<RTMData_t, RTMDevPtr_t> *P0Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP0Grid,
                             RTMCube<RTMData_t, RTMDevPtr_t> *P1Grid, RTMCube<RTMData_t, RTMDevPtr_t> *PP1Grid,
                             RTMStencil<RTMData_t, RTMDevPtr_t> *stencil,
                             const RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);
};

#endif