#ifndef _RTMPLATFORM_H
#define _RTMPLATFORM_H

#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMParam.hpp>
#include <RTMGridPartition.hpp>

using namespace std;
using namespace rtmparam;

/**
 * @brief Class RTMPlatform
 * @see RTMKernel.cpp
 */
template <typename GridData_type, typename DevPtr_type>
class RTMPlatform
{
protected:
    RTMParam *rtmParam = nullptr;
    RTMProcessLimits *pLimits = nullptr;
public:
     virtual void initRTMPlatform()=0;
     virtual void destroyRTMPlatform()=0;
    /**
     * @brief Function taperAllBorders
     * @details Attenuates all border energy values
     * @param grid reference to  RTMCube<GridData_type, DevPtr_type>
     */
    virtual void rtmTaperAllBorders( RTMCube<GridData_type, DevPtr_type> *grid, RTMTaperFunction<GridData_type,DevPtr_type> *rtmTaper) = 0;

    /**
     * @brief Function taperUpperBorders
     * @details Attenuates upper border energy values
     * @param grid reference to  RTMCube<GridData_type, DevPtr_type>
     */
    virtual void rtmTaperUpperBorders( RTMCube<GridData_type, DevPtr_type> *grid, RTMTaperFunction<GridData_type,DevPtr_type> *rtmTaper) = 0;

    /**
     * @brief Function rtmStep
     * @param PGrid reference to  RTMCube<GridData_type, DevPtr_type>
     * @param PPGrid reference to  RTMCube<GridData_type, DevPtr_type>
     * @param velmodel reference to RTMVelocityModel<GridData_type,DevPtr_type>
     */
    virtual void rtmStep( RTMCube<GridData_type, DevPtr_type> *PGrid,  RTMCube<GridData_type, DevPtr_type> *PPGrid, RTMStencil<GridData_type,DevPtr_type> *stencil,
                         const RTMVelocityModel<GridData_type,DevPtr_type> &v2dt2Grid) = 0;

     virtual void rtmStepMultipleWave( RTMCube<GridData_type, DevPtr_type> *P0Grid,  RTMCube<GridData_type, DevPtr_type> *PP0Grid,
                        RTMCube<GridData_type, DevPtr_type> *P1Grid,  RTMCube<GridData_type, DevPtr_type> *PP1Grid, 
                         RTMStencil<GridData_type,DevPtr_type> *stencil,
                         const RTMVelocityModel<GridData_type,DevPtr_type> &v2dt2Grid) = 0;

    /**
     * @brief Function rtmImageCondition
     * @details Cross-correlation image condition
     * @param imgGrid reference to  RTMCube<GridData_type, DevPtr_type>
     * @param PGrid reference to  RTMCube<GridData_type, DevPtr_type>
     * @param PRGrid reference to  RTMCube<GridData_type, DevPtr_type>
     */
    virtual void rtmImageCondition( RTMCube<GridData_type, DevPtr_type> *imgGrid,  RTMCube<GridData_type, DevPtr_type> *PSGrid,  RTMCube<GridData_type, DevPtr_type> *PRGrid) = 0;

    /**
     * rtmAddSourceEnegery
     * */
    virtual void rtmApplySource( RTMCube<GridData_type, DevPtr_type> *PGrid, RTMSeismicSource<GridData_type,DevPtr_type> *srcGrid, uint32_t it) = 0;

    /**
     * rtmRestoreReceiverData
     * */
    virtual void rtmRestoreReceiverData( RTMCube<GridData_type, DevPtr_type> *PGrid, RTMReceiverGrid<GridData_type,DevPtr_type> *rcvGrid, uint32_t it) = 0;

    /**
     * rtmSaveReceiverData
     * */
    virtual void rtmSaveReceiverData( RTMCube<GridData_type, DevPtr_type> *PGrid, RTMReceiverGrid<GridData_type,DevPtr_type> *rcvGrid, uint32_t it) = 0;

    /**
     * rtmSaveUpperBorder
     * */
    virtual void rtmSaveUpperBorder( RTMCube<GridData_type, DevPtr_type> *PGrid, RTMGridCollection<GridData_type,DevPtr_type> *upbGrid, uint32_t it) = 0;

    /**
     * rtmRestoreUpperBorder
     * */
    virtual void rtmRestoreUpperBorder( RTMCube<GridData_type, DevPtr_type> *PGrid, RTMGridCollection<GridData_type,DevPtr_type> *upbGrid, uint32_t it) = 0;
    
    
    virtual void rtmUpdateFreqContributions(int it, int iw, int lw, 
                               RTMCube<GridData_type, DevPtr_type> *PSGrid,  RTMCube<GridData_type, DevPtr_type> *PRGrid,
                              RTMGridCollection<GridData_type,DevPtr_type> *PSReGrid, RTMGridCollection<GridData_type,DevPtr_type> *PSImGrid,
                              RTMGridCollection<GridData_type,DevPtr_type> *PRReGrid, RTMGridCollection<GridData_type,DevPtr_type> *PRImGrid,
                              RTMPlane<GridData_type,DevPtr_type> * kernelRe, RTMPlane<GridData_type,DevPtr_type> * kernelIm) = 0;

     virtual void rtmFreqDomainImageCondition(int iw, int lw, RTMVector<GridData_type,DevPtr_type> * wList,
                              RTMCube<GridData_type, DevPtr_type> *imgGrid,
                              RTMGridCollection<GridData_type,DevPtr_type> *PSReGrid, RTMGridCollection<GridData_type,DevPtr_type> *PSImGrid,
                              RTMGridCollection<GridData_type,DevPtr_type> *PRReGrid, RTMGridCollection<GridData_type,DevPtr_type> *PRImGrid) = 0;
    
    
    /**
     * @brief Construct a new RTMPlatform object
     * @param _param reference to RTMParam
     */
    RTMPlatform(RTMParam &_rtmParam, RTMProcessLimits &_pLimits)
    {
        rtmParam = &_rtmParam;
        pLimits = &_pLimits;
    }
};
#endif