#ifndef RTMKERNEL_H
#define RTMKERNEL_H

#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMParam.hpp>
#include <RTMGridPartition.hpp>
#include <RTMGPUPlatform.hpp>
#include <RTMFPGAPlatform.hpp>

using namespace std;
using namespace rtmparam;

enum class RTMKernelState
{
    RTM_KERNEL_STATE_IDLE,
    RTM_KERNEL_STATE_FORWARD,
    RTM_KERNEL_STATE_BACKWARD
};
/**
 * @brief Class RTMKernelReport Inheritance of class RTMReport
 * @see RTMParam.hpp
 */
class RTMKernelReport : public rtmparam::RTMReport
{

public:
    /**
     * @brief Construct a new RTMKernelReport object
     * @param fname 
     */
    RTMKernelReport(string fname)
        : RTMReport{fname}
    {
        addParam("forward_time", jsontype::number_float);
        addParam("backward_time", jsontype::number_float);
        addParam("modeling_time", jsontype::number_float);
        addParam("migration_time", jsontype::number_float);
        addParam("propagfunc_time", jsontype::number_float);
        addParam("propagfunc_percentual", jsontype::number_float);
        addParam("mpi_time", jsontype::number_float);
        addParam("mpi_percentual", jsontype::number_float);
    }

    float rtmModelingTime = 0.0;
    float rtmMigrationTime = 0.0;
    float rtmForwardTime = 0.0;
    float rtmBackwardTime = 0.0;
    float propagFuncTime = 0.0;
    float propagFuncPercentual = 0.0;
    float mpiFuncTime = 0.0;
    float mpiFuncPercentual = 0.0; 

    float rtmModelingAVGTime = 0.0;
    float rtmMigrationAVGTime = 0.0;
    float rtmForwardAVGTime = 0.0;
    float rtmBackwardAVGTime = 0.0;
    float propagFuncAVGTime = 0.0;
    float mpiFuncAVGTime = 0.0; 

    unsigned long rtmMigrationCounter = 0;
    unsigned long rtmModelingCounter = 0;
    unsigned long rtmForwardCounter = 0;
    unsigned long rtmBackwardCounter = 0;
    unsigned long propagFuncCounter = 0;
    unsigned long mpiFuncCounter = 0;

    string &toString();
    /**
     * @brief Function calcAVGTime
     * @see RTMKernel.cpp
     */
    void calcAVGTime();
    /**
     * @brief Function saveToFile
     * @see RTMKernel.cpp
     */
    void saveToFile();
};

/**
 * @brief Class RTMMPIKernel
 * 
 */
class RTMMPIKernel
{
protected:
    int nProcesses; ///< total number of processes
    int processRank; ///< current process ID
public:
    /**
     * @brief Construct a new RTMMPIKernel object
     */
    RTMMPIKernel()
    {
#ifdef RTM_MPI   
        MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);///< number of processes
        MPI_Comm_rank(MPI_COMM_WORLD, &processRank); ///< the rank of the process
#else
        nProcesses = 1;
        processRank = 0;
#endif
    }
    /**
     * @brief Destroy the RTMMPIKernel object
     */
    ~RTMMPIKernel(){
#ifdef RTM_MPI 
        MPI_Finalize();
#endif
    }
};

/**
 * @brief Class RTMKernel
 * @see RTMKernel.cpp
 */
class RTMKernel
{
protected:
    RTMParam *rtmParam; ///< pointer to RTMParam
    RTMKernelReport *report; ///< pointer to RTMKernelReport
    RTMKernelState KERNEL_STATE; ///< 0 = IDLE; 1 = FWD; 2=BwD
    RTMPlatform<RTMData_t, RTMDevPtr_t> * defaultPlatform;
    RTMAccPlatform * accPlatform;
    RTMCPUPlatform * cpuPlatform;
    bool usingAcc;
public:
    /**
     * @brief Construct a new RTMKernel object
     * @param _param reference to RTMParam
     */
    RTMKernel(RTMParam &_param)
    {
        rtmParam = &_param;
#ifdef RTM_ACC
        usingAcc = true;
#else
        usingAcc = false;
#endif
    }

    void setAccPlatform(){
        usingAcc = true;
        defaultPlatform = accPlatform;
    }
    void setCpuPlatform(){
        usingAcc = false;
        defaultPlatform = cpuPlatform;
    }
    bool isUsingAcc(){
        return usingAcc;
    }

    /**
     * @brief Pure Virtual Function initKernel()
     * @details Initialization RTM Kernel
     */
    virtual void initKernel() = 0;
    /**
     * @brief Pure Virtual Function destroyKernel()
     */
    virtual void destroyKernel() = 0;

    /**
     * @brief Pure Virtual Function rtmMigrate
     * @details Performs an RTM Migration for given shot descriptor.
     *          Shot Image is stored into descriptor's ' RTMCube<T>'
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    virtual void rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                            RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel) = 0;

    /**
     * @brief Pure Virtual Function rtmModel
     * @details Performs an RTM Modeling for a given shot, considering
     *          its current list of receivers positions.
     *          Receiver data is stored into descriptor's RTMSeismicReceiver list.
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    virtual void rtmModel(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                          RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel) = 0;

    /**
     * @brief Function printKernelProgress
     * @param kname pointer to char
     * @param sx Source coordinate - X 
     * @param sy Source coordinate - Y 
     * @param sz Source coordinate - Z 
     * @param it time step
     * @param nt number of time steps
     * @param time_s Time in seconds
     */
    void printKernelProgress(char * kname, int sx, int sy, int sz, int it, int nt, 
                            float time_s, int printStep=100);
   
    /**
     * @brief Get the Report object
     * @return RTMKernelReport& 
     */
    RTMKernelReport &getReport()
    {
        return *report;
    }

};

/**
 * @brief Class RTMDistributedGridKernel Inheritance of class RTMMPIKernel
 */
class RTMDistributedGridKernel : public RTMKernel, public RTMMPIKernel
{
protected:
    RTMProcessLimits processLimits;
    RTMProcessLimits imagingLimits;
public:
    RTMDistributedGridKernel(RTMParam &_param, RTMProcessLimits & _pLimits)
    : RTMKernel{_param},RTMMPIKernel{}
    {
        processLimits = _pLimits;
        imagingLimits = _pLimits;
    }

    /**
     * @brief Function joinDistributedReceivers
     * @param sDesc reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     */
    virtual void joinDistributedReceivers(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &sDesc);

       /**
     * @brief Virtual fucntion rtmSynchDistributedGrids
     * @details Distributed implementations of this class need to merge vectors
     *          spread accross multiple processes, so they must overload this method.
     *          For single-process implementations, however, this method is just empty
     */
    virtual void rtmSynchDistributedGrid(RTMCube<RTMData_t, RTMDevPtr_t> *pGrid);

private:
    /**
     * @brief Function requestGhostZone
     * @param gZone reference to RTMGhostZone
     * @param tag Integer value 
     */
    void requestGhostZone(RTMGhostZone &gZone, RTMCube<RTMData_t, RTMDevPtr_t> *pGrid);

    /**
     * @brief Function answerGhostZoneRequest
     */
    void answerGhostZoneRequest(RTMCube<RTMData_t, RTMDevPtr_t> *pGrid);

    /**
     * @brief Function recvGridGhostZone
     * @param gZone reference to RTMGhostZone
     * @param pGrid pointer to   RTMCube<RTMData_t, RTMDevPtr_t>
     * @param MPI_GRID_TAG Integer value
     */
    void recvGridGhostZone(RTMGhostZone &gZone, RTMCube<RTMData_t, RTMDevPtr_t> *pGrid, int MPI_GRID_TAG);

    /**
     * @brief Function sendGridGhostZone
     * @param gZone reference to RTMGhostZone
     * @param reqRank Integer value
     * @param pGrid  pointer to   RTMCube<RTMData_t, RTMDevPtr_t>
     * @param MPI_GRID_TAG Integer value
     */
    void sendGridGhostZone(RTMGhostZone &gZone, int reqRank, RTMCube<RTMData_t, RTMDevPtr_t> *pGrid, int MPI_GRID_TAG);
    

    RTMProcessLimits & getImagingLimits(){return imagingLimits;}
    RTMProcessLimits & getProcessLimits(){return processLimits;}
};

/**
 * @brief Class RTMFiniteDifferencesKernel Inheritance of classes RTMKernel and RTMDistributedGridKernel
 * @see RTMFiniteDifferencesKernel.cpp and RTMDistributedGridKernel.cpp
 */
class RTMFiniteDifferencesKernel : public RTMDistributedGridKernel
{
private:
    bool finiteDifferencesKernelInitialized = false;
    vector<RTMData_t> *derivativesVec;
protected:
    RTMCube<RTMData_t, RTMDevPtr_t> *pSrcGrid   = nullptr;  ///< P(t) pressure grid
    RTMCube<RTMData_t, RTMDevPtr_t> *ppSrcGrid  = nullptr; /// P(t+1) pressure grid
    RTMCube<RTMData_t, RTMDevPtr_t> *pRcvGrid   = nullptr; ///< PR(nt) backpropagation grid
    RTMCube<RTMData_t, RTMDevPtr_t> *ppRcvGrid  = nullptr; ///< PPR(nt-1) backpropagation grid
    RTMCube<RTMData_t, RTMDevPtr_t> *imgGrid    = nullptr;
    RTMCube<RTMData_t, RTMDevPtr_t> *snap0      = nullptr;  ///< P(NT-2) pressure grid from FWD prop.
    RTMCube<RTMData_t, RTMDevPtr_t> *snap1      = nullptr; ///< P(NT-1) pressure grid from FWD prop.

    RTMStencil<RTMData_t, RTMDevPtr_t> *stencil = nullptr;

    RTMTaperFunction<RTMData_t, RTMDevPtr_t> * rtmTaper; ///< pointer to RTMTaperFunction<RTMData_t, RTMDevPtr_t>
    
public:
    /**
     * @brief Construct a new RTMFiniteDifferencesKernel object
     * @param _param reference to RTMParam
     * @param _pLimits reference to RTMProcessLimits
     */
    RTMFiniteDifferencesKernel(RTMParam &_param, RTMProcessLimits & _pLimits);

    void initKernel();
    void destroyKernel();

    virtual void rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                            RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel) = 0;

    virtual void rtmModel(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                          RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);


    /**
     * @brief Get the Stencil Kernel object
     * @return RTMStencil<RTMData_t, RTMDevPtr_t>& 
     */
    RTMStencil<RTMData_t, RTMDevPtr_t> &getStencilKernel()
    {
        return *stencil;
    }
protected:
    
    void createRTMGrid( RTMCube<RTMData_t, RTMDevPtr_t> ** pGrid, bool createOnAcc);

    void joinDistributedPSGrids();
    void joinDistributedPRGrids();
    void rtmAcousticFiniteDiffModeling(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                          RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);

    void rtmAcousticFiniteDiffModeling_RemoveDirectWave(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                          RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);
};

/**
 * @brief Class RTMHBCKernel Inheritance of class RTMFiniteDifferencesKernel
 * @details RTM kernel for Hybrid Boundary Condition (HBC)
 * @see RTMKernel_HBC.cpp
 */
class RTMHBCKernel : public RTMFiniteDifferencesKernel
{

private:
    bool hbcKernelInitialized = false;
    RTMGridCollection<RTMData_t, RTMDevPtr_t> *upbGrid = nullptr; ///< upper-border grid
public:
    /**
     * @brief Construct a new RTMHBCKernel object
     * 
     * @param _param reference to RTMParam
     * @param _pLimits reference to RTMProcessLimits
     */
    RTMHBCKernel(RTMParam &_param, RTMProcessLimits & _pLimits)
        : RTMFiniteDifferencesKernel{_param, _pLimits}
    {
    }

    virtual void initKernel();
    virtual void destroyKernel();

    /** 
     * Performs an RTM Migration for given shot descriptor.
     * Shot Image is stored into descriptor's ' RTMCube<T>'
     *  
    */
    virtual void rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                            RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);

private:
    /**
     * @brief Function rtmHBCForward
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    void rtmHBCForward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,  RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);
    /**
     * @brief Function rtmHBCBackward
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    void rtmHBCBackward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor, RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);
};


/**
 * @brief Class RTMRBCKernel Inheritance of class RTMFiniteDifferencesKernel
 * @details RTM kernel for Random Boundary Condition (RBC)
 * @see RTMKernel_RBC.cpp
 */
class RTMRBCKernel : public RTMFiniteDifferencesKernel
{

public:
    /**
     * @brief Construct a new RTMRBCKernel object
     * @param _param 
     * @param _pLimits 
     */
    RTMRBCKernel(RTMParam &_param, RTMProcessLimits & _pLimits)
        : RTMFiniteDifferencesKernel{_param, _pLimits}
    {
    }

    virtual void initKernel();
    virtual void destroyKernel();

    /** 
     * Performs an RTM Migration for given shot descriptor.
     * Shot Image is stored into descriptor's ' RTMCube<T>' 
     */
    virtual void rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                            RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);

private:
    /**
     * @brief Function rtmRBCForward
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    void rtmRBCForward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,  RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);
    /**
     * @brief Function rtmRBCBackward
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    void rtmRBCBackward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor, RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);
};


/**
 * @brief Class RTMFreqDomainKernel Inheritance of class RTMFiniteDifferencesKernel
 * @details RTM kernel for Hybrid Boundary Condition (FreDomain)
 * @see RTMKernel_FreDomain.cpp
 */
class RTMFreqDomainKernel : public RTMFiniteDifferencesKernel
{

protected:
    bool    freqDomainKernelInitialized = false;
    float   freq_min;
    float   freq_max;
    float   df;
    size_t  nw;
    size_t  nwstep;
    
    // real and imaginary grids for source and receiver waves
    RTMGridCollection<RTMData_t, RTMDevPtr_t>    * srcReGrid = nullptr;
    RTMGridCollection<RTMData_t, RTMDevPtr_t>    * srcImGrid = nullptr;
    RTMGridCollection<RTMData_t, RTMDevPtr_t>    * rcvReGrid = nullptr;
    RTMGridCollection<RTMData_t, RTMDevPtr_t>    * rcvImGrid = nullptr;

    // real and imaginary frequency kernels
    RTMPlane<RTMData_t, RTMDevPtr_t>  * kernelRe = nullptr;
    RTMPlane<RTMData_t, RTMDevPtr_t>  * kernelIm = nullptr;

    // frequency vector
    RTMVector <RTMData_t, RTMDevPtr_t>   * wList = nullptr;
public:
    /**
     * @brief Construct a new RTMFreqDomainKernel object
     * 
     * @param _param reference to RTMParam
     * @param _pLimits reference to RTMProcessLimits
     */
    RTMFreqDomainKernel(RTMParam &_param, RTMProcessLimits & _pLimits)
        : RTMFiniteDifferencesKernel{_param, _pLimits}
    {
    }

    virtual void initKernel();
    virtual void destroyKernel();

    /** 
     * Performs an RTM Migration for given shot descriptor.
     * Shot Image is stored into descriptor's ' RTMCube<T>'
     *  
    */
    virtual void rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                            RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);

private:
    /**
     * @brief Function rtmFreDomainForward
     * @param shotDescriptor reference to RTMShotDescriptor<RTMData_t, RTMDevPtr_t>
     * @param velmodel reference to RTMVelocityModel<RTMData_t, RTMDevPtr_t>
     */
    void rtmFreqDomainPropagation(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,  RTMVelocityModel<RTMData_t, RTMDevPtr_t> &velmodel);

    // auxiliary functions
    void rtmInitWList(RTMVector<RTMData_t, RTMDevPtr_t> * wList);
    void rtmInitFreqKernels(RTMPlane<RTMData_t, RTMDevPtr_t> * kRe, RTMPlane<RTMData_t, RTMDevPtr_t> * kIm);
    void rtmUpdatePressureGrids(RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid);

    void setDistributedImagingLimits();
    void stitchDistributedOutputImage( RTMCube<RTMData_t, RTMDevPtr_t> * outImage);
    void createFrequencyComponentsGrids();
};

#endif