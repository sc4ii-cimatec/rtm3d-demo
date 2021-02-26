#ifndef RTMCONTROLLER_H
#define RTMCONTROLLER_H


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <RTMBase.hpp>
#include <RTMParam.hpp>
#include <RTMGrid.hpp>
#include <RTMKernel.hpp>
#include <RTMSeismic.hpp>


using namespace rtmparam;
/**
 * @brief Class RTMController
 * @see RTMController.cpp
 */
class RTMController
{
protected:
    RTMParam *rtmParam = nullptr;
    RTMVelocityModel<RTMData_t, RTMDevPtr_t> *rtmVelModel = nullptr;
    RTMCube<RTMData_t, RTMDevPtr_t> *rtmOutputImage = nullptr;
    vector<RTMShotDescriptor<RTMData_t, RTMDevPtr_t>> *rtmShotDescriptors = nullptr;
    RTMKernel *rtmKernel = nullptr;
    RTMVelocityModel<RTMData_t, RTMDevPtr_t> * v2dt2SubGrid = nullptr;
    RTMProcessLimits processLimits;
public:
    /**
     * @brief Construct a new RTMController object
     * @param _param reference to RTMParam
     * @param _processLimits reference to RTMProcessLimits
     */
    RTMController(RTMParam &_param, RTMProcessLimits & _processLimits)
    {
        rtmParam = &_param;
        processLimits = _processLimits;
    }
    /**
     * @brief Copy constructor RTMController
     */
    RTMController()
    {
    }
    /**
     * @brief Destroy the RTMController object
     */
    ~RTMController()
    {
    }
    /**
     * @brief setRTMKernel
     * @param _kernel reference to RTMKernel
     */
    void setRTMKernel(RTMKernel &_kernel)
    {
        rtmKernel = &_kernel;
    }
    /**
     * @brief setRTMParam
     * @param _p reference to RTMParam
     */
    void setRTMParam(RTMParam &_p)
    {
        rtmParam = &_p;
    }

    /**
     * @brief Pure Virtual Function initRTMController()
     * @details Initialize RTM Serial Controler
     */
    virtual void initRTMController();
    /**
     * @brief Pure Virtual Function runMigrationProcess()
     * @details RTM Serial Migration Process
     */
    virtual void runMigrationProcess() = 0;
    /**
     * @brief Pure Virtual Function runModelingProcess()
     * @details RTM Serial Modeling Process
     */
    virtual void runModelingProcess() = 0;
    /**
     * @brief Function loadRTMVelocitiesModel()
     * @details Load the Velocity model
     */
    virtual void loadRTMVelocitiesModel();
    
    virtual void stackDistributedOutputImage( RTMCube<RTMData_t, RTMDevPtr_t> * outImage);

    virtual void stitchDistributedOutputImage( RTMCube<RTMData_t, RTMDevPtr_t> * outImage);

    virtual void joinDistributedOutputImage( RTMCube<RTMData_t, RTMDevPtr_t> * outImage, bool stack=true);

protected:
    /**
     * @brief Function generateRTMShotDescriptors
     * @details Generate RTM Shot Descriptors
     */
    void generateRTMShotDescriptors();

    /**
     * These functions are necessary in order to avoid allocating
     * all source and receivers memory at the begining of the process.
     * By using these functions the RTMController can allocate/free shot's memory
     * only when it is going to process it. 
     */
    void loadShotDescriptorData(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> & desc);
    void unloadShotDescriptorData(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> & desc);
};

/**
 * @brief Class RTMSerialController Inheritance of class RTMController
 * @details This controller is designed to migrate shots
 *          sequentially on a single workstation (NO_MPI).
 * @see RTMSerialController.cpp
 */
class RTMAcousticController : public RTMController
{
public:
    /**
     * @brief Construct a new RTMSerialController object
     * @param _param 
     * @param _processLimits 
     */
    RTMAcousticController(RTMParam & _param, RTMProcessLimits & _processLimits)
    :RTMController{_param, _processLimits}
    {
    }
    /**
     * @brief Construct a new RTMSerialController object
     * @param _param 
     * @param _kernel 
     * @param _processLimits 
     */
    RTMAcousticController(RTMParam & _param, RTMKernel & _kernel, RTMProcessLimits & _processLimits)
    :RTMController{_param, _processLimits}
    {
        setRTMKernel(_kernel);
    }

    virtual void    runMigrationProcess();
    virtual void    runModelingProcess();
protected:
    void rtmSerialMigration();
    void rtmDistShotModeling();
    void rtmDistGridModeling();
    void rtmDistShotMigration();
    void rtmDistGridMigration();
    void rtmCreateDeviceBuffers();
    void rtmRemoveDeviceBuffers();
};

/**
 * @brief Class RTMFPGAController inheritance of class RTMController
 */
class RTMFPGAController : public RTMController
{
    protected:

    public:
        /**
         * @brief Construct a new RTMFPGAController object
         * @param _param 
         * @param _processLimits 
         */
        RTMFPGAController(RTMParam & _param, RTMProcessLimits & _processLimits)
        :RTMController{_param, _processLimits}
        {
        }
        /**
         * @brief Construct a new RTMFPGAController object
         * @param _param 
         * @param _kernel 
         * @param _processLimits 
         */
        RTMFPGAController(RTMParam & _param, RTMKernel & _kernel, RTMProcessLimits & _processLimits)
        :RTMController{_param, _processLimits}
        {
            setRTMKernel(_kernel);
        }

        virtual void    initRTMController()=0;
        virtual void    runMigrationProcess()=0;
        virtual void    runModelingProcess()=0;
};

#endif

