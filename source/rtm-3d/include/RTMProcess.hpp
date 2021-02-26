#ifndef _RTMPROCESS_H
#define _RTMPROCESS_H


#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMParam.hpp>
#include <RTMGridPartition.hpp>
#include <RTMController.hpp>

/**
 * @brief Class RTMProcess
 */
class RTMProcess
{
protected:
    int nProcesses; ///< total number of processes
    int pGlobalRank; ///< global process ID
    int pLocalRank; /// < local process ID

    string              inputJSONFile;
    RTMParam *          rtmInputParam = nullptr; ///< pointer to RTMParam
    RTMController *     rtmController = nullptr; ///< pointer to RTMController
    RTMKernel *         rtmKernel = nullptr; ///< pointer to RTMKernel.
    RTMProcessLimits    pLimits; /// < process limits

public:

    RTMProcess()
    {
        pGlobalRank= 0;
        pLocalRank = 0;
        nProcesses=1;
    }
    /**
     * @brief Construct a new RTMProcess object
     * @param inputFile 
     */
    RTMProcess(string & inputFile)
    {
        inputJSONFile = inputFile;
        pGlobalRank= 0;
        pLocalRank = 0;
        nProcesses=1;
    }
    /// Class destructor
    ~RTMProcess(){
#ifdef RTM_MPI 
        MPI_Finalize();
#endif
    }
    /**
     * Loads RTM Param and get MPI processRank information
     * @return true if success
    */
    bool initRTMProcess();

    /// Pure Virtual Function rtm()
    virtual void rtm() = 0;
    
    /**
     * @brief _enterSequentialZone
     */
    static void _enterSequentialZone();
    /**
     * @brief _leaveSequentialZone
     */
    static void _leaveSequentialZone();

protected:
    /**
     * Every RTM process has a global and
     * a local MPI rank. Local rank
     * is the rank of a process inside
     * work node. Global rank is a unique
     * identificator accross all nodes
    */
    int getProcessLocalRank();

    /**
     * @brief Set the Process Limits
     * @param _param pointer to RTMParam
     */
    void setProcessLimits(RTMParam *_param);

    /**
     * @brief setDistributedGridLimits
     * @param _param pointer to RTMParam
     */
    void setDistributedGridLimits(RTMParam *rtmParam);

    /**
     * @brief Set the Ghost Zones object
     * @param pLim reference to RTMProcessLimits
     * @param nxe integer value, represent the total of samples in X  with border
     * @param nye integer value, represent the total of samples in Y with border
     * @param nze integer value, represent the total of samples in Z  with border
     * @param xSections integer value 
     */
    void setGhostZones(RTMProcessLimits &pLim, int nxe, int nye, int nze, int xSections);
};

class RTMAcousticProcess : public RTMProcess
{
public: 
    RTMAcousticProcess(string & inputFile):RTMProcess{inputFile}{};

    RTMAcousticProcess(RTMParam & _rtmParam):RTMProcess{}
    {
        rtmInputParam = &_rtmParam;
        inputJSONFile = _rtmParam.getJsonFilePath();

    };

    /// Pure Virtual Function rtm()
    virtual void rtm();
};




#endif