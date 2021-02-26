#ifndef RTMDISTKERNEL_H
#define RTMDISTKERNEL_H

#include <fstream>
#include <vector>
#include <Grid.hpp>
#include <RTMBase.hpp>
#include <RTMParam.hpp>
#include <RTMGrid.hpp>
#include <RTMSeismic.hpp>
#include <RTMException.hpp>
#include <RTMUtil.hpp>

using namespace std;

using namespace rtmparam;

#define RTM_MPI_MAX_GZONES 4
/**
 * @brief Get MPI error string
 * @param error Integer value
 * @return string 
 */
string getMPIErrorString(int error);
/**
 * @brief Function PRINT_MPI_STATUS()
 * @param pRank 
 * @param callPoint 
 * @param status 
 */
void  PRINT_MPI_STATUS(int pRank, string callPoint, MPI_Status & status);
/**
 * @brief Struct RTMRect
 */
typedef struct RTMRect{
    int xStart;
    int xEnd;
    int yStart;
    int yEnd;
    /**
     * @brief Assignment operator 
     * @param srce 
     * @return struct RTMRect& 
     */
    struct RTMRect& operator=(const RTMRect& srce) { 
        xStart = srce.xStart;
        yStart = srce.yStart;
        xEnd = srce.xEnd;
        yEnd = srce.yEnd;
        return *this; 
    }
}RTMRect;
/**
 * @brief Struct RTMGhostZone
 */
typedef struct RTMGhostZone{
    RTMRect gLimits;
    int remoteOwnerProcess;
    unsigned int gLength; // number of points
    bool remote;
    /**
     * @brief Assignment operator
     * @param srce 
     * @return struct RTMGhostZone& 
     */
    struct RTMGhostZone& operator=(const RTMGhostZone& srce) { 
        gLimits = srce.gLimits;
        remoteOwnerProcess = srce.remoteOwnerProcess;
        remote = srce.remote;
        gLength = srce.gLength;
        return *this; 
    }

}RTMGhostZone;

/**
 * @brief Struct RTMProcessLimits
 */
typedef struct RTMProcessLimits{
    int pRank;
    int lRank;
    int nProcesses;
    int shadowLength;
    int validZones;
    bool validProcessArea;
    RTMRect  processArea;

    /**
     * Each process may have up to 4 ghost zones 
     * around it
     *    x
     *    -------> 
     *  
     *      (0)AAAAAAAA
     *   (3) ............ (1) 
     *    D  ............  B
     *    D  ............  B
     *    D  ............  B
     *    D  ............  B
     *    D  ............  B
     *      (2)CCCCCCCCC
     * */
    RTMGhostZone gzone[RTM_MPI_MAX_GZONES];
    /**
     * @brief Assignment operator
     * @param srce 
     * @return struct RTMProcessLimits& 
     */
    struct RTMProcessLimits& operator=(const RTMProcessLimits& srce) { 
        processArea = srce.processArea;
        pRank = srce.pRank;
        lRank = srce.lRank;
        nProcesses = srce.nProcesses;
        validZones = srce.validZones;
        validProcessArea = srce.validProcessArea;
        shadowLength = srce.shadowLength;
        gzone[0] = srce.gzone[0];
        gzone[1] = srce.gzone[1];
        gzone[2] = srce.gzone[2];
        gzone[3] = srce.gzone[3];
        return *this; 
    }
}RTMProcessLimits;

#endif // RTMDISTKERNEL_H