#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include <Misc.hpp>
#include <RTMGrid.hpp>
#include <RTMController.hpp>
#include <RTMGridPartition.hpp>
#include <RTMKernel.hpp>
#include <RTMBase.hpp>
#include <RTM.hpp>

#ifdef RTM_MPI  
#include <mpi.h>
#endif

#ifdef RTM_MPI
MPI_Comm RTMMPIComm::rtm_mpi_comm_valid;
MPI_Comm RTMMPIComm::rtm_mpi_comm_invalid;
#endif

string getMPIErrorString(int error)
{
    switch (error)
    {
        case(0):
            return "MPI_SUCCESS";
        case(1):
            return "MPI_ERR_BUFFER";
        case(2):
            return "MPI_ERR_COUNT";
        case(3):
            return "MPI_ERR_TYPE";
        case(4):
            return "MPI_ERR_TAG";
        case(5):
            return "MPI_ERR_COMM";
        case(6):
            return "MPI_ERR_RANK";
        case(7):
            return "MPI_ERR_ROOT";
        case(8):
            return "MPI_ERR_GROUP";
        case(9):
            return "MPI_ERR_OP";
        case(10):
            return "MPI_ERR_TOPOLOGY";
        case(11):
            return "MPI_ERR_DIMS";
        case(12):
            return "MPI_ERR_ARG";
        case(13):
            return "MPI_ERR_UNKNOWN";
        case(14):
            return "MPI_ERR_TRUNCATE";
        case(15):
            return "MPI_ERR_OTHER";
        case(16):
            return "MPI_ERR_INTERN";
        case(17):
            return "MPI_ERR_IN_STATUS";
        case(18):
            return "MPI_ERR_PENDING";
        case(19):
            return "MPI_ERR_REQUEST";
        case(20):
            return "MPI_ERR_ACCESS";
        case(21):
            return "MPI_ERR_AMODE";
        case(22):
            return "MPI_ERR_BAD_FILE";
        case(23):
            return "MPI_ERR_CONVERSION";
        case(24):
            return "MPI_ERR_DUP_DATAREP";
        case(25):
            return "MPI_ERR_FILE_EXISTS";
        case(26):
            return "MPI_ERR_FILE_IN_USE";
        case(27):
            return "MPI_ERR_FILE";
        case(28):
            return "MPI_ERR_INFO";
        case(29):
            return "MPI_ERR_INFO_KEY";
        case(30):
            return "MPI_ERR_INFO_VALUE";
        case(31):
            return "MPI_ERR_INFO_NOKEY";
        case(32):
            return "MPI_ERR_IO";
        case(33):
            return "MPI_ERR_NAME";
        case(34):
            return "MPI_ERR_NO_MEM";
        case(35):
            return "MPI_ERR_NOT_SAME";
        case(36):
            return "MPI_ERR_NO_SPACE";
        case(37):
            return "MPI_ERR_NO_SUCH_FILE";
        case(38):
            return "MPI_ERR_PORT";
        case(39):
            return "MPI_ERR_QUOTA";
        case(40):
            return "MPI_ERR_READ_ONLY";
        case(41):
            return "MPI_ERR_SERVICE";
        case(42):
            return "MPI_ERR_SPAWN";
        case(43):
            return "MPI_ERR_UNSUPPORTED_DATAREP";
        case(44):
            return "MPI_ERR_UNSUPPORTED_OPERATION";
        case(45):
            return "MPI_ERR_WIN";
        case(46):
            return "MPI_ERR_BASE";
        case(47):
            return "MPI_ERR_LOCKTYPE";
        case(48):
            return "MPI_ERR_KEYVAL";
        case(49):
            return "MPI_ERR_RMA_CONFLICT";
        case(50):
            return "MPI_ERR_RMA_SYNC";
        case(51):
            return "MPI_ERR_SIZE";
        case(52):
            return "MPI_ERR_DISP";
        case(53):
            return "MPI_ERR_ASSERT";
        case(54):
            return "MPI_ERR_LASTCODE";
        default:
            return "MPI_ERR_LASTCODE";
    }
}

void PRINT_MPI_STATUS(int pRank, string callPoint, MPI_Status & status)
{
#ifdef RTM_MPI
    char *msg = new char[512];
    sprintf(msg, "> P[%d] - %s - MPI_STATUS: .MPI_SOURCE=%d .MPI_TAG=0x%08X .MPI_ERROR=%d (%s) ",
    pRank, callPoint.c_str(), status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR, getMPIErrorString(status.MPI_ERROR).c_str());
    cout << (string(msg)) << endl << flush;

    delete msg;
#endif
}

void  RTMDistributedGridKernel::recvGridGhostZone(RTMGhostZone &gZone, RTMCube<RTMData_t, RTMDevPtr_t> *pGrid, int MPI_GRID_TAG)
{
#ifdef RTM_MPI
    int ix, iy, iz, ik;
    MPI_Status mpiStatus;
    RTMData_t * gData = new RTMData_t[gZone.gLength];
    MPI_Recv(gData, gZone.gLength, 
        MPI_FLOAT, gZone.remoteOwnerProcess, MPI_GRID_TAG, MPI_COMM_WORLD,&mpiStatus);

    ik=0;
    if(mpiStatus.MPI_ERROR==0){
        for(ix=gZone.gLimits.xStart; ix<gZone.gLimits.xEnd; ix++){
        for(iy=gZone.gLimits.yStart; iy<gZone.gLimits.yEnd; iy++){
            for(iz=0; iz<pGrid->getNZ(); iz++){
                pGrid->setWithGlobalCoordinate(ix,iy,iz, gData[ik]);
                ik++;
            }
        }
    }
    }else{
        PRINT_MPI_STATUS(processRank, __func__, mpiStatus);
        exit(-1);
    }
    delete gData;
#endif
}

void  RTMDistributedGridKernel::sendGridGhostZone(RTMGhostZone &gZone, int reqRank, 
RTMCube<RTMData_t, RTMDevPtr_t> *pGrid, int MPI_GRID_TAG){
#ifdef RTM_MPI
    int ix, iy, iz, ik;
    MPI_Status mpiStatus;
    int len = (gZone.gLimits.xEnd-gZone.gLimits.xStart)*(gZone.gLimits.yEnd-gZone.gLimits.yStart)*pGrid->getNZ();
    if (len<=0){
        string s("[ GhostZone Limits Error ]");
        RTMException ex(s);
        throw ex;
    }
    ik=0;
    RTMData_t * gData = new RTMData_t[gZone.gLength];
    for(ix=gZone.gLimits.xStart; ix<gZone.gLimits.xEnd; ix++){
        for(iy=gZone.gLimits.yStart; iy<gZone.gLimits.yEnd; iy++){
            for(iz=0; iz<pGrid->getNZ(); iz++){
                gData[ik] = pGrid->getWithGlobalCoordinate(ix,iy,iz);
                ik++;
            }
        }
    }
    MPI_Send(gData, gZone.gLength, MPI_FLOAT, reqRank, MPI_GRID_TAG, MPI_COMM_WORLD);
    delete gData;
#endif
}

void  RTMDistributedGridKernel::requestGhostZone(RTMGhostZone &gZone, RTMCube<RTMData_t, RTMDevPtr_t> *pGrid)
{
#ifdef RTM_MPI
    int vec[4] = {gZone.gLimits.xStart, gZone.gLimits.xEnd, gZone.gLimits.yStart, gZone.gLimits.yEnd};
    MPI_Send(vec, 4, MPI_INT, gZone.remoteOwnerProcess, RTM_MPI_REQ_GZONE, MPI_COMM_WORLD);
    recvGridGhostZone(gZone, pGrid, RTM_MPI_RCV_GZONE);
#endif
}

void RTMDistributedGridKernel::answerGhostZoneRequest(RTMCube<RTMData_t, RTMDevPtr_t> *pGrid)
{
#ifdef RTM_MPI
    MPI_Status mpiStatus;
    RTMGhostZone gZone;
    int vec[4];
    MPI_Recv(vec, 4, MPI_INT, MPI_ANY_SOURCE, RTM_MPI_REQ_GZONE, MPI_COMM_WORLD, &mpiStatus);
    if (mpiStatus.MPI_ERROR){
        PRINT_MPI_STATUS(processRank, __func__ , mpiStatus);
        exit(-1);
    }
    gZone.gLimits.xStart = vec[0];
    gZone.gLimits.xEnd = vec[1];
    gZone.gLimits.yStart = vec[2];
    gZone.gLimits.yEnd = vec[3];
    gZone.gLength = (gZone.gLimits.xEnd - gZone.gLimits.xStart) * (gZone.gLimits.yEnd - gZone.gLimits.yStart) * pGrid->getNZ();
    int requesterRank = mpiStatus.MPI_SOURCE;

    sendGridGhostZone(gZone, requesterRank, pGrid, RTM_MPI_RCV_GZONE);

#endif
}

void RTMDistributedGridKernel::rtmSynchDistributedGrid(RTMCube<RTMData_t, RTMDevPtr_t> *pGrid)
{
#ifdef RTM_MPI
    if(nProcesses==1){
        return;
    }
    timepoint t0=tic();    
    MPI_Status mpiStatus;
    int validZones = processLimits.validZones; 
    if(validZones<=0)
    {
        return;
    }
    #pragma omp parallel num_threads(2)
    {
        int tid = omp_get_thread_num();
        if (tid%2==0){
            // even thread handles requests...
            for (int i0 = 0; i0 < RTM_MPI_MAX_GZONES; i0++)
            {
                if (processLimits.gzone[i0].remote)
                {
                    // request remote ghost zone
                    //printf("> P[%d.%d]-->P[%d] (Z%d) \n", processRank, tid, processLimits.gzone[i0].remoteOwnerProcess, i0); fflush(stdout);
                    requestGhostZone(processLimits.gzone[i0], pGrid);
                    //printf("< P[%d.%d]<--P[%d] (Z%d) \n", processRank, tid, processLimits.gzone[i0].remoteOwnerProcess, i0); fflush(stdout);
                }
            }
        }else{
            // odd thread handles responses...
            for (int j0 = 0; j0 < validZones; j0++)
            {
                // waits for remoze gzone request
                //printf(">> P[%d.%d]: wait=%d/%d \n", processRank, tid, j0, validZones); fflush(stdout);
                answerGhostZoneRequest(pGrid);
            }
            //printf(">> P[%d.%d]: done=%d/%d \n", processRank, tid, validZones, validZones); fflush(stdout);
        }
        #pragma omp barrier
        tid=0;
    }
    MPI_Barrier(RTM_COMM_VALID);
    report->mpiFuncTime += elapsed_s(t0,toc());
    report->mpiFuncCounter++;
#endif
}

void RTMDistributedGridKernel::joinDistributedReceivers(RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc){
#ifdef RTM_MPI
    if(nProcesses==1 || processLimits.validZones==0){
        return;
    }
    RTMReceiverGrid<RTMData_t,RTMDevPtr_t>* rcvGrid = sDesc.getReceiverGrid();
    
    int pk, ix, iy, it;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    size_t vec[5];
    size_t startX =vec[0]= processLimits.processArea.xStart;
    size_t endX   =vec[1]= processLimits.processArea.xEnd;
    size_t startY =vec[2]= processLimits.processArea.yStart;
    size_t endY   =vec[3]= processLimits.processArea.yEnd;
    size_t rcvCnt=0;

    // calc how many receivers in the process area
    for (ix=0; ix<rcvGrid->getNX(); ix++){
        for (iy=0; iy<rcvGrid->getNY(); iy++){
            int rx = (rcvGrid->getOffsetX() + ix*rcvGrid->getDistanceX())+blen;
            int ry = (rcvGrid->getOffsetY() + iy*rcvGrid->getDistanceY())+blen;
            int rz = rcvGrid->getOffsetZ()+blen;
            if((rx >= startX && rx<endX)&&
                (ry >= startY && ry<endY)){
                rcvCnt+=ntstep;
            }
        }
    }
    vec[4] = rcvCnt;

    MPI_Status mpiStatus;
    int processRank, nProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);// number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank); //the rank of the process
    
    if(processRank==0){
        for(pk=1; pk<nProcesses; pk++){
            // request process area
            int pTag = RTM_MPI_REQ_PAREA;
            MPI_Send(&pTag, 1, MPI_INT, pk, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);
            MPI_Recv(vec, 5, MPI_UINT64_T, pk, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD, &mpiStatus);

            if (mpiStatus.MPI_ERROR)
            {
                PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
            }
            startX = vec[0];
            endX = vec[1];
            startY = vec[2];
            endY = vec[3];
            size_t gSize = vec[4];
            
            
            RTMData_t *gData = new RTMData_t[gSize];
            pTag = RTM_MPI_SEISMRCV_TAG;
            MPI_Send(&pTag, 1, MPI_INT, pk, RTM_MPI_SEISMRCV_TAG, MPI_COMM_WORLD);
            MPI_Recv(gData, gSize, MPI_FLOAT, pk, RTM_MPI_SEISMRCV_TAG, MPI_COMM_WORLD, &mpiStatus);
            if (mpiStatus.MPI_ERROR)
            {
                PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
            }
            //printf("P[%d] %d-%d | %d - %d : length=%d \n", pk, startX, endX, startY, endY, length); fflush(stdout);

            size_t k0=0;
            for (ix=0; ix<rcvGrid->getNX(); ix++){
                for (iy=0; iy<rcvGrid->getNY(); iy++){
                    int rx = (rcvGrid->getOffsetX() + ix*rcvGrid->getDistanceX())+blen;
                    int ry = (rcvGrid->getOffsetY() + iy*rcvGrid->getDistanceY())+blen;
                    int rz = rcvGrid->getOffsetZ()+blen;
                    if((rx >= startX && rx<endX)&&
                        (ry >= startY && ry<endY)){
                        for (it=0; it<ntstep; it++){
                            rcvGrid->set(gData[k0],ix,iy,it);
                            k0++;
                        }
                    }
                }
            }
            delete gData;
        }
    }else{
        RTMData_t *gData = new RTMData_t[rcvCnt];
        size_t k0=0;
        for (ix=0; ix<rcvGrid->getNX(); ix++){
            for (iy=0; iy<rcvGrid->getNY(); iy++){
                int rx = (rcvGrid->getOffsetX() + ix*rcvGrid->getDistanceX())+blen;
                int ry = (rcvGrid->getOffsetY() + iy*rcvGrid->getDistanceY())+blen;
                int rz = rcvGrid->getOffsetZ()+blen;
                if((rx >= startX && rx<endX)&&
                    (ry >= startY && ry<endY)){
                    for (it=0; it<ntstep; it++){
                        gData[k0] = rcvGrid->get(ix,iy,it);
                        k0++;
                    }
                }
            }
        }
        // send process area
        int pTag = 0;
        MPI_Recv(&pTag, 1, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD, &mpiStatus);
        if (mpiStatus.MPI_ERROR)
        {
            PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
        }
        MPI_Send(vec, 5, MPI_UINT64_T, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);
        
        MPI_Recv(&pTag, 1, MPI_INT, 0, RTM_MPI_SEISMRCV_TAG, MPI_COMM_WORLD, &mpiStatus);
        if (mpiStatus.MPI_ERROR)
        {
            PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
        }
        MPI_Send(gData, rcvCnt, MPI_FLOAT, 0, RTM_MPI_SEISMRCV_TAG, MPI_COMM_WORLD);
        delete gData;
    }
#endif
}
