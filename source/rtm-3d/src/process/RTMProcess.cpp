#include <cstdlib>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <omp.h>
#ifdef RTM_MPI
#include <mpi.h>
#endif
#include <RTM.hpp>

bool RTMProcess::initRTMProcess()
{
#ifdef RTM_MPI   
        //MPI_Init(NULL, NULL);      ///< initialize MPI environment
        int provided;
        MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcesses); ///< number of processes
        MPI_Comm_rank(MPI_COMM_WORLD, &pGlobalRank); ///< the rank of the process
        pLocalRank = getProcessLocalRank();

        char pName[MPI_MAX_PROCESSOR_NAME];
        char pMsg[2048];
        int pNameLen;
        MPI_Get_processor_name( pName, &pNameLen );
        sprintf(pMsg, "RTMProcess Init: pRank=%d; lRank=%d; nProcesses=%d host=%s ",
            pGlobalRank, pLocalRank, nProcesses, pName);
        RTM_PRINT(string(pMsg), true);
#endif
    if (rtmInputParam==nullptr){
        /* Create global instance of RTMParam*/
        rtmInputParam = new RTMParam(inputJSONFile);
        /** Loading input parameters from JSON file */
        rtmInputParam->loadRTMParam();
    }

    /** set process grid boundaries and ghost zoens */
    setProcessLimits(rtmInputParam);

    // only pRank=0 saves snapshots files
    rtmInputParam->save_snapshots =(rtmInputParam->save_snapshots && (nProcesses==1));
    // only pRank=0 prints
    rtmInputParam->verbose =(rtmInputParam->verbose && (pGlobalRank==0));
    
    omp_set_num_threads(rtmInputParam->nthreads);
    /* If it gets to this point, everything is OK*/
    return true;
}

int RTMProcess::getProcessLocalRank(){
    int lRank = 0;
#ifdef RTM_MPI 
    {
        char nodeName[ 1024 ];
        int len;
        // Get the host name
        assert(MPI_Get_processor_name( nodeName, &len )==MPI_SUCCESS);

        // Hash "hostName" into a number (hostId) to ease the communication
        size_t hostId = std::hash<std::string>{}( nodeName );

        // Allocate storage on each process to hold the hostId from all MPI processes
        size_t * hostIds = new size_t[ nProcesses ];

        // Communicate all individual hostId
        assert(MPI_Allgather( &hostId, 1, MPI_UINT64_T, hostIds, 1, 
        MPI_UINT64_T, MPI_COMM_WORLD )==MPI_SUCCESS );

        // Determine which processes are on this node
        int lSize = 0;
        int * localProcesses = new int [ nProcesses ];
        for ( int i = 0; i < nProcesses; i++ )
        {
            if ( hostId == hostIds[ i ] )
            {
                localProcesses[ lSize ] = i;

                if ( pGlobalRank == i )
                    lRank = lSize;

                lSize++;
            }
        }

    }
#endif
    return lRank;
}

void RTMProcess::setProcessLimits(RTMParam *rtmParam)
{
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int half_order = rtmParam->stencil_order/2;
    int nxe = nx + 2 * blen;
    int nye = ny + 2 * blen;
    int nze = nz + 2 * blen;

    if (rtmParam->distributed_grid){
        setDistributedGridLimits(rtmParam);
    }else{
        pLimits.nProcesses = nProcesses;
        pLimits.pRank = pGlobalRank;
        pLimits.lRank = pLocalRank;
        pLimits.validZones = 0;
        pLimits.validProcessArea = true;
        pLimits.processArea.xStart = 0;
        pLimits.processArea.xEnd = nxe;
        pLimits.processArea.yStart=0;
        pLimits.processArea.yEnd=nye;
        pLimits.gzone[0].remote = false;
        pLimits.gzone[1].remote = false;
        pLimits.gzone[2].remote = false;
        pLimits.gzone[3].remote = false;
    }
#ifdef RTM_MPI
    // group all processes with valid area into a single
    // communicator. This is necessary to create
    // synch barriers only between valid processes.
    if(pLimits.validProcessArea){
        MPI_Comm_split(MPI_COMM_WORLD, RTM_COMM_COLOR_VALID, pGlobalRank, &RTM_COMM_VALID);
    }else{
        MPI_Comm_split(MPI_COMM_WORLD, RTM_COMM_COLOR_INVALID, pGlobalRank, &RTM_COMM_INVALID);
    }
#endif
    
}

void RTMProcess::setDistributedGridLimits(RTMParam *rtmParam)
{
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int half_order = rtmParam->stencil_order/2;
    int nxe = nx + 2 * blen;
    int nye = ny + 2 * blen;
    int nze = nz + 2 * blen;

    int ySections = 1;
    int xSections = 1;
    int res = 0;
    int procs = nProcesses;

    do{
        res = procs%2;
        procs/=2;
        if(res==0)ySections *=2;
    }while(res==0);
        
    xSections = nProcesses/ySections;

    if (xSections==1 && ySections > 2 && (ySections%2==0)){
        xSections = 2;
        ySections = ySections/2;
    }

    int x, y;
    int xLen = ((nxe-st_order)/xSections)+1;
    int yLen = ((nye-st_order)/ySections)+1;

    int xStart = half_order;
    int xEnd = xStart+xLen;
    int yStart = half_order;
    int yEnd = yStart + yLen;
    int pCount = 0;
    RTMProcessLimits procLimit;
    for (y=0; y < ySections; y++){
        for (x=0; x < xSections; x++){
            procLimit.shadowLength = half_order;
            procLimit.pRank = pCount;
            procLimit.lRank = pCount;
            procLimit.nProcesses = nProcesses;
            procLimit.processArea.xStart = xStart-half_order;
            procLimit.processArea.yStart = yStart-half_order;
            procLimit.processArea.xEnd = xEnd+half_order;
            procLimit.processArea.yEnd = yEnd+half_order;
            procLimit.validProcessArea = true;
            if(procLimit.processArea.yStart >= nye){
                procLimit.processArea.yStart = nye;
                procLimit.validProcessArea = false;
            }
            if(procLimit.processArea.xStart >= nxe){
                procLimit.processArea.xStart = nxe;
                procLimit.validProcessArea = false;
            }
            if(procLimit.processArea.yEnd > nye) procLimit.processArea.yEnd = nye;
            if(procLimit.processArea.xEnd > nxe) procLimit.processArea.xEnd = nxe;
            
            setGhostZones(procLimit, nxe, nye, nze, xSections);
            if(pCount==pGlobalRank){
                pLimits = procLimit;
                pLimits.lRank = pLocalRank;
                if (rtmParam->verbose){
                    /*****        START OF SEQUENTIAL ZONE     *****/
                    RTMProcess::_enterSequentialZone();
                    printf(">************************************************* \n");
                    printf(">* [nxe=%d nye=%d] xSections = %d xLen = %d | ySections = %d yLen = %d \n", 
                    nxe, nye, xSections, xLen, ySections, yLen);
                    printf(">* Proc=%02d: valid=%d;  valid_zones=%d; x[%03d-%03d];y[%03d-%03d] \n", 
                    pCount,procLimit.validProcessArea,procLimit.validZones,
                    procLimit.processArea.xStart, procLimit.processArea.xEnd, 
                    procLimit.processArea.yStart, procLimit.processArea.yEnd);
                    printf(">*      GZA.xSt=%03d  .xEnd=%03d .ySt=%03d .yEnd=%03d  remote=%d remoteOwner=%d \n", 
                        procLimit.gzone[0].gLimits.xStart, procLimit.gzone[0].gLimits.xEnd,
                        procLimit.gzone[0].gLimits.yStart, procLimit.gzone[0].gLimits.yEnd, procLimit.gzone[0].remote,
                        procLimit.gzone[0].remoteOwnerProcess);
                    printf(">*      GZB.xSt=%03d  .xEnd=%03d .ySt=%03d .yEnd=%03d  remote=%d remoteOwner=%d \n", 
                    procLimit.gzone[1].gLimits.xStart, procLimit.gzone[1].gLimits.xEnd,
                    procLimit.gzone[1].gLimits.yStart, procLimit.gzone[1].gLimits.yEnd, procLimit.gzone[1].remote,
                    procLimit.gzone[1].remoteOwnerProcess);
                    printf(">*      GZC.xSt=%03d  .xEnd=%03d .ySt=%03d .yEnd=%03d  remote=%d remoteOwner=%d \n", 
                    procLimit.gzone[2].gLimits.xStart, procLimit.gzone[2].gLimits.xEnd,
                    procLimit.gzone[2].gLimits.yStart, procLimit.gzone[2].gLimits.yEnd, procLimit.gzone[2].remote,
                    procLimit.gzone[2].remoteOwnerProcess);
                    printf(">*      GZD.xSt=%03d  .xEnd=%03d .ySt=%03d .yEnd=%03d  remote=%d remoteOwner=%d \n", 
                    procLimit.gzone[3].gLimits.xStart, procLimit.gzone[3].gLimits.xEnd,
                    procLimit.gzone[3].gLimits.yStart, procLimit.gzone[3].gLimits.yEnd, procLimit.gzone[3].remote,
                    procLimit.gzone[3].remoteOwnerProcess);
                    printf(">************************************************* \n");
                    
                    RTMProcess::_leaveSequentialZone();
                    /*****        END OF SEQUENTIAL ZONE      *****/
                }

            }            
            xStart+=xLen;
            xEnd+=xLen;
            pCount++;
        }
        yStart+=yLen;
        yEnd+=yLen;
        xStart = half_order;
        xEnd = xStart+xLen;
    }
}


void RTMProcess::setGhostZones(RTMProcessLimits &pLim, int nxe, int nye, int nze, int xSections)
{   // ghost zone A
    int xStart = pLim.processArea.xStart;
    int yStart = pLim.processArea.yStart;
    int xEnd =  pLim.processArea.xEnd;
    int yEnd =  pLim.processArea.yEnd;
    pLim.validZones=0;

    /**
     *  Each process may have up to 4 ghost zones around it
     *    x
     *    -------> 
     *  
     *       (0)AAAAAAAA
     *   (3) ............ (1) 
     *    D  ............  B
     *    D  ............  B
     *    D  ............  B
     *    D  ............  B
     *    D  ............  B
     *      (2)CCCCCCCCC
     * */
    /*******************************************************/
    pLim.gzone[0].gLimits.xStart = xStart+pLim.shadowLength;
    pLim.gzone[0].gLimits.xEnd = xEnd-pLim.shadowLength;
    pLim.gzone[0].gLimits.yStart = yStart;
    pLim.gzone[0].gLimits.yEnd = yStart+pLim.shadowLength;
    if (pLim.gzone[0].gLimits.yEnd <= pLim.shadowLength){
        pLim.gzone[0].remoteOwnerProcess = pLim.pRank;
        pLim.gzone[0].remote = false;
    }else{
        pLim.gzone[0].remoteOwnerProcess=pLim.pRank-xSections;
        pLim.gzone[0].remote=true;
        pLim.validZones++;
    }
    
    
    /*******************************************************/
    pLim.gzone[1].gLimits.xStart = xEnd - pLim.shadowLength;
    pLim.gzone[1].gLimits.xEnd = xEnd;
    pLim.gzone[1].gLimits.yStart = yStart+pLim.shadowLength;
    pLim.gzone[1].gLimits.yEnd = yEnd-pLim.shadowLength;
    if (pLim.gzone[1].gLimits.xStart >= nxe-pLim.shadowLength){
        pLim.gzone[1].remoteOwnerProcess = pLim.pRank;
        pLim.gzone[1].remote = false;
    }else{
        pLim.gzone[1].remoteOwnerProcess=pLim.pRank+1;
        pLim.gzone[1].remote=true;
        pLim.validZones++;
    }

    /*******************************************************/
    pLim.gzone[2].gLimits.xStart = xStart+pLim.shadowLength;
    pLim.gzone[2].gLimits.xEnd = xEnd-pLim.shadowLength;
    pLim.gzone[2].gLimits.yStart = yEnd-pLim.shadowLength;
    pLim.gzone[2].gLimits.yEnd = yEnd;
    if (pLim.gzone[2].gLimits.yStart >= nye-pLim.shadowLength){
        pLim.gzone[2].remoteOwnerProcess = pLim.pRank;
        pLim.gzone[2].remote = false;
    }else{
        pLim.gzone[2].remoteOwnerProcess=pLim.pRank+xSections;
        pLim.gzone[2].remote=true;
        pLim.validZones++;
    }

    /*******************************************************/
    pLim.gzone[3].gLimits.xStart = xStart;
    pLim.gzone[3].gLimits.xEnd = xStart+pLim.shadowLength;
    pLim.gzone[3].gLimits.yStart = yStart+pLim.shadowLength;
    pLim.gzone[3].gLimits.yEnd = yEnd-pLim.shadowLength;
    if (pLim.gzone[3].gLimits.xStart <= pLim.shadowLength){
        pLim.gzone[3].remoteOwnerProcess = pLim.pRank;
        pLim.gzone[3].remote = false;
    }else{
        pLim.gzone[3].remoteOwnerProcess=pLim.pRank-1;
        pLim.gzone[3].remote=true;
        pLim.validZones++;
    }

    int k0=0;
    for (k0=0; k0<4; k0++){
        pLim.gzone[k0].gLength= (pLim.gzone[k0].gLimits.xEnd-pLim.gzone[k0].gLimits.xStart)*
            (pLim.gzone[k0].gLimits.yEnd-pLim.gzone[k0].gLimits.yStart)*nze;
    }
}

void RTMProcess::_enterSequentialZone(){
#ifdef RTM_MPI
    int processRank= 0;
    int nProcesses=1; 
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses); ///< number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank); ///< the rank of the process
    MPI_Status mpiStatus;
    int token = 0xf5;
    if (processRank){
         MPI_Recv(&token, 1, MPI_INT, processRank-1, 0, MPI_COMM_WORLD,&mpiStatus);
         if (mpiStatus.MPI_ERROR){
             PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
         }
    }
#endif
}

void RTMProcess::_leaveSequentialZone(){
#ifdef RTM_MPI
    int processRank= 0;
    int nProcesses=1; 
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses); ///< number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank); ///< the rank of the process    
    MPI_Status mpiStatus;
    int token = 0xf5;
    if(processRank < (nProcesses-1)){
        MPI_Send(&token, 1, MPI_INT, processRank+1, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}