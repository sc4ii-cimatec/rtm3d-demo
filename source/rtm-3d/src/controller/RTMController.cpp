#include <iostream>
#include <vector>
#include <Misc.hpp>
#include <RTM.hpp>


using namespace std;
using namespace rtmparam;

void RTMController::initRTMController()
{
    try
    {
        RTM_PRINT("Initializing RTM Controller...", rtmParam->verbose);
        loadRTMVelocitiesModel();
        generateRTMShotDescriptors();
    }catch (exception &e)
    {
        cout << "> Error:\n"
             << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }
}

/* RTMController::loadRTMVelocitiesModel */
void RTMController::loadRTMVelocitiesModel(){

    try
    {
        if (!rtmParam->isValidated()){
            string s("Invalid input parameters.");
            RTMParamException ex(s);
            throw ex;
        }
        // load velocity model
        if (rtmParam->save_vpe_file){
            RTM_PRINT("Loading Velocity Model '"+ rtmParam->mname + "' at '"+ rtmParam->vpfile + "'." , rtmParam->verbose);
        }else{
           
        }
        rtmVelModel = new RTMVelocityModel<RTMData_t,RTMDevPtr_t>(rtmParam->nx, rtmParam->ny, rtmParam->nz,
                                                        rtmParam->mname, rtmParam->vpfile);

        // extend velocities grid
        RTM_PRINT("Extending '"+ rtmParam->mname +"' borders (blen=" 
        + to_string(rtmParam->blen) + ")", rtmParam->verbose);
        rtmVelModel->extendBorders(rtmParam->blen);

        RTMBoundaryCondition bc;
        if (rtmParam->modeling)
        {
            RTM_PRINT("Initializing velocity model with 'ABC' borders... ", rtmParam->verbose);
            bc = RTMBoundaryCondition::ABC; // modeling always uses ABC
            rtmVelModel->initBorders(bc);
        }
        else
        {
            bc = getRTMBoundaryCondition(rtmParam->boundary_condition);
            string bcName = rtmParam->boundary_condition;
            TO_UPPER(bcName);
            if(rtmParam->load_vpe_from_file){
                 RTM_PRINT("Loading Extended Velocity Model '"+ rtmParam->mname
                  + "' from '"+ rtmParam->vpefile + "'." , rtmParam->verbose);
                rtmVelModel->loadFromFile(rtmParam->vpefile);
            }else{
                RTM_PRINT("Initializing velocity model with '"+ bcName +"' borders... ", rtmParam->verbose);
                rtmVelModel->initBorders(bc);
            }
        }
        if (rtmParam->save_vpe_file)
        {   
            RTM_PRINT("Saving Extended Velocity Model '"+ rtmParam->mname
                  + "' from '"+ rtmParam->vpefile + "'." , rtmParam->verbose);
            rtmVelModel->saveToFile(rtmParam->vpefile);
        }
        
        /** 
         * Vel model must be v2dt2 to avoid extra multiplications
         * during RTM propagation **/
        
        int stX = processLimits.processArea.xStart;
        int eX = processLimits.processArea.xEnd;
        int stY = processLimits.processArea.yStart;
        int eY = processLimits.processArea.yEnd;
        int stZ = 0;
        int eZ = rtmParam->nz+2*rtmParam->blen;
        int nxe = rtmParam->nx+2*rtmParam->blen;
        int nye = rtmParam->ny+2*rtmParam->blen;
        
        if((eX-stX)==nxe&&(eY-stY)==nye){
            /**
             * In this case, process limits
             * are equal to model's extended dimensions.
             * So we can use memory already allocated
             * for 'rtmVelModel'. This is a way
             * to avoid bad_alloc exceptions
             * with larger models when not running
             * distributed grid controllers.
             * */
            v2dt2SubGrid = rtmVelModel;
        }else if(processLimits.validProcessArea){
            RTM_PRINT("Creating V2DT2 subgrid..." , rtmParam->verbose);
            RTMCube<RTMData_t, RTMDevPtr_t> * velSubGrid = rtmVelModel->getSubGrid(stX, eX, stY, eY, stZ, eZ);
            v2dt2SubGrid =static_cast<RTMVelocityModel<RTMData_t,RTMDevPtr_t>*>(velSubGrid);
        }else{
            // just to keep the pointer safe
            v2dt2SubGrid = new RTMVelocityModel<RTMData_t,RTMDevPtr_t>(1,1,1); 
        }
        v2dt2SubGrid->power2();
        v2dt2SubGrid->multiplyBy(rtmParam->dt*rtmParam->dt);
    }
    catch (exception &e)
    {
        cout << "> Error: " << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }
}

/* RTMSerialController::generateRTMShotDescriptors */
void RTMController::generateRTMShotDescriptors()
{

    int nshots = rtmParam->source_count_x * rtmParam->source_count_y;
    RTM_PRINT("Generating RTM Shot Descriptors (nShots="+to_string(nshots)+")... ", rtmParam->verbose);

    rtmShotDescriptors = new vector<RTMShotDescriptor<RTMData_t,RTMDevPtr_t>>(nshots);

    RTMShotDescriptor<RTMData_t,RTMDevPtr_t> *desc;
    for (int i0 = 0; i0 < nshots; i0++)
    {
        desc = new RTMShotDescriptor<RTMData_t,RTMDevPtr_t>();
        desc->setNT(rtmParam->nt);
        (*rtmShotDescriptors)[i0] = *desc;
    }
    int descCnt = 0;
    RTM_PRINT("Setting shot positions... ", rtmParam->verbose);
    for (int ix = 0; ix < rtmParam->source_count_x; ix++)
    {
        for (int iy = 0; iy < rtmParam->source_count_y; iy++)
        {
            desc = &(*rtmShotDescriptors)[descCnt];
            /* Creates source coordinate*/
            int sx = ix * rtmParam->source_distance_x + rtmParam->source_start_x;
            int sy = iy * rtmParam->source_distance_y + rtmParam->source_start_y;
            int sz = rtmParam->source_depth_z;
            /* Checks whether source position is within the model dimensions*/
            if (sx >= rtmParam->nx || sy >= rtmParam->ny || sz >= rtmParam->nz)
            {
                RTMGridCoordinate c(sx, sy, sz);
                string s("[ Invalid Source Position: " + c.toString() + " ]");
                RTMException ex(s);
                throw ex;
            }
            /**
             * Creates RTMSeismicSource object.
             * */
            RTMSeismicSource<RTMData_t,RTMDevPtr_t> *src = new RTMSeismicSource<RTMData_t,RTMDevPtr_t>(sx, sy, sz, rtmParam->nt);
            src->setDt(rtmParam->dt);
            src->setFpeak(rtmParam->fpeak);
            (*rtmShotDescriptors)[descCnt].setSource(src);
            descCnt++;
        }
    }
}


void RTMController::loadShotDescriptorData(RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc){
    // load source samples vector
    sDesc.loadSourceSamples();

    // init seismi receiver vector and load its data
    sDesc.loadReceiverGrid(rtmParam);
#ifdef RTM_ACC
    RTMReceiverGrid<RTMData_t,RTMDevPtr_t> * rcvGrid = sDesc.getReceiverGrid();
    if (rcvGrid!=nullptr)
        rcvGrid->createDeviceBuffer(); // rcvGrid on ACC must be created here
#endif
    
}

void RTMController::unloadShotDescriptorData(RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc){
    /* unload receivers samples */
#ifdef RTM_ACC
    RTMReceiverGrid<RTMData_t,RTMDevPtr_t> *rcvGrid = sDesc.getReceiverGrid();
    if (rcvGrid!=nullptr)
        rcvGrid->removeDeviceBuffer(); // rcvGrid on ACC must be removed here
#endif
    sDesc.unloadReceiverGrid();
    // unload source samples vector
    sDesc.unloadSource();
}

void RTMController::stackDistributedOutputImage(RTMCube<RTMData_t, RTMDevPtr_t> * outImage){
    joinDistributedOutputImage(outImage, true);
}
void RTMController::stitchDistributedOutputImage(RTMCube<RTMData_t, RTMDevPtr_t> * outImage){
    joinDistributedOutputImage(outImage, false);
}

void RTMController::joinDistributedOutputImage(RTMCube<RTMData_t, RTMDevPtr_t> * outImage,  bool stack){
#ifdef RTM_MPI
    int pk, ix, iy, iz;
    int nx = outImage->getNX();
    int ny = outImage->getNY();
    int nz = outImage->getNZ();
    int half_order = rtmParam->stencil_order/2;
    int blen = rtmParam->blen;
    MPI_Status mpiStatus;
    int processRank, nProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);// number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank); //the rank of the process
    
    if(processRank==0){
        for(pk=1; pk<nProcesses; pk++){
            //printf("---> P%d requested IMG from R%d \n", processRank, pk); fflush(stdout);
            // request process area
            int pTag = RTM_MPI_REQ_PAREA;
            MPI_Send(&pTag, 1, MPI_INT, pk, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);

            int vec[4];
            MPI_Recv(vec, 4, MPI_INT, pk, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD, &mpiStatus);

            if (mpiStatus.MPI_ERROR)
            {
                PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
            }
            int startX = vec[0];
            int endX = vec[1];
            int startY = vec[2];
            int endY = vec[3];
            if (startX==-1 || startY==-1 || endX==-1 || endY==-1){
                //invalid area
                continue;
            }else{
                unsigned int length = (endX-startX)*(endY-startY)*nz;
                RTMData_t *gData = new RTMData_t[length];
                pTag = RTM_MPI_OUTIMG_TAG;
                MPI_Send(&pTag, 1, MPI_INT, pk, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD);
                MPI_Recv(gData, length, MPI_FLOAT, pk, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD, &mpiStatus);
                if (mpiStatus.MPI_ERROR)
                {
                    PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
                }

                int k0=0;
                for (ix = startX; ix < endX; ix++){
                    for (iy = startY; iy < endY; iy++)
                    {
                        for (iz = 0; iz < nz; iz++)
                        {
                            RTMData_t remoteVal = gData[k0++];
                            if(stack){
                                RTMData_t localVal = outImage->get(ix, iy, iz);
                                outImage->set((remoteVal + localVal), ix, iy, iz);
                            }else{
                                outImage->set((remoteVal), ix, iy, iz);
                            }
                        }
                    }
                }
                delete gData;
            }
            //printf("---> P%d received IMG from R%d \n", processRank, pk); fflush(stdout);
        }
    }else{
        int vec[4];
        int startX=0, endX=0, startY=0, endY=0;
        RTMData_t *gData=NULL;
        unsigned int length=0;
        int pTag = 0;
        MPI_Recv(&pTag, 1, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD, &mpiStatus);
        if (mpiStatus.MPI_ERROR)
        {
            PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
        }
        //printf("<--- P%d received request IMG from R%d \n", processRank, 0); fflush(stdout);
        if(processLimits.validProcessArea){
            startX =vec[0]= processLimits.processArea.xStart;
            endX   =vec[1]= processLimits.processArea.xEnd;
            startY =vec[2]= processLimits.processArea.yStart;
            endY   =vec[3]= processLimits.processArea.yEnd;
            length = (endX-startX)*(endY-startY)*nz;
            gData = new RTMData_t[length];

            int k0=0;
            for (ix = startX; ix < endX; ix++){
                for (iy = startY; iy < endY; iy++)
                {
                    for (iz = 0; iz < nz; iz++)
                    {
                        gData[k0++] = outImage->get(ix, iy, iz);
                    }
                }
            }
            // send process area
            MPI_Send(vec, 4, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);
            
            MPI_Recv(&pTag, 1, MPI_INT, 0, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD, &mpiStatus);
            if (mpiStatus.MPI_ERROR)
            {
                PRINT_MPI_STATUS(processRank,__func__, mpiStatus);
            }
            MPI_Send(gData, length, MPI_FLOAT, 0, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD);
            delete gData;
            //printf("<--- P%d delivered IMG to R%d \n", processRank, 0); fflush(stdout);
        }else{
            // send invalid process area
            vec[0]=-1;
            vec[1]=-1;
            vec[2]=-1;
            vec[3]=-1;
            MPI_Send(vec, 4, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);
            //printf("<--- P%d delivered invalid img to R%d \n", processRank, 0); fflush(stdout);
        }
    }
#endif
}