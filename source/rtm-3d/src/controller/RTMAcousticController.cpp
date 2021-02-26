#include <cstdlib>
#include <iostream>
#include <fstream>
#include <RTM.hpp>
#include <RTMController.hpp>

using namespace std;

void RTMAcousticController::runMigrationProcess()
{
    if (rtmShotDescriptors == NULL)
    {
        cout << "> Error: Shot Descriptors List is not initialized\n"
             << endl;
        cout << "> Aborting!" << endl;
        system("exit");
    }
#ifdef RTM_ACC
    rtmCreateDeviceBuffers();
#endif
    if(rtmParam->distributed_grid){
        rtmDistGridMigration();
    }else{
        if(rtmParam->fmig && rtmParam->fmig_distributed_imaging){
            rtmSerialMigration();
        }else{
            rtmDistShotMigration();
        }
    }
#ifdef RTM_ACC
    rtmRemoveDeviceBuffers();
#endif   
}

void RTMAcousticController::runModelingProcess()
{
   if (rtmShotDescriptors == NULL)
    {
        cout << "> Error: Shot Descriptors List is not initialized\n"
             << endl;
        cout << "> Aborting!" << endl;
        system("exit");
    }
#ifdef RTM_ACC
    rtmCreateDeviceBuffers();
#endif
    if(rtmParam->distributed_grid){
        rtmDistGridModeling();
    }else{
        rtmDistShotModeling();
    }
#ifdef RTM_ACC
    rtmRemoveDeviceBuffers();
#endif
}

void RTMAcousticController::rtmSerialMigration()
{

    int nxe = rtmParam->nx + 2 * rtmParam->blen;
    int nye = rtmParam->ny + 2 * rtmParam->blen;
    int nze = rtmParam->nz + 2 * rtmParam->blen;
    int pxe = processLimits.processArea.xEnd - processLimits.processArea.xStart;
    int pye = processLimits.processArea.yEnd - processLimits.processArea.yStart;
    rtmOutputImage = new RTMCube<RTMData_t, RTMDevPtr_t>(nxe, nye, nze);
    rtmOutputImage->setBorderLength(rtmParam->blen);
    rtmOutputImage->fill(0.0);
    RTM_PRINT("Running Serial RTM Migration... ", rtmParam->verbose);
    rtmKernel->initKernel();
    for (int s0 = 0; s0 < rtmShotDescriptors->size(); s0++)
    {
        if(processLimits.validProcessArea){
            RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc = (*rtmShotDescriptors)[s0];
            
            int sx = sDesc.getSource()->getX();
            int sy = sDesc.getSource()->getY();
            int sz = sDesc.getSource()->getZ();
            /* create shot image grid within the limits of the process */
            sDesc.createShotImage(pxe, pye, nze);
            sDesc.getShotImage().setBorderLength(rtmParam->blen);
            sDesc.getShotImage().fill(0.0);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            loadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            if (rtmParam->verbose)
            {
                char *MSG = new char[100];
                sprintf(MSG, "Shot %d/%d (%d,%d,%d):", s0+1, 
                        rtmShotDescriptors->size(), sx, sy, sz);
                RTM_PRINT(string(MSG), rtmParam->verbose);
            }

            /* run shot migration */
            rtmKernel->rtmMigrate(sDesc, *v2dt2SubGrid);

            /* stack shot image */
            rtmOutputImage->stackRegion(sDesc.getShotImage(), processLimits.processArea.xStart,
            processLimits.processArea.xEnd, processLimits.processArea.yStart, 
            processLimits.processArea.yEnd, 0, nze);
                    
            /* destroy shot image grid */
            sDesc.destroyShotImage();

            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            unloadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
        }
    }
    rtmKernel->destroyKernel();
    rtmOutputImage->removeBorders(rtmParam->blen);

    // save shot output image
    string imgFile;
    // filter first
    RTM_PRINT("Applying Laplacian Filter to output image... ", rtmParam->verbose);
    RTMData_t derivativesVec[] = {rtmParam->dx,rtmParam->dy,rtmParam->dz};
    RTMStencil<RTMData_t,RTMDevPtr_t, RTM_NDIM_3D> fkernel(RTM_LAPFILTER_ORDER,derivativesVec);
    
    rtmOutputImage->filter(fkernel);
    RTM_MIGIMG_NAME(imgFile, rtmParam->outdir, rtmParam->mname,
                    rtmParam->source_count_x*rtmParam->source_count_y, rtmParam->nx,
                    rtmParam->ny, rtmParam->nz, rtmParam->nt, processLimits.nProcesses);
    RTM_PRINT("Saving rtm output image at: "+imgFile, rtmParam->verbose);
    if(processLimits.pRank==0){
        rtmOutputImage->saveToFile(imgFile);
    }
    delete rtmOutputImage;
}

void RTMAcousticController::rtmDistShotModeling()
{
    int nProcesses = 1;
    int processRank = 0;
#ifdef RTM_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
#endif

    RTM_PRINT("Running Distributed Shot RTM Modeling Process... ", rtmParam->verbose);
    int nt = rtmParam->nt;
    char *MSG = new char[100];
    for (int sh0 = 0; sh0 < rtmShotDescriptors->size(); sh0+=nProcesses)
    {
        
        int shotNumber = sh0 + processRank;
        if(shotNumber < rtmShotDescriptors->size()){
            RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc = (*rtmShotDescriptors)[shotNumber];
            int sx = sDesc.getSource()->getX();
            int sy = sDesc.getSource()->getY();
            int sz = sDesc.getSource()->getZ();
            
            sprintf(MSG, "Shot %d/%d (%d,%d,%d):", shotNumber+1, 
                    rtmShotDescriptors->size(), sx, sy, sz);
            RTM_PRINT(string(MSG), rtmParam->verbose);
            
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            RTM_PRINT("Allocating Shot receivers memory...", rtmParam->verbose);
            loadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/

            // run modeling process for a single shot
            rtmKernel->rtmModel(sDesc, *v2dt2SubGrid);
            
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            unloadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
        }
    }
    delete MSG;

}

void RTMAcousticController::rtmDistShotMigration()
{
    /**
     * RTMDistShotController Migration Process
     * */
    int nProcesses = 1;
    int processRank = 0;
#ifdef RTM_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
#endif
    int nxe = rtmParam->nx + 2 * rtmParam->blen;
    int nye = rtmParam->ny + 2 * rtmParam->blen;
    int nze = rtmParam->nz + 2 * rtmParam->blen;
    int pxe = processLimits.processArea.xEnd - processLimits.processArea.xStart;
    int pye = processLimits.processArea.yEnd - processLimits.processArea.yStart;
    rtmOutputImage = new RTMCube<RTMData_t, RTMDevPtr_t>(nxe, nye, nze);
    rtmOutputImage->setBorderLength(rtmParam->blen);
    rtmOutputImage->fill(0.0);
    RTM_PRINT("Running Distributed Shot RTM Migration Process... ", rtmParam->verbose);
    rtmKernel->initKernel();
    for (int sh0 = 0; sh0 < rtmShotDescriptors->size(); sh0+=nProcesses)
    {
        int shotNumber = sh0 + processRank;
        if(shotNumber < rtmShotDescriptors->size()){
            RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc = (*rtmShotDescriptors)[shotNumber];
            int sx = sDesc.getSource()->getX();
            int sy = sDesc.getSource()->getY();
            int sz = sDesc.getSource()->getZ();
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            RTM_PRINT("Loading Shot receivers data...", rtmParam->verbose);
            loadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/

            /* create shot image grid within the limits of the process */
            sDesc.createShotImage(pxe, pye, nze);
            sDesc.getShotImage().setBorderLength(rtmParam->blen);
            sDesc.getShotImage().fill(0.0);
            
            if (rtmParam->verbose)
            {
                char *MSG = new char[100];
                sprintf(MSG, "Shot %d/%d (%d,%d,%d):", shotNumber+1, 
                    rtmShotDescriptors->size(), sx, sy, sz);
                RTM_PRINT(string(MSG), rtmParam->verbose);
            }

            /* run shot migration */     
            rtmKernel->rtmMigrate(sDesc, *v2dt2SubGrid);

            /* stack shot image */
            rtmOutputImage->stackRegion(sDesc.getShotImage(), processLimits.processArea.xStart,
            processLimits.processArea.xEnd, processLimits.processArea.yStart, 
            processLimits.processArea.yEnd, 0, nze);
                    
            /* destroy shot image grid */
            sDesc.destroyShotImage();

            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            unloadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
        }
    }
    // now save
#ifdef RTM_MPI
    timepoint t0=tic();
    /* statck instead of joining*/
    stackDistributedOutputImage(rtmOutputImage);
    rtmKernel->getReport().mpiFuncTime += elapsed_s(t0,tic());
    rtmKernel->getReport().mpiFuncCounter++;
#endif
    rtmKernel->destroyKernel();
    rtmOutputImage->removeBorders(rtmParam->blen);

    // save shot output image
    std::string imgFile;
#ifdef RTM_MPI
    imgFile += ".NP"+to_string(nProcesses);
#endif 
    // filter first
    RTM_PRINT("Applying Laplacian Filter to output image... ", rtmParam->verbose);
    RTMData_t derivativesVec[] = {rtmParam->dx,rtmParam->dy,rtmParam->dz};
    RTMStencil<RTMData_t,RTMDevPtr_t, RTM_NDIM_3D> fkernel(RTM_LAPFILTER_ORDER,derivativesVec);
    rtmOutputImage->filter(fkernel);

    RTM_MIGIMG_NAME(imgFile, rtmParam->outdir, rtmParam->mname,
                    rtmParam->source_count_x*rtmParam->source_count_y, rtmParam->nx,
                    rtmParam->ny, rtmParam->nz, rtmParam->nt, nProcesses);
    RTM_PRINT("Saving rtm output image at: "+imgFile, rtmParam->verbose);
#ifdef RTM_MPI
    int pRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
    if(pRank==0) // only pRank=0 saves outputimage
#endif
    {
        rtmOutputImage->saveToFile(imgFile);
    }
    delete rtmOutputImage;
}

void RTMAcousticController::rtmDistGridModeling()
{
    RTM_PRINT("Running RTM Modeling with Distributed Grid Process... ", rtmParam->verbose);
    int nt = rtmParam->nt;
    char *MSG = new char[100];
    rtmKernel->initKernel();
    for (int s0 = 0; s0 < rtmShotDescriptors->size(); s0++)
    {
        if(processLimits.validProcessArea){
            RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc = (*rtmShotDescriptors)[s0];
            int sx = sDesc.getSource()->getX();
            int sy = sDesc.getSource()->getY();
            int sz = sDesc.getSource()->getZ();

            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            loadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/

            if (rtmParam->verbose)
            {
                sprintf(MSG, "Shot %d/%d (%d,%d,%d):", s0+1, 
                        rtmShotDescriptors->size(), sx, sy, sz);
                RTM_PRINT(string(MSG), rtmParam->verbose);
            }

            // run modeling process for a single shot
            rtmKernel->rtmModel(sDesc, *v2dt2SubGrid);

            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            unloadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
        }
    }
    rtmKernel->destroyKernel();
    delete MSG;

}

void RTMAcousticController::rtmDistGridMigration()
{
    /**
     * rtmDistGridMigration Process
     * */
    int nxe = rtmParam->nx + 2 * rtmParam->blen;
    int nye = rtmParam->ny + 2 * rtmParam->blen;
    int nze = rtmParam->nz + 2 * rtmParam->blen;
    int pxe = processLimits.processArea.xEnd - processLimits.processArea.xStart;
    int pye = processLimits.processArea.yEnd - processLimits.processArea.yStart;
    rtmOutputImage = new RTMCube<RTMData_t, RTMDevPtr_t>(nxe, nye, nze);
    rtmOutputImage->setBorderLength(rtmParam->blen);
    rtmOutputImage->fill(0.0);
    RTM_PRINT("Running Distributed Grid RTM Migration... ", rtmParam->verbose);
    rtmKernel->initKernel();
    for (int s0 = 0; s0 < rtmShotDescriptors->size(); s0++)
    {
        if(processLimits.validProcessArea){
            RTMShotDescriptor<RTMData_t,RTMDevPtr_t> &sDesc = (*rtmShotDescriptors)[s0];
            
            int sx = sDesc.getSource()->getX();
            int sy = sDesc.getSource()->getY();
            int sz = sDesc.getSource()->getZ();
            /**
             * create shot image grid within the limits
             * of the process
             *  */
            sDesc.createShotImage(pxe, pye, nze);
            sDesc.getShotImage().setBorderLength(rtmParam->blen);
            sDesc.getShotImage().fill(0.0);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            loadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            if (rtmParam->verbose)
            {
                char *MSG = new char[100];
                sprintf(MSG, "Shot %d/%d (%d,%d,%d):", s0+1, 
                        rtmShotDescriptors->size(), sx, sy, sz);
                RTM_PRINT(string(MSG), rtmParam->verbose);
            }

            /* run shot migration */     
            rtmKernel->rtmMigrate(sDesc, *v2dt2SubGrid);

            /* stack shot image */
            rtmOutputImage->stackRegion(sDesc.getShotImage(), processLimits.processArea.xStart,
            processLimits.processArea.xEnd, processLimits.processArea.yStart, 
            processLimits.processArea.yEnd, 0, nze);
                    
            /* destroy shot image grid */
            sDesc.destroyShotImage();

            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
            unloadShotDescriptorData(sDesc);
            /** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
        }
    }
    // now save
#ifdef RTM_MPI
    timepoint t0=tic();
    stitchDistributedOutputImage(rtmOutputImage);
    rtmKernel->getReport().mpiFuncTime += elapsed_s(t0,tic());
    rtmKernel->getReport().mpiFuncCounter++;
#endif
    rtmKernel->destroyKernel(); 
    rtmOutputImage->removeBorders(rtmParam->blen);

    // save shot output image
    string imgFile;
    int nProcesses = 1;
#ifdef RTM_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    imgFile += ".NP"+to_string(nProcesses);
#endif 
    RTM_PRINT("Applying Laplacian Filter to output image... ", rtmParam->verbose);
    // filter first
    RTMData_t derivativesVec[] = {rtmParam->dx,rtmParam->dy,rtmParam->dz};
    RTMStencil<RTMData_t,RTMDevPtr_t, RTM_NDIM_3D> fkernel(RTM_LAPFILTER_ORDER,derivativesVec);
    rtmOutputImage->filter(fkernel);

    RTM_MIGIMG_NAME(imgFile, rtmParam->outdir, rtmParam->mname,
                    rtmParam->source_count_x*rtmParam->source_count_y, rtmParam->nx,
                    rtmParam->ny, rtmParam->nz, rtmParam->nt, nProcesses);
    RTM_PRINT("Saving rtm output image at: "+imgFile, rtmParam->verbose);
#ifdef RTM_MPI
    int pRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
    if(pRank==0) // only pRank=0 saves outputimage
#endif
    {
        rtmOutputImage->saveToFile(imgFile);
    }
    delete rtmOutputImage;
}

void RTMAcousticController::rtmCreateDeviceBuffers()
{
#ifdef RTM_ACC
    if (v2dt2SubGrid != NULL){
        v2dt2SubGrid->createDeviceBuffer();
        v2dt2SubGrid->moveToDevice();
    }
#endif
}

void RTMAcousticController::rtmRemoveDeviceBuffers()
{
#ifdef RTM_ACC
    if (v2dt2SubGrid != nullptr){
        v2dt2SubGrid->removeDeviceBuffer();
    }
#endif
}