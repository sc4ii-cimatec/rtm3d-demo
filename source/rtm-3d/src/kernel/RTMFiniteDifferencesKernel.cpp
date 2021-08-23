#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <RTM.hpp>

RTMFiniteDifferencesKernel::RTMFiniteDifferencesKernel(RTMParam &_param, RTMProcessLimits &_pLimits)
    : RTMDistributedGridKernel{_param, _pLimits}
{ // create kernel defaultPlatform instance
        accPlatform = new RTMAccPlatform(*rtmParam, processLimits);
        accPlatform->initRTMPlatform();
        // create additional CPU defaultPlatform for verification
        cpuPlatform = new RTMCPUPlatform(*rtmParam, processLimits);
        cpuPlatform->initRTMPlatform();
#ifdef RTM_ACC
        usingAcc = true;
        defaultPlatform = accPlatform;
#else
        usingAcc = false;
        defaultPlatform = cpuPlatform;
#endif
}

void RTMFiniteDifferencesKernel::initKernel()
{
    if (!finiteDifferencesKernelInitialized)
    {
        RTM_PRINT("Initializing RTM Finite Differences Kernel... ", rtmParam->verbose);
        int ndims = RTM_NDIM_3D; // 3 comes from 3D
        int nx = rtmParam->nx;
        int ny = rtmParam->ny;
        int nz = rtmParam->nz;
        int nt = rtmParam->nt;
        int ntstep = rtmParam->ntstep;
        int blen = rtmParam->blen;
        int st_order = rtmParam->stencil_order;
        int nxe = nx + 2 * blen;
        int nye = ny + 2 * blen;
        int nze = nz + 2 * blen;

        // create kernel RTMGrids
        createRTMGrid(&pRcvGrid, true);
        createRTMGrid(&ppRcvGrid, true);
        createRTMGrid(&pSrcGrid, true);
        createRTMGrid(&ppSrcGrid, true);
        createRTMGrid(&imgGrid, true);
        createRTMGrid(&snap0, false);
        createRTMGrid(&snap1, false);

        // create stencil kernel grid
        RTMData_t derivativesVec[] = {rtmParam->dx, rtmParam->dy, rtmParam->dz};
        stencil = new RTMStencil<RTMData_t, RTMDevPtr_t, RTM_NDIM_3D>(st_order, derivativesVec);

        // create taper function grid
        rtmTaper = new RTMTaperFunction<RTMData_t, RTMDevPtr_t>(rtmParam->blen, rtmParam->taper_factor);
        if (isUsingAcc())
        {
            // create stencil and taper buffers on device
            stencil->createDeviceBuffer();
            rtmTaper->createDeviceBuffer();
            stencil->moveToDevice();
            rtmTaper->moveToDevice();
        }

        // kernel initialized
        finiteDifferencesKernelInitialized = true;
    }
}
void RTMFiniteDifferencesKernel::destroyKernel()
{
    if (finiteDifferencesKernelInitialized)
    {
        if (snap0 != nullptr)
            delete snap0;
        if (snap1 != nullptr)
            delete snap1;
        if (pSrcGrid != nullptr)
            delete pSrcGrid;
        if (ppSrcGrid != nullptr)
            delete ppSrcGrid;
        if (pRcvGrid != nullptr)
            delete pRcvGrid;
        if (ppRcvGrid != nullptr)
            delete ppRcvGrid;
        if (imgGrid != nullptr)
            delete imgGrid;
        if (stencil != nullptr)
            delete stencil;
        if (rtmTaper != nullptr)
            delete rtmTaper;
        finiteDifferencesKernelInitialized = false;
    }
}

/**         
 * Performs an RTM Modeling for a given shot, considering
 * its current list of receivers positions. 
 * Receiver data is stored into descriptor's RTMSeismicReceiver list.
 * 
 * OBS: Finite Differences Modeling process always uses ABC borders, thus
 * it makes more sense to be implemented in base class only.
 * */
void RTMFiniteDifferencesKernel::rtmModel(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                          RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    try
    {
        timepoint t0 = tic(), t1;
#ifdef RTM_ACC_FPGA
        RTMFPGAPlatform * fpgaPlatform = dynamic_cast<RTMFPGAPlatform *>(accPlatform);
        fpgaPlatform->rtmSeismicModeling(&shotDescriptor,stencil,rtmTaper,v2dt2Grid);
#else
        rtmAcousticFiniteDiffModeling(shotDescriptor, v2dt2Grid);
        //rtmAcousticFiniteDiffModeling_RemoveDirectWave(shotDescriptor, v2dt2Grid);
#endif

        report->rtmModelingTime += elapsed_s(t0, toc());
        report->rtmModelingCounter++;
    }
    catch (RTMException &e)
    {
        cout << "> Error:\n"
             << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }
}

void RTMFiniteDifferencesKernel::rtmAcousticFiniteDiffModeling(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                                               RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    RTM_PRINT("RTM Modeling Forward Propagation:", rtmParam->verbose);
    timepoint t0 = tic(), t1;
    int itstep = 0, lt = 0, it = 0, ix = 0, iy = 0, iz = 0;
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_FORWARD;
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int nxe = nx + 2 * blen;
    int nye = ny + 2 * blen;
    int nze = nz + 2 * blen;
    RTMData_t dt2 = rtmParam->dt * rtmParam->dt;
    int sx = shotDescriptor.getSource()->getX();
    int sy = shotDescriptor.getSource()->getY();
    int sz = shotDescriptor.getSource()->getZ();

    string snapshotsFile;
    if (rtmParam->save_snapshots)
    {
        int snap_step = rtmParam->snapshot_step;
        snapshotsFile += rtmParam->tmpdir + "/RTMSNAPSHOT_S" + to_string(sx + blen) + "x" + to_string(sy + blen) +
                         "_" + to_string(nt / snap_step) + "x" + to_string(nye) +
                         "x" + to_string(nze) + ".snapshot.RTM_MOD";
        /* touch and create a new file */
        ofstream ofs(snapshotsFile, ios::out | ios::binary);
        ofs.close();
        RTM_PRINT("Saving RTM_MOD snapshots file at " + snapshotsFile, rtmParam->verbose);
    }
    // contour values are zero
    pSrcGrid->reset();  // fill method updates device buffers as well
    ppSrcGrid->reset(); // fill method updates device buffers as well

    RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid = shotDescriptor.getReceiverGrid();

    for (itstep = 0; itstep < nt; itstep += ntstep)
    {
        // foward propagation
        for (lt = 0; lt < ntstep; lt++)
        {
            // update global time step
            it = itstep + lt;
            // swap pointers
            RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);

            // attenuate borders
            defaultPlatform->rtmTaperAllBorders(pSrcGrid, rtmTaper);
            defaultPlatform->rtmTaperAllBorders(ppSrcGrid, rtmTaper);

            // propagate wave for one timestep
            timepoint t2 = tic();
            defaultPlatform->rtmStep(pSrcGrid, ppSrcGrid, stencil, v2dt2Grid);

            // update propagation function time report
            report->propagFuncTime += elapsed_s(t2, toc());
            report->propagFuncCounter++;

            // apply source
            defaultPlatform->rtmApplySource(ppSrcGrid, shotDescriptor.getSource(), it);

            // join all distributed PS grids, if any. 
            // !!! MUST BE BEFORE SAVING RECEIVER DATA
            joinDistributedPSGrids();

            // save reciever grid
            defaultPlatform->rtmSaveReceiverData(ppSrcGrid, rcvGrid, lt);

            // update snapshots file
            if (rtmParam->save_snapshots && (it % rtmParam->snapshot_step) == 0 && processRank==0)
            {
                /* For now, snapshots are always on X dimension*/
                if (isUsingAcc())
                {
                    ppSrcGrid->moveFromDevice();
                }
                RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = ppSrcGrid->get2DSection(sx + blen, RTMDim::Xdim);
                secGrid->appendTofile(snapshotsFile);
                delete secGrid;
            }
            
            // print kernel progress
            printKernelProgress("+ RTM_MOD", sx, sy, sz, it, nt, elapsed_s(t0, toc()));
        } // for (lt...

        // join all distributed RCV grids, if any.
        // in case an acc is used, get the latest version from device
        if (isUsingAcc())
        {
            rcvGrid->moveFromDevice();
        }
        joinDistributedReceivers(shotDescriptor);

        string seismFile;
        RTM_SEISMOGRAM_NAME(seismFile, rtmParam->datdir, sx, sy, sz,
                            rtmParam->receiver_count_x, rtmParam->receiver_count_y,
                            itstep, rtmParam->ntstep);
        // save shot seismic traces
        RTM_PRINT("", rtmParam->verbose);
        RTM_PRINT("Saving shot seismic traces to '" + seismFile + "' ...", rtmParam->verbose);
        if (rtmParam->distributed_grid)
        {
            // with disGrid only pRank=0 saves .seism file
            if (processRank == 0)
            {
                rcvGrid->saveToFile(seismFile);
            }
        }
        else
        {
            rcvGrid->saveToFile(seismFile);
        }
    } // for (itstep ...
}

void RTMFiniteDifferencesKernel::joinDistributedPSGrids()
{
    if (nProcesses > 1 && processLimits.validZones > 0)
    {
        if (isUsingAcc())
        {
            // int deviceID;
            // deviceID=processLimits.lRank%4;
            // CUDACHECK(cudaSetDevice(deviceID));

            // processLimits.pRank,deviceID, processLimits.lRank, pSrcGrid->getDevPtr(),
            // pSrcGrid->getDevPtr()); fflush(stdout);
            pSrcGrid->moveFromDevice();
            ppSrcGrid->moveFromDevice();
        }

        rtmSynchDistributedGrid(ppSrcGrid);
        rtmSynchDistributedGrid(pSrcGrid);
        if (isUsingAcc())
        {
            pSrcGrid->moveToDevice();
            ppSrcGrid->moveToDevice();
        }
    }
}

void RTMFiniteDifferencesKernel::joinDistributedPRGrids()
{
    if (nProcesses > 1 && processLimits.validZones > 0)
    {
        if (isUsingAcc())
        {
            pSrcGrid->moveFromDevice();
            ppSrcGrid->moveFromDevice();
            pRcvGrid->moveFromDevice();
            ppRcvGrid->moveFromDevice();
        }
        rtmSynchDistributedGrid(pSrcGrid);
        rtmSynchDistributedGrid(ppSrcGrid);
        rtmSynchDistributedGrid(pRcvGrid);
        rtmSynchDistributedGrid(ppRcvGrid);
        if (isUsingAcc())
        {
            pSrcGrid->moveToDevice();
            ppSrcGrid->moveToDevice();
            pRcvGrid->moveToDevice();
            ppRcvGrid->moveToDevice();
        }
    }
}

void RTMFiniteDifferencesKernel::rtmAcousticFiniteDiffModeling_RemoveDirectWave(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                                                                RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    try
    {
        RTM_PRINT("RTM Modeling Forward Propagation:", rtmParam->verbose);
        timepoint t0 = tic(), t1;
        KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_FORWARD;
        size_t itstep = 0, lt = 0, it = 0, ix = 0, iy = 0, iz = 0;
        size_t nx = rtmParam->nx;
        size_t ny = rtmParam->ny;
        size_t nz = rtmParam->nz;
        size_t nt = rtmParam->nt;
        size_t ntstep = rtmParam->ntstep;
        size_t blen = rtmParam->blen;
        size_t st_order = rtmParam->stencil_order;
        size_t nxe = nx + 2 * blen;
        size_t nye = ny + 2 * blen;
        size_t nze = nz + 2 * blen;
        RTMData_t dt2 = rtmParam->dt * rtmParam->dt;
        size_t sx = shotDescriptor.getSource()->getX();
        size_t sy = shotDescriptor.getSource()->getY();
        size_t sz = shotDescriptor.getSource()->getZ();

        string snapshotsFile;
        if (rtmParam->save_snapshots)
        {
            int snap_step = rtmParam->snapshot_step;
            snapshotsFile += rtmParam->tmpdir + "/RTMSNAPSHOT_S" + to_string(sx + blen) + "x" + to_string(sy + blen) +
                             "_" + to_string(nt / snap_step) + "x" + to_string(nye) +
                             "x" + to_string(nze) + ".snapshot.RTM_MOD";
            /* touch and create a new file */
            ofstream ofs(snapshotsFile, ios::out | ios::binary);
            ofs.close();
            RTM_PRINT("Saving RTM_MOD snapshots file at " + snapshotsFile, rtmParam->verbose);
        }
        // contour values are zero
        pSrcGrid->reset();  // fill method updates device buffers as well
        ppSrcGrid->reset(); // fill method updates device buffers as well
        pRcvGrid->reset();
        ppRcvGrid->reset();

        // create shadow velgrid
        RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid = shotDescriptor.getReceiverGrid();
        RTMVelocityModel<RTMData_t, RTMDevPtr_t> *filterV2dt2Grid = v2dt2Grid.createDirectWaveFilterGrid(blen);
        RTMReceiverGrid<RTMData_t, RTMDevPtr_t> &filterRcvGrid = *(new RTMReceiverGrid<RTMData_t, RTMDevPtr_t>(rcvGrid->getNX(), rcvGrid->getNY(), rcvGrid->getNZ()));
        filterRcvGrid.setDistanceX(rcvGrid->getDistanceX());
        filterRcvGrid.setDistanceY(rcvGrid->getDistanceY());
        filterRcvGrid.setOffsetX(rcvGrid->getOffsetX());
        filterRcvGrid.setOffsetY(rcvGrid->getOffsetY());
        filterRcvGrid.setOffsetZ(rcvGrid->getOffsetZ());
        if (isUsingAcc())
        {
            filterV2dt2Grid->createDeviceBuffer();
            filterV2dt2Grid->moveToDevice();
            filterRcvGrid.createDeviceBuffer();
            filterRcvGrid.reset();
        }
        for (itstep = 0; itstep < nt; itstep += ntstep)
        {
            // foward propagation
            for (lt = 0; lt < ntstep; lt++)
            {
                // update global time step
                it = itstep + lt;

                // swap pointers
                RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);
                RTMGRID_SWAP(&pRcvGrid, &ppRcvGrid);

                // attenuate borders
                defaultPlatform->rtmTaperAllBorders(pSrcGrid, rtmTaper);
                defaultPlatform->rtmTaperAllBorders(ppSrcGrid, rtmTaper);
                defaultPlatform->rtmTaperAllBorders(pRcvGrid, rtmTaper);
                defaultPlatform->rtmTaperAllBorders(ppRcvGrid, rtmTaper);

                // propagate wave for one timestep
                timepoint t2 = tic();
                defaultPlatform->rtmStep(pSrcGrid, ppSrcGrid, stencil, v2dt2Grid);
                //defaultPlatform->rtmStep(pRcvGrid, ppRcvGrid, stencil, *filterV2dt2Grid);

                // update propagation function time report
                report->propagFuncTime += elapsed_s(t2, toc());
                report->propagFuncCounter++;

                // apply source
                defaultPlatform->rtmApplySource(ppSrcGrid, shotDescriptor.getSource(), it);
                defaultPlatform->rtmApplySource(ppRcvGrid, shotDescriptor.getSource(), it);

                // join all distributed PS grids, if any. 
                // !!! MUST BE BEFORE SAVING RECEIVER DATA
                joinDistributedPSGrids();

                // save reciever grid
                defaultPlatform->rtmSaveReceiverData(ppSrcGrid, rcvGrid, lt);
                defaultPlatform->rtmSaveReceiverData(ppRcvGrid, &filterRcvGrid, lt);

                // update snapshots file
                if (rtmParam->save_snapshots && (it % rtmParam->snapshot_step) == 0)
                {
                    /* For now, snapshots are always on X dimension*/
                    RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = ppSrcGrid->get2DSection(sx + blen, RTMDim::Xdim);
                    secGrid->appendTofile(snapshotsFile);
                    delete secGrid;
                }

                // print kernel progress
                printKernelProgress("+ RTM_MOD", sx, sy, sz, it, nt, elapsed_s(t0, toc()));
            }
            // join all distributed RCV grids, if any.
            // in case an acc is used, get the latest version from device
            if (isUsingAcc())
            {
                rcvGrid->moveFromDevice();
                filterRcvGrid.moveFromDevice();
            }
            joinDistributedReceivers(shotDescriptor);

            // sub direct wave
            rcvGrid->subtractBy(filterRcvGrid);

            string seismFile;
            RTM_SEISMOGRAM_NAME(seismFile, rtmParam->datdir, sx, sy, sz,
                                rtmParam->receiver_count_x, rtmParam->receiver_count_y,
                                itstep, rtmParam->ntstep);
            // save shot seismic traces
            RTM_PRINT("", rtmParam->verbose);
            RTM_PRINT("Saving shot seismic traces to '" + seismFile + "' ...", rtmParam->verbose);
            if (rtmParam->distributed_grid)
            {
                // with disGrid only pRank=0 saves .seism file
                if (processRank == 0)
                {
                    rcvGrid->saveToFile(seismFile);
                }
            }
            else
            {
                rcvGrid->saveToFile(seismFile);
            }
        } // for (int it = 0; it<nt; it++)
        if (isUsingAcc())
        {
            filterRcvGrid.removeDeviceBuffer();
            filterV2dt2Grid->removeDeviceBuffer();
        }
        filterV2dt2Grid->destroyGrid();
        filterRcvGrid.destroyGrid();
    }
    catch (RTMException &e)
    {
        throw e;
    }
}

void RTMFiniteDifferencesKernel::createRTMGrid(RTMCube<RTMData_t, RTMDevPtr_t> **pGrid, bool createOnAcc)
{
    /**
     * RTMGrids are created respecting the process
     * area dimensions
     * */
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int blen = rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int nxe = processLimits.processArea.xEnd - processLimits.processArea.xStart;
    int nye = processLimits.processArea.yEnd - processLimits.processArea.yStart;
    int nze = nz + 2 * blen;

    *pGrid = new RTMCube<RTMData_t, RTMDevPtr_t>(nxe, nye, nze);
    (*pGrid)->setBorderLength(blen);
    (*pGrid)->setStartX(processLimits.processArea.xStart);
    (*pGrid)->setStartY(processLimits.processArea.yStart);
    (*pGrid)->setEndX(processLimits.processArea.xEnd);
    (*pGrid)->setEndY(processLimits.processArea.yEnd);

    if (createOnAcc)
    {
        if (isUsingAcc())
        {
            (*pGrid)->createDeviceBuffer();
        }
    }
    (*pGrid)->reset();
}