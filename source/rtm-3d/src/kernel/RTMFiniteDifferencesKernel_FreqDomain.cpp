#include <cstdlib>
#include <cassert>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <Misc.hpp>
#include <RTMGrid.hpp>
#include <RTM.hpp>
#include <RTMController.hpp>
#include <RTMKernel.hpp>

using namespace std;

void RTMFreqDomainKernel::rtmInitWList(RTMVector<RTMData_t, RTMDevPtr_t> *wList)
{
    if (wList == NULL)
    {
        string s("[ ERROR: WList not allocated! ]");
        RTMException ex(s);
        throw ex;
    }
    int iw;
    float pi = acosf(-1.); //M_PI;
    for (iw = 0; iw < wList->size(); iw++)
    {
        float lw = (float)iw * 1.0;
        wList->set((2. * pi * ((lw * df) + freq_min) * rtmParam->dt), iw);
        //printf("***** Wlist[%d] = %.04f Hz, lw=%.2f \n", iw, wList->get(iw), lw);
    }
}

void RTMFreqDomainKernel::rtmInitFreqKernels(RTMPlane<RTMData_t, RTMDevPtr_t> *kRe, RTMPlane<RTMData_t, RTMDevPtr_t> *kIm)
{
    if (wList == NULL || kRe == NULL || kIm == NULL)
    {
        string s("[ ERROR: Freq Kernels not allocated! ]");
        RTMException ex(s);
        throw ex;
    }
    int it = 0, iw = 0;
    for (iw = 0; iw < nw; iw++)
    {
        for (it = 0; it < rtmParam->nt; it++)
        {
            kRe->set(cosf((*wList)[iw] * (it * 1.0)), iw, it);
            kIm->set(sinf((*wList)[iw] * (it * 1.0)), iw, it);
        }
    }
}

void RTMFreqDomainKernel::createFrequencyComponentsGrids()
{
    int inx = imagingLimits.processArea.xEnd - imagingLimits.processArea.xStart;
    int iny = imagingLimits.processArea.yEnd - imagingLimits.processArea.yStart;
    int inz = rtmParam->nz;
    int startX = imagingLimits.processArea.xStart;
    int startY = imagingLimits.processArea.yStart;
    int startZ = 0;
    int endX = imagingLimits.processArea.xEnd;
    int endY = imagingLimits.processArea.yEnd;
    int endZ = rtmParam->nz;

    srcReGrid = new RTMGridCollection<RTMData_t, RTMDevPtr_t>(nwstep, inx, iny, inz);
    srcImGrid = new RTMGridCollection<RTMData_t, RTMDevPtr_t>(nwstep, inx, iny, inz);
    rcvReGrid = new RTMGridCollection<RTMData_t, RTMDevPtr_t>(nwstep, inx, iny, inz);
    rcvImGrid = new RTMGridCollection<RTMData_t, RTMDevPtr_t>(nwstep, inx, iny, inz);

    srcReGrid->setStartX(startX);
    srcReGrid->setEndX(endX);
    srcReGrid->setStartY(startY);
    srcReGrid->setEndY(endY);
    srcReGrid->setStartZ(startZ);
    srcReGrid->setEndZ(endZ);
    srcImGrid->setStartX(startX);
    srcImGrid->setEndX(endX);
    srcImGrid->setStartY(startY);
    srcImGrid->setEndY(endY);
    srcImGrid->setStartZ(startZ);
    srcImGrid->setEndZ(endZ);

    rcvReGrid->setStartX(startX);
    rcvReGrid->setEndX(endX);
    rcvReGrid->setStartY(startY);
    rcvReGrid->setEndY(endY);
    rcvReGrid->setStartZ(startZ);
    rcvReGrid->setEndZ(endZ);
    rcvImGrid->setStartX(startX);
    rcvImGrid->setEndX(endX);
    rcvImGrid->setStartY(startY);
    rcvImGrid->setEndY(endY);
    rcvImGrid->setStartZ(startZ);
    rcvImGrid->setEndZ(endZ);

    if (isUsingAcc())
    {
        // create stencil and taper buffers on device
        srcReGrid->createDeviceBuffer();
        srcImGrid->createDeviceBuffer();
        rcvReGrid->createDeviceBuffer();
        rcvImGrid->createDeviceBuffer();
    }
}

void RTMFreqDomainKernel::initKernel()
{
    RTMFiniteDifferencesKernel::initKernel();
    if (!freqDomainKernelInitialized)
    {
        RTM_PRINT("Initializing RTM Freqquency Domain Kernel environment... ",
                  rtmParam->verbose);
        setDistributedImagingLimits();

        int nx = rtmParam->nx;
        int ny = rtmParam->ny;
        int nz = rtmParam->nz;
        int nt = rtmParam->nt;
        int ntstep = rtmParam->ntstep;
        int blen = rtmParam->blen;
        int st_order = rtmParam->stencil_order;
        float dt = rtmParam->dt;

        freq_max = 1.0 * rtmParam->fmig_max_freq;
        freq_min = 1.0 * rtmParam->fmig_min_freq;

        // calc df and nw
        df = 1. / ((nt * 1.0) * dt);
        float nf = ((freq_max - freq_min) / df);
        nw = static_cast<size_t>(nf);

        if (isUsingAcc())
        {
            nwstep = rtmParam->fmig_nwstep;
        }
        else
        {
            // CPU builds don't need to partition
            // nwstep parameter.
            nwstep = nw;
        }
        // define freq values;
        wList = new RTMVector<RTMData_t, RTMDevPtr_t>(nw);
        rtmInitWList(wList);
        // create and init freq kernels
        kernelIm = new RTMPlane<RTMData_t, RTMDevPtr_t>(nw, nt);
        kernelRe = new RTMPlane<RTMData_t, RTMDevPtr_t>(nw, nt);
        rtmInitFreqKernels(kernelRe, kernelIm);

        // square wList
        wList->power2();

        // create freq 4D grids and allocate then on acc, if any...
        createFrequencyComponentsGrids();

        if (isUsingAcc())
        {
            // create stencil and taper buffers on device
            wList->createDeviceBuffer();
            wList->moveToDevice();
            kernelRe->createDeviceBuffer();
            kernelRe->moveToDevice();
            kernelIm->createDeviceBuffer();
            kernelIm->moveToDevice();
        }

        report = new RTMKernelReport(rtmParam->outdir + "/rtmHBCKernelReport.txt");
        report->mname = rtmParam->mname;
        report->nx = rtmParam->nx;
        report->ny = rtmParam->ny;
        report->nz = rtmParam->nz;
        report->nt = rtmParam->nt;
        report->nw = nw;
        freqDomainKernelInitialized = true;
    }
}

void RTMFreqDomainKernel::destroyKernel()
{
    if (freqDomainKernelInitialized)
    {
        RTMFiniteDifferencesKernel::destroyKernel();
        defaultPlatform->destroyRTMPlatform();
        if (srcReGrid != nullptr)
            delete srcReGrid;
        if (srcImGrid != nullptr)
            delete srcImGrid;
        if (rcvReGrid != nullptr)
            delete rcvReGrid;
        if (rcvImGrid != nullptr)
            delete rcvImGrid;
        if (kernelRe != nullptr)
            delete kernelRe;
        if (kernelIm != nullptr)
            delete kernelIm;
        if (wList != nullptr)
            delete wList;
        freqDomainKernelInitialized = false;
    }
}
/** 
 * Performs an RTM Migration for given shot descriptor.
 * Shot Image is stored into descriptor's 'RTMGrid <T>'
 *  
*/
void RTMFreqDomainKernel::rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                     RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{

    timepoint t0 = tic();
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_FORWARD;
    rtmFreqDomainPropagation(shotDescriptor, v2dt2Grid);
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_IDLE;
    timepoint t1 = tic();
    report->rtmMigrationTime += elapsed_s(t0, t1);
    report->rtmMigrationCounter++;
}

void RTMFreqDomainKernel::rtmFreqDomainPropagation(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                                   RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    timepoint t0 = tic(), t1, t2;
    int itstep = 0, lt = 0, st = 0, rt = 0, gw = 0, lw = 0, iw = 0, ix = 0, iy = 0, iz = 0;
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int hf_order = rtmParam->stencil_order / 2;
    int nxe = nx + 2 * blen;
    int nye = ny + 2 * blen;
    int nze = nz + 2 * blen;
    /* source position does not consider the extended grid */
    int sx = shotDescriptor.getSource()->getX();
    int sy = shotDescriptor.getSource()->getY();
    int sz = shotDescriptor.getSource()->getZ();

    char dbgMsg[256];
    sprintf(dbgMsg, "\tRTM Frequency Domain Propagation (NW=%d; DF=%.5f):", nw, df);
    RTM_PRINT(dbgMsg, rtmParam->verbose);

    string snapshotsFile;
    if (rtmParam->save_snapshots)
    {
        int snap_step = rtmParam->snapshot_step;
        snapshotsFile += rtmParam->tmpdir + "/RTMSNAPSHOT_S" + to_string(sx) + "x" + to_string(sy) +
                         "_" + to_string(nt / snap_step) + "x" + to_string(nye) +
                         "x" + to_string(nze) + ".snapshot.FREQIMG";
        /* touch and create a new file */
        ofstream ofs(snapshotsFile, ios::out | ios::binary);
        ofs.close();

        RTM_PRINT("Saving RTM_FWD snapshots file at " + snapshotsFile, rtmParam->verbose);
        /* touch and create a new file */
    }
    string seismFile;
    for (gw = 0; gw < nw; gw += nwstep)
    {
        // reset all grids
        pSrcGrid->reset();
        ppSrcGrid->reset();
        pRcvGrid->reset();
        ppRcvGrid->reset();
        srcReGrid->reset();
        srcImGrid->reset();
        rcvReGrid->reset();
        rcvImGrid->reset();
        /* rtm propagation */
        t1 = tic();
        sprintf(dbgMsg, "\t+ Freq[%d/%d (nwstep=%d)] = %.2f Hz (%.2f rad/s);",
                gw + 1, nw, nwstep, freq_min + gw * df, wList->get(gw) / (wList->get(gw) > 0 ? wList->get(gw) : 1));
        RTM_PRINT(dbgMsg, rtmParam->verbose);
        for (itstep = nt - ntstep; itstep >= 0; itstep -= ntstep)
        {
            // loads rcv grid from file
            RTM_SEISMOGRAM_NAME(seismFile, rtmParam->datdir, sx, sy, sz,
                                rtmParam->receiver_count_x, rtmParam->receiver_count_y,
                                itstep, ntstep);
            RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid = shotDescriptor.getReceiverGrid();
            rcvGrid->loadFromFile(seismFile);
            if (isUsingAcc())
            {
                rcvGrid->moveToDevice();
            }
            for (lt = ntstep - 1; lt >= 0; lt--)
            {
                // update global time step
                rt = itstep + lt;
                st = nt - 1 - rt;
                timepoint t2 = tic();
                rtmUpdatePressureGrids(v2dt2Grid);
                // update propagation function time report
                report->propagFuncTime += elapsed_s(t2, toc());
                report->propagFuncCounter++;

                // apply source
                defaultPlatform->rtmApplySource(ppSrcGrid, shotDescriptor.getSource(), st);

                // restore receiver energy
                defaultPlatform->rtmRestoreReceiverData(ppRcvGrid, rcvGrid, lt);

                // update all freqs. contribution
                for (lw = 0, iw = gw; (lw < nwstep && iw < nw); lw++, iw++)
                {
                    t2 = tic();
                    defaultPlatform->rtmUpdateFreqContributions(st, iw, lw, ppSrcGrid, ppRcvGrid,
                                                         srcReGrid, srcImGrid, rcvReGrid, rcvImGrid, kernelRe, kernelIm);
                    // freq domain image condition
                    defaultPlatform->rtmFreqDomainImageCondition(iw, lw, wList,
                                                          imgGrid, srcReGrid, srcImGrid, rcvReGrid, rcvImGrid);
                    //printf(">> Lw[%d]: %.5f ms \n", lw, elapsed_ms(t2, toc()) );
                }

                // update snapshots file
                if (rtmParam->save_snapshots && (st % rtmParam->snapshot_step) == 0)
                {
                    /// For now, snapshots are always on X dimension
                    if (isUsingAcc())
                    {
                        imgGrid->moveFromDevice();
                    }
                    RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = imgGrid->get2DSection(sx + blen, RTMDim::Xdim);
                    secGrid->appendTofile(snapshotsFile);
                    delete secGrid;
                }

                // swap pointers
                RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);
                RTMGRID_SWAP(&pRcvGrid, &ppRcvGrid);

                // print kernel progress...
                printKernelProgress("+ FMIG", sx, sy, sz, st, nt, elapsed_s(t1, toc()));
            } // for (lt...
        }     // for (itstep...
    }
    if (isUsingAcc())
    {
        imgGrid->moveFromDevice();
    }
    shotDescriptor.getShotImage().copyData(*imgGrid);

    if (rtmParam->fmig_distributed_imaging)
    {
        t2 = tic();
        stitchDistributedOutputImage(&shotDescriptor.getShotImage());
        report->mpiFuncTime += elapsed_s(t2, toc());
        report->mpiFuncCounter++;
    }

    // update rtmForwardTime report
    t1 = tic();
    report->rtmForwardTime += elapsed_s(t0, t1);
    report->rtmForwardCounter++;
}

void RTMFreqDomainKernel::rtmUpdatePressureGrids(RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    // propagate wave for one timestep (t+dt)
    defaultPlatform->rtmStepMultipleWave(pSrcGrid, ppSrcGrid, pRcvGrid, ppRcvGrid, stencil, v2dt2Grid);

    // attenuate borders
    defaultPlatform->rtmTaperAllBorders(pSrcGrid, rtmTaper);
    defaultPlatform->rtmTaperAllBorders(ppSrcGrid, rtmTaper);
    defaultPlatform->rtmTaperAllBorders(pRcvGrid, rtmTaper);
    defaultPlatform->rtmTaperAllBorders(ppRcvGrid, rtmTaper);
}

void RTMFreqDomainKernel::setDistributedImagingLimits()
{

    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int nxe = nx;
    int nye = ny;
    int nze = nz;

    if (rtmParam->fmig_distributed_imaging)
    {
        int ySections = 1;
        int xSections = 1;
        int res = 0;
        int procs = nProcesses;

        do
        {
            res = procs % 2;
            procs /= 2;
            if (res == 0)
                ySections *= 2;
        } while (res == 0);

        xSections = nProcesses / ySections;
        if (xSections == 1 && ySections > 2 && (ySections % 2 == 0))
        {
            xSections = 2;
            ySections = ySections / 2;
        }

        int x, y;
        int xLen = ((nxe) / xSections) + 1;
        int yLen = ((nye) / ySections) + 1;

        int xStart = 0;
        int xEnd = xStart + xLen;
        int yStart = 0;
        int yEnd = yStart + yLen;
        int pCount = 0;
        RTMProcessLimits procLimit;
        for (y = 0; y < ySections; y++)
        {
            for (x = 0; x < xSections; x++)
            {
                procLimit.shadowLength = 0;
                procLimit.pRank = pCount;
                procLimit.lRank = pCount;
                procLimit.nProcesses = nProcesses;
                procLimit.processArea.xStart = xStart;
                procLimit.processArea.yStart = yStart;
                procLimit.processArea.xEnd = xEnd;
                procLimit.processArea.yEnd = yEnd;
                procLimit.validProcessArea = true;
                imagingLimits.validZones = 0;
                if (procLimit.processArea.yStart >= nye)
                {
                    procLimit.processArea.yStart = nye;
                    procLimit.validProcessArea = false;
                }
                if (procLimit.processArea.xStart >= nxe)
                {
                    procLimit.processArea.xStart = nxe;
                    procLimit.validProcessArea = false;
                }
                if (procLimit.processArea.yEnd > nye)
                    procLimit.processArea.yEnd = nye;
                if (procLimit.processArea.xEnd > nxe)
                    procLimit.processArea.xEnd = nxe;

                if (pCount == processLimits.pRank)
                {
                    imagingLimits = procLimit;
                    imagingLimits.lRank = processLimits.lRank;
                    /*****        START OF SEQUENTIAL ZONE     *****/
                    /**/
                    // RTMProcess::_enterSequentialZone();

                    printf(">************************************************* \n");
                    printf(">* [nxe=%d nye=%d] xSections = %d xLen = %d | ySections = %d yLen = %d \n",
                           nxe, nye, xSections, xLen, ySections, yLen);
                    printf(">* Proc=%02d xSt=%03d xEnd=%03d      ySt=%03d yEnd=%03d \n", pCount,
                           procLimit.processArea.xStart, procLimit.processArea.xEnd,
                           procLimit.processArea.yStart, procLimit.processArea.yEnd);
                    printf(">************************************************* \n");

                    // RTMProcess::_leaveSequentialZone();

                    /*****        END OF SEQUENTIAL ZONE      *****/
                }
                xStart += xLen;
                xEnd += xLen;
                pCount++;
            }
            yStart += yLen;
            yEnd += yLen;
            xStart = 0;
            xEnd = xStart + xLen;
        }
    }
    else
    {
        imagingLimits.shadowLength = 0;
        imagingLimits.nProcesses = nProcesses;
        imagingLimits.processArea.xStart = 0;
        imagingLimits.processArea.yStart = 0;
        imagingLimits.processArea.xEnd = nx;
        imagingLimits.processArea.yEnd = ny;
        imagingLimits.validProcessArea = true;
        imagingLimits.validZones = 0;
        imagingLimits.pRank = processLimits.pRank;
        imagingLimits.lRank = processLimits.lRank;
    }
}

void RTMFreqDomainKernel::stitchDistributedOutputImage(RTMCube<RTMData_t, RTMDevPtr_t> *outImage)
{
#ifdef RTM_MPI
    int coordinatesVec[6];
    int startX, endX, startY, endY, startZ, endZ;
    int pk, ix, iy, iz;
    int half_order = rtmParam->stencil_order / 2;
    int blen = rtmParam->blen;
    MPI_Status mpiStatus;
    int processRank, nProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);  // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank); //the rank of the process

    if (processRank == 0)
    {
        for (pk = 1; pk < nProcesses; pk++)
        {
            //printf("---> P%d requested IMG from R%d \n", processRank, pk); fflush(stdout);
            // request process area
            int pTag = RTM_MPI_REQ_PAREA;
            MPI_Send(&pTag, 1, MPI_INT, pk, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);

            MPI_Recv(coordinatesVec, 6, MPI_INT, pk, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD, &mpiStatus);

            if (mpiStatus.MPI_ERROR)
            {
                PRINT_MPI_STATUS(processRank, __func__, mpiStatus);
            }
            startX = coordinatesVec[0];
            endX = coordinatesVec[1];
            startY = coordinatesVec[2];
            endY = coordinatesVec[3];
            startZ = coordinatesVec[4];
            endZ = coordinatesVec[5];
            if (startX == -1 || startY == -1 || startZ == -1 ||
                endX == -1 || endY == -1 || endZ == -1)
            {
                //invalid area
                //printf("---> P%d received invalid area from R%d \n", processRank, pk); fflush(stdout);
                continue;
            }
            else
            {
                unsigned int length = (endX - startX) * (endY - startY) * (endZ - startZ);
                RTMData_t *gData = new RTMData_t[length];
                pTag = RTM_MPI_OUTIMG_TAG;
                MPI_Send(&pTag, 1, MPI_INT, pk, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD);
                MPI_Recv(gData, length, MPI_FLOAT, pk, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD, &mpiStatus);
                if (mpiStatus.MPI_ERROR)
                {
                    PRINT_MPI_STATUS(processRank, __func__, mpiStatus);
                }

                int k0 = 0;
                for (ix = startX; ix < endX; ix++)
                {
                    for (iy = startY; iy < endY; iy++)
                    {
                        for (iz = startZ; iz < endZ; iz++)
                        {
                            RTMData_t remoteVal = gData[k0++];
                            outImage->set((remoteVal), ix + blen, iy + blen, iz + blen);
                        }
                    }
                }
                delete gData;
            }
            //printf("---> P%d received IMG from R%d \n", processRank, pk); fflush(stdout);
        }
    }
    else
    {
        int startX = 0, endX = 0, startY = 0, endY = 0;
        RTMData_t *gData = NULL;
        unsigned int length = 0;
        int pTag = 0;
        MPI_Recv(&pTag, 1, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD, &mpiStatus);
        if (mpiStatus.MPI_ERROR)
        {
            PRINT_MPI_STATUS(processRank, __func__, mpiStatus);
        }
        //printf("<--- P%d received request IMG \n", processRank); fflush(stdout);
        if (imagingLimits.validProcessArea)
        {
            startX = coordinatesVec[0] = imagingLimits.processArea.xStart;
            endX = coordinatesVec[1] = imagingLimits.processArea.xEnd;
            startY = coordinatesVec[2] = imagingLimits.processArea.yStart;
            endY = coordinatesVec[3] = imagingLimits.processArea.yEnd;
            startZ = coordinatesVec[4] = 0;
            endZ = coordinatesVec[5] = rtmParam->nz;
            length = (endX - startX) * (endY - startY) * (endZ - startZ);
            gData = new RTMData_t[length];
            int k0 = 0;
            for (ix = startX; ix < endX; ix++)
            {
                for (iy = startY; iy < endY; iy++)
                {
                    for (iz = startZ; iz < endZ; iz++)
                    {
                        gData[k0++] = outImage->get(ix + blen, iy + blen, iz + blen);
                    }
                }
            }
            // send process area
            MPI_Send(coordinatesVec, 6, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);

            MPI_Recv(&pTag, 1, MPI_INT, 0, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD, &mpiStatus);
            if (mpiStatus.MPI_ERROR)
            {
                PRINT_MPI_STATUS(processRank, __func__, mpiStatus);
            }
            MPI_Send(gData, length, MPI_FLOAT, 0, RTM_MPI_OUTIMG_TAG, MPI_COMM_WORLD);
            delete gData;
            //printf("<--- P%d delivered IMG \n", processRank); fflush(stdout);
        }
        else
        {
            // send invalid process area
            coordinatesVec[0] = -1;
            coordinatesVec[1] = -1;
            coordinatesVec[2] = -1;
            coordinatesVec[3] = -1;
            coordinatesVec[4] = -1;
            coordinatesVec[5] = -1;
            MPI_Send(coordinatesVec, 6, MPI_INT, 0, RTM_MPI_REQ_PAREA, MPI_COMM_WORLD);
            //printf("<--- P%d sent invalid area\n", processRank); fflush(stdout);
        }
    }
#endif
}