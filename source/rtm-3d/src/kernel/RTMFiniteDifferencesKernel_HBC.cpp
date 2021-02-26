#include <cstdlib>
#include <cassert>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <Misc.hpp>
#include <RTMGrid.hpp>
#include <RTM.hpp>
#include <RTMController.hpp>
#include <RTMKernel.hpp>

using namespace std;

void RTMHBCKernel::initKernel()
{
    RTMFiniteDifferencesKernel::initKernel();
    RTM_PRINT("Initializing RTM HBC Kernel environment... ", rtmParam->verbose);

    if (!hbcKernelInitialized)
    {
        int nx = rtmParam->nx;
        int ny = rtmParam->ny;
        int nz = rtmParam->nz;
        int nt = rtmParam->nt;
        int ntstep = rtmParam->ntstep;
        int blen = rtmParam->blen;
        int st_order = rtmParam->stencil_order;
        int nxe = processLimits.processArea.xEnd - processLimits.processArea.xStart;
        int nye = processLimits.processArea.yEnd - processLimits.processArea.yStart;
        int nze = nz + 2 * blen;

        // upper-border grid is tricky...
        upbGrid = new RTMGridCollection<RTMData_t, RTMDevPtr_t>(ntstep, nxe, nye, (st_order / 2));
        upbGrid->setBorderLength(blen);
        upbGrid->reset();
        if (isUsingAcc()){
            upbGrid->createDeviceBuffer();
        }
        report = new RTMKernelReport(rtmParam->outdir + "/rtmHBCKernelReport.txt");
        report->mname = rtmParam->mname;
        report->nx = rtmParam->nx;
        report->ny = rtmParam->ny;
        report->nz = rtmParam->nz;
        report->nt = rtmParam->nt;
        hbcKernelInitialized = true;
    }
}

void RTMHBCKernel::destroyKernel()
{
    if(hbcKernelInitialized){
        RTMFiniteDifferencesKernel::destroyKernel();
        if(upbGrid!=nullptr)delete upbGrid;
        defaultPlatform->destroyRTMPlatform();
        hbcKernelInitialized = false;
    }
}
/** 
 * Performs an RTM Migration for given shot descriptor.
 * Shot Image is stored into descriptor's 'RTMGrid <T>'
 *  
*/
void RTMHBCKernel::rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                              RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{

    timepoint t0 = tic();
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_FORWARD;
#ifdef RTM_ACC_FPGA
    timepoint fwd_st=tic();
    RTMFPGAPlatform * fpgaPlatform = dynamic_cast<RTMFPGAPlatform *>(accPlatform);
    fpgaPlatform->rtmForwardPropagation(&shotDescriptor, stencil, rtmTaper, 
    v2dt2Grid, snap0, snap1, upbGrid);
    report->rtmForwardTime += elapsed_s(fwd_st, toc());
    report->rtmForwardCounter++;
#else
    rtmHBCForward(shotDescriptor, v2dt2Grid);
#endif

    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_BACKWARD;

    /**!!!! Set CPU Platform for Backward Propagation !!!!**/
    setCpuPlatform();
    /**!!!! Set CPU Platform for Backward Propagation !!!!**/

    rtmHBCBackward(shotDescriptor, v2dt2Grid);
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_IDLE;
    timepoint t1 = tic();
    report->rtmMigrationTime += elapsed_s(t0, t1);
    report->rtmMigrationCounter++;
}

void RTMHBCKernel::rtmHBCForward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                 RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    RTM_PRINT("RTM HBC Forward Propagation:", rtmParam->verbose);
    timepoint t0 = tic(), t1;
    char pmsg[256];
    int itstep = 0, lt = 0, it = 0, ix = 0, iy = 0, iz = 0;
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
    RTMData_t dt2 = rtmParam->dt * rtmParam->dt;
    /* source position does not consider the extended grid */
    int sx = shotDescriptor.getSource()->getX();
    int sy = shotDescriptor.getSource()->getY();
    int sz = shotDescriptor.getSource()->getZ();

    string snapshotsFile;
    if (rtmParam->save_snapshots)
    {
        int snap_step = rtmParam->snapshot_step;
        snapshotsFile += rtmParam->tmpdir + "/RTMSNAPSHOT_S" + to_string(sx) + "x" + to_string(sy) +
                         "_" + to_string(nt / snap_step) + "x" + to_string(nye) +
                         "x" + to_string(nze) + ".snapshot.HBC_FWD";
        /* touch and create a new file */
        ofstream ofs(snapshotsFile, ios::out | ios::binary);
        ofs.close();

        RTM_PRINT("Saving RTM_FWD snapshots file at " + snapshotsFile, rtmParam->verbose);
        /* touch and create a new file */
    }

    pSrcGrid->reset();
    ppSrcGrid->reset();

    string upbFile;
    /* Forward propagation */
    for (itstep = 0; itstep < nt; itstep += ntstep)
    {
        for (lt = 0; lt < ntstep; lt++)
        {
            // update global time step
            it = itstep + lt;

            // // swap pointers
            RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);
            // attenuate borders
            defaultPlatform->rtmTaperUpperBorders(pSrcGrid, rtmTaper);
            defaultPlatform->rtmTaperUpperBorders(ppSrcGrid, rtmTaper);

            // propagate wave for one timestep (t+dt)
            timepoint t2 = tic();
            defaultPlatform->rtmStep(pSrcGrid, ppSrcGrid, stencil, v2dt2Grid);

            // update propagation function time report
            report->propagFuncTime += elapsed_s(t2, toc());
            report->propagFuncCounter++;

            // apply source
            defaultPlatform->rtmApplySource(ppSrcGrid, shotDescriptor.getSource(), it);

            // save upper border
            defaultPlatform->rtmSaveUpperBorder(ppSrcGrid, upbGrid, lt);

            // update snapshots file
            if (rtmParam->save_snapshots && (it % rtmParam->snapshot_step) == 0)
            {
                /* For now, snapshots are always on X dimension*/
                if (isUsingAcc()){
                    ppSrcGrid->moveFromDevice();
                }
                RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = ppSrcGrid->get2DSection(sx + blen, RTMDim::Xdim);
                secGrid->appendTofile(snapshotsFile);
                delete secGrid;
            }
            // join all distributed grids, if any.
            joinDistributedPSGrids();

            // print kernel progress...
            sprintf(pmsg, "+[P%d] HBC_FWD", processLimits.pRank);
            printKernelProgress(pmsg, sx, sy, sz, it, nt, elapsed_s(t0, toc()));
            //printKernelProgress("+ HBC_FWD", sx, sy, sz, it, nt, elapsed_s(t0, toc()));
        } // for (lt...
        if (isUsingAcc()){
            upbGrid->moveFromDevice();
        }
        // saves intermediary large upper-boundary files
        int unxe = upbGrid->getNX();
        int unye = upbGrid->getNY();
        RTM_HBCUPB_NAME(upbFile, rtmParam->tmpdir, sx, sy, sz, unxe, unye,
                        hf_order, ntstep, itstep, processRank, nProcesses);
        upbGrid->saveToFile(upbFile);
    } // for (itstep...
    if (isUsingAcc()){
        // move from acc
        pSrcGrid->moveFromDevice();
        ppSrcGrid->moveFromDevice();
    }
    // saves the last two pressure grids for backward propagation
    snap0->copyData(*pSrcGrid);
    snap1->copyData(*ppSrcGrid);

    // update rtmForwardTime report
    t1 = tic();
    report->rtmForwardTime += elapsed_s(t0, t1);
    report->rtmForwardCounter++;
}

void RTMHBCKernel::rtmHBCBackward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                                  RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    RTM_PRINT("RTM HBC Backward Propagation:", rtmParam->verbose);
    timepoint t0 = tic(), t1;
    char pmsg[256];
    int itstep = 0, lt = 0, it = 0, ix = 0, iy = 0, iz = 0;
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int receiver_depth = rtmParam->receiver_depth_z + rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int hf_order = rtmParam->stencil_order / 2;
    int nxe = nx + 2 * blen;
    int nye = ny + 2 * blen;
    int nze = nz + 2 * blen;
    RTMData_t dt2 = rtmParam->dt * rtmParam->dt;

    // source position does not consider the extended grid
    int sx = shotDescriptor.getSource()->getX();
    int sy = shotDescriptor.getSource()->getY();
    int sz = shotDescriptor.getSource()->getZ();

    string snapshotsFile;
    if (rtmParam->save_snapshots)
    {
        int snap_step = rtmParam->snapshot_step;
        snapshotsFile += rtmParam->tmpdir + "/RTMSNAPSHOT_S" + to_string(sx) + "x" + to_string(sy) +
                         "_" + to_string(nt / snap_step) + "x" + to_string(nye) +
                         "x" + to_string(nze) + ".snapshot.HBC_BWD";
        /* touch and create a new file */
        ofstream ofs(snapshotsFile, ios::out | ios::binary);
        ofs.close();

        RTM_PRINT("Saving RTM_BWD snapshots file at " + snapshotsFile, rtmParam->verbose);
        /* touch and create a new file */
    }
    // reset pressure grids
    pSrcGrid->reset();
    ppSrcGrid->reset();
    pRcvGrid->reset();
    ppRcvGrid->reset();
    imgGrid->reset();

    string seismFile;
    string upbFile;
    for (itstep = nt - ntstep; itstep >= 0; itstep -= ntstep)
    {
        // loads upb grid from file
        int unxe = upbGrid->getNX();
        int unye = upbGrid->getNY();
        RTM_HBCUPB_NAME(upbFile, rtmParam->tmpdir, sx, sy, sz, unxe, unye, hf_order,
                        ntstep, itstep, processRank, nProcesses);
        upbGrid->loadFromFile(upbFile);

        // loads rcv grid from file
        RTM_SEISMOGRAM_NAME(seismFile, rtmParam->datdir, sx, sy, sz,
                            rtmParam->receiver_count_x, rtmParam->receiver_count_y,
                            itstep, ntstep);
        RTMReceiverGrid<RTMData_t, RTMDevPtr_t> *rcvGrid = shotDescriptor.getReceiverGrid();
        rcvGrid->loadFromFile(seismFile);
        if (isUsingAcc()){
            upbGrid->moveToDevice();
            rcvGrid->moveToDevice();
        }
        for (lt = ntstep - 1; lt >= 0; lt--)
        {
            // update global time step
            it = itstep + lt;
            // Reconstruct source wavefield
            if (it == (nt - 1)){
                ppSrcGrid->copyData(*snap1);
                if (isUsingAcc()){
                    ppSrcGrid->moveToDevice();
                }
            }else if (it == (nt - 2)){
                ppSrcGrid->copyData(*snap0);
                if (isUsingAcc()){
                    ppSrcGrid->moveToDevice();
                }
            }else
            {
                // propagate source wave for one timestep (t+dt)
                t1 = tic();
                defaultPlatform->rtmStep(pSrcGrid, ppSrcGrid, stencil, v2dt2Grid);

                // update propagation function time report
                report->propagFuncTime += elapsed_s(t1, toc());
                report->propagFuncCounter++;

                // restore upper border upb to src grid
                defaultPlatform->rtmRestoreUpperBorder(ppSrcGrid, upbGrid, lt);
            }
            // swap pointers
            RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);

            // attenuate receiver grid borders
            defaultPlatform->rtmTaperUpperBorders(pRcvGrid, rtmTaper);
            defaultPlatform->rtmTaperUpperBorders(ppRcvGrid, rtmTaper);

            // propagate receiver wave for one timestep (t-dt)
            t1 = tic();
            defaultPlatform->rtmStep(pRcvGrid, ppRcvGrid, stencil, v2dt2Grid);

            // update propagation function time report
            report->propagFuncTime += elapsed_s(t1, toc());
            report->propagFuncCounter++;

            // restore receiver energy
            defaultPlatform->rtmRestoreReceiverData(ppRcvGrid, rcvGrid, lt);

            // apply image condition */
            defaultPlatform->rtmImageCondition(imgGrid, pSrcGrid, ppRcvGrid);

            // update snapshots file
            if (rtmParam->save_snapshots && (it % rtmParam->snapshot_step) == 0)
            {
                if (isUsingAcc()){
                    ppRcvGrid->moveFromDevice();
                }
                /* For now, snapshots are always on X dimension*/
                RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = ppRcvGrid->get2DSection(sx + blen, RTMDim::Xdim);
                secGrid->appendTofile(snapshotsFile);
                delete secGrid;
            }
            // swap receiver pointers
            RTMGRID_SWAP(&pRcvGrid, &ppRcvGrid);
            
            // join distributed PS and PR grids, if any.
            joinDistributedPSGrids();
            joinDistributedPRGrids();

            // print kernel progress;
            sprintf(pmsg, "+[P%d] HBC_BWD", processLimits.pRank);
            printKernelProgress(pmsg, sx, sy, sz, nt - it, nt, elapsed_s(t0, toc()));
            //printKernelProgress("+ HBC_BWD", sx, sy, sz, nt - it, nt, elapsed_s(t0, toc()));
        }
        // remove border file from FS
        remove(upbFile.c_str());
    }
    // stack migrated image onto shot img grid
    if (isUsingAcc()){
        imgGrid->moveFromDevice();
    }
    shotDescriptor.getShotImage().copyData(*imgGrid);

    // update rtmBackwardTime report
    t1 = tic();
    report->rtmBackwardTime += elapsed_s(t0, t1);
    report->rtmBackwardCounter++;
}