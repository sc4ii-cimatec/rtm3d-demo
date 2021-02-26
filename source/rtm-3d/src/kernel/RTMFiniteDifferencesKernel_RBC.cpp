#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <Misc.hpp>
#include <RTMGrid.hpp>
#include <RTM.hpp>
#include <RTMController.hpp>
#include <RTMKernel.hpp>

using namespace std;

void RTMRBCKernel::initKernel()
{
    RTMFiniteDifferencesKernel::initKernel();
    RTM_PRINT("Initializing RTM RBC Kernel environment... ", rtmParam->verbose);

    report = new RTMKernelReport(rtmParam->outdir + "/rtmRBCKernelReport.txt");
    report->mname = rtmParam->mname;
    report->nx = rtmParam->nx;
    report->ny = rtmParam->ny;
    report->nz = rtmParam->nz;
    report->nt = rtmParam->nt;
}

void RTMRBCKernel::destroyKernel()
{
    RTMFiniteDifferencesKernel::destroyKernel();
    defaultPlatform->destroyRTMPlatform();
}

/** 
 * Performs an RTM Migration for given shot descriptor.
 * Shot Image is stored into descriptor's 'RTMGrid <T>'
 * */
void RTMRBCKernel::rtmMigrate(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor,
                              RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    timepoint t0 = tic();
#ifdef RTM_ACC_FPGA
    setCpuPlatform();
#endif
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_FORWARD;
    rtmRBCForward(shotDescriptor, v2dt2Grid);
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_BACKWARD;
    rtmRBCBackward(shotDescriptor, v2dt2Grid);
    KERNEL_STATE = RTMKernelState::RTM_KERNEL_STATE_IDLE;
    timepoint t1 = tic();
    report->rtmMigrationTime += elapsed_s(t0, t1);
    report->rtmMigrationCounter++;
}


void RTMRBCKernel::rtmRBCForward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor, 
                                RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    RTM_PRINT("RTM RBC Forward Propagation:", rtmParam->verbose);
    timepoint t0 = tic(), t1;
    int itstep=0, lt=0, it=0, ix=0, iy=0, iz=0;
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int hf_order = rtmParam->stencil_order/2;
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
                         "x" + to_string(nze) + ".snapshot.RBC_FWD";
        /* touch and create a new file */
        ofstream ofs(snapshotsFile, ios::out | ios::binary);
        ofs.close();

        RTM_PRINT("Saving RTM_FWD snapshots file at " + snapshotsFile, rtmParam->verbose);
        /* touch and create a new file */
    } 
    // contour values are zero
	pSrcGrid->reset();
    ppSrcGrid->reset();
    string upbFile;
	for (itstep = 0; itstep < nt; itstep += ntstep)
	{
	    // Forward propagation 
		for (lt = 0; lt < ntstep; lt++)
		{
            // update global time step
			it = itstep + lt; 
            
            // swap pointers
            RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);

            // propagate wave for one timestep (t+dt)
            timepoint t2 = tic();
            defaultPlatform->rtmStep(pSrcGrid, ppSrcGrid, stencil, v2dt2Grid);
            
            // update propagation function time report
            report->propagFuncTime += elapsed_s(t2, toc());
            report->propagFuncCounter++;

            // apply source
            defaultPlatform->rtmApplySource(ppSrcGrid, shotDescriptor.getSource(), it);
            
            // update snapshots file
            if (rtmParam->save_snapshots && (it % rtmParam->snapshot_step) == 0)
            {
                /* For now, snapshots are always on X dimension*/
                RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = ppSrcGrid->get2DSection(sx+blen, RTMDim::Xdim);
                secGrid->appendTofile(snapshotsFile);
                delete secGrid;
            }
            // join all distributed grids, if any.
            joinDistributedPSGrids();

            // print kernel progress...
            printKernelProgress("+ RBC_FWD", sx, sy, sz, it, nt, elapsed_s(t0,toc()));
		}// for (lt...
	}// for (itstep...
    if (isUsingAcc()){
        pSrcGrid->moveFromDevice();
        ppSrcGrid->moveFromDevice();
    }
    // saves the last two pressure grids for backward propagation
    snap0->copyData(*pSrcGrid);
    snap1->copyData(*ppSrcGrid);
    
    // update rtmForwardTime report 
    t1 = tic();
    report->rtmForwardTime += elapsed_s(t0,t1);
    report->rtmForwardCounter++;
}

void RTMRBCKernel::rtmRBCBackward(RTMShotDescriptor<RTMData_t, RTMDevPtr_t> &shotDescriptor, 
                                RTMVelocityModel<RTMData_t, RTMDevPtr_t> &v2dt2Grid)
{
    RTM_PRINT("RTM RBC Backward Propagation:", rtmParam->verbose);
    timepoint t0 = tic(), t1;
    int itstep=0, lt=0, it=0, ix=0, iy=0, iz=0;
    int nx = rtmParam->nx;
    int ny = rtmParam->ny;
    int nz = rtmParam->nz;
    int nt = rtmParam->nt;
    int ntstep = rtmParam->ntstep;
    int blen = rtmParam->blen;
    int receiver_depth = rtmParam->receiver_depth_z+rtmParam->blen;
    int st_order = rtmParam->stencil_order;
    int hf_order = rtmParam->stencil_order/2;
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
                         "x" + to_string(nze) + ".snapshot.RBC_BWD";
        /* touch and create a new file */
        ofstream ofs(snapshotsFile, ios::out | ios::binary);
        ofs.close();

        RTM_PRINT("Saving RTM_FWD snapshots file at " + snapshotsFile, rtmParam->verbose);
        /* touch and create a new file */
    }
    // reset pressure grids
    pSrcGrid->reset();
    ppSrcGrid->reset();
	pRcvGrid->reset();
    ppRcvGrid->reset();
    imgGrid->reset();
    string seismFile;
    for (itstep = nt-ntstep; itstep >=0; itstep -= ntstep)
	{        
        // loads rcv grid from file
        RTMReceiverGrid<RTMData_t, RTMDevPtr_t>* rcvGrid = shotDescriptor.getReceiverGrid();
        RTM_SEISMOGRAM_NAME(seismFile, rtmParam->datdir, sx, sy, sz,
                                rtmParam->receiver_count_x, rtmParam->receiver_count_y, 
                                itstep, ntstep);
        rcvGrid->loadFromFile(seismFile);
        if (isUsingAcc()){
            rcvGrid->moveToDevice();
        }
		for (lt = ntstep-1; lt >=0; lt--)
		{
            // update global time step
			it = itstep + lt;

            // Reconstruct source wavefield 
            if(it==(nt-1)){
                ppSrcGrid->copyData(*snap1);
                if (isUsingAcc()){
                    ppSrcGrid->moveToDevice();
                }
            }else if(it==(nt-2)){
                ppSrcGrid->copyData(*snap0);
                if (isUsingAcc()){
                    ppSrcGrid->moveToDevice();
                }
            }else{
                // propagate source wave for one timestep (t+dt)
                t1 = tic();
                defaultPlatform->rtmStep(pSrcGrid, ppSrcGrid, stencil, v2dt2Grid);
            
                // update propagation function time report
                report->propagFuncTime += elapsed_s(t1, toc());
                report->propagFuncCounter++;
            }
			// swap pointers
            RTMGRID_SWAP(&pSrcGrid, &ppSrcGrid);
            
            // propagate receiver wave for one timestep (t-dt)
            t1 = tic();
            defaultPlatform->rtmStep(pRcvGrid, ppRcvGrid, stencil, v2dt2Grid);

            // update propagation function time report
            report->propagFuncTime += elapsed_s(t1, toc());
            report->propagFuncCounter++;
		
            // restore receiver energy
            defaultPlatform->rtmRestoreReceiverData(ppRcvGrid, rcvGrid, lt);
			
            // apply image condition */
            defaultPlatform->rtmImageCondition(imgGrid,pSrcGrid,ppRcvGrid);

            // update snapshots file
            if (rtmParam->save_snapshots && (it % rtmParam->snapshot_step) == 0)
            {
                /* For now, snapshots are always on X dimension*/
                ppRcvGrid->moveFromDevice();
                RTMPlane<RTMData_t, RTMDevPtr_t> *secGrid = ppRcvGrid->get2DSection(sx+blen, RTMDim::Xdim);
                secGrid->appendTofile(snapshotsFile);
                delete secGrid;
            }
            // swap receiver pointers
			RTMGRID_SWAP(&pRcvGrid, &ppRcvGrid);
            
            // join distributed PS and PR grids, if any.
            joinDistributedPSGrids();
            joinDistributedPRGrids();
			
            // print kernel progress
            printKernelProgress("+ RBC_BWD", sx, sy, sz, nt-it, nt, elapsed_s(t0,toc()));
		}
    }
    if (isUsingAcc()){
        imgGrid->moveFromDevice();
    }
    // stack migrated image onto shot img grid
    shotDescriptor.getShotImage().copyData(*imgGrid);

    // update rtmBackwardTime report
    t1 = tic();
    report->rtmBackwardTime += elapsed_s(t0,t1);
    report->rtmBackwardCounter++;
}