#include <cstdlib>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <RTM.hpp>


void RTMAcousticProcess::rtm(){
    try{
        RTM_PRINT("Initializing RTM Acoustic Process...", true);
        // init process limits and load input parameters
        initRTMProcess();
        // print input parameters
        rtmInputParam->printRTMParam();
        // create RTMKernel
        if(rtmInputParam->fmig){
            rtmKernel = dynamic_cast<RTMKernel*>(new RTMFreqDomainKernel(*rtmInputParam, pLimits)); 
        }else{
            if(getRTMBoundaryCondition(rtmInputParam->boundary_condition)==RTMBoundaryCondition::RBC){
                rtmKernel = dynamic_cast<RTMKernel*>(new RTMRBCKernel(*rtmInputParam, pLimits)); 
            }else{
                rtmKernel = dynamic_cast<RTMKernel*>(new RTMHBCKernel(*rtmInputParam, pLimits));
            }
        }
        // init Kernel and RTM Controller...
        rtmKernel->initKernel();
        
        // init RTMController
        rtmController = dynamic_cast<RTMController*> (new RTMAcousticController (*rtmInputParam, pLimits));
        rtmController->initRTMController();
        rtmController->setRTMKernel(*rtmKernel);
        
        // Running RTM Process... 
        omp_set_num_threads(rtmInputParam->nthreads);
        if(rtmInputParam->modeling){
            rtmController->runModelingProcess();
        }else if (rtmInputParam->migration){
            rtmController->runMigrationProcess();
        }

        // print RTMReport
        rtmKernel->getReport().calcAVGTime();
        if(rtmInputParam->verbose){
            cout << "> Time Report: " << endl;
            cout << rtmKernel->getReport().toString() << endl;
        }
    }catch(RTMException &e){
        cout << "> Error:\n" << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }catch(exception &e){
        cout << "> Error:\n" << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }catch(...){
        cout << endl;
        cout << "> Error: Rougue Exception Flying Around. Aborting!" << endl;
        exit(EXIT_FAILURE);
    }
}

