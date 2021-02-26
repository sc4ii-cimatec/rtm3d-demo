#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <Misc.hpp>
#include <RTMGrid.hpp>
#include <RTM.hpp>
#include <RTMController.hpp>
#include <RTMKernel.hpp>

string &RTMKernelReport::toString()
{
    string *str = new string("{");
    *str += "\n  \"mname\": " + (mname);
    *str += "\n  \"nx\":" + to_string(nx);
    *str += "\n  \"ny\":" + to_string(ny);
    *str += "\n  \"nz\":" + to_string(nz);
    *str += "\n  \"nt\":" + to_string(nt);
    *str += "\n  \"modeling_time\":" + to_string(rtmModelingTime);
    *str += "\n";
    *str += "\n  \"migration_time\":" + to_string(rtmMigrationTime);
    *str += "\n  \"forward_time\":" + to_string(rtmForwardTime);
    *str += "\n  \"backward_time\":" + to_string(rtmBackwardTime);
    *str += "\n  \"propagfunc_time\":" + to_string(propagFuncTime);
    *str += "\n  \"propagfunc_percentual\":" + to_string(propagFuncPercentual);
    *str += "\n";
    *str += "\n  \"mpi_time: \":" + to_string(mpiFuncTime);
    *str += "\n  \"mpi_percentual: \":" + to_string(mpiFuncPercentual);
    *str += "\n";
    *str += "\n  \"forward_avg_time\":" + to_string(rtmForwardAVGTime);
    *str += "\n  \"backward_avg_time\":" + to_string(rtmBackwardAVGTime);
    *str += "\n  \"migration_avg_time\":" + to_string(rtmMigrationAVGTime);
    *str += "\n  \"modeling_avg_time\":" + to_string(rtmModelingAVGTime);
    *str += "\n  \"propagfunc_avg_time\":" + to_string(propagFuncAVGTime);
    *str += "\n  \"mpi_avg_time: \":" + to_string(mpiFuncAVGTime);
    
    *str += "\n}\n";
    return *str;
}

void RTMKernelReport::calcAVGTime()
{

    rtmModelingCounter  = rtmModelingCounter != 0 ? rtmModelingCounter : 1;
    rtmMigrationCounter  = rtmMigrationCounter != 0 ? rtmMigrationCounter : 1;
    rtmForwardCounter   = rtmForwardCounter != 0 ? rtmForwardCounter : 1;
    rtmBackwardCounter  = rtmBackwardCounter != 0 ? rtmBackwardCounter : 1;
    propagFuncCounter   = propagFuncCounter != 0 ? propagFuncCounter : 1;
    mpiFuncCounter      = mpiFuncCounter != 0 ? mpiFuncCounter : 1;

    rtmModelingAVGTime  = (rtmModelingTime) / (rtmModelingCounter * 1.0);
    rtmMigrationAVGTime = (rtmMigrationTime) / (rtmMigrationCounter * 1.0);
    rtmForwardAVGTime   = (rtmForwardTime) / (rtmForwardCounter * 1.0);
    rtmBackwardAVGTime  = (rtmBackwardTime) / (rtmBackwardCounter * 1.0);
    propagFuncAVGTime   = (propagFuncTime) / (propagFuncCounter * 1.0);
    mpiFuncAVGTime      = (mpiFuncTime)/(mpiFuncCounter*1.0);

    if (rtmMigrationTime > 0.0){
        propagFuncPercentual = 100.0*(propagFuncTime/rtmMigrationTime);
        mpiFuncPercentual = 100.0*(mpiFuncTime/rtmMigrationTime);
    }else if(rtmModelingTime>0.0){
        propagFuncPercentual = 100.0*(propagFuncTime/rtmModelingTime);
        mpiFuncPercentual = 100.0*(mpiFuncTime/rtmModelingTime);
    }else{
        propagFuncPercentual = 0.0;
        mpiFuncPercentual = 0.0;
    }
     
}

void RTMKernelReport::saveToFile()
{
    ofstream ofs(jsonFilePath);
    if (ofs.fail())
    {
        return;
    }
    string str = toString();
    ofs.write((str.c_str()), str.size() * sizeof(char));
}

void RTMKernel::printKernelProgress(char * kname, int sx, int sy, int sz, int it, int nt, 
    float time_s, int printStep){
    
    bool newl = (it==(nt-1) || it==nt);
    bool printCondition = (((it) % printStep == 0) || newl) && rtmParam->verbose;
    if (printCondition){
        char *progressMsg = new char[250];
        if(time_s < 600){// ten minutes
            sprintf(progressMsg, "\r>\t%s(%d,%d,%d) = %d/%d (%.1f%%) (%.2f s)   ",
                    kname, (sx), (sy), (sz),it, nt, 100.0*(it*1.0)/(nt*1.0), time_s);
        }else{
            sprintf(progressMsg, "\r>\t%s(%d,%d,%d) = %d/%d (%.1f%%) (%.2f min)   ",
                    kname, (sx), (sy), (sz),it, nt, 100.0*(it*1.0)/(nt*1.0), time_s/60.0);
        }
        cout << (string(progressMsg)) << flush;
        if (newl)
            cout << endl;
        delete progressMsg;
        fflush(stdout);
    }
}