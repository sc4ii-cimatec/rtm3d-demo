#ifndef _RTMDEVICEGRID_H
#define _RTMDEVICEGRID_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <Misc.hpp>

using namespace std;


template<typename HostPtr_type, typename DevPtr_type>
class DeviceGrid
{
protected:
    size_t          BUFFER_LENGTH; // in bytes
    int             deviceID;
    DevPtr_type     * DEV_PTR  = nullptr;
    HostPtr_type    * HOST_PTR = nullptr;
public:
    DeviceGrid(/* args */){}

    int getDeviceID(){
        return deviceID;
    }
    void setDeviceID(int _devID){
        deviceID = _devID;
    }
    virtual void createDeviceBuffer(uint32_t bufferFlags=0)=0;
    virtual void removeDeviceBuffer()=0;
    virtual void moveFromDevice()=0;
    virtual void moveToDevice()=0;
    virtual void devMemSet(HostPtr_type val) = 0;

    void hostMemSet(HostPtr_type val){
        if(HOST_PTR!=nullptr){
            void * ptr = reinterpret_cast<void*>(HOST_PTR);
            memset(ptr, val, BUFFER_LENGTH);
        }
    }

    DevPtr_type* getDevPtr() const{
        return DEV_PTR;
    }
    HostPtr_type* getHostPtr() const{
        return HOST_PTR;
    }

    void setDevPtr(DevPtr_type** _devPtr){
        DEV_PTR = *_devPtr;
    }
    void setHostPtr(HostPtr_type** _hostPtr){
        HOST_PTR = *_hostPtr;
    }
};

template<typename HostPtr_type, typename DevPtr_type=HostPtr_type>
class CPUDeviceGrid : public DeviceGrid<HostPtr_type,DevPtr_type>
{
public:
    CPUDeviceGrid(/* args */):DeviceGrid<HostPtr_type,DevPtr_type>{}{}

    virtual void createDeviceBuffer(uint32_t bufferFlags=0){}
    virtual void removeDeviceBuffer(){}
    virtual void moveFromDevice(){}
    virtual void moveToDevice(){}
    virtual void devMemSet(HostPtr_type val){}
};


#endif
