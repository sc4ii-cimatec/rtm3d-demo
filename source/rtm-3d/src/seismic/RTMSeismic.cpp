#include <RTM.hpp>

template<>
void RTMShotDescriptor<RTMData_t,RTMDevPtr_t>::loadReceiverGrid(RTMParam * rtmParam)
{

    RTMReceiverGrid<RTMData_t,RTMDevPtr_t> *rcv;
    int sx = getSource()->getX();
    int sy = getSource()->getY();
    int sz = getSource()->getZ();
    int rx = rtmParam->receiver_count_x;
    int ry = rtmParam->receiver_count_y;
    int rt = rtmParam->nt;
    /* Checks whether receiver position is within the model dimensions*/
    if (rx > rtmParam->nx || ry > rtmParam->ny)
    {
        RTMGridCoordinate c(rx, ry, rt);
        string msg("[ Invalid Receiver Grid Size: " + c.toString() + " ]");
        throw *(new RTMException(msg));
    }else
    {
        rcv = new RTMReceiverGrid<RTMData_t,RTMDevPtr_t>(rx, ry, rtmParam->ntstep);
    }
    rcv->setDistanceX(rtmParam->receiver_distance_x);
    rcv->setDistanceY(rtmParam->receiver_distance_y);
    rcv->setOffsetX(rtmParam->receiver_start_x);
    rcv->setOffsetY(rtmParam->receiver_start_y);
    rcv->setOffsetZ(rtmParam->receiver_depth_z);
    recvGrid = rcv;
}

