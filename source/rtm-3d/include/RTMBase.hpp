#ifndef RTMBASE_H
#define RTMBASE_H

#ifdef RTM_MPI
#include <mpi.h>
#else
#define MPI_Status int
#endif

#define RTM_NDIM_3D                 3
#define RTM_NDIM_2D                 2
#define RTM_NDIM_1D                 1

#define RTM_LAPFILTER_ORDER         2

#define rtmMpiTag_t                 int
#define RTM_MPI_REQ_PAREA           0x4000
#define RTM_MPI_REQ_GZONE           0x1000
#define RTM_MPI_RCV_GZONE           0x2000
#define RTM_MPI_PGRID_TAG           (RTM_MPI_REQ_GZONE | 0x0001)
#define RTM_MPI_PPGRID_TAG          (RTM_MPI_REQ_GZONE | 0x0002)
#define RTM_MPI_PRGRID_TAG          (RTM_MPI_REQ_GZONE | 0x0004)
#define RTM_MPI_PPRGRID_TAG         (RTM_MPI_REQ_GZONE | 0x0008)
#define RTM_MPI_UPBGRID_TAG         (RTM_MPI_REQ_GZONE | 0x0010)
#define RTM_MPI_OUTIMG_TAG          (RTM_MPI_REQ_GZONE | 0x0020)
#define RTM_MPI_SEISMRCV_TAG        (RTM_MPI_REQ_GZONE | 0x0040)

#ifdef RTM_MPI
class RTMMPIComm{ 
    public: 
        static MPI_Comm rtm_mpi_comm_valid;
        static MPI_Comm rtm_mpi_comm_invalid;
};
#define RTM_COMM_VALID RTMMPIComm::rtm_mpi_comm_valid
#define RTM_COMM_INVALID RTMMPIComm::rtm_mpi_comm_invalid
#endif

#define RTM_COMM_COLOR_VALID        0x0001
#define RTM_COMM_COLOR_INVALID      0x0004


#include <RTMParam.hpp>
#include <RTMException.hpp>
#include <RTMAcc.hpp>

using namespace std;

using rtmparam::RTMParam;
using rtmparam::RTMParamException;

#define RTMData_t                  float
#if defined (RTM_ACC_GPU)
#define RTMDevPtr_t                 RTMData_t
#define DefaultDeviceGrid           GPUDeviceGrid
#define RTMAccPlatform              RTMGPUPlatform
#elif defined (RTM_ACC_FPGA)
#define DefaultDeviceGrid           FPGADeviceGrid
#define RTMDevPtr_t                 CLBuffer_t
#define RTMAccPlatform              RTMFPGAPlatform
#else
#define DefaultDeviceGrid           CPUDeviceGrid
#define RTMDevPtr_t                 RTMData_t
#define RTMAccPlatform              RTMCPUPlatform
#endif


/**
 * @brief Function static inline RTM_MIGIMG_NAME()
 * @param str  image file
 * @param datdir  output directory
 * @param mname modeling name
 * @param ns  shot window in x * y
 * @param nx  number of samples in x
 * @param ny  number of samples in y
 * @param nz  number of samples in z
 * @param nt  number total of time steps
 * @param np  number of process
 */
static inline void RTM_MIGIMG_NAME(string &str, string &datdir, string &mname,
                                       int ns, int nx, int ny, int nz, int nt, int np)
{
    str.clear();

    std::for_each(mname.begin(), mname.end(), [](char & c){
	    c = ::toupper(c);
    });


    str += datdir + "/" + "RTMIMG_"+ mname +"_NS" + to_string(ns) + "_NXNYNZ_" + to_string(nx) + "x" 
           + to_string(ny) + "x" + to_string(nz) + "_NT" + to_string(nt) + ".bin.NP"+ to_string(np);
}

/**
 * @brief Function static inline RTM_HBCUPB_NAME()
 * @details Saves intermediary large upper-boundary files to avoid memory blow-up
 * @param str up border file
 * @param datdir directory for data
 * @param sx source position in x
 * @param sy source position in y
 * @param sz source position in z
 * @param nx number of samples in x
 * @param ny number of samples in y
 * @param nz number of samples in z
 * @param ntstep number of time steps (ntstep<nt)
 * @param it time steps
 * @param np number of process
 */
static inline void RTM_HBCUPB_NAME(string &str, string &datdir,
                                       int sx, int sy, int sz, int nx, int ny, int nz, 
                                       int ntstep, int it, 
                                       int pRank, int nP)
{
    str.clear();
    str += datdir + "/" + "RTMHBC_S" + to_string(sx) + "x" + to_string(sy) + "x" +
           to_string(sz) + "_NXNYNZ_" + to_string(nx) + "x" 
           + to_string(ny) + "x" + to_string(nz) + "_NTSTEP" + to_string(ntstep) + "_IT" + to_string(it) + ".upb.rank"
           +to_string(pRank) + "x" + to_string(nP);
}

/**
 * @brief Function static inline RTM_SEISMOGRAM_NAME()
 * @param str Input Seismic File
 * @param datdir synthetic seismogram is saved into 'datdir' directory.
 * @param sx source position in x
 * @param sy source position in y
 * @param sz source position in z
 * @param nx number of samples in x
 * @param ny number of samples in y
 * @param nt number total of time steps
 */
static inline void RTM_SEISMOGRAM_NAME(string &str, string &datdir,
                                       int sx, int sy, int sz, int nx, int ny, int startnt, int ntstep)
{
    str.clear();
    str += datdir + "/" + "RTMSEISMOGRAM_S" + to_string(sx) + "x" + to_string(sy) + "x" +
           to_string(sz) + "_NXNY_" + to_string(nx) + "x" + to_string(ny) + 
           "_NTSTART" + to_string(startnt) +"_NTSTEP" +to_string(ntstep) + ".seism";
}

/**
 * @brief Function static inline RTM_SEISMTRACE_NAME()
 * @param str SeismicTraces file
 * @param datdir 'datdir' directory
 * @param sx source position in x
 * @param sy source position in y
 * @param sz source position in z
 * @param rx receiver position in x
 * @param ry receiver position in y
 * @param rz receiver position in z
 * @param nt number total of time step
 */
static inline void RTM_SEISMTRACE_NAME(string &str, string &datdir, int sx, int sy, int sz, int rx, int ry, int rz, int nt)
{
    str.clear();
    str += datdir + "/" + "RTMSEISM_S" + to_string(sx) + "x" + to_string(sy) + "x" +
           to_string(sz) + "_R" + to_string(rx) + "x" + to_string(ry) + "x" +
           to_string(rz) + "_NT" + to_string(nt) + ".trace";
}
/**
 * @brief Function static inline RTM_PRINT()
 * @param msg 
 * @param b 
 */
static inline void RTM_PRINT(const string &msg, bool b)
{
    if (b)
    {
        std::cout << "> " + msg << endl << flush;
    }
}

#endif // RTMBASE_H