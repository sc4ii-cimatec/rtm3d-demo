#ifndef RTMPARAM_H
#define RTMPARAM_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <tuple>
#include <nlohmann/json.hpp>

namespace rtmparam
{

using namespace std;
using json = nlohmann::json;
using jsontype = json::value_t;

void parameterType(json::value_t &t, string &name);

/**
 * @brief Class RTMParamException Inheritance of class exception
 */
class RTMParamException : public exception
{
private:
    string what_msg; ///< Message for exception

public:
    /**
     * @brief Construct a new RTMParamException object
     * @param _what_msg 
     */
    RTMParamException(string &_what_msg)
    {
        what_msg = _what_msg;
    }
    RTMParamException(RTMParamException &e)
    {
        what_msg = e.getWhatMsg();
    }
    ~RTMParamException(){

    }
    RTMParamException& operator=(RTMParamException &e){
        what_msg = e.getWhatMsg();
        return *this;
    }
    string & getWhatMsg(){
        return what_msg;
    }

    /**
     * @brief Exception specifications
     * @return const char* string with explanatory information
     */    
    virtual const char *what() const throw()
    {
        return what_msg.c_str();
    }
};

/**
 * @brief Class RTMJsonParser
 */
class RTMJsonParser
{

protected:
    json jsonParser;
    string jsonFilePath;

    bool validated = false;

    vector<tuple<string, json::value_t>> paramList;

public:
    virtual bool validateParameters() = 0;
    virtual string &toString() =0;
    string& getJsonFilePath(){return jsonFilePath;}
};

/**
 * @brief Class RTMParam Inheritance of class RTMJsonParser
 * @details Initializes, loads, validates RTM parameters
 * @see RTMParam.cpp
 */
class RTMParam : public RTMJsonParser
{
private:
    void loadValues();

public:
    /**
     * @brief Construct a new RTMParam object
     */
    RTMParam()
    {
        initParamList();
    }
    /**
     * @brief Construct a new RTMParam object
     * @param iFilePath 
     */
    RTMParam(string iFilePath)
    {
        jsonFilePath = iFilePath;
    }

    /**
     * @brief Function initParamList
     * @details Parameters List defines input parameters
     *          and their types
     */
    void initParamList()
    {
        vector<tuple<string, json::value_t>> pList =
        {
            {"mname", json::value_t::string},
            {"vpfile", json::value_t::string},
            {"vpefile", json::value_t::string},
            {"save_vpe_file", json::value_t::boolean},
            {"outdir", json::value_t::string},
            {"tmpdir", json::value_t::string},
            {"datdir", json::value_t::string},
            {"nx", json::value_t::number_unsigned},
            {"ny", json::value_t::number_unsigned},
            {"nz", json::value_t::number_unsigned},
            {"blen", json::value_t::number_unsigned}, 
            {"width_m", json::value_t::number_unsigned},
            {"height_m", json::value_t::number_unsigned},
            {"depth_m", json::value_t::number_unsigned},
            {"nt", json::value_t::number_unsigned},
            {"dz", json::value_t::number_float},
            {"dy", json::value_t::number_float},
            {"dx", json::value_t::number_float},
            {"dt", json::value_t::number_float},
            {"fpeak", json::value_t::number_float},
            {"fmax", json::value_t::number_float},
            {"taper_factor", json::value_t::number_float},
            {"stencil_order", json::value_t::number_unsigned},
            {"modeling", json::value_t::boolean},
            {"migration", json::value_t::boolean},
            {"rtm_type", json::value_t::string},
            {"boundary_condition", json::value_t::string},
            {"source_total", json::value_t::number_unsigned},
            {"source_start_x", json::value_t::number_unsigned},
            {"source_start_y", json::value_t::number_unsigned},
            {"source_depth_z", json::value_t::number_unsigned},
            {"source_distance_x", json::value_t::number_unsigned},
            {"source_distance_y", json::value_t::number_unsigned},
            {"source_count_x", json::value_t::number_unsigned},
            {"source_count_y", json::value_t::number_unsigned},
            {"receiver_total", json::value_t::number_unsigned},
            {"receiver_start_x", json::value_t::number_unsigned},
            {"receiver_start_y", json::value_t::number_unsigned},
            {"receiver_depth_z", json::value_t::number_unsigned},
            {"receiver_distance_x", json::value_t::number_unsigned},
            {"receiver_distance_y", json::value_t::number_unsigned},
            {"receiver_count_x", json::value_t::number_unsigned},
            {"receiver_count_y", json::value_t::number_unsigned},
            {"ntstep", json::value_t::number_unsigned},
            {"nthreads", json::value_t::number_unsigned},
            {"save_vpe_file", json::value_t::boolean},
            {"load_vpe_from_file", json::value_t::boolean},
            {"save_snapshots", json::value_t::boolean},
            {"distributed_grid", json::value_t::boolean},
            {"snapshot_step", json::value_t::number_unsigned},

            {"fmig", json::value_t::boolean},
            {"fmig_max_freq", json::value_t::number_float},
            {"fmig_min_freq", json::value_t::number_float},
            {"fmig_nwstep", json::value_t::number_unsigned},
            {"fmig_distributed_imaging", json::value_t::boolean},

            
            {"fpga_xclbin", json::value_t::string},
            {"fpga_fwd_kernel", json::value_t::string},
            {"fpga_bwd_kernel", json::value_t::string}
        };

        for (int k0 = 0; k0 < pList.size(); k0++)
        {
            paramList.push_back(pList[k0]);
        }
    }

public:
    string mname; ///< modeling name
    string vpfile; ///< velocity model
    string vpefile; ///< extended velocity model
    string outdir; ///< output directory
    string tmpdir; ///< directory for data
    string rtm_type; ///< type fo RTM
    string datdir; ///< synthetic seismogram is saved into 'datdir' directory.
    string boundary_condition; ///< boundary condition
    string fpga_fwd_kernel; ///< 
    string fpga_bwd_kernel; ///< 
    string fpga_xclbin; ///<
    unsigned int nx; ///< number of samples in x
    unsigned int ny; ///< number of samples in y
    unsigned int nz; ///< number of samples in z
    unsigned int blen; ///< border length
    unsigned int width_m; ///<
    unsigned int height_m; ///<
    unsigned int depth_m; ///<
    unsigned int nt; ///< number total of time steps
    unsigned int stencil_order; ///< FD order
    unsigned int source_total; ///< number of sources
    unsigned int source_start_x; ///< source position in x
    unsigned int source_start_y; ///< source position in y
    unsigned int source_depth_z; ///< source depth in z
    unsigned int source_distance_x; ///< source distance in x
    unsigned int source_distance_y; ///< source distance in y
    unsigned int source_count_x; ///< shot window in x
    unsigned int source_count_y; ///< shot window in y
    unsigned int receiver_total; ///< number of receiver
    unsigned int receiver_start_x; ///< receiver position in x
    unsigned int receiver_start_y; ///< receiver position in y
    unsigned int receiver_depth_z; ///< receiver depth in z
    unsigned int receiver_distance_x; ///< receiver distance in x
    unsigned int receiver_distance_y; ///< receiver distance in y
    unsigned int receiver_count_x; ///< receiver count in x
    unsigned int receiver_count_y; ///< receiver count in y
    unsigned int ntstep; ///< number of time steps (ntstep<nt)
    unsigned int snapshot_step; ///< snapshot step
    unsigned int nthreads;
    float dz; ///< sampling interval in z
    float dy; ///< sampling interval in y
    float dx; ///< sampling interval in x
    float dt; ///< sampling interval in t
    float fpeak; ///< source peak frequency
    float fmax; ///< maximum frequency
    float taper_factor; ///< tapper factor
    
    bool modeling; ///< souce peak frequency
    bool migration; ///< migration
    bool save_vpe_file; ///< save extended velocity model
    bool load_vpe_from_file;
    bool verbose; ///< verbose
    bool save_snapshots; ///< save snapshots
    bool distributed_grid; /// < distributed_grid

    // frequency domain migration parameters
    bool    fmig;
    bool    fmig_distributed_imaging;
    float   fmig_max_freq;
    float   fmig_min_freq;
    unsigned int fmig_nwstep;

    /**
     * @brief Function loadRTMParam
     * @details Load RTM parameters
     */
    void loadRTMParam();
    /**
     * @brief Function printRTMParam
     */
    void printRTMParam();
    /**
     * @brief function validateParameters
     * @return true 
     * @return false 
     */
    virtual bool validateParameters();
    /**
     * @brief function toString
     * @return string& 
     */
    virtual string &toString();


    template <typename T>
    T &get(string name)
    {
        return jsonParser[name];
    }

    bool isValidated()
    {
        return validated;
    }
};

/**
 * @brief Class RTMReport Inheritance of class RTMJsonParser
 */
class RTMReport : public RTMJsonParser
{

public:
    string mname;
    string rtm_type;
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    unsigned int nt;
    unsigned int nw;

    /**
     * @brief Construct a new RTMReport object
     */
    RTMReport()
    {
        initParamList();
    }

    /**
     * @brief Construct a new RTMReport object
     * @param iFilePath 
     */
    RTMReport(string iFilePath)
    {
        jsonFilePath = iFilePath;
        initParamList();
    }

    /**
     * @brief Function initParamList
     */
    void initParamList()
    {
        vector<tuple<string, json::value_t>> pList =
        {
            {"mname", json::value_t::string},
            {"rtm_type", json::value_t::string},
            {"nx", json::value_t::number_unsigned},
            {"ny", json::value_t::number_unsigned},
            {"nz", json::value_t::number_unsigned},
            {"nt", json::value_t::number_unsigned},
            {"nw", json::value_t::number_unsigned}
        };

        for (int k0 = 0; k0 < pList.size(); k0++)
        {
            paramList.push_back(pList[k0]);
        }
    }

    template <typename T>
    T &get(string name)
    {
        return jsonParser[name];
    }

    template <typename T>
    void set(string name, T val)
    {
        jsonParser[name] = val;
    }

    /**
     * @brief Function addParam
     * @param name 
     * @param ptype 
     */
    void addParam(string name, json::value_t ptype)
    {
        tuple<string, json::value_t> *tp = new tuple<string, json::value_t>(name, ptype);
        paramList.push_back(*tp);
    }

    /**
     * @brief Function validateParameters
     * @return true 
     * @return false 
     */
    bool validateParameters(){
        return true;
    }

    /**
     * @brief Pure virtual function saveToFile()
     */
    virtual void saveToFile() = 0;
};

} // namespace rtmparam end

#endif