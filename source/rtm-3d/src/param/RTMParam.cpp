#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <tuple>
#include <nlohmann/json.hpp>
#include <RTM.hpp>
#include <RTMParam.hpp>

using namespace std;

using namespace rtmparam;
using json = nlohmann::json;


const string & RTM_HEADER_MSG = "3D Acoustic RTM Program";

/* RTMParam::parameterType(json::value_t &t, string &name) */
void rtmparam::parameterType(json::value_t &t, string &name)
{

    switch (t)
    {
    case json::value_t::null:
        name = "Null";
        break;
    case json::value_t::string:
        name = "String";
        break;
    case json::value_t::object:
        name = "Object";
        break;
    case json::value_t::array:
        name = "Array";
        break;
    case json::value_t::boolean:
        name = "Boolean";
        break;
    case json::value_t::number_integer:
        name = "Signed_Integer";
        break;
    case json::value_t::number_unsigned:
        name = "Unsigned_Integer";
        break;
    case json::value_t::number_float:
        name = "Float";
        break;
    case json::value_t::discarded:
        name = "Discarded";
        break;
    default:
        name = "String";
        break;
    }
}

/* RTMParam::validateParameters() */
bool RTMParam::validateParameters()
{

    ifstream ifs(jsonFilePath);

    if (ifs.fail())
    {
        std::string msg = "[ Input File '" + jsonFilePath + "' doesn't exists. Please check input file path. ]";
        RTMParamException ex(msg);
        throw ex;
        msg = false;
    }
    jsonParser = json::parse(ifs);

    vector<tuple<string, json::value_t>>::iterator ptr;

    // Displaying vector elements using begin() and end()
    for (ptr = paramList.begin(); ptr < paramList.end(); ptr++)
    {
        tuple<string, json::value_t> &tp = (*ptr);
        string pname = std::get<0>(tp);
        json::value_t valtype = std::get<1>(tp);
        
        json::value_t exptype = jsonParser[pname].type();
        if (jsonParser[pname].type() != valtype)
        {
            string expTypeName;
            string valTypeName;
            parameterType(exptype, expTypeName);
            parameterType(valtype, valTypeName);
            std::string msg = "[ Input parameter '" + pname + "' has wrong type! Expecpted: '" + valTypeName + "' got: '" + expTypeName + "' ]";
            RTMParamException ex(msg);
            throw ex;
        }
    }
    validated = true;
    return true;
}

/* RTMParam::loadValues() */
void RTMParam::loadValues(){

    mname                       = jsonParser["mname"];
    vpfile                      = jsonParser["vpfile"];
    vpefile                     = jsonParser["vpefile"];
    outdir                      = jsonParser["outdir"];
    tmpdir                      = jsonParser["tmpdir"];
    rtm_type                    = jsonParser["rtm_type"];
    boundary_condition          = jsonParser["boundary_condition"];
    datdir                      = jsonParser["datdir"];
    nx                          = jsonParser["nx"];
    ny                          = jsonParser["ny"];
    nz                          = jsonParser["nz"];
    blen                        = jsonParser["blen"];
    width_m                     = jsonParser["width_m"];
    height_m                    = jsonParser["height_m"];
    depth_m                     = jsonParser["depth_m"];
    nt                          = jsonParser["nt"];
    stencil_order               = jsonParser["stencil_order"];
    source_total                = jsonParser["source_total"];
    source_start_x              = jsonParser["source_start_x"];
    source_start_y              = jsonParser["source_start_y"];
    source_depth_z              = jsonParser["source_depth_z"];
    source_distance_x           = jsonParser["source_distance_x"];
    source_distance_y           = jsonParser["source_distance_y"];
    source_count_x              = jsonParser["source_count_x"];
    source_count_y              = jsonParser["source_count_y"];
    receiver_total              = jsonParser["receiver_total"];
    receiver_start_x            = jsonParser["receiver_start_x"];
    receiver_start_y            = jsonParser["receiver_start_y"];
    receiver_depth_z            = jsonParser["receiver_depth_z"];
    receiver_distance_x         = jsonParser["receiver_distance_x"];
    receiver_distance_y         = jsonParser["receiver_distance_y"];
    receiver_count_x            = jsonParser["receiver_count_x"];
    receiver_count_y            = jsonParser["receiver_count_y"];
    ntstep                      = jsonParser["ntstep"];
    nthreads                    = jsonParser["nthreads"];
    dz                          = jsonParser["dz"];
    dy                          = jsonParser["dy"];
    dx                          = jsonParser["dx"];
    dt                          = jsonParser["dt"];
    fpeak                       = jsonParser["fpeak"];
    fmax                        = jsonParser["fmax"];
    taper_factor                = jsonParser["taper_factor"];
    modeling                    = jsonParser["modeling"];
    migration                   = jsonParser["migration"];
    save_vpe_file               = jsonParser["save_vpe_file"];
    load_vpe_from_file          = jsonParser["load_vpe_from_file"];
    verbose                     = jsonParser["verbose"];
    save_snapshots              = jsonParser["save_snapshots"];
    snapshot_step               = jsonParser["snapshot_step"];
    distributed_grid            = jsonParser["distributed_grid"];

    fmig                        = jsonParser["fmig"];
    fmig_max_freq               = jsonParser["fmig_max_freq"];
    fmig_min_freq               = jsonParser["fmig_min_freq"];
    fmig_nwstep                 = jsonParser["fmig_nwstep"];
    fmig_distributed_imaging    = jsonParser["fmig_distributed_imaging"];

    fpga_xclbin                = jsonParser["fpga_xclbin"];
    fpga_fwd_kernel            = jsonParser["fpga_fwd_kernel"];
    fpga_bwd_kernel            = jsonParser["fpga_bwd_kernel"];
}

string & RTMParam::toString(){

    string * str = new string("#####################################################################");

    *str += "\n  # mname                         = "+  ( mname);
    *str += "\n  # vpfile                        = "+ ( vpfile);
    *str += "\n  # vpefile                       = "+ ( vpefile);
    *str += "\n  # outdir                        = "+ ( outdir);
    *str += "\n  # tmpdir                        = "+ ( tmpdir);
    *str += "\n  # rtm_type                      = "+ ( rtm_type);
    *str += "\n  # nthreads                      = "+ to_string(nthreads);
    *str += "\n  # boundary_condition            = "+ boundary_condition;
    *str += "\n  # datdir                        = "+ ( datdir);
    *str += "\n  # nx                            = "+ to_string(nx);
    *str += "\n  # ny                            = "+ to_string(ny);
    *str += "\n  # nz                            = "+ to_string(nz);
    *str += "\n  # blen                          = "+ to_string(blen);
    *str += "\n  # nt                            = "+ to_string(nt);
    *str += "\n  # ntstep                        = "+ to_string(ntstep);
    *str += "\n  # dz                            = "+ to_string(dz);
    *str += "\n  # dy                            = "+ to_string(dy);
    *str += "\n  # dx                            = "+ to_string(dx);
    *str += "\n  # dt                            = "+ to_string(dt);
    *str += "\n  # fpeak                         = "+ to_string(fpeak);
    *str += "\n  # fmax                          = "+ to_string(fmax);
    *str += "\n  # taper_factor                  = "+ to_string(taper_factor);
    *str += "\n  # stencil_order                 = "+ to_string(stencil_order);
    *str += "\n  # source_total                  = "+ to_string(source_total);
    *str += "\n  # source_start_x                = "+ to_string(source_start_x);
    *str += "\n  # source_start_y                = "+ to_string(source_start_y);
    *str += "\n  # source_depth_z                = "+ to_string(source_depth_z);
    *str += "\n  # source_distance_x             = "+ to_string(source_distance_x);
    *str += "\n  # source_distance_y             = "+ to_string(source_distance_y);
    *str += "\n  # source_count_x                = "+ to_string(source_count_x);
    *str += "\n  # source_count_y                = "+ to_string(source_count_y);
    *str += "\n  # receiver_total                = "+ to_string(receiver_total);
    *str += "\n  # receiver_start_x              = "+ to_string(receiver_start_x);
    *str += "\n  # receiver_start_y              = "+ to_string(receiver_start_y);
    *str += "\n  # receiver_depth_z              = "+ to_string(receiver_depth_z);
    *str += "\n  # receiver_distance_x           = "+ to_string(receiver_distance_x);
    *str += "\n  # receiver_distance_y           = "+ to_string(receiver_distance_y);
    *str += "\n  # receiver_count_x              = "+ to_string(receiver_count_x);
    *str += "\n  # receiver_count_y              = "+ to_string(receiver_count_y);
    *str += "\n  # modeling                      = "+ to_string(modeling);
    *str += "\n  # migration                     = "+ to_string(migration);
    *str += "\n  # distributed_grid              = "+ to_string(distributed_grid);
    *str += "\n  # save_vpe_file                 = "+ to_string(save_vpe_file);
    *str += "\n  # load_vpe_from_file            = "+ to_string(load_vpe_from_file);
    *str += "\n  # verbose                       = "+ to_string(verbose);
    *str += "\n  # save_snapshots                = "+ to_string(save_snapshots);
    *str += "\n  # snapshots_step                = "+ to_string(snapshot_step);
    *str += "\n  # ";
    *str += "\n  # Freq. Domain Migration: ";
    *str += "\n  #   fmig                        = "+ (to_string(fmig));
    *str += "\n  #   fmig_max_freq               = "+ (to_string(fmig_max_freq));
    *str += "\n  #   fmig_min_freq               = "+ (to_string(fmig_min_freq));
    *str += "\n  #   fmig_nwstep                 = "+ (to_string(fmig_nwstep));
    *str += "\n  #   fmig_distributed_imaging    = "+ (to_string(fmig_distributed_imaging));
    *str += "\n  # ";
    *str += "\n  # HW Acceleration: ";
    *str += "\n  #   fpga_xclbin                 = "+ (fpga_xclbin);
    *str += "\n  #   fpga_fwd_kernel             = "+ (fpga_fwd_kernel);
    *str += "\n  #   fpga_bwd_kernel             = "+ (fpga_bwd_kernel);

    *str += "\n> #####################################################################";
    return *str;
}

/* RTMParam::loadRTMParam */
void RTMParam::loadRTMParam(){

    try
    {
        validateParameters();
        // load all parameter values
        loadValues();

        // check if nt is smaller than nt step
        if(nt < ntstep){
            string msg("[ RTMParam Error: NT smaller than NTSTEP ]");
		    RTMParamException ex(msg);
            throw ex;
	    }
        if (fmig){
            if(fmig_max_freq <= fmig_min_freq){
                string msg("[ RTMParam Error: freq_max smaller than freq_min ]");
		        RTMParamException ex(msg);
                throw ex;
            }
            if(distributed_grid){
                string msg("[ RTMParam Error: Can't use 'distributed_grid' on fmig migrations. Set it to false on input json file ]");
		        RTMParamException ex(msg);
                throw ex;
            }
        }
    }catch (RTMParamException &e)
    {
        cout << "> Error:\n" << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }
}

void RTMParam::printRTMParam(){
    RTM_PRINT(RTM_HEADER_MSG, verbose);
    RTM_PRINT(toString(), verbose);
}