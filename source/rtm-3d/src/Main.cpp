#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <RTM.hpp>


using namespace std;

int main(int argc, char const *argv[])
{
    if (argc<2){
        cout << "> You forgot the input.json file!" << endl;
        return -1;
    }
    try{
        string inputFile(argv[1]);

        RTMParam rtmParam(inputFile);
        /** Loading input parameters from JSON file */
        rtmParam.loadRTMParam();

        /* fpga kernel path passed as input*/
        if (argc==3){
            string xclbinFile(argv[2]);
            rtmParam.fpga_xclbin = xclbinFile;
        }

        // creates RTMProcess instance
        RTMAcousticProcess rtmProcess(rtmParam);
        rtmProcess.rtm();
    }catch(RTMException &e){
        cout << "> Error:\n" << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }catch(RTMParamException &e){
        cout << "> Error:\n" << e.what() << endl;
        cout << "> Aborting!" << endl;
        exit(EXIT_FAILURE);
    }
    catch(...){
       cout << endl;
       cout << "> Error: Rougue Exception Flying Around. Aborting!" << endl;
       exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
