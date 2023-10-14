#pragma once
#include <iostream>

#include "functions_library.cuh"
#include "parameter_load.h"
#include <cstdlib>

#include <curand.h>
#include <curand_kernel.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

#include <sstream>

#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <string>
#include <map>

#include <string>
#include <vector>
#include <queue>

#include <thread>
#include <mutex>
#include <shared_mutex>

using namespace std;

class simulator_Master
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    string output_Folder_location;
    string intermediate_Folder_location;

    int CUDA_device_number;
    int gpu_Limit;
    int tot_Blocks;
    int tot_ThreadsperBlock;

    int CPU_cores;
    string output_Folder;
    string output_Folder_Sequences;
    string intermediate_Folders;
    string intermediate_sequence_Store;

    string node_Master_location;

    string multi_Read;

    string network_Model = "NA";
    int Total_number_of_Nodes = 0;

    string connection_Model = "FIXED";

    int BA_FIXED = 1;

    int BA_NB_sucesses = 0;
    float BA_NB_probability = 0;

    float BA_Poisson_mean = 0;

    int SCM_number_of_caves = 0;
    int SCM_number_of_nodes_per_cave = 0;
    // int Total_num_Nodes_SCM = 0;

    int DCM_number_of_caves = 0;
    int *per_cave_Stride;

    int DC_ND_sucesses = 0;
    float DC_ND_probability = 0;

    float DC_Poisson_mean = 0;

    float DC_percent_Neighbouring = 0;
    float DC_percent_Global_freedom = 0;

    float generation_Time = 0;

    float shape_days_in_Host = 0;
    float scale_days_in_Host = 0;

    int trials_Sampling = -1;
    float sampling_trials = 0;
    float sampling_probability = 0;
    string sampling_effect= "No";

    string progeny_distribution_Model = "NA";

    int progeny_NB_sucesses = 0;
    float progeny_NB_probability = 0;

    float progeny_Poisson_mean = 0;

    float progeny_Gamma_shape = 0;
    float progeny_Gamma_scale = 0;

    int number_of_node_Profiles = 0;
    string node_Profile_folder_Location = "";
    int *profile_Distributions;

    string output_Network_location = "";
    string network_File_location = "";

public:
    simulator_Master(string parameter_Master_Location);

    void configure_Network_Profile(string network_Profile_File, parameter_load &Parameters);

    void ingress();

    void network_Manager(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions);

    void BA_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions);
    void SCM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions);
    void DCM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions);
    void RANDOM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions);

    void node_Master_Manager(functions_library &functions);
};