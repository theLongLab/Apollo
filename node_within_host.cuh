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

class node_within_host
{
private:
    int host_Index, cave_ID, host_ID, profile_ID;
    int num_Generation;

    int infectious_Load;
    int terminal_Load;
    float sampling_Effect;

    int num_Tissues;
    int *cell_Limit;

    string status = "Susceptible";

    int current_Generation = -1;
    int *current_Viral_load_per_Tissue;

    // Remember to clear after getting the indexes in th current tissue;
    vector<set<int>> removed_by_Transfer_Indexes;

public:
    node_within_host();

    void setHost(int host_Index, int cave_ID, int host_ID, int profile_ID, int num_Tissues);
    void setNum_Generation(int num_Generation);
    void setInfectious_Load(int infectious_Load);
    void setTerminal_Load(int terminal_Load);
    void setSampling_Effect(float sampling_Effect);
    void setCell_Limit(vector<int> cell_Limit_vec);

    int get_host_Index();
    string get_Name();
    string get_Status();
    int get_Profile();
    int get_Generation();
    int *get_current_Viral_load_per_Tissue();

    int get_Load(int &num_tissues_Calc, int *tissue_array);

    int infectious_status(int &num_tissues_Calc, int *tissue_array);
    int terminal_status(int &num_tissues_Calc, int *tissue_array);

    void print_All();

    void begin_Infection(functions_library &functions, string &intermediary_Sequence_location,
                         int entry_tissues, int *entry_array, int &max_sequences_per_File);
    void transfer_Infection(functions_library &functions, string &intermediary_Sequence_location, string &source_Target_file_Location,
                            int &source_Index, int &source_Generation, string &source_Name, int *source_current_Viral_load_per_Tissue,
                            int num_viruses_to_transfer,
                            int &entry_tissues, int *entry_array, int exit_Load, int &exit_tissues, int *exit_array,
                            int &max_sequences_per_File,
                            vector<vector<pair<int, int>>> &indexed_Source_Folders,
                            mt19937 &gen);
    void run_Generation();

    void intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions);
};