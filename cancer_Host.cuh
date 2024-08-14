#pragma once
#include <iostream>

#include "functions_library.cuh"
#include "parameter_load.h"
#include "node_within_host.cuh"

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

class cancer_Host
{
private:
    int *current_cell_load_per_Tissue;

    int *dead_Particle_count;

    int *parents_Prev_generation;

    vector<pair<string, string>> to_write_Sequence_Store;
    vector<string> converted_Sequences;

    int num_Tissues = 0;

public:
    cancer_Host();

    vector<set<int>> removed_by_Transfer_Indexes;

    void initialize(functions_library &functions,
                    vector<string> &tissue_Names,
                    string &intermediary_Sequence_location, string &first_Infection,
                    int &current_Generation,
                    string &output_Node_location,
                    int &max_sequences_per_File);

    void intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions, int &current_Generation);

    void sequence_Write_Configurator(vector<string> &sequence_Write_Store_All, vector<string> &sequence_Write_Store,
                                     int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                     vector<char> &seq_Status,
                                     string sequence_Profiles_Location, string tissue, int current_Generation);

    void partial_Write_Check(vector<string> &sequence_Write_Store_All,
                             const string &folder_Location, int &last_seq_Num,
                             vector<char> &seq_Status,
                             string sequence_Profiles_Location, string tissue, int current_Generation);

    void simulate_Generations(functions_library &functions,
                              int &overall_Generations, float &date_Increment,
                              int &stop_Type,
                              int &stop_gen_Mode,
                              int &stop_generations_Count, float &decimal_Date, float &stop_Date,
                              vector<string> &tissue_Names,
                              int &terminal_tissues, int *terminal_array,
                              string source_sequence_Data_folder,
                              string &enable_Folder_management, string &enable_Compression,
                              int &terminal_Load,
                              string &output_Node_location,
                              vector<vector<float>> &time_Ratios_per_Tissue, vector<vector<string>> &phase_Type_per_tissue, vector<vector<pair<float, float>>> &phase_paramaters_per_Tissue);

    int terminal_status(int &num_tissues, int *tissue_array, string &source_sequence_Data_folder,
                        string &enable_Folder_management, string &enable_Compression, int &terminal_Load);

    int get_Load(int &num_tissues_Calc, int *tissue_array);

    void compress_Folder(string path, string &enable_Compression);

    string get_generation_Phase(int &overall_Generations,
                                vector<float> time_Generation,
                                vector<string> phase_Type,
                                vector<pair<float, float>> phase_paramaters_per_Tissue,
                                float &variable_1, float &variable_2);
};