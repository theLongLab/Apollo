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
    shared_mutex g_mutex;
    mutex gi_mutex;

    int *current_cell_load_per_Tissue;

    int *dead_Particle_count;

    int *parents_Prev_generation;

    vector<pair<string, string>> to_write_Sequence_Store;
    vector<string> converted_Sequences;

    int num_Tissues = 0;

    string generational_Summary;
    string sequence_Profiles;
    // string cells_of_parents_location;
    // string cells_of_progeny_location;
    string sequence_parent_Progeny_relationships;

    int CPU_cores = 0;

    set<int> found_Tissue_Folder_Indexes;

    int genome_Length;

public:
    cancer_Host();

    vector<set<int>> removed_by_Transfer_Indexes;

    void initialize(functions_library &functions,
                    vector<string> &tissue_Names,
                    string &intermediary_Sequence_location, string &first_Infection,
                    int &current_Generation,
                    string &output_Node_location,
                    int &max_sequences_per_File);

    void intialize_Tissues(string &host_Folder, vector<vector<pair<string, string>>> &tissue_Sequences, functions_library &functions, int &current_Generation);

    void sequence_Write_Configurator(vector<pair<string, string>> &sequence_Write_Store_All, vector<pair<string, string>> &sequence_Write_Store,
                                     int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                     string sequence_Profiles_Location, string tissue, int current_Generation);

    void partial_Write_Check(vector<pair<string, string>> &sequence_Write_Store_All,
                             const string &folder_Location, int &last_seq_Num,
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
                              vector<vector<float>> &time_Ratios_per_Tissue, vector<vector<string>> &phase_Type_per_tissue, vector<vector<pair<float, float>>> &phase_paramaters_per_Tissue,
                              int &max_Cells_at_a_time,
                              string &multi_Read, int &CPU_cores, int &num_Cuda_devices,
                              int &genome_Length);

    int terminal_status(int &num_tissues, int *tissue_array, string &source_sequence_Data_folder,
                        string &enable_Folder_management, string &enable_Compression, int &terminal_Load);

    int get_Load(int &num_tissues_Calc, int *tissue_array);

    void compress_Folder(string path, string &enable_Compression);

    string get_generation_Phase(int &overall_Generations,
                                vector<float> time_Generation,
                                vector<string> phase_Type,
                                vector<pair<float, float>> phase_paramaters_per_Tissue,
                                float &variable_1, float &variable_2);

    vector<pair<int, int>> get_Rounds(int &total_Count, int &gpu_Max_Limit);

    void simulate_cell_Round(functions_library &functions, string &multi_Read, int &num_Cuda_devices,
                             int &num_of_Cells, int &start, int &stop,
                             int *parents_in_Tissue, int &tissue, string tissue_Name,
                             vector<pair<int, int>> &indexed_Tissue_Folder,
                             string source_sequence_Data_folder,
                             int &overall_Generations,
                             int &last_index_Seq_Written,
                             mt19937 &gen);

    void replication_Generation_thread(int core_ID,
                                       vector<int> &parent_IDs, vector<string> &collected_Sequences,
                                       int &start, int &stop, int cell_Count);

    vector<string> find_Sequences_Master(int &offset, int &tissue, string &tissue_Name, functions_library &functions, string &folder_Path, int *parents_in_Tissue, int &num_Sequences, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation, vector<int> &parent_IDs, int &last_index_Seq_Written, mt19937 &gen);
    void thread_find_Files(int offset, int start, int stop, int *parents_in_Tissue, vector<pair<int, int>> &indexed_Tissue_Folder);
};