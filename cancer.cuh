#pragma once
#include <iostream>

#include "functions_library.cuh"
#include "parameter_load.h"
#include "cancer_Host.cuh"
#include "simulator_Master.cuh"

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

class cancer
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    string output_Folder_location;
    string intermediate_Folder_location;
    string output_Node_location;

    string enable_Folder_management = "NO";
    string enable_Compression = "NO";

    int max_sequences_per_File = 0;
    int max_Cells_at_a_time = 0;

    string start_Date;
    // string stop_after_generations = "NO";
    string first_Infection = "Random";

    string node_Master_location;
    string sequence_Master_location;

    int stop_gen_Mode = 0;
    int stop_generations_Count = 0;
    float stop_Date = 0;

    float generation_Time = 0;

    int CPU_cores;
    string multi_Read;

    int *CUDA_device_IDs;
    int num_Cuda_devices;
    int gpu_Limit;
    int *tot_Blocks;
    int *tot_ThreadsperBlock;

    string parent_Sequence_Folder;
    float *Reference_fitness_survivability_proof_reading;
    int *mutation_proof_Reading_availability;

    int mutation_Hotspots = 0;

    float **A_0_mutation;
    float **T_1_mutation;
    float **G_2_mutation;
    float **C_3_mutation;

    float **mutation_hotspot_parameters;

    int num_tissues_per_Node = 0;
    vector<string> tissue_Names;

    int terminal_tissues = 0;
    int *terminal_array;

    string viral_Migration = "No";
    float **viral_Migration_Values;
    int *migration_start_Generation;

    string profile_Name = "";

    int terminal_Load = -1;

    // per tissue number of cells;
    // Inactive = -1;
    int *profile_tissue_Limits;

    vector<int> replication_phases_tissues;
    vector<vector<float>> time_Ratios_per_Tissue;
    vector<vector<string>> phase_Type_per_tissue;
    vector<vector<pair<float, float>>> phase_paramaters_per_Tissue;

    // columns
    // 0 = position-1;
    // 1 = A
    // 2 = T
    // 3 = G
    // 4 = C
    float **sequence_Fitness_changes;
    float **sequence_Survivability_changes;
    float **sequence_Proof_reading_changes;

    float **sequence_replication_factor_changes;
    float **sequence_mutation_rate_changes;
    float **sequence_generation_death_changes;
    float **sequence_replication_prob_changes;
    float **sequence_metastatic_prob_changes;

    // 0 = Reference replication factor
    // 1 = Reference mutation rate factor
    // 2 = Reference generational death
    // 3 = Reference replication probability
    // 4 = Reference metastatic probability
    float *Reference_cancer_parameters;

    // 0 = fitness
    // 1 = survivability
    // 2 = proof reading
    int *num_effect_Segregating_sites;

    // 0 = Reference replication factor
    // 1 = Reference mutation rate factor
    // 2 = Reference generational death
    // 3 = Reference replication probability
    // 4 = Reference metastatic probability
    int *num_effect_Segregating_sites_Cancer;

    int genome_Length = 0;

    string intermediary_Sequence_location;
    string intermediary_Index_location;
    string reference_Sequences;

    vector<pair<string, string>> all_sequences_String;

    int count_tajima_Regions = 0;
    int **tajima_regions_Start_Stop;
    string reference_Genome_location = "";

public:
    cancer(string parameter_Master_Location);

    void ingress();

    void node_Master_Manager(functions_library &functions);
    void sequence_Master_Manager(functions_library &functions);

    int **process_Reference_Sequences(functions_library &functions, vector<string> collect_Sequences, int &genome_Length, int &num_of_Sequences_current, float *replication_probs, float *gen_Death_probs,
                                      float *replication_Factor, float *metastatic_Prob, float *Survivability);

    vector<string> read_Reference_Sequences(vector<int> &tissue_Sequence_Count);
    void write_Reference_Sequences(vector<string> collect_Sequences, vector<int> &tissue_Sequence_Count, functions_library &functions);

    vector<pair<string, string>> convert_Sequences_Master(int **sequences, int &genome_Length, int &num_of_Sequences_current, float *replication_probs, float *gen_Death_probs,
                                                          float *replication_Factor, float *metastatic_Prob, float *Survivability);
    void sequence_to_string_Threads(int start, int stop, int **sequences, int genome_Length, float *replication_probs, float *gen_Death_probs,
                                    float *replication_Factor, float *metastatic_Prob, float *Survivability);

    void sequence_Write_Configurator(vector<pair<string, string>> &sequence_Write_Store_All, vector<pair<string, string>> sequence_Write_Store,
                                     int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num);

    void partial_Write_Check(vector<pair<string, string>> &sequence_Write_Store_All,
                             const string &folder_Location, int &last_seq_Num);
};