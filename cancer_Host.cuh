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

    // vector<pair<string, string>> to_write_Sequence_Store;
    // vector<string> converted_Sequences;

    int num_Tissues = 0;

    string generational_Summary;
    string sequence_Profiles;
    // string cells_of_parents_location;
    // string cells_of_progeny_location;
    string sequence_parent_Progeny_relationships;

    int CPU_cores = 0;

    set<int> found_Tissue_Folder_Indexes;

    int genome_Length;

    vector<pair<string, string>> to_write_Sequence_Store_NEXT_Gen;
    // vector<pair<string, string>> to_write_Sequence_Store_THIS_Gen;

    vector<vector<vector<pair<string, string>>>> to_write_Sequence_Store_OTHER_Gens;
    vector<vector<int>> last_index_Seq_Written_OTHERs;

public:
    cancer_Host();

    vector<set<int>> removed_by_Transfer_Indexes;

    void initialize(functions_library &functions,
                    vector<string> &tissue_Names,
                    string &intermediary_Sequence_location, string &first_Infection,
                    int &current_Generation,
                    string &output_Node_location,
                    int &max_sequences_per_File);

    void intialize_Tissues(string &host_Folder, vector<vector<pair<string, string>>> &tissue_Sequences, vector<vector<string>> &profile_Lines_Tissues,
                           functions_library &functions, int &current_Generation);

    void sequence_Write_Configurator(vector<pair<string, string>> &sequence_Write_Store_All, vector<pair<string, string>> &sequence_Write_Store, vector<string> &profile_Lines,
                                     int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                     string sequence_Profiles_Location, string tissue, int current_Generation);

    void partial_Write_Check(vector<pair<string, string>> &sequence_Write_Store_All,
                             const string &folder_Location, int &last_seq_Num,
                             string sequence_Profiles_Location, string tissue, int current_Generation);

    void full_Write_Sequences_NEXT_Generation(int &max_sequences_per_File, string next_Generation_location,
                                              functions_library &functions, vector<pair<string, string>> &to_write_Sequence_Store);

    void remainder_Write_Sequences_NEXT_Generation(string next_Generation_location,
                                                   functions_library &functions, vector<pair<string, string>> &to_write_Sequence_Store);

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
                              string &multi_Read, int &CPU_cores, int &num_Cuda_devices, int *CUDA_device_IDs,
                              int &genome_Length,
                              float *Reference_fitness_survivability_proof_reading, float *Reference_cancer_parameters,
                              float **A_0_mutation,
                              float **T_1_mutation,
                              float **G_2_mutation,
                              float **C_3_mutation,
                              int &mutation_Hotspots,
                              float **mutation_hotspot_parameters,
                              int *num_effect_Segregating_sites,
                              float **sequence_Survivability_changes,
                              float **sequence_Proof_reading_changes,
                              int *num_effect_Segregating_sites_Cancer,
                              float **sequence_replication_factor_changes,
                              float **sequence_mutation_rate_changes,
                              float **sequence_generation_death_changes,
                              float **sequence_replication_prob_changes,
                              float **sequence_metastatic_prob_changes,
                              int &max_sequences_per_File,
                              float **viral_Migration_Values, int *migration_start_Generation,
                              int &count_tajima_Regions, int **tajima_regions_Start_Stop,
                              string &reference_Genome_location,
                              int *tissue_selection_Position_Count,
                              int *Survivability_Positions,
                              int *Proof_Positions,
                              int *Replication_factor_Positions,
                              int *Mutation_rate_factor_Positions,
                              int *Generation_death_Positions,
                              int *Replication_prob_Positions,
                              int *Metastatic_Positions,
                              float **tissues_ATGC_positions_Survivability,
                              float **tissues_ATGC_positions_Proof,
                              float **tissues_ATGC_positions_Replication_factor,
                              float **tissues_ATGC_positions_Mutation_rate_factor,
                              float **tissues_ATGC_positions_Generation_death,
                              float **tissues_ATGC_positions_Replication_prob,
                              float **tissues_ATGC_positions_Metastatic);

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

    void simulate_cell_Round(functions_library &functions, string &multi_Read, int &num_Cuda_devices, int *CUDA_device_IDs,
                             int &num_of_Cells, int &start, int &stop,
                             int *parents_in_Tissue, int &tissue, string tissue_Name,
                             vector<pair<int, int>> &indexed_Tissue_Folder,
                             string this_Gen_intermediary_Sequences,
                             int &overall_Generations,
                             int &last_index_Seq_Written,
                             mt19937 &gen,
                             float *Reference_fitness_survivability_proof_reading, float *Reference_cancer_parameters,
                             float **A_0_mutation,
                             float **T_1_mutation,
                             float **G_2_mutation,
                             float **C_3_mutation,
                             int &mutation_Hotspots,
                             float **mutation_hotspot_parameters,
                             int *num_effect_Segregating_sites,
                             float **sequence_Survivability_changes,
                             float **sequence_Proof_reading_changes,
                             int *num_effect_Segregating_sites_Cancer,
                             float **sequence_replication_factor_changes,
                             float **sequence_mutation_rate_changes,
                             float **sequence_generation_death_changes,
                             float **sequence_replication_prob_changes,
                             float **sequence_metastatic_prob_changes,
                             int &max_sequences_per_File, string &intermediary_Tissue_folder, string &source_sequence_Data_folder,
                             int &last_Progeny_written_this_Gen, string rapid_Progeny_Location,
                             int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List,
                             int *tissue_selection_Position_Count,
                             int *Survivability_Positions,
                             int *Proof_Positions,
                             int *Replication_factor_Positions,
                             int *Mutation_rate_factor_Positions,
                             int *Generation_death_Positions,
                             int *Replication_prob_Positions,
                             int *Metastatic_Positions,
                             float **tissues_ATGC_positions_Survivability,
                             float **tissues_ATGC_positions_Proof,
                             float **tissues_ATGC_positions_Replication_factor,
                             float **tissues_ATGC_positions_Mutation_rate_factor,
                             float **tissues_ATGC_positions_Generation_death,
                             float **tissues_ATGC_positions_Replication_prob,
                             float **tissues_ATGC_positions_Metastatic);

    // void replication_Generation_thread(int gpu, cudaStream_t *streams,
    //                                    char *cuda_full_Char, char *full_Char,
    //                                    int **cuda_progeny_Sequences,
    //                                    int &start, int &stop, int cell_Count,
    //                                    int *CUDA_device_IDs);

    string find_Sequences_Master(int &offset, int &tissue, string &tissue_Name, functions_library &functions, string &folder_Path, int *parents_in_Tissue, int &num_Sequences, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation, vector<int> &parent_IDs, float *parents_Elapsed, int &last_index_Seq_Written, mt19937 &gen,
                                 int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List);
    void thread_find_Files(int offset, int start, int stop, int *parents_in_Tissue, vector<pair<int, int>> &indexed_Tissue_Folder);

    void thread_Sequence_to_String_Cancer(int start, int stop, int **progeny_Sequences);

    vector<pair<int, int>> compile_Progeny(functions_library &functions,
                                           string &intermediary_Tissue_folder, string &rapid_Progeny_Location,
                                           int &parent_Cells_Found,
                                           float *progeny_Elapsed, float **progeny_Configuration_Cancer,
                                           int &last_index_Seq_Written, int &overall_Generations, string &tissue_Name,
                                           int **progeny_Sequences, int &tissue_Index,
                                           mt19937 &gen, vector<int> &parent_IDs,
                                           string &source_sequence_Data_folder,
                                           int &last_Progeny_written_this_Gen,
                                           int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List);

    // void rerun_Progeny_THIS_gen(functions_library &functions, vector<pair<int, int>> &rerun_Progeny, vector<pair<int, int>> &start_stop_Per_GPU, int &num_Cuda_devices, int &tot_Parents,
    //                             int **parent_sequences_INT, float *parents_Elapsed,
    //                             float *cuda_Reference_fitness_survivability_proof_reading[],
    //                             float *cuda_Reference_cancer_parameters[],
    //                             float **cuda_sequence_replication_factor_changes[],
    //                             int mutation_Hotspots, float **cuda_mutation_hotspot_parameters[],
    //                             float **cuda_A_0_mutation[],
    //                             float **cuda_T_1_mutation[],
    //                             float **cuda_G_2_mutation[],
    //                             float **cuda_C_3_mutation[],
    //                             int *cuda_num_effect_Segregating_sites[],
    //                             int *cuda_num_effect_Segregating_sites_Cancer[],
    //                             float **cuda_sequence_Survivability_changes[],
    //                             float **cuda_sequence_Proof_reading_changes[],
    //                             float **cuda_sequence_mutation_rate_changes[],
    //                             float **cuda_sequence_generation_death_changes[],
    //                             float **cuda_sequence_replication_prob_changes[],
    //                             float **cuda_sequence_metastatic_prob_changes[],
    //                             int *CUDA_device_IDs,
    //                             string &intermediary_Tissue_folder, string &rapid_Progeny_Location,
    //                             int &last_index_Seq_Written, int &overall_Generations, string &tissue_Name,
    //                             int &tissue_Index, mt19937 &gen, string &source_sequence_Data_folder, int &last_Progeny_written_this_Gen,
    //                             int &max_sequences_per_File,
    //                             int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List);

    void cell_Migration_set(int &max_Limit, multiset<pair<float, int>> &migration_cell_List, pair<float, int> candidate_Cell);

    void migration_of_Cells(string &source_sequence_Data_folder, vector<string> &tissue_Names,
                            int source, vector<pair<int, int>> tissue_migration_Targets_amount, multiset<pair<float, int>> &migration_cell_List,
                            int &overall_Generations,
                            functions_library &functions);

    void thread_find_Files_2(int start, int stop, vector<int> &cell_Indexes, vector<pair<int, int>> &indexed_Tissue_Folder);

    void calculate_Tajima(functions_library &functions,
                          int &num_Regions, int **tajima_regions_Start_Stop,
                          int &overall_Generations, string &tissue_Name, int &tissue_Index,
                          int *CUDA_device_IDs,
                          string sequence_Tissue_Folder, int &max_Cells_at_a_time,
                          string output_Tajima_File,
                          string reference_Genome_location);

    void process_Tajima_String(string &all_Sequences, int &count_Track, int &num_Regions, int **cuda_tajima_regions_Start_Stop,
                               char *cuda_Reference_Genome, int *cuda_per_Region,
                               functions_library &functions);

    long int combos_N(int count);
    long int fact_half(int count);
};