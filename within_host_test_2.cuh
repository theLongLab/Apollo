#pragma once
#include <iostream>

#include "functions_library.cuh"
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

class within_host_test_2
{
private:
    mt19937 gen;

    shared_mutex g_mutex;
    mutex gi_mutex;

    string references_Folder;

    string replication_profile_file;
    string sequence_profile_file;

    int CUDA_device_number;
    int gpu_Limit;
    int at_a_Time_cells = 1;
    int tot_Blocks;
    int tot_ThreadsperBlock;

    int CPU_cores;
    string output_Folder;
    string intermediate_Folders;
    string intermediate_sequence_Store;
    string intermediate_profile_Store;
    string multi_READ = "NO";

    string write_Progeny_parent = "NO";
    string progeny_Parent_folder;
    string progeny_File;
    string progeny_Recombination_File;
    string sequence_Profiles;
    string write_Sequences = "NO";
    string progeny_sequences_Folder;

    int genome_SIZE;

    float cell_scale, cell_shape = 0;
    float cell_r, cell_prob = 0;

    int cell_Limit = -1;
    float total_Cells_shape = 0;
    float total_Cells_scale = 0;

   // int cells_Available = -1;

    // float progeny_scale, progeny_shape = 0;
    // //string progeny_distribution_Type;
    // float progeny_mean, progeny_dispersion = 0;

    int mutation_Activate_parent = 0;
    int recombination_Activate_parent = 0;
    int proof_reading_Activate_parent = 0;

    int **sequence_Mutation_tracker;

    int mutation_hotspots = -1;
    float **mutation_rates_Hotspot_generation;
    // vector<pair<int, int>> mutation_Regions;
    int **mutation_Regions_start_stop;

    // 0 = Growth
    // 1 = Stationary
    // 2 = Depriciation
    int *generation_modes;
    int *phase_Modes;
    float **phase_parameters;

    vector<pair<float, int>> index_Fitness;

    // A = 0
    // T = 1
    // G = 2
    // C = 3

    float **A_0_mutation;
    float **T_1_mutation;
    float **G_2_mutation;
    float **C_3_mutation;

    int fitness_points = -1;
    //! might not need
    int *fitness_point_Locations;

    float **A_0_fitness;
    float **T_1_fitness;
    float **G_2_fitness;
    float **C_3_fitness;

    int survivability_points = -1;
    int *survivability_point_Locations;

    float **A_0_survivability;
    float **T_1_survivability;
    float **G_2_survivability;
    float **C_3_survivability;

    int recombination_hotspots = -1;
    int **recombination_hotspots_start_stop;
    // float *recombination_probability;

    int **recomb_effects_prob_selectivity_fitness;

    vector<vector<string>> effect_type_All_hotspots;
    vector<vector<int>> effect_positions_All_hotspots;
    vector<vector<vector<string>>> effect_base_mutations_All_hotspots;

    int rows_Prob = 0;
    int rows_Selectivity = 0;
    int rows_Fitness = 0;
    int rows_Survivability = 0;

    float **A_0_fitness_Recombination;
    float **T_1_fitness_Recombination;
    float **G_2_fitness_Recombination;
    float **C_3_fitness_Recombination;

    float **A_0_selectivity_Recombination;
    float **T_1_selectivity_Recombination;
    float **G_2_selectivity_Recombination;
    float **C_3_selectivity_Recombination;

    float **A_0_probability_Recombination;
    float **T_1_probability_Recombination;
    float **G_2_probability_Recombination;
    float **C_3_probability_Recombination;

    float **A_0_survivability_Recombination;
    float **T_1_survivability_Recombination;
    float **G_2_survivability_Recombination;
    float **C_3_survivability_Recombination;

    int proof_Reading_mutations = -1;

    //! might not need
    int *position_Proof_reading_mutations;

    float **A_0_probability_Proof_reading;
    float **T_1_probability_Proof_reading;
    float **G_2_probability_Proof_reading;
    float **C_3_probability_Proof_reading;

    float *sequences_Proof_reading_probability;
    float **sequences_Survivability;

    float **current_gen_Parent_data;

public:
    within_host_test_2(int CUDA_device_number, int CPU_cores, int gpu_Limit,
                       string output_Folder,
                       string reference_genome_Location, string replication_profile_file,
                       mt19937 gen,
                       string sequence_Profile_file,
                       string write_Progeny_parents, string write_Sequences,
                       string intermediate_Folders,
                       string multi_READ, int at_a_Time_cells);

    void ingress(float rep_time, float host_days, string mode);

    float *parent_Fitness(int &parents_in_current_generation, float **CUDA_current_gen_Parent_data);
    // int **standard_Progeny_Numbers(int &parents_in_current_generation, float *cuda_Parental_Fitness, functions_library &function, float **CUDA_current_gen_Parent_data);

    void allocate_Phases(int num_Generations);
    void configure_Mutation_Profiles(int num_Generations);
    void fitness_Profiles();
    void survivability_Profiles();
    void configure_Recombination_Profiles();
    void configure_proof_reading_Profiles();
    void configure_sequence_Profile(int parents_in_current_generation, vector<string> parent_headers);

    void get_Cell_distribution(string &distribution_Type);
    void get_Progeny_distribution(functions_library &function);

    void read_Profiles_multi_Thread(int start, int stop,
                                    int generation_Current, string parent_Profiles_Store, int start_Index, vector<int> surviving_Progeny);

    struct CustomComparator
    {
        bool operator()(const pair<float, int> &a, const pair<float, int> &b) const
        {
            return a.first > b.first;
        }
    };
};