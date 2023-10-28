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
    string sequence_Master_location;

    string multi_Read;

    string network_Model = "NA";
    int Total_number_of_Nodes = 0;

    vector<vector<pair<int, int>>> each_Nodes_Connection;
    vector<pair<int, int>> all_node_IDs;
    vector<int> node_Profiles;

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

    string progeny_distribution_Model = "NA";

    int progeny_NB_sucesses = 0;
    float progeny_NB_probability = 0;

    float progeny_Poisson_mean = 0;

    float progeny_Gamma_shape = 0;
    float progeny_Gamma_scale = 0;

    int number_of_node_Profiles = 0;
    string node_Profile_folder_Location = "";
    float *node_profile_Distributions;
    vector<string> profile_names;

    // 0 = No Change
    // 1 = Removed
    //-1 = Less infectious
    float **node_sampling_effect;

    int num_tissues_per_Node = 0;
    vector<string> tissue_Names;

    int entry_tissues = 0;
    int *entry_array;

    int infectious_tissues = 0;
    int *infectious_array;

    int terminal_tissues = 0;
    int *terminal_array;

    string viral_Migration = "No";
    float **viral_Migration_Values;

    float **profile_tissue_Limits;

    int *num_replication_phases;
    float **tissue_replication_data;

    string parent_Sequence_Folder;
    // row 0 = mutations
    // row 1 = recombination
    // row 2 = proof reading
    // 0 = inactive 1 = activated
    int *mutation_recombination_proof_Reading_availability;

    float *Reference_fitness_survivability_proof_reading;

    // columns
    // 0 = position-1;
    // 1 = A
    // 2 = T
    // 3 = G
    // 4 = C
    float **sequence_Fitness_changes;
    float **sequence_Survivability_changes;
    float **sequence_Proof_reading_changes;

    // 0 = fitness
    // 1 = survivability
    // 2 = proof reading
    int *num_effect_Segregating_sites;

    int mutation_Hotspots = 0;
    int recombination_Hotspots = 0;

    float **A_0_mutation;
    float **T_1_mutation;
    float **G_2_mutation;
    float **C_3_mutation;

    float **mutation_hotspot_parameters;
    float **recombination_hotspot_parameters;

    int *tot_prob_selectivity;

    int *recombination_prob_Stride;
    int *recombination_select_Stride;

    float **recombination_Prob_matrix;
    float **recombination_Select_matrix;

    string output_Network_location = "";
    string network_File_location = "";

public:
    simulator_Master(string parameter_Master_Location);

    void configure_Network_Profile(string network_Profile_File, parameter_load &Parameters);

    void ingress();

    void network_Manager(functions_library &functions);

    void BA_Model_Engine(functions_library &functions);
    void SCM_Model_Engine(functions_library &functions);
    void DCM_Model_Engine(functions_library &functions);
    void RANDOM_Model_Engine(functions_library &functions);

    void node_Master_Manager(functions_library &functions);

    void sequence_Master_Manager(functions_library &functions);

    void node_Profile_assignment(functions_library &functions);
};