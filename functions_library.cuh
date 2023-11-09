#pragma once
#include <iostream>
#include <filesystem>
#include <fstream>

#include "cuda_runtime.h"
#include <curand.h>
#include "device_launch_parameters.h"
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
#include <sstream>

#include <curand_kernel.h>

#include <algorithm>
#include <random>
#include <chrono>

#include <string>
#include <vector>

#include <thread>
#include <mutex>
#include <shared_mutex>

#include <set>

using namespace std;

class functions_library
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;
    // vector<T> print;

    int tot_Blocks;
    int tot_ThreadsperBlock;

    int *CUDA_device_IDs;
    int num_Cuda_devices;
    int *tot_Blocks_array;
    int *tot_ThreadsperBlock_array;

    int gpu_Limit;
    int CPU_cores;

    set<int> parent_Indexes_Unique;
    vector<string> sequences;

public:
    string mode = "PARENT";

    int genome_Size;

    string cells_of_parents;
    string cells_of_progeny;

    string progeny_distribution_Type;
    float progeny_scale, progeny_shape = 0;
    float progeny_mean, progeny_dispersion = 0;

    float progeny_prob;
    int progeny_r;

    int recombination_hotspots = 0;
    int proof_reading_Activate_parent = 0;

    int mutation_hotspots = 0;
    float **CUDA_mutation_rates_Hotspot_generation;
    int **CUDA_mutation_Regions_start_stop;

    float **CUDA_A_0_mutation;
    float **CUDA_T_1_mutation;
    float **CUDA_G_2_mutation;
    float **CUDA_C_3_mutation;

    int *CUDA_stride_Array;
    int **CUDA_recombination_hotspots_start_stop;

    float **CUDA_A_0_Recombination;
    float **CUDA_T_1_Recombination;
    float **CUDA_G_2_Recombination;
    float **CUDA_C_3_Recombination;

    float **CUDA_A_0_fitness;
    float **CUDA_T_1_fitness;
    float **CUDA_G_2_fitness;
    float **CUDA_C_3_fitness;

    float **CUDA_A_0_survivability;
    float **CUDA_T_1_survivability;
    float **CUDA_G_2_survivability;
    float **CUDA_C_3_survivability;

    float **CUDA_A_0_probability_Proof_reading;
    float **CUDA_T_1_probability_Proof_reading;
    float **CUDA_G_2_probability_Proof_reading;
    float **CUDA_C_3_probability_Proof_reading;

    int **CUDA_sequence_Mutation_tracker;

    functions_library();
    functions_library(int tot_Blocks, int tot_ThreadsperBlock, int gpu_Limit, int CPU_cores);
    functions_library(int *tot_Blocks_array, int *tot_ThreadsperBlock_array, int *CUDA_device_IDs, int num_Cuda_devices, int gpu_Limit, int CPU_cores);

    void print_Cuda_device(int cuda_ID, int &tot_Blocks, int &tot_ThreadsperBlock);
    void print_Cuda_devices(vector<string> cuda_IDs, int *CUDA_device_IDs, int num_Cuda_devices, int *tot_Blocks, int *tot_ThreadsperBlock);

    void config_Folder(string location, string type_Folder);
    void create_File(string location, string headers);
    void create_File(string location);
    void config_File_crash(string location, string headers);
    void config_File_delete_create(string location, string headers);
    void config_File_delete_create(string location);

    void config_File_progeny(string location);

    string read_Reference(string file_location, string &header, int &genome_Size);
    string read_Reference(string file_location, int &genome_Size);

    void split(vector<string> &line_Data, string line, char delim);

    string to_Upper_Case(const string &text);

    vector<string> get_Files(string folder, string extension);

    int **create_INT_2D_arrays(int rows, int columns);
    float **create_FLOAT_2D_arrays(int rows, int columns);

    int **create_Fill_2D_array(int rows, int columns, int fill_Value);
    float **create_Fill_2D_array_FLOAT(int rows, int columns, float fill_Value);

    int **Fill_2D_array_CUDA(int rows, int columns, int fill_Value, int **cuda_Progeny_numbers);

    int **create_Progeny_Array(int parents_in_current_generation, int *stride_Progeny_Index, int total_Progeny, int num_Hotspots, int **cuda_Progeny_numbers, int fill_Value);
    void progeny_Recombination_parents_array(int **progeny_recom_Index_Cuda, int total_Progeny, int num_Hotspots,
                                             int *parent_and_their_cell_CUDA, int **cell_and_their_viruses_CUDA, int *per_Cell_max_viruses_CUDA,
                                             float **CUDA_current_gen_Parent_data, int max_Count, int num_Unique_cells);
    void progeny_Shuffle(int **progeny_recom_Index_Cuda, int num_Hotspots, int parents_in_current_generation, int *stride_Progeny_Index_CUDA);

    float **create_current_Progeny_data(int total_Progeny, int num_Hotspots, int **progeny_recom_Index_Cuda, float **CUDA_current_gen_Parent_data, float *CUDA_parent_Proof_reading_probability, float *CUDA_progeny_Proof_reading_probability);

    void array_Copy_Float(float **target_2D, float *D1_array, int row, int num_Values);

    float *normal_distribution_CUDA(int num_of_values, float mean, float st_deviation);
    int *poisson_distribution_CUDA(int num_of_values, float mean);
    float *gamma_distribution_CUDA(int num_of_values, float shape, float scale);
    int *binomial_distribution_CUDA(int num_of_values, float prob, int progeny_Number);
    int *negative_binomial_distribution_CUDA(int num_of_values, float mean, float dispersion_Parameter);

    float beta_Distribution(float &alpha, float &beta, mt19937 &gen);

    void get_base_mutation(string query, char &base, int &mutation);

    float **float_2D_Array_load_to_CUDA(float **host_Array, int rows, int columns);
    int **int_2D_Array_load_to_CUDA(int **host_Array, int rows, int columns);

    int *copy_1D_to_CUDA_INT(int *host_Array, int num_Values);
    float *copy_1D_to_CUDA_FLOAT(float *host_Array, int num_Values);

    int **load_to_Host(int **cuda_2D_Array, int rows, int columns);
    float **load_to_Host_FLOAT(float **cuda_2D_Array, int rows, int columns);

    int **create_CUDA_2D_int(int rows, int columns);
    float **create_CUDA_2D_FLOAT(int rows, int columns);

    void clear_Array_INT(int **CUDA_2D_array, int rows);
    void clear_Array_FLOAT(float **CUDA_2D_array, int rows);

    void clear_Array_float_CPU(float **array_cpu, int rows);
    void clear_Array_int_CPU(int **array_cpu, int rows);

    // ! CONFIGURE (distribution_Type == "Negative binomial")
    int **progeny_distribution_CUDA(string &distribution_Type,
                                    int &num_of_parents,
                                    float &shape, float &scale,
                                    float &mean, float &dispersion_Parameter,
                                    float *cuda__Parent_finess,
                                    float **CUDA_current_gen_Parent_data, int recombination_hotspots);

    int sum_CUDA(int **cuda_Array_input, int num_Elements);

    void mutate_Sequences(int **cuda_progeny_Sequences, float **CUDA_current_gen_Progeny_data, float *CUDA_progeny_Proof_reading_probability,
                          int current_Generation, int num_Mutation_hotspots,
                          float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop, int **CUDA_sequence_Mutation_tracker,
                          int total_Progeny, int proof_Reading_Activate);

    int **create_progeny_Sequences(int **cuda_parent_Sequences, int **progeny_recom_Index_Cuda, int num_Hotspots, int total_Progeny, int genome_SIZE, int **CUDA_recombination_hotspots_start_stop);

    void find_Unique_values(int **CUDA_Array_2D, int total_Elements, int num_Recom_hotspots,
                            int start_Index,
                            string multi_READ, vector<string> &parent_IDs, string &parent_Sequences_Store,
                            int &genome_Size,
                            string &progeny_Sequences_Store,
                            int **CUDA_recombination_hotspots_start_stop,
                            int mutation_Activate,
                            float **CUDA_current_gen_Progeny_data, float *CUDA_progeny_Proof_reading_probability, int current_Generation, float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop, int **CUDA_sequence_Mutation_tracker, int proof_Reading_Activate,
                            vector<string> &progeny_IDs);

    void hard_Load_progeny(int generation, string &parent_Sequences_Store, int genome_Size,
                           vector<string> &parent_IDs,
                           int total_Progeny,
                           int **progeny_recom_Index_Cuda, int num_Hotspots, int **CUDA_recombination_hotspots_start_stop,
                           string multi_READ,
                           string &progeny_Sequences_Store, int mutation_Activate,
                           float **CUDA_current_gen_Progeny_data, float *CUDA_progeny_Proof_reading_probability, float **CUDA_mutation_rates_Hotspot_generation, int **CUDA_mutation_Regions_start_stop, int **CUDA_sequence_Mutation_tracker, int proof_Reading_Activate,
                           vector<string> &progeny_IDs);

    string read_nFASTA(string file_location, int genome_Size);
    void read_nFASTA_multi_Read(string parent_Sequences_Store, int genome_Size, int start, int stop, int *parent_Indexes, int generation_Current);
    void read_nFASTA_multi_Read_CELLS(string parent_Sequences_Store, int genome_Size, int start, int stop, int **parent_Indexes_parents_and_their_Cells, int generation_Current);

    void write_nFASTA_multi_WRITE(string progeny_Sequences_Store, int sum_Progeny_in_Generation, int genome_Size, int start, int stop, int **progeny_Sequences, int generation_Current);

    void unique_Collect_Threads(int *unique_values_Array, int start, int stop);

    void process_Cells(string &multi_READ, int &generation_Current, int &sum_Progeny_in_Generation,
                       int &num_of_Cells, int &start, int &stop,
                       vector<int> &cells_Start_Stop,
                       string &parent_Profiles_Store, int *parents,
                       string &parent_Sequences_Store,
                       string &progeny_File, string &progeny_Recombination_File, string &sequence_Profiles,
                       string &progeny_Sequences_Store, string progeny_Profile_Store,
                       int &processed_Cells, vector<int> &surviving_Progeny);

    void read_Profiles_Single(string file_Names, float **current_gen_Parent_data, float *parent_Proof_reading_probability, int &parent_Index);
    void read_Profiles_multi_Thread(int start, int stop, int cell_Start, int *parents,
                                    float **current_gen_Parent_data, float *parent_Proof_reading_probability,
                                    int generation_Current, string parent_Profiles_Store,
                                    int *parent_Indexes);

    void read_Profiles_multi_Thread_CELL(int start, int stop, int **parent_Indexes_parents_and_their_Cells,
                                         float **current_gen_Parent_data, float *parent_Proof_reading_probability,
                                         int generation_Current, string parent_Profiles_Store,
                                         float **parent_survivability_Probabilities);

    // void split(vector<string> &line_Data, string line, char delim);

    int get_base_Index(string base);

    int binary_Search(vector<int> &values, int value);

    string clean_Line(string &line);

    void process_Reference_Sequences(vector<string> collect_Sequences, int &genome_Length);
};