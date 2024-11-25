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

class extract_seq
{
private:
    mt19937 gen;

    shared_mutex g_mutex;
    mutex gi_mutex;

    int CUDA_device_number;
    int gpu_Limit;
    int tot_Blocks;
    int tot_ThreadsperBlock;

    int CPU_cores;
    string output_Folder;
    string output_Folder_Sequences;
    string intermediate_Folders;
    string intermediate_sequence_Store;

    string multi_READ = "NO";

    vector<string> N_sequences_Store;

    vector<pair<string, string>> unique_Sequences;
    vector<int> sequence_Count;

public:
    extract_seq(int CUDA_device_number, int CPU_cores, int gpu_Limit,
                string intermediate_Folders,
                string output_Folder,
                string multi_READ);

    void ingress();

    void read_N_fasta(int start, int stop, int start_Index, int generation, string folder);
};