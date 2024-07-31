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

class cancer
{
private:
    shared_mutex g_mutex;
    mutex gi_mutex;

    string output_Folder_location;
    string intermediate_Folder_location;

    string enable_Folder_management = "NO";
    string enable_Compression = "NO";

    int max_sequences_per_File = 0;
    int max_Cells_at_a_time = 0;

    string start_Date;
    string stop_after_generations = "NO";
    string first_Infection = "Random";

    string node_Master_location;
    string sequence_Master_location;

    int stop_gen_Mode = 0;
    int stop_generations_Count = 0;
    float stop_Date = 0;

    int CPU_cores;
    string multi_Read;

    int *CUDA_device_IDs;
    int num_Cuda_devices;
    int gpu_Limit;
    int *tot_Blocks;
    int *tot_ThreadsperBlock;

public:
    cancer(string parameter_Master_Location);

    void ingress();
};