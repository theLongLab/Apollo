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
#include <unordered_map>

#include <string>
#include <vector>
#include <queue>

#include <thread>
#include <mutex>
#include <shared_mutex>

using namespace std;

class get_Variants
{
private:
    string reference_Sequence;
    vector<string> variant_Lines;

    int num_Threads;

    int lines_To_process;

    string output_File_location;

public:
    get_Variants(string reference_Sequence_location, string variant_File_location, string output_File_location);

    void ingress();

    void process_Lines(int thread_ID, int start, int stop, functions_library &functions);
};