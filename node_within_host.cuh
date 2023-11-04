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

class node_within_host
{
private:
    int host_Index, cave_ID, host_ID, profile_ID;
    int num_Generation;

    int infectious_Load;
    int terminal_Load;
    float sampling_Effect;

    int num_Tissues;
    int *cell_Limit;

    string status = "Susceptible";
    int current_Generation = 0;

public:
    node_within_host();

    void setHost(int host_Index, int cave_ID, int host_ID, int profile_ID);
    void setNum_Generation(int num_Generation);
    void setInfectious_Load(int infectious_Load);
    void setTerminal_Load(int terminal_Load);
    void setSampling_Effect(float sampling_Effect);
    void setCell_Limit(vector<int> cell_Limit_vec);

    void print_All();

    void begin_Infection();
    void run_Generation();
};