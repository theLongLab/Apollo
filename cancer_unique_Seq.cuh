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

class cancer_unique_Seq
{
private:
    string intermediate_Folder_location;
    string output_Folder_location;

    struct sequence_Details
    {
        int count = 0;
        string Replication_Prob;
        string Gen_Death_prob;
        string Replication_Factor;
        string metastatic_Prob;
        string survivability;
    };

public:
    cancer_unique_Seq(string parameter_Master_Location);

    void ingress();
};