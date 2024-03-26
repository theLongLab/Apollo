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

class pedigree2r
{
private:
    string pedigree_Folder_location = "";
    vector<pair<string, string>> files_to_Process;

public:
    pedigree2r(string parameter_Master_Location);

    void ingress();
};