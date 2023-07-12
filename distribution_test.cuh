#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <curand.h>
#include <curand_kernel.h>

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
#include <sstream>

#include <iomanip>
#include <string>
#include <map>
#include <random>

using namespace std;

class distribution_test
{

private:
    int tot_Blocks;
    int tot_ThreadsperBlock;
    int gpu_Limit;

   int CPU_cores;

public:
    distribution_test();

    void ingress();
};