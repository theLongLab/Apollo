#pragma once
#include <iostream>

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

using namespace std;

class within_host_test
{
private:
    string full_Reference_genome;
    string reference_Genome_header;
    int genome_Size;

    /**
     * @param tot_Blocks defines number of GPU blocks that are available
     **/
    int tot_Blocks;
    /**
     * @param tot_ThreadsperBlock defines number of threads that are available per GPU block.
     **/
    int tot_ThreadsperBlock;

public:
    within_host_test(int cuda_Device_ID, string reference_Genome_location);
    void ingress();
};