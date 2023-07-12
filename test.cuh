#pragma once
#include <iostream>

#include "cuda_runtime.h"
#include <curand.h>
#include "device_launch_parameters.h"
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
#include <sstream>

#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

class test
{
private:
public:
    test();
};