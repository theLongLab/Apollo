#include "distribution_test.cuh"
#include "functions_library.cuh"

distribution_test::distribution_test()
{
    functions_library function = functions_library();
    cout << "Listing all CUDA capable devices:" << endl;
    // int nDevices;
    // cudaGetDeviceCount(&nDevices);
    // for (int i = 0; i < nDevices; i++)
    // {

    //     cout << endl;
    //     cudaDeviceProp prop;
    //     cudaGetDeviceProperties(&prop, i);
    //     cout << "GPU number\t: " << i << endl;
    //     cout << "GPU name\t: " << prop.name << endl;
    //     size_t l_free = 0;
    //     size_t l_Total = 0;
    //     cudaError_t error_id = cudaMemGetInfo(&l_free, &l_Total);
    //     cout << "GPU memory (GB)\t: " << l_Total / (1000 * 1000 * 1000) << endl;
    //     cout << "GPU number of multiprocessor(s)\t: " << prop.multiProcessorCount << endl;
    //     cout << "GPU block(s) per multiprocessor\t: " << prop.maxBlocksPerMultiProcessor << endl;
    //     cout << "GPU thread(s) per block\t: " << prop.maxThreadsPerBlock << endl;
    //     cout << endl;

    //     cudaError_t err = cudaGetLastError();

    //     if (err != cudaSuccess)
    //     {
    //         printf("CUDA Error: %s\n", cudaGetErrorString(err));
    //     }
    // }

    function.print_Cuda_device(0, this->tot_Blocks, this->tot_ThreadsperBlock);
}

__global__ void setup(curandState_t *d_states)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < (100 * 100))
    {
        curand_init(tid + 7, tid, 0, &d_states[tid]);
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void normal(curandState_t *d_states, double *values)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < (100 * 100))
    {
        curandState_t localState = d_states[tid];
        double z = (double)(curand_poisson(&localState, 2));
        values[tid] = z;
        tid += blockDim.x * gridDim.x;
    }
}

void distribution_test::ingress()
{
    this->gpu_Limit = 1000;
    this->CPU_cores = 5;

    functions_library function = functions_library(this->tot_Blocks, this->tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    random_device rd;
    mt19937 gen(rd());

    cout << "DT testing" << endl
         << endl;

    int number = 200005;

    // float *normal_distribution = function.normal_distribution_CUDA(number, 100, 2);

    // int *poisson_distribution = function.poisson_distribution_CUDA(number, 0.039);

    int *neg_binomial_distribution = function.negative_binomial_distribution_CUDA(number, 3, 5);

    fstream results;
    results.open("result_test_dis.csv", ios::out);

    for (size_t i = 0; i < number; i++)
    {
        // cout << neg_binomial_distribution[i] << "\n";
        results << neg_binomial_distribution[i] << "\n";
    }

    results.close();

    free(neg_binomial_distribution);

    cout << endl;
}