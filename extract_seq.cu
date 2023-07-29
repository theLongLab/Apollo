#include "extract_seq.cuh"
#include "functions_library.cuh"
#include "parameter_load.h"

extract_seq::extract_seq(int CUDA_device_number, int CPU_cores, int gpu_Limit,
                         string intermediate_Folders, string output_Folder,
                         string multi_READ)
{
    functions_library function = functions_library();
    cout << "Configuring sequence extraction\n\n";

    this->CPU_cores = CPU_cores;
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    this->multi_READ = multi_READ;
    cout << "Multiple read and write: " << this->multi_READ << endl
         << endl;

    this->CUDA_device_number = CUDA_device_number;
    function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

    this->gpu_Limit = gpu_Limit;

    cout << "Per round GPU max unit: " << this->gpu_Limit << endl
         << endl;

    this->output_Folder = output_Folder;
    this->output_Folder_Sequences = output_Folder + "/within_host_Sequences";
    function.config_Folder(this->output_Folder_Sequences, "Within host sequences");
    this->intermediate_Folders = intermediate_Folders;
    this->intermediate_sequence_Store = this->intermediate_Folders + "/sequences";
}

void extract_seq::ingress()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    int num_Generations;
}