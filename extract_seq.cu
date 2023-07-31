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

    // vector<string> sequence_Folders;

    int num_Generations = 0;

    cout << "Intermediate sequences folder: " << this->intermediate_sequence_Store << endl;

    for (const auto &entry : filesystem::directory_iterator(intermediate_sequence_Store))
    {
        string file_Name = entry.path().filename();
        // cout << file_Name.substr(0, 11) << endl;
        if (file_Name.substr(0, 11) == "generation_")
        {
            // cout << entry.path() << '\n';
            // sequence_Folders.push_back(entry.path());
            num_Generations++;
        }
    }

    if (num_Generations != 0)
    {
        cout << endl
             << num_Generations << " generation(s) worth of sequence(s) were found.\n\n";

        for (int folder_Index = 0; folder_Index < num_Generations; folder_Index++)
        {
            string folder = intermediate_sequence_Store + "/generation_" + to_string(folder_Index);
            string command_Tar = "tar -xzf " + folder + ".tar.gz";

            int result = system(command_Tar.c_str());

            if (result == 0)
            {
                cout << "Extraction successful: " << folder + ".tar.gz" << endl;
                cout << "Processing folder: " << folder << endl
                     << endl;

                int sequence_Files_count = 0;

                for (const auto &entry : filesystem::directory_iterator(folder))
                {
                    string file_Name = entry.path().filename();

                    sequence_Files_count++;

                    // cout << "Processing sequence: " << file_Name << endl;
                }
                cout << sequence_Files_count << " sequence(s) found." << endl;
            }
            else
            {
                cout << "Failed to extract the sequence folder: " << folder << endl;
                exit(-1);
            }
            break;
        }
    }
    else
    {
        cout << "No generational sequence data was found.\n";
    }
}