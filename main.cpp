#include <iostream>
#include <filesystem>

#include <random>

#include <algorithm>

#include <unistd.h>
#include <bits/stdc++.h>
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "distribution_test.cuh"
#include "test.cuh"
#include "within_host_test.cuh"

#include "within_host_test_2.cuh"
#include "network.cuh"

#include "parameter_load.h"

#include "functions_library.cuh"
#include "extract_seq.cuh"

#include "simulator_Master.cuh"

using namespace std;

int main(int argc, char *argv[])
{
     // upload CHECK NOW
     // Hello

     functions_library function_main = functions_library();

     cout << "Simulator\n\n";

     // distribution_test dt = distribution_test();
     // dt.ingress();

     // exit(-1);

     string function(argv[1]);

     /**
      * * Functions are converted to lowercase formats so that they will not be case sensetive.
      **/
     transform(function.begin(), function.end(), function.begin(), ::tolower);

     string parameter_MASTER_file(argv[2]);

     cout << "Selected function: " << function << endl;
     cout << "Master parameter file location: " << parameter_MASTER_file << endl
          << endl;

     // int seed_value = rd();

     parameter_load Parameters = parameter_load();

     if (function == "--simulator")
     {
          cout << "Simulator for within host viral replication\n\n";

          simulator_Master simulator = simulator_Master(parameter_MASTER_file);
          simulator.ingress();
     }
     else if (function == "--extract")
     {
          cout << "Extraction\n";

          vector<string> parameters_List = {
              "\"CUDA Device ID\"",
              "\"CPU cores\"",
              "\"GPU max units\"",
              "\"Intermediate folders\"",
              "\"Output folders\"",
              "\"Multi read\""};

          vector<string> found_Parameters = Parameters.get_parameters(parameter_MASTER_file, parameters_List);

          string output_Folders = Parameters.get_STRING(found_Parameters[4]);
          function_main.config_Folder(output_Folders, "Output");

          //     extract_seq(int CUDA_device_number, int CPU_cores, int gpu_Limit,
          //                 string intermediate_Folders,
          //                 string output_Folder,
          //                 string multi_READ);

          extract_seq e_seq = extract_seq(Parameters.get_INT(found_Parameters[0]), Parameters.get_INT(found_Parameters[1]), Parameters.get_INT(found_Parameters[2]), Parameters.get_STRING(found_Parameters[3]), output_Folders, Parameters.get_STRING(found_Parameters[5]));
          e_seq.ingress();

          exit(-1);
     }
     else if (function == "--network")
     {
          vector<string> parameters_List = {
              "\"CUDA Device ID\"",
              "\"CPU cores\"",
              "\"GPU max units\"",
              "\"Multi read\""};
          // int CUDA_device_number, int CPU_cores, int gpu_Limit
          cout << "Network test\n";

          vector<string> found_Parameters = Parameters.get_parameters(parameter_MASTER_file, parameters_List);

          network net = network(Parameters.get_INT(found_Parameters[0]), Parameters.get_INT(found_Parameters[1]), Parameters.get_INT(found_Parameters[2]), Parameters.get_STRING(found_Parameters[3]), 1);
          net.ncbi_gene_Count();
          exit(-1);
     }

     exit(-1);

     vector<string> parameters_List = {
         "\"CUDA Device ID\"",
         "\"Seed\"",
         "\"Parent sequences folder\"",
         "\"Replication profile folder\"",
         "\"Sequence profile folder\"",
         "\"CPU cores\"",
         "\"GPU max units\"",
         "\"Output folders\"",
         "\"Write progeny parent\"",
         "\"Write progeny sequences\"",
         "\"Intermediate folders\"",
         "\"Multi read\"",
         "\"Number of cells\"",
         "\"Mode\""};

     // cout << parameters_List.size();

     vector<string> found_Parameters = Parameters.get_parameters(parameter_MASTER_file, parameters_List);

     // for (string p : found_Parameters)
     // {
     //      cout << p << endl;
     // }

     int CUDA_Device = Parameters.get_INT(found_Parameters[0]);
     int CPU_cores = Parameters.get_INT(found_Parameters[5]);
     int GPU_max_units = Parameters.get_INT(found_Parameters[6]);

     string output_Folders = Parameters.get_STRING(found_Parameters[7]);
     string intermediate_Folders = Parameters.get_STRING(found_Parameters[10]);

     string multi_READ = Parameters.get_STRING(found_Parameters[11]);
     int at_a_time_Cells_simulate = Parameters.get_INT(found_Parameters[12]);
     string mode = Parameters.get_STRING(found_Parameters[13]);

     function_main.config_Folder(output_Folders, "Output");
     function_main.config_Folder(intermediate_Folders, "Intermediate");

     int seed_value;

     if (found_Parameters[1] == "\"x\"")
     {
          random_device rd;
          seed_value = rd();
     }
     else
     {
          seed_value = stoi(found_Parameters[1]);
     }

     cout << "Seed value: " << seed_value << endl
          << endl;

     string parent_SEQ_folder = Parameters.get_STRING(found_Parameters[2]);
     string replication_Profile_folder = Parameters.get_STRING(found_Parameters[3]);
     string sequence_Profile_folder = Parameters.get_STRING(found_Parameters[4]);
     string write_Progeny_parents = Parameters.get_STRING(found_Parameters[8]);
     string write_Sequences = Parameters.get_STRING(found_Parameters[9]);

     transform(write_Progeny_parents.begin(), write_Progeny_parents.end(), write_Progeny_parents.begin(), ::toupper);
     transform(write_Sequences.begin(), write_Sequences.end(), write_Sequences.begin(), ::toupper);
     transform(multi_READ.begin(), multi_READ.end(), multi_READ.begin(), ::toupper);
     transform(mode.begin(), mode.end(), mode.begin(), ::toupper);

     parameters_List.clear();
     found_Parameters.clear();

     cout << "Collecting replication profiles from: " << replication_Profile_folder << endl;
     vector<string> replication_profile_Files = function_main.get_Files(replication_Profile_folder, "json");

     parameters_List = {"\"Shape replication time\"", "\"Scale replication time\"",
                        "\"Shape days in host\"", "\"Scale days in host\""};

     found_Parameters = Parameters.get_parameters(replication_profile_Files[0], parameters_List);

     // cout << parent_SEQ_folder << endl
     //      << replication_Profile_folder << endl;

     mt19937 gen(seed_value);

     float shape_rep_time = Parameters.get_FLOAT(found_Parameters[0]);
     float scale_rep_time = Parameters.get_FLOAT(found_Parameters[1]);

     float shape_days_host = Parameters.get_FLOAT(found_Parameters[2]);
     float scale_deviation_host_time = Parameters.get_FLOAT(found_Parameters[3]);

     // parameter_load Parameters = parameter_load(parameter_file);
     // Parameters.get_parameters(CUDA_Device, parent_SEQ_folder,
     //                           mean_rep_time, standard_deviation_rep_time,
     //                           mean_days_host, standard_deviation_host_time);

     cout << endl;

     cout << "Shape replication time: " << shape_rep_time << endl;
     cout << "Scale replication time: " << scale_rep_time << endl
          << endl;

     cout << "Shape days in host: " << shape_days_host << endl;
     cout << "Scale days in host: " << scale_deviation_host_time << endl
          << endl;

     gamma_distribution<float> distribution_rep(shape_rep_time, scale_rep_time);
     gamma_distribution<float> distribution_host(shape_days_host, scale_deviation_host_time);

     cout << "Collecting sequence profiles from: " << sequence_Profile_folder << endl;
     vector<string> sequence_profile_Files = function_main.get_Files(sequence_Profile_folder, "json");
     cout << endl;
     // cout << distribution(gen) << endl;

     // within_host_test_2(int CUDA_device_number,int CPU_cores,int gpu_Limit, string output_Folder, string reference_genome_Location, string replication_profile_file, mt19937 gen, string sequence_Profile_file);

     within_host_test_2 wht_2 = within_host_test_2(CUDA_Device, CPU_cores, GPU_max_units, output_Folders, parent_SEQ_folder, replication_profile_Files[0], gen, sequence_profile_Files[0], write_Progeny_parents, write_Sequences, intermediate_Folders, multi_READ, at_a_time_Cells_simulate);
     wht_2.ingress(distribution_rep(gen), distribution_host(gen), mode);

     // distribution_test distribution = distribution_test();
     // distribution.ingress();

     // test testing = test();

     // string reference_Genome_location = "reference_Genomes/reference_Sequence.fasta";
     // int cuda_Device_ID = 0;

     // // within_host_test(int cuda_Device, string reference_Genome_location)
     // within_host_test wht = within_host_test(cuda_Device_ID, reference_Genome_location);
     // wht.ingress();

     cout << "\nSimulator concluded sucessfully" << endl;

     return 0;
}