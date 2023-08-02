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

    fstream sequences_files_Summary;
    sequences_files_Summary.open(output_Folder_Sequences + "/sequences_Summary.csv", ios::out);
    sequences_files_Summary << "First_occurrence\tSequence\tCount\n";
    sequences_files_Summary.close();

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
                cout << sequence_Files_count << " sequence(s) found.\n"
                     << endl;

                if (multi_READ == "YES")
                {
                    cout << "Intiating multi read of N_sequences.\n";

                    int full_Rounds = sequence_Files_count / this->gpu_Limit;
                    int partial_Rounds = sequence_Files_count % this->gpu_Limit;

                    vector<pair<int, int>> start_stops;

                    for (int full = 0; full < full_Rounds; full++)
                    {
                        int start = full * this->gpu_Limit;
                        int stop = start + this->gpu_Limit;
                        start_stops.push_back(make_pair(start, stop));
                    }

                    if (partial_Rounds != 0)
                    {
                        int start = sequence_Files_count - partial_Rounds;
                        start_stops.push_back(make_pair(start, sequence_Files_count));
                    }

                    for (size_t i = 0; i < start_stops.size(); i++)
                    {
                        vector<thread> threads_vec;
                        int num_of_values_current = start_stops[i].second - start_stops[i].first;
                        int num_per_Core = num_of_values_current / this->CPU_cores;
                        int remainder = num_of_values_current % this->CPU_cores;

                        for (int sequence = 0; sequence < num_of_values_current; sequence++)
                        {
                            N_sequences_Store.push_back("");
                        }

                        for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
                        {
                            int start_Cell = core_ID * num_per_Core;
                            int stop_Cell = start_Cell + num_per_Core;

                            threads_vec.push_back(thread{&extract_seq::read_N_fasta, this, start_Cell, stop_Cell, start_stops[i].first, folder_Index, folder});
                        }

                        if (remainder != 0)
                        {
                            int start_Cell = num_of_values_current - remainder;
                            int stop_Cell = num_of_values_current;

                            threads_vec.push_back(thread{&extract_seq::read_N_fasta, this, start_Cell, stop_Cell, start_stops[i].first, folder_Index, folder});
                        }

                        for (thread &t : threads_vec)
                        {
                            if (t.joinable())
                            {
                                t.join();
                            }
                        }

                        threads_vec.clear();

                        cout << endl;

                        // for (size_t i = 0; i < N_sequences_Store.size(); i++)
                        // {
                        //     cout << N_sequences_Store[i] << endl
                        //          << endl;
                        // }

                        for (int query = 0; query < N_sequences_Store.size(); query++)
                        {
                            int found = 0;

                            for (int check = 0; check < unique_Sequences.size(); check++)
                            {
                                if (N_sequences_Store[query] == unique_Sequences[check].second)
                                {
                                    sequence_Count[check] = sequence_Count[check] + 1;
                                    found = 1;
                                    break;
                                }
                            }

                            if (found == 0)
                            {
                                string sequence_ID = to_string(folder_Index) + "_" + to_string(query + start_stops[i].first);

                                unique_Sequences.push_back(make_pair(sequence_ID, N_sequences_Store[query]));
                                sequence_Count.push_back(1);
                            }
                        }
                        N_sequences_Store.clear();
                    }
                }
                else
                {
                    cout << "Intiating single read of N_sequences.\n";
                    for (int sequence = 0; sequence < sequence_Files_count; sequence++)
                    {
                        string sequence_Name = to_string(folder_Index) + "_" + to_string(sequence) + ".nfasta";
                        cout << "Processing sequence : " << folder << "/" << sequence_Name;

                        fstream sequence_File;
                        sequence_File.open(folder + "/" + sequence_Name, ios::in);

                        string sequence_Full = "";

                        if (sequence_File.is_open())
                        {
                            string line;

                            while (getline(sequence_File, line))
                            {
                                sequence_Full.append(line);
                            }
                            sequence_File.close();
                            N_sequences_Store.push_back(sequence_Full);
                        }
                    }
                }
            }
            else
            {
                cout << "Failed to extract the sequence folder: " << folder << endl;
                exit(-1);
            }

            string command_remove = "rm -r " + folder;

            int result_remove = system(command_remove.c_str());

            if (result_remove == 0)
            {
                cout << "Purged the sequence folder: " << folder << endl
                     << endl;
            }
            else
            {
                cout << "Failed to purge the sequence folder: " << folder << endl
                     << endl;
                exit(-1);
            }
            // if (folder_Index == 1)
            // {
            //     break;
            // }
        }

        cout << "Writing results of sequence summary\n";

        sequences_files_Summary.open(output_Folder_Sequences + "/sequences_Summary.csv", ios::app);

        for (int unique = 0; unique < unique_Sequences.size(); unique++)
        {
            sequences_files_Summary << unique_Sequences[unique].first << "\t";

            for (int base = 0; base < unique_Sequences[unique].second.size(); base++)
            {
                if (unique_Sequences[unique].second.at(base) == '0')
                {
                    sequences_files_Summary << "A";
                }
                else if (unique_Sequences[unique].second.at(base) == '1')
                {
                    sequences_files_Summary << "T";
                }
                else if (unique_Sequences[unique].second.at(base) == '2')
                {
                    sequences_files_Summary << "G";
                }
                else if (unique_Sequences[unique].second.at(base) == '3')
                {
                    sequences_files_Summary << "C";
                }
            }

            sequences_files_Summary << "\t";

            sequences_files_Summary << sequence_Count[unique] << "\n";
        }

        sequences_files_Summary.close();

        cout << "Completed extraction\n";
    }
    else
    {
        cout << "No generational sequence data was found.\n";
    }
}

void extract_seq::read_N_fasta(int start, int stop, int start_Index, int generation, string folder)
{
    vector<string> N_sequences_Store_partial;

    // cout << "check\n";

    for (int i = start; i < stop; i++)
    {
        int sequence = i + start_Index;

        string sequence_Name = to_string(generation) + "_" + to_string(sequence) + ".nfasta";

        cout << "Processing sequence : " << folder << "/" << sequence_Name << endl;

        fstream sequence_File;
        sequence_File.open(folder + "/" + sequence_Name, ios::in);

        string sequence_Full = "";

        if (sequence_File.is_open())
        {
            string line;

            while (getline(sequence_File, line))
            {
                sequence_Full.append(line);
            }
            sequence_File.close();
            N_sequences_Store_partial.push_back(sequence_Full);
        }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (int i = start; i < stop; i++)
    {
        N_sequences_Store[i] = N_sequences_Store_partial[i - start];
    }
}