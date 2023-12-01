#include "node_within_host.cuh"

node_within_host::node_within_host()
{
    cout << "Intializing host: ";
}

void node_within_host::setHost(int host_Index, int cave_ID, int host_ID, int profile_ID, int num_Tissues)
{
    this->host_Index = host_Index;
    this->cave_ID = cave_ID;
    this->host_ID = host_ID;
    this->profile_ID = profile_ID;
    this->num_Tissues = num_Tissues;
    cout << this->cave_ID << "_" << host_ID << endl;
}

void node_within_host::setNum_Generation(int num_Generation)
{
    this->num_Generation = num_Generation;
}

void node_within_host::setInfectious_Load(int infectious_Load)
{
    this->infectious_Load = infectious_Load;
}

void node_within_host::setTerminal_Load(int terminal_Load)
{
    this->terminal_Load = terminal_Load;
}

void node_within_host::setSampling_Effect(float sampling_Effect)
{
    this->sampling_Effect = sampling_Effect;
}

void node_within_host::setCell_Limit(vector<int> cell_Limit_vec)
{
    num_Tissues = cell_Limit_vec.size();
    this->cell_Limit = (int *)malloc(sizeof(int) * num_Tissues);

    for (size_t i = 0; i < num_Tissues; i++)
    {
        cell_Limit[i] = cell_Limit_vec[i];
    }
}

void node_within_host::print_All()
{
    cout << host_Index << "\t"
         << cave_ID << "_" << host_ID << "\t"
         << profile_ID << "\t"
         << num_Generation << "\t"
         << infectious_Load << "\t"
         << terminal_Load << "\t"
         << sampling_Effect;

    for (size_t i = 0; i < num_Tissues; i++)
    {
        cout << "\t" << cell_Limit[i];
    }
    cout << endl;
}

void node_within_host::begin_Infection(functions_library &functions, string &intermediary_Sequence_location,
                                       int entry_tissues, int *entry_array, int &max_sequences_per_File,
                                       string &output_Node_location,
                                       vector<string> &tissue_Names)
{
    // FIRST NODE OF INFECTION IN THE HOST

    string host_Folder = intermediary_Sequence_location + "/" + to_string(host_Index);
    functions.config_Folder(host_Folder, to_string(cave_ID) + "_" + to_string(host_ID));

    vector<vector<string>> tissue_Sequences;
    intialize_Tissues(host_Folder, tissue_Sequences, functions);

    string reference_Sequences = intermediary_Sequence_location + "/reference_Sequences";

    vector<string> files;

    for (const auto &entry : filesystem::directory_iterator(reference_Sequences))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".nfasta")
        {
            files.push_back(entry.path().string());
        }
    }

    vector<string> Sequences;
    vector<string> Sequence_IDs;

    cout << endl;
    vector<char> seq_Status;
    for (int file = 0; file < files.size(); file++)
    {
        cout << "Reading file: " << files[file] << endl;
        fstream nfasta;
        nfasta.open(files[file], ios::in);

        if (nfasta.is_open())
        {
            string line;
            string sequence = "";

            while (getline(nfasta, line))
            {
                if (line.at(0) != '>')
                {
                    sequence.append(line);
                }
                else
                {
                    Sequence_IDs.push_back(line);
                    if (sequence != "")
                    {
                        Sequences.push_back(sequence);
                        sequence = "";
                    }
                }
            }

            if (sequence != "")
            {
                Sequences.push_back(sequence);
                sequence = "";
            }

            random_device rd; // Will be used to obtain a seed for the random number engine
            mt19937 gen(rd());
            uniform_int_distribution<int> entry_Tissue_select(0, entry_tissues - 1);

            cout << endl;

            for (int sequence = 0; sequence < Sequences.size(); sequence++)
            {
                int tissue_Index = entry_array[entry_Tissue_select(gen)];
                cout << "Sequence " << sequence + 1 << " infects tissue: " << tissue_Index << endl;
                tissue_Sequences[tissue_Index].push_back(Sequences[sequence]);
            }

            for (int tissue = 0; tissue < entry_tissues; tissue++)
            {
                if (tissue_Sequences[entry_array[tissue]].size() > 0)
                {
                    current_Viral_load_per_Tissue[entry_array[tissue]] = tissue_Sequences[entry_array[tissue]].size();
                    functions.config_Folder(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(entry_array[tissue]) + " Generation 0");

                    if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                    {
                        functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                        functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                        functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                    }

                    vector<string> sequence_Write_Store_All;
                    int last_seq_Num = 0;
                    functions.sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[entry_array[tissue]],
                                                          max_sequences_per_File, host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num, seq_Status,
                                                          output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation);
                    functions.partial_Write_Check(sequence_Write_Store_All,
                                                  host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num, seq_Status,
                                                  output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation);
                }
            }
        }
        else
        {
            cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << files[file] << endl;
            exit(-1);
        }
    }
    set_Infected();
    // for (int tissue = 0; tissue < num_Tissues; tissue++)
    // {
    //     cout << current_Viral_load_per_Tissue[tissue] << endl;
    // }
    // exit(-1);
}

string node_within_host::transfer_Infection(functions_library &functions, string &intermediary_Sequence_location, string &source_Target_file_Location,
                                            int &source_Index, int &source_Generation, string &source_Name, int *source_current_Viral_load_per_Tissue,
                                            int num_viruses_to_transfer,
                                            int &entry_tissues, int *entry_array, int exit_Load, int &exit_tissues, int *exit_array,
                                            vector<set<int>> &source_removed_by_Transfer_Indexes,
                                            int &max_sequences_per_File,
                                            vector<vector<pair<int, int>>> &indexed_Source_Folders,
                                            string &Host_source_target_network_location,
                                            string &output_Node_location,
                                            vector<string> &tissue_Names,
                                            mt19937 &gen)
{
    if (exit_Load > 0)
    {
        cout << "\nNode " << this->cave_ID << "_" << this->host_ID << " is being infected by " << source_Name << endl;

        if (num_viruses_to_transfer > exit_Load)
        {
            num_viruses_to_transfer = exit_Load;
        }

        if (num_viruses_to_transfer > 0)
        {
            string host_Folder = intermediary_Sequence_location + "/" + to_string(host_Index);

            if (current_Generation == -1)
            {
                functions.config_Folder(host_Folder, to_string(cave_ID) + "_" + to_string(host_ID));
                vector<vector<string>> tissue_Sequences;
                intialize_Tissues(host_Folder, tissue_Sequences, functions);
            }

            cout << "Attempting to transfer " << num_viruses_to_transfer << " viral particle(s)\n";

            uniform_int_distribution<> distribution_exit_Tissue(0, exit_tissues - 1);

            vector<set<int>> unique_indexes_to_Remove_Tissues;

            for (int tissue = 0; tissue < exit_tissues; tissue++)
            {
                set<int> init_Set;
                unique_indexes_to_Remove_Tissues.push_back(init_Set);
            }

            // cout << "Attempting to transfer " << num_viruses_to_transfer << " viral particle(s)\n";

            for (int particle = 0; particle < num_viruses_to_transfer; particle++)
            {
                int exit_tissue_Index = distribution_exit_Tissue(gen);
                if (source_current_Viral_load_per_Tissue[exit_array[exit_tissue_Index]] > 0)
                {
                    uniform_int_distribution<> distribution_particle(0, source_current_Viral_load_per_Tissue[exit_array[exit_tissue_Index]] - 1);
                    unique_indexes_to_Remove_Tissues[exit_tissue_Index].insert(distribution_particle(gen));
                }
            }

            cout << "Viral particle(s) and their exit tissue(s) have been indentifed\n";

            // vector<vector<int>> indexes_to_Remove;

            vector<vector<string>> seq_to_Write;
            vector<vector<string>> source_Seq_Data;
            for (int init = 0; init < entry_tissues; init++)
            {
                vector<string> initialize;
                seq_to_Write.push_back(initialize);
                source_Seq_Data.push_back(initialize);
            }

            uniform_int_distribution<> entry_Select(0, entry_tissues - 1);

            for (int tissue = 0; tissue < exit_tissues; tissue++)
            {
                vector<int> init_Tissue(unique_indexes_to_Remove_Tissues[tissue].begin(), unique_indexes_to_Remove_Tissues[tissue].end());

                if (init_Tissue.size() > 0)
                {
                    cout << "Exit tissue: " << exit_array[tissue] + 1 << endl;

                    vector<int> indexes_of_Seq_write;

                    for (int transfer_Cell = 0; transfer_Cell < init_Tissue.size(); transfer_Cell++)
                    {
                        auto it = source_removed_by_Transfer_Indexes[exit_array[tissue]].find(init_Tissue[transfer_Cell]);

                        if (it == source_removed_by_Transfer_Indexes[exit_array[tissue]].end())
                        {
                            // not present
                            indexes_of_Seq_write.push_back(init_Tissue[transfer_Cell]);
                            source_removed_by_Transfer_Indexes[exit_array[tissue]].insert(init_Tissue[transfer_Cell]);
                        }
                    }
                    if (indexes_of_Seq_write.size() > 0)
                    {
                        int valid_Sequences = 0;
                        // cout << "Collecting " << indexes_of_Seq_write.size() << " sequence(s)\n";
                        vector<string> collected_Sequences = functions.find_Sequences_Master(source_Target_file_Location, indexes_of_Seq_write, exit_array[tissue], indexed_Source_Folders[exit_array[tissue]], source_Generation, valid_Sequences);
                        cout << "Assinging sequence(s) to entry tissue(s)\n";

                        if (valid_Sequences != 0)
                        {
                            if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                            {
                                functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                                functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                                functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                            }
                        }

                        for (int check_Seq = 0; check_Seq < collected_Sequences.size(); check_Seq++)
                        {
                            if (collected_Sequences[check_Seq] != "")
                            {
                                int entry_Tissue_index = entry_Select(gen);
                                seq_to_Write[entry_Tissue_index].push_back(collected_Sequences[check_Seq]);
                                // sequence_Profile << host << "_" << tissue << "_" << last_seq_Num << "\t" << host << tissue<<endl;
                                source_Seq_Data[entry_Tissue_index].push_back(source_Name + "_" + tissue_Names[exit_array[tissue]] + "_" + to_string(source_Generation) + "_" + to_string(indexes_of_Seq_write[check_Seq]));
                            }
                        }
                    }
                }
            }
            vector<char> seq_Status;
            cout << "Writing sequence(s) to entry tissue(s)\n";
            int infected_Check = 0;
            for (int tissue = 0; tissue < entry_tissues; tissue++)
            {
                if (seq_to_Write[tissue].size() > 0)
                {
                    if (!filesystem::exists(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation)))
                    {
                        functions.config_Folder(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(entry_array[tissue]) + " Generation 0");
                    }

                    vector<string> sequence_Write_Store_All;
                    vector<int> indexes_Written;
                    functions.sequence_Write_Configurator_transfer(sequence_Write_Store_All, seq_to_Write[tissue],
                                                                   max_sequences_per_File, host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), current_Viral_load_per_Tissue[entry_array[tissue]], seq_Status,
                                                                   output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation,
                                                                   indexes_Written);
                    functions.partial_Write_Check_transfer(sequence_Write_Store_All,
                                                           host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), current_Viral_load_per_Tissue[entry_array[tissue]], seq_Status,
                                                           output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", get_Name(), tissue_Names[entry_array[tissue]], current_Generation, indexes_Written);
                    infected_Check = 1;

                    //(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                    fstream source;
                    fstream target;
                    source.open(output_Node_location + "/" + source_Name + "/sequence_parent_Progeny_relationships.csv", ios::app);
                    target.open(output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv", ios::app);
                    // functions.create_File(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", "Sequence_ID\tHost\tTissue");
                    fstream target_Profiles;
                    target_Profiles.open(output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv", ios::app);
                    fstream source_Profiles;
                    source_Profiles.open(output_Node_location + "/" + source_Name + "/sequence_Profiles.csv", ios::app);
                    for (int transfers = 0; transfers < indexes_Written.size(); transfers++)
                    {
                        source << source_Seq_Data[tissue][transfers] << "\t" << get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers]) << "\tTransmission" << endl;
                        target << source_Seq_Data[tissue][transfers] << "\t" << get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers]) << "\tTransmission" << endl;
                        target_Profiles << source_Seq_Data[tissue][transfers] << "\t" << source_Name << "\t" << tissue_Names[entry_array[tissue]] << endl;
                        source_Profiles << get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers]) << "\t" << get_Name() << "\t" << tissue_Names[entry_array[tissue]] << endl;
                    }
                    source.close();
                    target.close();
                    target_Profiles.close();
                }
            }
            if (infected_Check == 1)
            {
                if (status == "Susceptible")
                {
                    set_Infected();
                }

                fstream write_source_Target;
                write_source_Target.open(Host_source_target_network_location, ios::app);

                if (write_source_Target.is_open())
                {
                    cout << "Writing host's source target relationship\n";
                    write_source_Target << source_Name << "\t" << this->get_Name() << endl;
                    write_source_Target.close();
                }
                else
                {
                    cout << "ERROR: UNABLE TO OPEN SOURCE TARGET FILE: " << Host_source_target_network_location << "\n";
                    exit(-1);
                }
            }
        }
    }
    else
    {
        cout << source_Name << " has no viral particles in the exit tissues\n";
    }

    return status;
}

int node_within_host::get_generation_Phase(int generation, int *num_replication_phases, float **tissue_replication_data, int *tissue_param_profile_Stride, int &tissue,
                                           float &variable_1, float &variable_2)
{
    cout << "Getting generation phase\n";
    int gen_Phase = -1;

    int num_Phases = num_replication_phases[(profile_ID * num_Tissues) + tissue];

    int num_phases_per_tissue = 0;
    int tissue_Check = 0;

    for (int param_Index = tissue_param_profile_Stride[profile_ID]; param_Index < tissue_param_profile_Stride[profile_ID + 1]; param_Index++)
    {
        if (tissue_Check == tissue)
        {
            float time_Check = 0;
            float current_Generation_Ratio = (float)generation / (float)num_Generation;
            for (int phases = param_Index; phases < (param_Index + num_Phases); phases++)
            {
                time_Check = time_Check + tissue_replication_data[phases][0];
                if (current_Generation_Ratio < time_Check)
                {
                    variable_1 = tissue_replication_data[phases][2];
                    variable_2 = tissue_replication_data[phases][3];
                    gen_Phase = tissue_replication_data[phases][1];
                    return gen_Phase;
                    // break;
                }
            }

            //   break;
        }

        num_phases_per_tissue++;

        if (num_phases_per_tissue == num_replication_phases[(profile_ID * num_Tissues) + tissue_Check])
        {
            num_phases_per_tissue = 0;
            tissue_Check++;
        }
    }

    return gen_Phase;
}

void node_within_host::run_Generation(functions_library &functions, string &multi_Read, int &max_Cells_at_a_time, int &gpu_Limit, int *CUDA_device_IDs, int &num_Cuda_devices, int &genome_Length,
                                      string source_sequence_Data_folder, string &output_Node_location,
                                      vector<string> &tissue_Names,
                                      int *num_replication_phases, float **tissue_replication_data, int *tissue_param_profile_Stride,
                                      int terminal_tissues, int *terminal_array,
                                      int **cell_Distribution_Type, vector<pair<float, float>> &viral_distribution_per_Tissue_param,
                                      float *Reference_fitness_survivability_proof_reading,
                                      int *mutation_recombination_proof_Reading_availability,
                                      int *num_effect_Segregating_sites,
                                      float **sequence_Fitness_changes,
                                      float **sequence_Survivability_changes,
                                      float **sequence_Proof_reading_changes,
                                      int &mutation_Hotspots,
                                      float **mutation_hotspot_parameters,
                                      float **A_0_mutation,
                                      float **T_1_mutation,
                                      float **G_2_mutation,
                                      float **C_3_mutation,
                                      int &recombination_Hotspots,
                                      float **recombination_hotspot_parameters,
                                      int *tot_prob_selectivity,
                                      int *recombination_prob_Stride,
                                      int *recombination_select_Stride,
                                      float **recombination_Prob_matrix,
                                      float **recombination_Select_matrix,
                                      float *progeny_distribution_parameters_Array,
                                      mt19937 &gen)
{
    cout << "\nSimulating generation " << current_Generation << " of " << num_Generation << " for " << get_Name() << endl
         << endl;

    if (current_Generation < num_Generation)
    {
        cout << "Calculating actual particles in each tissue: \n";
        int *real_Particle_count_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
        int sum_Check = 0;
        for (int tissue = 0; tissue < num_Tissues; tissue++)
        {
            real_Particle_count_per_Tissue[tissue] = current_Viral_load_per_Tissue[tissue] - removed_by_Transfer_Indexes[tissue].size() - dead_Particle_count[tissue];
            cout << tissue_Names[tissue] << " tissue: " << real_Particle_count_per_Tissue[tissue] << endl;
            sum_Check = sum_Check + real_Particle_count_per_Tissue[tissue];
        }

        if (sum_Check > 0)
        {
            if (terminal_status(terminal_tissues, terminal_array) != 1)
            {
                cout << "\nIntiating simulation\n";
                // cout << profile_ID << endl;

                vector<vector<pair<int, int>>> indexed_Source_Folders = functions.index_sequence_Folders(source_sequence_Data_folder, num_Tissues, current_Generation, multi_Read);

                for (int tissue = 0; tissue < num_Tissues; tissue++)
                {
                    if (real_Particle_count_per_Tissue[tissue] > 0)
                    {
                        // real_Particle_count_per_Tissue[tissue] = 123;

                        cout << "\nSimulating " << real_Particle_count_per_Tissue[tissue] << " particle(s) for " << tissue_Names[tissue] << " tissue\n"
                             << endl;

                        // cout << profile_ID << endl;

                        // for (int generation = current_Generation; generation < num_Generation; generation++)
                        // {
                        //     float variable_1, variable_2;
                        //     int gen_Phase = get_generation_Phase(generation, num_replication_phases, tissue_replication_data, tissue_param_profile_Stride, tissue,
                        //                                          variable_1, variable_2);
                        //     cout << "Gen phase " << generation << ": " << gen_Phase << endl;
                        //     cout << variable_1 << "\t" << variable_2 << endl;
                        // }
                        // exit(-1);

                        cout << "Identifying indexes to remove\n";
                        set<int> check_to_Remove;
                        //// Account for dead file
                        if (dead_Particle_count[tissue] > 0)
                        {
                            cout << "\nIdentifying dead viral indexe(s)\n";
                            // indexes_of_Dead = (int *)malloc(sizeof(int) * dead_Particle_count[tissue]);

                            fstream dead_File;
                            dead_File.open(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation) + "/dead_List.txt");
                            if (dead_File.is_open())
                            {
                                string line;
                                // int index = 0;
                                while (getline(dead_File, line))
                                {
                                    check_to_Remove.insert(stoi(line));
                                    // index++;
                                }
                                dead_File.close();
                            }
                            else
                            {
                                cout << "ERROR: UNABLE TO OPEN DEAD LIST FILE: " << source_sequence_Data_folder << "/" << tissue << "/generation_" << current_Generation << "/dead_List.txt" << endl;
                                exit(-1);
                            }
                        }

                        if (removed_by_Transfer_Indexes[tissue].size() > 0)
                        {
                            cout << "Identifying transfered viral indexe(s)\n";
                            for (auto it = removed_by_Transfer_Indexes[tissue].begin(); it != removed_by_Transfer_Indexes[tissue].end(); ++it)
                            {
                                int value = *it; // Dereference the iterator to get the value
                                check_to_Remove.insert(value);
                            }
                        }

                        // int *parents_in_Tissue = (int *)malloc(sizeof(int) * real_Particle_count_per_Tissue[tissue]);
                        int **parents_in_Tissue = functions.create_INT_2D_arrays(2, real_Particle_count_per_Tissue[tissue]);

                        // test
                        // check_to_Remove.insert(0);
                        // check_to_Remove.insert(1);
                        // check_to_Remove.insert(5);
                        // check_to_Remove.insert(99);
                        // cout << profile_ID << endl;
                        // for (int generation = current_Generation; generation < num_Generation; generation++)
                        //{

                        float variable_1, variable_2;
                        int gen_Phase = get_generation_Phase(current_Generation, num_replication_phases, tissue_replication_data, tissue_param_profile_Stride, tissue,
                                                             variable_1, variable_2);

                        // cout << "Gen phase " << generation << ": " << gen_Phase << endl;
                        // cout << variable_1 << "\t" << variable_2 << endl;

                        // cout << "Gen phase " << generation << ": " << gen_Phase << endl;
                        // cout << variable_1 << "\t" << variable_2 << endl;

                        vector<int> start_Stop_cells = assign_Cells(functions, parents_in_Tissue, real_Particle_count_per_Tissue[tissue], tissue,
                                                                    cell_Distribution_Type[profile_ID][tissue], viral_distribution_per_Tissue_param[tissue].first, viral_distribution_per_Tissue_param[tissue].second,
                                                                    check_to_Remove,
                                                                    gen_Phase, variable_1, variable_2,
                                                                    gen);

                        check_to_Remove.clear();
                        removed_by_Transfer_Indexes[tissue].clear();
                        dead_Particle_count[tissue] = 0;

                        cout << "Total number of cell(s) infected: " << start_Stop_cells.size() - 1 << endl;

                        if (start_Stop_cells.size() - 1 > 0)
                        {
                            for (int i = 0; i < start_Stop_cells.size() - 1; i++)
                            {
                                // cout << start_Stop_cells[i] << " : \t" << start_Stop_cells[i + 1] << endl;
                                for (int particle = start_Stop_cells[i]; particle < start_Stop_cells[i + 1]; particle++)
                                {
                                    // cout << parents_in_Tissue[0][particle] << " :\t" << parents_in_Tissue[1][particle] << endl;
                                    cout << parents_in_Tissue[1][particle] << "_" << parents_in_Tissue[0][particle] << ", ";
                                }
                                cout << endl;
                            }
                        }
                        // }
                        // exit(-1);
                        // vector<int> start_Stop_cells;
                        if (start_Stop_cells.size() - 1 > 0)
                        {
                            for (int i = 0; i < start_Stop_cells.size() - 1; i++)
                            {
                                // cout << start_Stop_cells[i] << " : \t" << start_Stop_cells[i + 1] << endl;
                                for (int particle = start_Stop_cells[i]; particle < start_Stop_cells[i + 1]; particle++)
                                {
                                    // cout << parents_in_Tissue[0][particle] << " :\t" << parents_in_Tissue[1][particle] << endl;
                                    cout << parents_in_Tissue[1][particle] << "_" << parents_in_Tissue[0][particle] << ", ";
                                }
                                cout << endl;
                            }

                            // exit(-1);

                            vector<pair<int, int>> cells_Rounds_start_stop;

                            int full_Rounds = (start_Stop_cells.size() - 1) / max_Cells_at_a_time;
                            int partial_Rounds = (start_Stop_cells.size() - 1) % max_Cells_at_a_time;

                            for (int full = 0; full < full_Rounds; full++)
                            {
                                int start = full * max_Cells_at_a_time;
                                int stop = start + max_Cells_at_a_time;
                                cells_Rounds_start_stop.push_back(make_pair(start, stop));
                            }

                            if (partial_Rounds != 0)
                            {
                                int start = (start_Stop_cells.size() - 1) - partial_Rounds;
                                cells_Rounds_start_stop.push_back(make_pair(start, (start_Stop_cells.size() - 1)));
                            }

                            int sequence_Count = 0;
                            int index_Last_Written = 0;

                            string intermediary_Tissue_folder = source_sequence_Data_folder + "/" + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation + 1);
                            string dead_List = intermediary_Tissue_folder + "/dead_List.txt";

                            if (filesystem::exists(intermediary_Tissue_folder))
                            {
                                functions.config_Folder(intermediary_Tissue_folder, to_string(current_Generation + 1) + " generation Tissue " + tissue_Names[tissue] + " sequences");
                            }

                            if (!filesystem::exists(dead_List))
                            {
                                functions.create_File(dead_List);
                            }

                            string sequence_Profiles = output_Node_location + "/" + get_Name() + "/sequence_Profiles.csv";
                            string sequence_parent_Progeny_relationships = output_Node_location + "/" + get_Name() + "/sequence_parent_Progeny_relationships.csv";

                            if (!filesystem::exists(output_Node_location + "/" + get_Name()))
                            {
                                functions.config_Folder(output_Node_location + "/" + get_Name(), get_Name() + " node");
                                functions.create_File(sequence_Profiles, "Sequence_ID\tHost\tTissue");
                                functions.create_File(sequence_parent_Progeny_relationships, "Source\tTarget\tType");
                            }
                            // size_t arraySize = sizeof(parents_in_Tissue) / sizeof(parents_in_Tissue[0]);

                            for (int cell_Round = 0; cell_Round < cells_Rounds_start_stop.size(); cell_Round++)
                            {
                                int num_of_Cells = cells_Rounds_start_stop[cell_Round].second - cells_Rounds_start_stop[cell_Round].first;
                                cout << "\nProcessing round " << cell_Round + 1 << " of " << cells_Rounds_start_stop.size() << ": " << num_of_Cells << " cell(s)" << endl;

                                int seqeunces_to_Process = start_Stop_cells[cells_Rounds_start_stop[cell_Round].second] - start_Stop_cells[cells_Rounds_start_stop[cell_Round].first];
                                cout << "Processing " << seqeunces_to_Process << " sequence(s) in total\n";

                                simulate_Cell_replication(functions, multi_Read, gpu_Limit, CUDA_device_IDs, num_Cuda_devices, source_sequence_Data_folder, indexed_Source_Folders[tissue],
                                                          genome_Length,
                                                          tissue, parents_in_Tissue,
                                                          start_Stop_cells, cells_Rounds_start_stop[cell_Round].first, cells_Rounds_start_stop[cell_Round].second, num_of_Cells,
                                                          seqeunces_to_Process,
                                                          sequence_Count,
                                                          Reference_fitness_survivability_proof_reading,
                                                          mutation_recombination_proof_Reading_availability,
                                                          num_effect_Segregating_sites,
                                                          sequence_Fitness_changes,
                                                          sequence_Survivability_changes,
                                                          sequence_Proof_reading_changes,
                                                          mutation_Hotspots,
                                                          mutation_hotspot_parameters,
                                                          A_0_mutation,
                                                          T_1_mutation,
                                                          G_2_mutation,
                                                          C_3_mutation,
                                                          recombination_Hotspots,
                                                          recombination_hotspot_parameters,
                                                          tot_prob_selectivity,
                                                          recombination_prob_Stride,
                                                          recombination_select_Stride,
                                                          recombination_Prob_matrix,
                                                          recombination_Select_matrix,
                                                          progeny_distribution_parameters_Array,
                                                          gen);
                            }
                            // free(parents_in_Tissue);
                            exit(-1);
                        }
                        functions.clear_Array_int_CPU(parents_in_Tissue, 2);
                    }
                    // cout << "Cell Limit: " << cell_Limit[tissue] << endl;

                    // cout << "Distribution type: " << cell_Distribution_Type[profile_ID][tissue] << endl;
                    // cout << viral_distribution_per_Tissue_param[tissue].first << "\t" << viral_distribution_per_Tissue_param[tissue].second << endl;
                }
            }
        }
        else
        {
            set_Removed();
        }
    }
    else
    {
        set_Removed();
    }

    // get each tissues generational phase
}

void node_within_host::write_Full_Sequences_Progeny(functions_library &functions,
                                                    int &CPU_cores,
                                                    string &tissue_Name,
                                                    int &num_Progeny_being_Processed,
                                                    int &genome_Length, int &recombination_Hotspots,
                                                    int &sequence_Count,
                                                    int &index_Last_Written,
                                                    int **parent_IDs,
                                                    int **progeny_Sequences, int *Dead_or_Alive, int **progeny_Configuration_Filled,
                                                    string intermediary_Tissue_folder, string dead_List, string sequence_Profiles, string sequence_parent_Progeny_relationships,
                                                    int &max_sequences_per_File, int &last_index_Seq_Written)
{
    cout << "Writing sequence profiles and parent progeny_Relationships\n";

    fstream sequence_Profile_File;
    sequence_Profile_File.open(sequence_Profiles, ios::app);
    fstream parent_Progeny_Relationships_File;
    parent_Progeny_Relationships_File.open(sequence_parent_Progeny_relationships, ios::app);

    string viral_prefix_Progeny = get_Name() + "_" + tissue_Name + "_" + to_string(current_Generation + 1) + "_";
    string viral_prefix_Parent = get_Name() + "_" + tissue_Name + "_" + to_string(current_Generation) + "_";

    // Cells of parent and progeny

    if (sequence_Profile_File.is_open() && parent_Progeny_Relationships_File.is_open())
    {
        for (int progeny = 0; progeny < num_Progeny_being_Processed; progeny++)
        {
            converted_Sequences.push_back("");
            sequence_Profile_File << viral_prefix_Progeny << to_string(sequence_Count)
                                  << "\t" << get_Name()
                                  << "\t" << tissue_Name << endl;

            parent_Progeny_Relationships_File << viral_prefix_Parent << to_string(parent_IDs[0][progeny_Configuration_Filled[progeny][0]])
                                              << "\t" << viral_prefix_Progeny << to_string(sequence_Count)
                                              << "\tprimary_parent" << endl;

            for (int hotspot_parent = 1; hotspot_parent < (1 + recombination_Hotspots); hotspot_parent++)
            {
                if (progeny_Configuration_Filled[progeny][hotspot_parent] != -1)
                {
                    parent_Progeny_Relationships_File << viral_prefix_Parent << to_string(parent_IDs[0][progeny_Configuration_Filled[progeny][hotspot_parent]])
                                                      << "\t" << viral_prefix_Progeny << to_string(sequence_Count)
                                                      << "\trecombination_parent_" << to_string(hotspot_parent) << endl;
                }
            }

            sequence_Count++;
        }

        sequence_Profile_File.close();
        parent_Progeny_Relationships_File.close();
    }
    else
    {
        cout << "ERROR IN OPENING BOTH OR EITHER ONE OF THE FOLLOWING PROGENY SEQUENCE INORMATION FILES.\n";
        cout << "SEQUENCE PROFILE FILE:  " << sequence_Profiles << endl;
        cout << "SEQUENCE PROGENY PARENT RELATIONSHIPS FILE:  " << sequence_parent_Progeny_relationships << endl;
        exit(-1);
    }

    int num_per_Core = num_Progeny_being_Processed / CPU_cores;
    int remainder = num_Progeny_being_Processed % CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
    {
        int start_Node = core_ID * num_per_Core;
        int stop_Node = start_Node + num_per_Core;

        threads_vec.push_back(thread{&node_within_host::thread_Sequence_to_String, this, start_Node, stop_Node, progeny_Sequences, genome_Length});
    }

    if (remainder != 0)
    {
        int start_Node = num_Progeny_being_Processed - remainder;
        int stop_Node = num_Progeny_being_Processed;

        threads_vec.push_back(thread{&node_within_host::thread_Sequence_to_String, this, start_Node, stop_Node, progeny_Sequences, genome_Length});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    for (int progeny = 0; progeny < num_Progeny_being_Processed; progeny++)
    {
        to_write_Sequence_Store.push_back(make_pair(to_string(Dead_or_Alive[progeny]), converted_Sequences[progeny]));
    }

    converted_Sequences.clear();
    //  int **progeny_Sequences, int *Dead_or_Alive, int **progeny_Configuration_Filled,
    functions.clear_Array_int_CPU(progeny_Sequences, num_Progeny_being_Processed);
    free(Dead_or_Alive);
    functions.clear_Array_int_CPU(progeny_Configuration_Filled, num_Progeny_being_Processed);

    if (to_write_Sequence_Store.size() >= max_sequences_per_File)
    {
        int full_Write_Count = to_write_Sequence_Store.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {
            string fasta_file_Location = intermediary_Tissue_folder + "/" + to_string(last_index_Seq_Written) + "_" + to_string(last_index_Seq_Written + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);

            fstream dead_List_File;
            dead_List_File.open(dead_List, ios::app);

            if (fasta_File.is_open() && dead_List_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_index_Seq_Written << "_";
                    if (to_write_Sequence_Store[write_Seq].first == 0)
                    {
                        fasta_File << "D";
                        dead_List_File << last_index_Seq_Written << endl;
                    }
                    else
                    {
                        fasta_File << "A";
                    }
                    fasta_File << endl;
                    fasta_File << c[write_Seq].second << endl;

                    last_index_Seq_Written++;
                }
                fasta_File.close();
                dead_List_File.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                cout << "OR DEAD LIST FILE COULD NOT BE OPENED: " << dead_List << endl;
                exit(-1);
            }
        }

        vector<pair<string, string>> to_write_Sequence_Store_temp;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < sequence_Write_Store_All.size(); fill++)
        {
            to_write_Sequence_Store_temp.push_back(make_pair(to_write_Sequence_Store[fill].first, to_write_Sequence_Store[fill].second));
        }

        to_write_Sequence_Store.clear();
        to_write_Sequence_Store = to_write_Sequence_Store_temp;
    }

    // CREATE A FUNCTION FOR PARTIAL
}

void node_within_host::thread_Sequence_to_String(int start, int stop, int **progeny_Sequences, int genome_Length)
{
    vector<string> converted_Sequences_Store;

    for (int progeny = start; progeny < stop; progeny++)
    {
        string sequence = "";
        for (int base = 0; base < genome_Length; base++)
        {
            sequence.append(to_string(progeny_Sequences[progeny][base]));
        }
        converted_Sequences_Store.push_back(sequence);
    }

    unique_lock<shared_mutex> ul(g_mutex);
    int index = 0;
    for (int progeny = start; progeny < stop; progeny++)
    {
        converted_Sequences[progeny] = converted_Sequences_Store[index];
        index++;
    }
}

void node_within_host::simulate_Cell_replication(functions_library &functions, string &multi_Read, int &gpu_Limit, int *CUDA_device_IDs, int &num_Cuda_devices, string &source_sequence_Data_folder, vector<pair<int, int>> &indexed_Tissue_Folder,
                                                 int &genome_Length,
                                                 int &tissue, int **parents_in_Tissue,
                                                 vector<int> &start_Stop_cells, int &start_Cell, int &stop_Cell, int &num_Cells,
                                                 int &Total_seqeunces_to_Process,
                                                 int &sequence_Count,
                                                 float *Reference_fitness_survivability_proof_reading,
                                                 int *mutation_recombination_proof_Reading_availability,
                                                 int *num_effect_Segregating_sites,
                                                 float **sequence_Fitness_changes,
                                                 float **sequence_Survivability_changes,
                                                 float **sequence_Proof_reading_changes,
                                                 int &mutation_Hotspots,
                                                 float **mutation_hotspot_parameters,
                                                 float **A_0_mutation,
                                                 float **T_1_mutation,
                                                 float **G_2_mutation,
                                                 float **C_3_mutation,
                                                 int &recombination_Hotspots,
                                                 float **recombination_hotspot_parameters,
                                                 int *tot_prob_selectivity,
                                                 int *recombination_prob_Stride,
                                                 int *recombination_select_Stride,
                                                 float **recombination_Prob_matrix,
                                                 float **recombination_Select_matrix,
                                                 float *progeny_distribution_parameters_Array,
                                                 mt19937 &gen)
{
    // gpu_Limit = 5;

    int full_Rounds = Total_seqeunces_to_Process / gpu_Limit;
    int partial_Rounds = Total_seqeunces_to_Process % gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * gpu_Limit;
        int stop = start + gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = Total_seqeunces_to_Process - partial_Rounds;
        start_stops.push_back(make_pair(start, Total_seqeunces_to_Process));
    }

    cout << "Retrieving parent sequences and configuring their profiles\n";

    // ! clear 2d array
    int **parent_Sequences = functions.create_INT_2D_arrays(Total_seqeunces_to_Process, genome_Length);
    // ! clear 2d array
    float **sequence_Configuration_standard;
    // int columns = 3 + (2 * recombination_Hotspots);
    sequence_Configuration_standard = functions.create_FLOAT_2D_arrays(Total_seqeunces_to_Process, 2 + (2 * recombination_Hotspots));

    // float **sequence_Configuration_recombination;
    // if (mutation_recombination_proof_Reading_availability[1] == 1)
    // {
    //     sequence_Configuration_recombination = functions.create_FLOAT_2D_arrays(Total_seqeunces_to_Process, 2 * recombination_Hotspots);
    // }

    // ! clear 2d array
    // int *parent_IDs = (int *)malloc(sizeof(int) * Total_seqeunces_to_Process);
    int **parent_IDs = functions.create_INT_2D_arrays(2, Total_seqeunces_to_Process);

    // cout << progeny_distribution_parameters_Array[0] << endl;
    // cout << progeny_distribution_parameters_Array[1] << endl;
    // cout << progeny_distribution_parameters_Array[2] << endl;
    // exit(-1);

    int cell_ID = 0;
    int *cell_Index = (int *)malloc(sizeof(int) * (num_Cells + 1));
    cell_Index[0] = 0;

    int check_Cell = -1;

    for (int round = 0; round < start_stops.size(); round++)
    {
        cout << "\nParent sequence processing round " << round + 1 << " of " << start_stops.size() << endl;

        vector<int> sequence_List;
        for (int parent = start_stops[round].first; parent < start_stops[round].second; parent++)
        {
            sequence_List.push_back(parents_in_Tissue[0][start_Stop_cells[start_Cell] + parent]);
            parent_IDs[0][parent] = parents_in_Tissue[0][start_Stop_cells[start_Cell] + parent];
            int current_Cell = parents_in_Tissue[1][start_Stop_cells[start_Cell] + parent];
            // parent_IDs[1][parent] = parents_in_Tissue[1][start_Stop_cells[start_Cell] + parent];

            if (parent != 0)
            {
                if (check_Cell != current_Cell)
                {
                    cell_ID++;
                    cell_Index[cell_ID] = parent;
                    check_Cell = current_Cell;
                }
                parent_IDs[1][parent] = cell_ID;
            }
            else
            {
                parent_IDs[1][parent] = cell_ID;
                check_Cell = current_Cell;
            }

            if (parent + 1 == start_stops[round].second)
            {
                cell_ID++;
                cell_Index[cell_ID] = parent + 1;
            }
        }

        vector<string> collected_Sequences = functions.find_Sequences_Master(source_sequence_Data_folder, sequence_List, tissue, indexed_Tissue_Folder, current_Generation);

        if (collected_Sequences.size() == sequence_List.size())
        {
            sequence_List.clear();
            // for (int test = 0; test < collected_Sequences.size(); test++)
            // {
            //     cout << collected_Sequences[test] << endl;
            // }
            // cout << endl;

            process_Sequences_get_Configuration(functions,
                                                collected_Sequences, CUDA_device_IDs, num_Cuda_devices, genome_Length,
                                                Reference_fitness_survivability_proof_reading,
                                                mutation_recombination_proof_Reading_availability,
                                                num_effect_Segregating_sites,
                                                sequence_Fitness_changes,
                                                sequence_Proof_reading_changes,
                                                recombination_Hotspots,
                                                recombination_hotspot_parameters,
                                                tot_prob_selectivity,
                                                recombination_prob_Stride,
                                                recombination_select_Stride,
                                                recombination_Prob_matrix,
                                                recombination_Select_matrix,
                                                progeny_distribution_parameters_Array,
                                                parent_Sequences,
                                                sequence_Configuration_standard,
                                                start_stops[round].first);
        }
        else
        {
            cout << "ERROR: WAS UNABLE TO FIND ALL REQUIRED SEQUENCES.\n";
            exit(-1);
        }

        cout << endl;
    }

    // for (int i = 0; i < num_Cells; i++)
    // {
    //     for (int parent = cell_Index[i]; parent < cell_Index[i + 1]; parent++)
    //     {
    //         cout << parent_IDs[1][parent] << "_" << parent_IDs[0][parent] << ",";
    //     }
    //     cout << endl;
    // }

    // exit(-1);

    cout << endl;

    for (int row = 0; row < Total_seqeunces_to_Process; row++)
    {
        for (int col = 0; col < genome_Length; col++)
        {
            cout << parent_Sequences[row][col];
        }
        cout << endl;
    }

    cout << endl;

    for (int row = 0; row < Total_seqeunces_to_Process; row++)
    {
        for (int col = 0; col < (2 + (2 * recombination_Hotspots)); col++)
        {
            cout << sequence_Configuration_standard[row][col] << "\t";
        }
        cout << endl;
    }

    cout << "\nAll parent sequences configured\n";
    int total_Progeny = 0;
    // ! clear 1d array
    // float *totals_Progeny_Selectivity = (float *)malloc(sizeof(float) * (1 + recombination_Hotspots));
    float **totals_Progeny_Selectivity = functions.create_Fill_2D_array_FLOAT(num_Cells, recombination_Hotspots, 0);
    ////  clear 1d array
    int *progeny_Stride = (int *)malloc(sizeof(int) * (Total_seqeunces_to_Process + 1));
    progeny_Stride[0] = 0;

    // for (int fill = 0; fill < (1 + recombination_Hotspots); fill++)
    // {
    //     totals_Progeny_Selectivity[fill] = 0;
    // }

    cout << "\nDetermining total progeny and configuring recombination hotspots\n";

    for (int cell = 0; cell < num_Cells; cell++)
    {
        for (int parent = cell_Index[cell]; parent < cell_Index[cell + 1]; parent++)
        {
            // cout << parent_IDs[1][parent] << "_" << parent_IDs[0][parent] << ",";
            progeny_Stride[parent + 1] = progeny_Stride[parent] + sequence_Configuration_standard[parent][0];

            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                totals_Progeny_Selectivity[cell][hotspot] = totals_Progeny_Selectivity[cell][hotspot] + sequence_Configuration_standard[parent][(hotspot * 2) + 3];
            }
        }
        cout << endl;
    }

    // for (int row = 0; row < Total_seqeunces_to_Process; row++)
    // {
    //     // totals_Progeny_Selectivity[0] = totals_Progeny_Selectivity[0] + sequence_Configuration_standard[row][0];
    //     progeny_Stride[row + 1] = progeny_Stride[row] + sequence_Configuration_standard[row][0];

    //     for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
    //     {
    //         totals_Progeny_Selectivity[hotspot + 1] = totals_Progeny_Selectivity[hotspot + 1] + sequence_Configuration_standard[row][(hotspot * 2) + 3];
    //     }
    // }

    total_Progeny = progeny_Stride[Total_seqeunces_to_Process];
    cout << "Total progeny to be simulated: " << total_Progeny << endl;

    cout << endl;
    for (int i = 0; i < num_Cells; i++)
    {
        for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
        {
            cout << totals_Progeny_Selectivity[i][hotspot] << "\t";
        }
        cout << endl;
    }

    cout << "\nIntiating Progeny configurations\n";
    cout << "Intializing GPU memory structures\n";

    cudaSetDevice(CUDA_device_IDs[0]);

    // !clear array
    // int **cuda_progeny_Configuration = functions.create_CUDA_2D_int(total_Progeny, 1 + recombination_Hotspots);
    int **progeny_Configuration = functions.create_INT_2D_arrays(total_Progeny, 1 + recombination_Hotspots);

    //// clear array
    int *cuda_progeny_Stride;
    cudaMallocManaged(&cuda_progeny_Stride, (Total_seqeunces_to_Process + 1) * sizeof(int));
    cudaMemcpy(cuda_progeny_Stride, progeny_Stride, (Total_seqeunces_to_Process + 1) * sizeof(int), cudaMemcpyHostToDevice);
    free(progeny_Stride);

    float **cuda_sequence_Configuration_standard = functions.float_2D_Array_load_to_CUDA(sequence_Configuration_standard, Total_seqeunces_to_Process, 2 + (2 * recombination_Hotspots));

    for (int round = 0; round < start_stops.size(); round++)
    {
        cout << "\nParent sequence processing round " << round + 1 << " of " << start_stops.size() << endl;

        progeny_Configurator(functions,
                             cuda_sequence_Configuration_standard, recombination_Hotspots,
                             start_stops[round].first, start_stops[round].second - start_stops[round].first,
                             progeny_Configuration, cuda_progeny_Stride, cuda_progeny_Stride[start_stops[round].second] - cuda_progeny_Stride[start_stops[round].first], cuda_progeny_Stride[start_stops[round].first]);
    }

    cout << "\nCopying test\n";
    // int **progeny_Configuration = functions.load_to_Host(cuda_progeny_Configuration, total_Progeny, 1 + recombination_Hotspots);

    for (int row = 0; row < total_Progeny; row++)
    {
        for (int col = 0; col < (1 + recombination_Hotspots); col++)
        {
            cout << progeny_Configuration[row][col] << "\t";
        }
        cout << endl;
    }

    // functions.clear_Array_INT(cuda_progeny_Configuration, total_Progeny);
    cudaFree(cuda_progeny_Stride);

    // create and save sequence
    start_stops.clear();

    full_Rounds = total_Progeny / gpu_Limit;
    partial_Rounds = total_Progeny % gpu_Limit;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * gpu_Limit;
        int stop = start + gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Progeny - partial_Rounds;
        start_stops.push_back(make_pair(start, total_Progeny));
    }

    for (int round = 0; round < start_stops.size(); round++)
    {
        progeny_Completion(functions,
                           CUDA_device_IDs, num_Cuda_devices,
                           genome_Length, Reference_fitness_survivability_proof_reading, mutation_recombination_proof_Reading_availability,
                           num_effect_Segregating_sites,
                           sequence_Survivability_changes,
                           recombination_Hotspots,
                           recombination_hotspot_parameters,
                           tot_prob_selectivity,
                           mutation_Hotspots,
                           A_0_mutation,
                           T_1_mutation,
                           G_2_mutation,
                           C_3_mutation,
                           mutation_hotspot_parameters,
                           parent_Sequences, Total_seqeunces_to_Process, sequence_Configuration_standard, parent_IDs, num_Cells, cell_Index, start_Cell,
                           progeny_Configuration, start_stops[round].second - start_stops[round].first,
                           totals_Progeny_Selectivity,
                           start_stops[round].first, start_stops[round].second, sequence_Count, source_sequence_Data_folder);
    }
}

__global__ void cuda_Progeny_Complete_Configuration(int genome_Length,
                                                    float *cuda_Reference_fitness_survivability_proof_reading,
                                                    int *cuda_num_effect_Segregating_sites,
                                                    float **cuda_sequence_Survivability_changes,
                                                    int recombination_Hotspots,
                                                    float **cuda_recombination_hotspot_parameters,
                                                    int *cuda_tot_prob_selectivity,
                                                    int mutation_Hotspots,
                                                    float **cuda_A_0_mutation,
                                                    float **cuda_T_1_mutation,
                                                    float **cuda_G_2_mutation,
                                                    float **cuda_C_3_mutation,
                                                    float **cuda_mutation_hotspot_parameters,
                                                    int **cuda_parent_Sequences, int **cuda_parent_IDs,
                                                    float **cuda_sequence_Configuration_standard,
                                                    int *cuda_cell_Index, int num_Cells,
                                                    float **cuda_totals_Progeny_Selectivity,
                                                    int **cuda_progeny_Configuration,
                                                    int **cuda_progeny_Sequences,
                                                    int *cuda_Dead_or_Alive,
                                                    int per_gpu_Progeny)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < per_gpu_Progeny)
    {
        for (int base = 0; base < genome_Length; base++)
        {
            cuda_progeny_Sequences[tid][base] = cuda_parent_Sequences[cuda_progeny_Configuration[tid][0]][base];
        }

        curandState localState;
        curand_init(clock64(), tid, 0, &localState);

        if (recombination_Hotspots > 0)
        {
            int get_Cell = cuda_parent_IDs[1][cuda_progeny_Configuration[tid][0]];

            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                if (cuda_progeny_Configuration[tid][hotspot + 1] != -1)
                {
                    // TEST BLOCK 1
                    float rand_num = curand_uniform(&localState);
                    float cumulative_prob = 0.0f;

                    int recomb_parent = -1;

                    for (int check = cuda_cell_Index[get_Cell]; check < cuda_cell_Index[get_Cell + 1]; check++)
                    {
                        cumulative_prob += (cuda_sequence_Configuration_standard[check][(hotspot * 2) + 3] / cuda_totals_Progeny_Selectivity[get_Cell][hotspot]);
                        if (rand_num < cumulative_prob)
                        {
                            recomb_parent = check;
                            break;
                        }
                    }

                    cuda_progeny_Configuration[tid][hotspot + 1] = recomb_parent;

                    // TEST BLOCK 2
                    if (recomb_parent != cuda_progeny_Configuration[tid][0])
                    {
                        for (int base = (cuda_recombination_hotspot_parameters[hotspot][0] - 1); base < cuda_recombination_hotspot_parameters[hotspot][1]; base++)
                        {
                            cuda_progeny_Sequences[tid][base] = cuda_parent_Sequences[recomb_parent][base];
                        }
                    }
                }
            }
        }

        // TEST BLOCK 3
        if (mutation_Hotspots > 0)
        {
            for (int hotspot = 0; hotspot < mutation_Hotspots; hotspot++)
            {
                int num_Mutations = -1;

                if (cuda_mutation_hotspot_parameters[hotspot][2] == 0)
                {
                    // Poisson
                    num_Mutations = curand_poisson(&localState, cuda_mutation_hotspot_parameters[hotspot][3]);
                }
                else if (cuda_mutation_hotspot_parameters[hotspot][2] == 1)
                {
                    // neg binomial
                    int failures = 0;
                    int successes = 0;

                    while (successes < cuda_mutation_hotspot_parameters[hotspot][3])
                    {
                        float rand_num = curand_uniform(&localState);
                        if (rand_num < cuda_mutation_hotspot_parameters[hotspot][4])
                        {
                            successes++;
                        }
                        else
                        {
                            failures++;
                        }
                    }

                    num_Mutations = failures;
                }
                else
                {
                    // fixed or binomial distribution
                    int count = 0;

                    int bases_in_Region = cuda_mutation_hotspot_parameters[hotspot][1] - (cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                    for (int trial = 0; trial < bases_in_Region; trial++)
                    {
                        if (curand_uniform(&localState) < cuda_mutation_hotspot_parameters[hotspot][3])
                        {
                            count++;
                        }
                    }

                    num_Mutations = count;
                }

                if (num_Mutations > 0)
                {
                    if (cuda_sequence_Configuration_standard[cuda_progeny_Configuration[tid][0]][1] != -1)
                    {
                        int count = 0;

                        int bases_in_Region = cuda_mutation_hotspot_parameters[hotspot][1] - (cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                        for (int trial = 0; trial < num_Mutations; trial++)
                        {
                            if (curand_uniform(&localState) < cuda_sequence_Configuration_standard[cuda_progeny_Configuration[tid][0]][1])
                            {
                                count++;
                            }
                        }
                        num_Mutations = num_Mutations - count;
                    }

                    if (num_Mutations > 0)
                    {
                        for (int mutation = 0; mutation < num_Mutations; mutation++)
                        {
                            int position = (int)(curand_uniform(&localState) * ((cuda_mutation_hotspot_parameters[hotspot][1] - 1) - (cuda_mutation_hotspot_parameters[hotspot][0] - 1) + 1)) + (cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                            float rand_num = curand_uniform(&localState);
                            float cumulative_prob = 0.0f;

                            int original_BASE = cuda_progeny_Sequences[tid][position];
                            int new_Base = 0;

                            for (int base = 0; base < 4; base++)
                            {
                                cumulative_prob += (original_BASE == 0)   ? cuda_A_0_mutation[hotspot][base]
                                                   : (original_BASE == 1) ? cuda_T_1_mutation[hotspot][base]
                                                   : (original_BASE == 2) ? cuda_G_2_mutation[hotspot][base]
                                                   : (original_BASE == 3) ? cuda_C_3_mutation[hotspot][base]
                                                                          : 0.0f;

                                if (rand_num < cumulative_prob)
                                {
                                    new_Base = base;
                                    break;
                                }
                            }

                            cuda_progeny_Sequences[tid][position] = new_Base;
                        }
                    }
                }
            }
        }

        // TEST BLOCK 4
        // Determine survivability
        float survivability = cuda_Reference_fitness_survivability_proof_reading[1];

        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[1]; pos++)
        {
            if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 0)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][1];
            }
            else if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 1)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][2];
            }
            else if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 2)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][3];
            }
            else if (cuda_progeny_Sequences[tid][(int)cuda_sequence_Survivability_changes[pos][0] - 1] == 3)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][4];
            }
        }

        if (survivability >= 1)
        {
            cuda_Dead_or_Alive[tid] = 1;
        }
        else if (survivability <= 0)
        {
            cuda_Dead_or_Alive[tid] = 0;
        }
        else
        {
            float survivability_Check = curand_uniform(&localState);
            cuda_Dead_or_Alive[tid] = (survivability_Check < survivability) ? 1 : 0;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void node_within_host::progeny_Completion(functions_library &functions,
                                          int *CUDA_device_IDs, int &num_Cuda_devices,
                                          int genome_Length, float *Reference_fitness_survivability_proof_reading, int *mutation_recombination_proof_Reading_availability,
                                          int *num_effect_Segregating_sites,
                                          float **sequence_Survivability_changes,
                                          int recombination_Hotspots,
                                          float **recombination_hotspot_parameters,
                                          int *tot_prob_selectivity,
                                          int mutation_Hotspots,
                                          float **A_0_mutation,
                                          float **T_1_mutation,
                                          float **G_2_mutation,
                                          float **C_3_mutation,
                                          float **mutation_hotspot_parameters,
                                          int **parent_Sequences, int num_Parent_sequence, float **sequence_Configuration_standard, int **parent_IDs, int num_Cells, int *cell_Index, int &start_Cell,
                                          int **progeny_Configuration, int num_Progeny_being_Processed,
                                          float **totals_Progeny_Selectivity,
                                          int start_Progeny, int stop_Progeny, int &progeny_Count, string write_Progeny_Folder)
{
    int num_of_Sequences_current = num_Progeny_being_Processed;
    cout << "\nConfiguring multi gpu distribution of " << num_of_Sequences_current << " sequence(s)\n";

    int standard_num_per_GPU = num_of_Sequences_current / num_Cuda_devices;
    int remainder = num_of_Sequences_current % num_Cuda_devices;

    vector<pair<int, int>> start_stop_Per_GPU;

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        int start = gpu * standard_num_per_GPU;
        int stop = start + standard_num_per_GPU;

        start_stop_Per_GPU.push_back(make_pair(start, stop));
    }

    start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

    cudaStream_t streams[num_Cuda_devices];
    cudaDeviceProp deviceProp;

    float *cuda_Reference_fitness_survivability_proof_reading[num_Cuda_devices];

    int *cuda_num_effect_Segregating_sites[num_Cuda_devices];
    float **cuda_sequence_Survivability_changes[num_Cuda_devices];

    float **cuda_recombination_hotspot_parameters[num_Cuda_devices];
    int *cuda_tot_prob_selectivity[num_Cuda_devices];

    float **cuda_A_0_mutation[num_Cuda_devices];
    float **cuda_T_1_mutation[num_Cuda_devices];
    float **cuda_G_2_mutation[num_Cuda_devices];
    float **cuda_C_3_mutation[num_Cuda_devices];

    float **cuda_mutation_hotspot_parameters[num_Cuda_devices];

    int **cuda_parent_Sequences[num_Cuda_devices];
    float **cuda_sequence_Configuration_standard[num_Cuda_devices];
    int **cuda_parent_IDs[num_Cuda_devices];

    int *cuda_cell_Index[num_Cuda_devices];

    float **cuda_totals_Progeny_Selectivity[num_Cuda_devices];

    int **cuda_progeny_Configuration[num_Cuda_devices];

    int **cuda_progeny_Sequences[num_Cuda_devices];
    int *cuda_Dead_or_Alive[num_Cuda_devices];

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaGetDeviceProperties(&deviceProp, gpu);
        cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

        cudaMallocManaged(&cuda_progeny_Sequences[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_progeny_Sequences[gpu][row], genome_Length * sizeof(int));
        }

        cudaMallocManaged(&cuda_Dead_or_Alive[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int));

        cudaMallocManaged(&cuda_Reference_fitness_survivability_proof_reading[gpu], 3 * sizeof(float));
        cudaMemcpy(cuda_Reference_fitness_survivability_proof_reading[gpu], Reference_fitness_survivability_proof_reading, 3 * sizeof(float), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_num_effect_Segregating_sites[gpu], 3 * sizeof(int));
        cudaMemcpy(cuda_num_effect_Segregating_sites[gpu], num_effect_Segregating_sites, 3 * sizeof(int), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_sequence_Survivability_changes[gpu], num_effect_Segregating_sites[1] * sizeof(float *));
        for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        {
            cudaMalloc((void **)&cuda_sequence_Survivability_changes[gpu][row], 5 * sizeof(float));
            cudaMemcpy(cuda_sequence_Survivability_changes[gpu][row], sequence_Survivability_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        }

        cudaMallocManaged(&cuda_tot_prob_selectivity[gpu], 2 * sizeof(int));

        cudaMallocManaged(&cuda_recombination_hotspot_parameters[gpu], recombination_Hotspots * sizeof(float *));
        if (recombination_Hotspots > 0)
        {
            cudaMemcpy(cuda_tot_prob_selectivity[gpu], tot_prob_selectivity, 2 * sizeof(int), cudaMemcpyHostToDevice);

            for (int row = 0; row < recombination_Hotspots; row++)
            {
                cudaMalloc((void **)&cuda_recombination_hotspot_parameters[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_recombination_hotspot_parameters[gpu][row], recombination_hotspot_parameters[row], 4 * sizeof(float), cudaMemcpyHostToDevice);
            }
        }

        cudaMallocManaged(&cuda_A_0_mutation[gpu], mutation_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_T_1_mutation[gpu], mutation_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_G_2_mutation[gpu], mutation_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_C_3_mutation[gpu], mutation_Hotspots * sizeof(float *));

        cudaMallocManaged(&cuda_mutation_hotspot_parameters[gpu], mutation_Hotspots * sizeof(float *));

        if (mutation_Hotspots > 0)
        {
            for (int row = 0; row < mutation_Hotspots; row++)
            {
                cudaMalloc((void **)&cuda_A_0_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_A_0_mutation[gpu][row], A_0_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_T_1_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_T_1_mutation[gpu][row], T_1_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_G_2_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_G_2_mutation[gpu][row], G_2_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_C_3_mutation[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_C_3_mutation[gpu][row], C_3_mutation[row], 4 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMalloc((void **)&cuda_mutation_hotspot_parameters[gpu][row], 5 * sizeof(float));
                cudaMemcpy(cuda_mutation_hotspot_parameters[gpu][row], mutation_hotspot_parameters[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
            }
        }

        cudaMallocManaged(&cuda_parent_Sequences[gpu], num_Parent_sequence * sizeof(int *));
        cudaMallocManaged(&cuda_sequence_Configuration_standard[gpu], num_Parent_sequence * sizeof(float *));
        for (int row = 0; row < num_Parent_sequence; row++)
        {
            cudaMalloc((void **)&cuda_parent_Sequences[gpu][row], genome_Length * sizeof(int));
            cudaMemcpy(cuda_parent_Sequences[gpu][row], parent_Sequences[row], genome_Length * sizeof(int), cudaMemcpyHostToDevice);

            cudaMalloc((void **)&cuda_sequence_Configuration_standard[gpu][row], (2 + (2 * recombination_Hotspots)) * sizeof(float));
            cudaMemcpy(cuda_sequence_Configuration_standard[gpu][row], sequence_Configuration_standard[row], (2 + (2 * recombination_Hotspots)) * sizeof(float), cudaMemcpyHostToDevice);
        }

        cudaMallocManaged(&cuda_parent_IDs[gpu], 2 * sizeof(int *));
        for (int row = 0; row < 2; row++)
        {
            cudaMalloc((void **)&cuda_parent_IDs[gpu][row], num_Parent_sequence * sizeof(int));
            cudaMemcpy(cuda_parent_IDs[gpu][row], parent_IDs[row], num_Parent_sequence * sizeof(int), cudaMemcpyHostToDevice);
        }

        // cudaMalloc(&cuda_cell_Index[gpu], (num_Cells + 1) * sizeof(int));
        cudaMallocManaged(&cuda_cell_Index[gpu], (num_Cells + 1) * sizeof(int));
        cudaMemcpy(cuda_cell_Index[gpu], cell_Index, (num_Cells + 1) * sizeof(int), cudaMemcpyHostToDevice);

        // num_Cells, recombination_Hotspots,

        cudaMallocManaged(&cuda_totals_Progeny_Selectivity[gpu], num_Cells * sizeof(float *));
        for (int row = 0; row < num_Cells; row++)
        {
            cudaMalloc((void **)&cuda_totals_Progeny_Selectivity[gpu][row], recombination_Hotspots * sizeof(float));
            cudaMemcpy(cuda_totals_Progeny_Selectivity[gpu][row], totals_Progeny_Selectivity[row], recombination_Hotspots * sizeof(float), cudaMemcpyHostToDevice);
        }

        cudaMallocManaged(&cuda_progeny_Configuration[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_progeny_Configuration[gpu][row], (1 + recombination_Hotspots) * sizeof(int));
            cudaMemcpy(cuda_progeny_Configuration[gpu][row], progeny_Configuration[row + start_stop_Per_GPU[gpu].first + start_Progeny], (1 + recombination_Hotspots) * sizeof(int), cudaMemcpyHostToDevice);
        }

        cudaStreamCreate(&streams[gpu]);
    }

    cout << "Loaded " << num_Progeny_being_Processed << " sequence(s) and all pre-requisites to the GPU(s)\nInitiating GPU(s) execution\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cuda_Progeny_Complete_Configuration<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(genome_Length,
                                                                                                                                            cuda_Reference_fitness_survivability_proof_reading[gpu],
                                                                                                                                            cuda_num_effect_Segregating_sites[gpu],
                                                                                                                                            cuda_sequence_Survivability_changes[gpu],
                                                                                                                                            recombination_Hotspots,
                                                                                                                                            cuda_recombination_hotspot_parameters[gpu],
                                                                                                                                            cuda_tot_prob_selectivity[gpu],
                                                                                                                                            mutation_Hotspots,
                                                                                                                                            cuda_A_0_mutation[gpu],
                                                                                                                                            cuda_T_1_mutation[gpu],
                                                                                                                                            cuda_G_2_mutation[gpu],
                                                                                                                                            cuda_C_3_mutation[gpu],
                                                                                                                                            cuda_mutation_hotspot_parameters[gpu],
                                                                                                                                            cuda_parent_Sequences[gpu], cuda_parent_IDs[gpu],
                                                                                                                                            cuda_sequence_Configuration_standard[gpu],
                                                                                                                                            cuda_cell_Index[gpu], num_Cells,
                                                                                                                                            cuda_totals_Progeny_Selectivity[gpu],
                                                                                                                                            cuda_progeny_Configuration[gpu],
                                                                                                                                            cuda_progeny_Sequences[gpu],
                                                                                                                                            cuda_Dead_or_Alive[gpu],
                                                                                                                                            start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first);
    }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaStreamSynchronize(streams[gpu]);
    }

    cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

    int **progeny_Configuration_Filled = functions.create_INT_2D_arrays(num_Progeny_being_Processed, (1 + recombination_Hotspots));
    int **progeny_Sequences = functions.create_INT_2D_arrays(num_Progeny_being_Processed, genome_Length);
    int *Dead_or_Alive = (int *)malloc(sizeof(int) * num_Progeny_being_Processed);

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMemcpy(progeny_Configuration_Filled[start_stop_Per_GPU[gpu].first + row], cuda_progeny_Configuration[gpu][row], (1 + recombination_Hotspots) * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(progeny_Sequences[start_stop_Per_GPU[gpu].first + row], cuda_progeny_Sequences[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
        }

        cudaMemcpy(Dead_or_Alive + start_stop_Per_GPU[gpu].first, cuda_Dead_or_Alive[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int), cudaMemcpyDeviceToHost);
    }

    cout << "Data received by host\n";

    cout << "Terminating GPU streams: ";
    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        cudaFree(cuda_Reference_fitness_survivability_proof_reading);

        for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        {
            cudaFree(cuda_sequence_Survivability_changes[gpu][row]);
        }
        cudaFree(cuda_sequence_Survivability_changes[gpu]);

        for (int row = 0; row < recombination_Hotspots; row++)
        {
            cudaFree(cuda_recombination_hotspot_parameters[gpu][row]);
        }
        cudaFree(cuda_recombination_hotspot_parameters[gpu]);

        cudaFree(cuda_tot_prob_selectivity);

        for (int row = 0; row < mutation_Hotspots; row++)
        {
            cudaFree(cuda_A_0_mutation[gpu][row]);
            cudaFree(cuda_T_1_mutation[gpu][row]);
            cudaFree(cuda_G_2_mutation[gpu][row]);
            cudaFree(cuda_C_3_mutation[gpu][row]);

            cudaFree(cuda_mutation_hotspot_parameters[gpu][row]);
        }
        cudaFree(cuda_A_0_mutation[gpu]);
        cudaFree(cuda_T_1_mutation[gpu]);
        cudaFree(cuda_G_2_mutation[gpu]);
        cudaFree(cuda_C_3_mutation[gpu]);

        cudaFree(cuda_mutation_hotspot_parameters[gpu]);

        for (int row = 0; row < num_Parent_sequence; row++)
        {
            cudaFree(cuda_parent_Sequences[gpu][row]);
            cudaFree(cuda_sequence_Configuration_standard[gpu][row]);
        }
        cudaFree(cuda_parent_Sequences[gpu]);
        cudaFree(cuda_sequence_Configuration_standard[gpu]);

        for (int row = 0; row < 2; row++)
        {
            cudaFree(cuda_parent_IDs[gpu][row]);
        }
        cudaFree(cuda_parent_IDs[gpu]);

        cudaFree(cuda_cell_Index);

        for (int row = 0; row < num_Cells; row++)
        {
            cudaFree(cuda_totals_Progeny_Selectivity[gpu][row]);
        }
        cudaFree(cuda_totals_Progeny_Selectivity[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaFree(cuda_progeny_Configuration[gpu][row]);
            cudaFree(cuda_progeny_Sequences[gpu][row]);
        }
        cudaFree(cuda_progeny_Configuration);
        cudaFree(cuda_progeny_Sequences);

        cudaFree(cuda_Dead_or_Alive);

        cudaStreamDestroy(streams[gpu]);
    }
    cout << "Completed\n";

    for (int test = 0; test < num_Progeny_being_Processed; test++)
    {
        for (size_t i = 0; i < recombination_Hotspots + 1; i++)
        {
            cout << progeny_Configuration_Filled[test][i] << "\t";
        }
        cout << endl;
    }
    cout << endl;

    for (int test = 0; test < 1; test++)
    {
        for (size_t i = 0; i < genome_Length; i++)
        {
            cout << progeny_Sequences[test][i];
        }
        cout << endl;
    }
    cout << endl;
    for (int test = 0; test < num_Progeny_being_Processed; test++)
    {
        cout << Dead_or_Alive[test] << endl;
    }

    // test when more than one sequence in cell
    // get_Name() << "_" << tissue_Names[entry_array[tissue]] << "_" << current_Generation << "_" << to_string(indexes_Written[transfers])
}

__global__ void cuda_Progeny_Configurator(int num_Parents_to_Process, int start_Index,
                                          float **cuda_sequence_Configuration_standard, int recombination_Hotspots,
                                          int **cuda_progeny_Configuration, int *cuda_progeny_Stride, int remove_Back)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Parents_to_Process)
    {
        int parent_Index = tid + start_Index;
        int progeny_Fill_start = cuda_progeny_Stride[parent_Index] - remove_Back;
        int progeny_Fill_end = cuda_progeny_Stride[parent_Index + 1] - remove_Back;

        for (int progeny = 0; progeny < cuda_sequence_Configuration_standard[parent_Index][0]; progeny++)
        {
            cuda_progeny_Configuration[progeny_Fill_start + progeny][0] = parent_Index;

            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                if (progeny < cuda_sequence_Configuration_standard[parent_Index][(hotspot * 2) + 2])
                {
                    cuda_progeny_Configuration[progeny_Fill_start + progeny][hotspot + 1] = parent_Index;
                }
                else
                {
                    cuda_progeny_Configuration[progeny_Fill_start + progeny][hotspot + 1] = -1;
                }
            }
        }

        for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
        {
            if (cuda_sequence_Configuration_standard[parent_Index][(hotspot * 2) + 2] > 0)
            {
                if (cuda_sequence_Configuration_standard[parent_Index][(hotspot * 2) + 2] != cuda_sequence_Configuration_standard[parent_Index][0])
                {
                    curandState state;
                    curand_init(clock64(), tid, 0, &state);
                    for (int i = progeny_Fill_start; i < progeny_Fill_end - 1; i++)
                    {

                        int j = curand(&state) % (progeny_Fill_end - i) + i;

                        int temp = cuda_progeny_Configuration[i][hotspot + 1];
                        cuda_progeny_Configuration[i][hotspot + 1] = cuda_progeny_Configuration[j][hotspot + 1];
                        cuda_progeny_Configuration[j][hotspot + 1] = temp;
                    }
                }
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void node_within_host::progeny_Configurator(functions_library &functions,
                                            float **cuda_sequence_Configuration_standard, int recombination_Hotspots,
                                            int start_Index, int num_Parents_to_Process,
                                            int **progeny_Configuration, int *cuda_progeny_Stride, int progeny_Total, int remove_Back)
{
    cout << "Configuring " << num_Parents_to_Process << " parents' " << progeny_Total << " progeny\n";

    int **cuda_progeny_Configuration = functions.create_CUDA_2D_int(progeny_Total, 1 + recombination_Hotspots);

    cuda_Progeny_Configurator<<<functions.tot_Blocks_array[0], functions.tot_ThreadsperBlock_array[0]>>>(num_Parents_to_Process, start_Index,
                                                                                                         cuda_sequence_Configuration_standard, recombination_Hotspots,
                                                                                                         cuda_progeny_Configuration, cuda_progeny_Stride, remove_Back);
    cudaDeviceSynchronize();

    for (int row = 0; row < progeny_Total; row++)
    {
        cudaMemcpy(progeny_Configuration[row + remove_Back], cuda_progeny_Configuration[row], (recombination_Hotspots + 1) * sizeof(cuda_progeny_Configuration[0][0]), cudaMemcpyDeviceToHost);
    }

    functions.clear_Array_INT(cuda_progeny_Configuration, progeny_Total);
}

__device__ float generateExponential(curandState *state, float lambda)
{
    float u = curand_uniform(state);
    return -logf(u) / lambda;
}

__global__ void cuda_Parent_configuration(int num_Sequences, int **sequence_INT, int genome_Length, char *sites, float **cuda_sequence_Configuration_standard,
                                          float *cuda_Reference_fitness_survivability_proof_reading, int *cuda_num_effect_Segregating_sites,
                                          float **cuda_sequence_Fitness_changes, float **cuda_sequence_Proof_reading_changes,
                                          int recombination_Hotspots, float **cuda_recombination_hotspot_parameters,
                                          int *cuda_recombination_prob_Stride, float **cuda_recombination_Prob_matrix,
                                          int *cuda_recombination_select_Stride, float **cuda_recombination_Select_matrix,
                                          float *cuda_progeny_distribution_parameters_Array)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Sequences)
    {
        int site_Start = tid * genome_Length;
        int site_End = site_Start + genome_Length;

        int bp_Pos = 0;

        for (int site = site_Start; site < site_End; site++)
        {
            if (sites[site] == 'A' || sites[site] == 'a' || sites[site] == '0')
            {
                sequence_INT[tid][bp_Pos] = 0;
            }
            else if (sites[site] == 'T' || sites[site] == 't' || sites[site] == '1')
            {
                sequence_INT[tid][bp_Pos] = 1;
            }
            else if (sites[site] == 'G' || sites[site] == 'g' || sites[site] == '2')
            {
                sequence_INT[tid][bp_Pos] = 2;
            }
            else if (sites[site] == 'C' || sites[site] == 'c' || sites[site] == '3')
            {
                sequence_INT[tid][bp_Pos] = 3;
            }

            bp_Pos++;
        }

        // Fitness
        float fitness = cuda_Reference_fitness_survivability_proof_reading[0];
        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[0]; pos++)
        {
            if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 0)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][1];
            }
            else if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 1)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][2];
            }
            else if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 2)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][3];
            }
            else if (sequence_INT[tid][(int)cuda_sequence_Fitness_changes[pos][0] - 1] == 3)
            {
                fitness = fitness * cuda_sequence_Fitness_changes[pos][4];
            }
        }

        curandState localState;
        curand_init(clock64(), tid, 0, &localState);

        int progeny = 1;

        if (cuda_progeny_distribution_parameters_Array[0] == 0)
        {
            int failures = 0;
            int successes = 0;

            while (successes < cuda_progeny_distribution_parameters_Array[1])
            {
                float rand_num = curand_uniform(&localState);
                if (rand_num < cuda_progeny_distribution_parameters_Array[2])
                {
                    successes++;
                }
                else
                {
                    failures++;
                }
            }

            progeny = failures;
        }
        else if (cuda_progeny_distribution_parameters_Array[0] == 1)
        {
            // progeny = (int)rand_gamma_node(&localState, cuda_progeny_distribution_parameters_Array[1], cuda_progeny_distribution_parameters_Array[2]);

            float sum = 0.0f;
            for (int j = 0; j < cuda_progeny_distribution_parameters_Array[1]; ++j)
            {
                sum += generateExponential(&localState, 1.0f / cuda_progeny_distribution_parameters_Array[2]);
            }
            progeny = (int)sum;
        }
        else if (cuda_progeny_distribution_parameters_Array[0] == 2)
        {
            progeny = curand_poisson(&localState, cuda_progeny_distribution_parameters_Array[1]);
        }

        cuda_sequence_Configuration_standard[tid][0] = progeny * fitness;

        // proof reading
        if (cuda_Reference_fitness_survivability_proof_reading[2] != -1)
        {
            float proof_Reading = cuda_Reference_fitness_survivability_proof_reading[2];

            for (int pos = 0; pos < cuda_num_effect_Segregating_sites[2]; pos++)
            {
                if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 0)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][1];
                }
                else if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 1)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][2];
                }
                else if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 2)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][3];
                }
                else if (sequence_INT[tid][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] == 3)
                {
                    proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][4];
                }
            }
            if (proof_Reading > 1)
            {
                proof_Reading = 1;
            }
            else if (proof_Reading < 0)
            {
                proof_Reading = 0;
            }
            cuda_sequence_Configuration_standard[tid][1] = proof_Reading;
        }
        else
        {
            cuda_sequence_Configuration_standard[tid][1] = -1;
        }

        if (recombination_Hotspots > 0)
        {
            for (int hotspot = 0; hotspot < recombination_Hotspots; hotspot++)
            {
                int index_Progeny = (hotspot * 2) + 2;
                int index_Selectivity = index_Progeny + 1;

                float probability = cuda_recombination_hotspot_parameters[hotspot][2];
                float selectivity = cuda_recombination_hotspot_parameters[hotspot][3];

                for (int stride = cuda_recombination_prob_Stride[hotspot]; stride < cuda_recombination_prob_Stride[hotspot + 1]; stride++)
                {
                    if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 0)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][1];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 1)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][2];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 2)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][3];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Prob_matrix[stride][0] - 1] == 3)
                    {
                        probability = probability + cuda_recombination_Prob_matrix[stride][4];
                    }
                }

                for (int stride = cuda_recombination_select_Stride[hotspot]; stride < cuda_recombination_select_Stride[hotspot + 1]; stride++)
                {
                    if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 0)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][1];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 1)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][2];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 2)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][3];
                    }
                    else if (sequence_INT[tid][(int)cuda_recombination_Select_matrix[stride][0] - 1] == 3)
                    {
                        selectivity = selectivity * cuda_recombination_Select_matrix[stride][4];
                    }
                }

                if (probability > 1)
                {
                    probability = 1;
                }
                else if (probability < 0)
                {
                    probability = 0;
                }

                int hotspot_Progeny = 0;

                for (int trial = 0; trial < cuda_sequence_Configuration_standard[tid][0]; trial++)
                {
                    if (curand_uniform(&localState) < probability)
                    {
                        hotspot_Progeny++;
                    }
                }

                cuda_sequence_Configuration_standard[tid][index_Progeny] = hotspot_Progeny;

                cuda_sequence_Configuration_standard[tid][index_Selectivity] = selectivity;
            }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void node_within_host::process_Sequences_get_Configuration(functions_library &functions,
                                                           vector<string> &collected_Sequences, int *CUDA_device_IDs, int &num_Cuda_devices, int &genome_Length,
                                                           float *Reference_fitness_survivability_proof_reading,
                                                           int *mutation_recombination_proof_Reading_availability,
                                                           int *num_effect_Segregating_sites,
                                                           float **sequence_Fitness_changes,
                                                           float **sequence_Proof_reading_changes,
                                                           int recombination_Hotspots,
                                                           float **recombination_hotspot_parameters,
                                                           int *tot_prob_selectivity,
                                                           int *recombination_prob_Stride,
                                                           int *recombination_select_Stride,
                                                           float **recombination_Prob_matrix,
                                                           float **recombination_Select_matrix,
                                                           float *progeny_distribution_parameters_Array,
                                                           int **parent_Sequences,
                                                           float **sequence_Configuration_standard,
                                                           int &start_Index)
{
    int num_of_Sequences_current = collected_Sequences.size();
    cout << "\nConfiguring multi gpu distribution of " << num_of_Sequences_current << " sequence(s)\n";

    int standard_num_per_GPU = num_of_Sequences_current / num_Cuda_devices;
    int remainder = num_of_Sequences_current % num_Cuda_devices;

    vector<pair<int, int>> start_stop_Per_GPU;

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        int start = gpu * standard_num_per_GPU;
        int stop = start + standard_num_per_GPU;

        start_stop_Per_GPU.push_back(make_pair(start, stop));
    }

    start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

    string all_Sequences = "";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        for (int sequence = start_stop_Per_GPU[gpu].first; sequence < start_stop_Per_GPU[gpu].second; sequence++)
        {
            all_Sequences.append(collected_Sequences[sequence]);
        }
    }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        for (int sequence = start_stop_Per_GPU[gpu].first; sequence < start_stop_Per_GPU[gpu].second; sequence++)
        {
            all_Sequences.append(collected_Sequences[sequence]);
        }
    }

    char *full_Char;
    full_Char = (char *)malloc((all_Sequences.size() + 1) * sizeof(char));
    strcpy(full_Char, all_Sequences.c_str());

    cudaStream_t streams[num_Cuda_devices];
    cudaDeviceProp deviceProp;

    char *cuda_full_Char[num_Cuda_devices];
    int **cuda_Sequence[num_Cuda_devices];
    float **cuda_sequence_Configuration_standard[num_Cuda_devices];

    float *cuda_Reference_fitness_survivability_proof_reading[num_Cuda_devices];
    int *cuda_mutation_recombination_proof_Reading_availability[num_Cuda_devices];
    int *cuda_num_effect_Segregating_sites[num_Cuda_devices];

    float **cuda_sequence_Fitness_changes[num_Cuda_devices];
    float **cuda_sequence_Proof_reading_changes[num_Cuda_devices];

    float **cuda_recombination_hotspot_parameters[num_Cuda_devices];
    int *cuda_tot_prob_selectivity[num_Cuda_devices];
    int *cuda_recombination_prob_Stride[num_Cuda_devices];
    int *cuda_recombination_select_Stride[num_Cuda_devices];
    float **cuda_recombination_Prob_matrix[num_Cuda_devices];
    float **cuda_recombination_Select_matrix[num_Cuda_devices];

    float *cuda_progeny_distribution_parameters_Array[num_Cuda_devices];

    // cout << "Sites: " << num_effect_Segregating_sites[0] << endl;
    // for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
    // {
    //     for (size_t i = 0; i < 5; i++)
    //     {
    //         cout << sequence_Fitness_changes[row][i] << "\t";
    //     }
    //     cout << endl;
    // }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaGetDeviceProperties(&deviceProp, gpu);
        cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

        cudaMalloc(&cuda_full_Char[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(char));
        cudaMemcpy(cuda_full_Char[gpu], full_Char + (start_stop_Per_GPU[gpu].first * genome_Length), (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * genome_Length * sizeof(char), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_Sequence[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(int *));
        cudaMallocManaged(&cuda_sequence_Configuration_standard[gpu], (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first) * sizeof(float *));
        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMalloc((void **)&cuda_Sequence[gpu][row], genome_Length * sizeof(int));
            cudaMalloc((void **)&cuda_sequence_Configuration_standard[gpu][row], (3 + (2 * recombination_Hotspots)) * sizeof(float));
        }

        cudaMallocManaged(&cuda_Reference_fitness_survivability_proof_reading[gpu], 3 * sizeof(float));
        cudaMemcpy(cuda_Reference_fitness_survivability_proof_reading[gpu], Reference_fitness_survivability_proof_reading, 3 * sizeof(float), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_mutation_recombination_proof_Reading_availability[gpu], 3 * sizeof(int));
        cudaMemcpy(cuda_mutation_recombination_proof_Reading_availability[gpu], mutation_recombination_proof_Reading_availability, 3 * sizeof(int), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_num_effect_Segregating_sites[gpu], 3 * sizeof(int));
        cudaMemcpy(cuda_num_effect_Segregating_sites[gpu], num_effect_Segregating_sites, 3 * sizeof(int), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_progeny_distribution_parameters_Array[gpu], 3 * sizeof(float));
        cudaMemcpy(cuda_progeny_distribution_parameters_Array[gpu], progeny_distribution_parameters_Array, 3 * sizeof(float), cudaMemcpyHostToDevice);

        cudaMallocManaged(&cuda_sequence_Fitness_changes[gpu], num_effect_Segregating_sites[0] * sizeof(float *));
        for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
        {
            cudaMalloc((void **)&cuda_sequence_Fitness_changes[gpu][row], 5 * sizeof(float));
            cudaMemcpy(cuda_sequence_Fitness_changes[gpu][row], sequence_Fitness_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        }

        // cudaMallocManaged(&cuda_sequence_Survivability_changes[gpu], num_effect_Segregating_sites[1] * sizeof(float *));
        // for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        // {
        //     cudaMalloc((void **)&cuda_sequence_Survivability_changes[gpu][row], 5 * sizeof(float));
        //     cudaMemcpy(cuda_sequence_Survivability_changes[gpu][row], sequence_Survivability_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        // }

        cudaMallocManaged(&cuda_sequence_Proof_reading_changes[gpu], num_effect_Segregating_sites[2] * sizeof(float *));
        for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
        {
            cudaMalloc((void **)&cuda_sequence_Proof_reading_changes[gpu][row], 5 * sizeof(float));
            cudaMemcpy(cuda_sequence_Proof_reading_changes[gpu][row], sequence_Proof_reading_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
        }

        cudaMallocManaged(&cuda_tot_prob_selectivity[gpu], 2 * sizeof(int));
        cudaMallocManaged(&cuda_recombination_prob_Stride[gpu], (recombination_Hotspots + 1) * sizeof(int));
        cudaMallocManaged(&cuda_recombination_select_Stride[gpu], (recombination_Hotspots + 1) * sizeof(int));

        cudaMallocManaged(&cuda_recombination_hotspot_parameters[gpu], recombination_Hotspots * sizeof(float *));
        cudaMallocManaged(&cuda_recombination_Prob_matrix[gpu], tot_prob_selectivity[0] * sizeof(float *));
        cudaMallocManaged(&cuda_recombination_Select_matrix[gpu], tot_prob_selectivity[1] * sizeof(float *));

        if (recombination_Hotspots > 0)
        {
            cudaMemcpy(cuda_tot_prob_selectivity[gpu], tot_prob_selectivity, 2 * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_recombination_prob_Stride[gpu], recombination_prob_Stride, (recombination_Hotspots + 1) * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(cuda_recombination_select_Stride[gpu], recombination_select_Stride, (recombination_Hotspots + 1) * sizeof(int), cudaMemcpyHostToDevice);

            for (int row = 0; row < recombination_Hotspots; row++)
            {
                cudaMalloc((void **)&cuda_recombination_hotspot_parameters[gpu][row], 4 * sizeof(float));
                cudaMemcpy(cuda_recombination_hotspot_parameters[gpu][row], recombination_hotspot_parameters[row], 4 * sizeof(float), cudaMemcpyHostToDevice);
            }

            for (int row = 0; row < tot_prob_selectivity[0]; row++)
            {
                cudaMalloc((void **)&cuda_recombination_Prob_matrix[gpu][row], 5 * sizeof(float));
                cudaMemcpy(cuda_recombination_Prob_matrix[gpu][row], recombination_Prob_matrix[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
            }

            for (int row = 0; row < tot_prob_selectivity[1]; row++)
            {
                cudaMalloc((void **)&cuda_recombination_Select_matrix[gpu][row], 5 * sizeof(float));
                cudaMemcpy(cuda_recombination_Select_matrix[gpu][row], recombination_Select_matrix[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
            }
        }

        cudaStreamCreate(&streams[gpu]);
    }

    cout << "Loaded " << num_of_Sequences_current << " sequence(s) and all pre-requisites to the GPU(s)\nInitiating GPU(s) execution\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        // (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first, cuda_Sequence[gpu], genome_Length, cuda_full_Char[gpu]);
        cuda_Parent_configuration<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first, cuda_Sequence[gpu], genome_Length, cuda_full_Char[gpu], cuda_sequence_Configuration_standard[gpu],
                                                                                                                                  cuda_Reference_fitness_survivability_proof_reading[gpu], cuda_num_effect_Segregating_sites[gpu],
                                                                                                                                  cuda_sequence_Fitness_changes[gpu], cuda_sequence_Proof_reading_changes[gpu],
                                                                                                                                  recombination_Hotspots, cuda_recombination_hotspot_parameters[gpu],
                                                                                                                                  cuda_recombination_prob_Stride[gpu], cuda_recombination_Prob_matrix[gpu],
                                                                                                                                  cuda_recombination_select_Stride[gpu], cuda_recombination_Select_matrix[gpu],
                                                                                                                                  cuda_progeny_distribution_parameters_Array[gpu]);
    }

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);
        cudaStreamSynchronize(streams[gpu]);
    }

    cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaMemcpy(parent_Sequences[start_stop_Per_GPU[gpu].first + row + start_Index], cuda_Sequence[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(sequence_Configuration_standard[start_stop_Per_GPU[gpu].first + row + start_Index], cuda_sequence_Configuration_standard[gpu][row], (2 + (2 * recombination_Hotspots)) * sizeof(float), cudaMemcpyDeviceToHost);
        }
    }
    cout << "Data received by host\n";

    cout << "Terminating GPU streams: ";
    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
    {
        cudaSetDevice(CUDA_device_IDs[gpu]);

        cudaFree(cuda_full_Char);

        for (int row = 0; row < (start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first); row++)
        {
            cudaFree(cuda_Sequence[gpu][row]);
            cudaFree(cuda_sequence_Configuration_standard[gpu][row]);
        }
        cudaFree(cuda_Sequence[gpu]);
        cudaFree(cuda_sequence_Configuration_standard[gpu]);

        cudaFree(cuda_Reference_fitness_survivability_proof_reading);
        cudaFree(cuda_mutation_recombination_proof_Reading_availability);
        cudaFree(cuda_num_effect_Segregating_sites);

        for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
        {
            cudaFree(cuda_sequence_Fitness_changes[gpu][row]);
        }
        cudaFree(cuda_sequence_Fitness_changes[gpu]);

        for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
        {
            cudaFree(cuda_sequence_Proof_reading_changes[gpu][row]);
        }
        cudaFree(cuda_sequence_Proof_reading_changes[gpu]);

        cudaFree(cuda_tot_prob_selectivity);
        cudaFree(cuda_recombination_prob_Stride);
        cudaFree(cuda_recombination_select_Stride);
        if (recombination_Hotspots > 0)
        {
            for (int row = 0; row < recombination_Hotspots; row++)
            {
                cudaFree(cuda_recombination_hotspot_parameters[gpu][row]);
            }

            for (int row = 0; row < tot_prob_selectivity[0]; row++)
            {
                cudaFree(cuda_recombination_Prob_matrix[gpu][row]);
            }

            for (int row = 0; row < tot_prob_selectivity[1]; row++)
            {
                cudaFree(cuda_recombination_Select_matrix[gpu][row]);
            }
        }
        cudaFree(cuda_recombination_hotspot_parameters[gpu]);
        cudaFree(cuda_recombination_Prob_matrix[gpu]);
        cudaFree(cuda_recombination_Select_matrix[gpu]);

        cudaFree(cuda_progeny_distribution_parameters_Array);

        cudaStreamDestroy(streams[gpu]);
    }
    cout << "Completed\n";
}

vector<int> node_within_host::assign_Cells(functions_library &functions, int **parents_in_Tissue, int num_Viral_particles, int &tissue,
                                           int distribution_Type, float &parameter_1, float &parameter_2,
                                           set<int> &check_to_Remove,
                                           int &gen_Phase, float &variable_1, float &variable_2,
                                           mt19937 &gen)
{
    cout << "\nAssigning cell(s) their virulant particle(s)\n";

    // if (parents_Prev_generation != 0)
    // {
    //     num_Viral_particles = parents_Prev_generation;
    // }

    vector<int> start_Stop_cells;

    int cells_Assigned = 0;
    int particles_Assigned = 0;

    start_Stop_cells.push_back(0);

    // int cells_Full = 0;

    int index_Track_removed = 0;
    int particle_ID = 0;

    // test
    // cell_Limit[tissue] = 2;

    vector<int> removals(check_to_Remove.begin(), check_to_Remove.end());

    do
    {
        // if (cell_Limit[tissue] != -1 && cells_Assigned >= cell_Limit[tissue])
        // {
        //     cells_Full = 1;
        //     break;
        // }

        int num_Particles_in_Cell;

        if (distribution_Type == 0)
        {
            binomial_distribution<int> num_Particles(parameter_1, parameter_2);
            num_Particles_in_Cell = num_Particles(gen);
        }
        else
        {
            gamma_distribution<float> num_Particles(parameter_1, parameter_2);
            num_Particles_in_Cell = (int)num_Particles(gen);
        }

        if (num_Particles_in_Cell > 0)
        {
            for (int cell = 0; cell < num_Particles_in_Cell; cell++)
            {
                parents_in_Tissue[1][particles_Assigned] = cells_Assigned;
                //// Account for dead
                if (index_Track_removed < removals.size())
                {
                    while (particle_ID == removals[index_Track_removed])
                    {
                        particle_ID++;
                        index_Track_removed++;
                    }
                }

                parents_in_Tissue[0][particles_Assigned] = particle_ID;

                particle_ID++;
                particles_Assigned++;

                if (particles_Assigned >= num_Viral_particles)
                {
                    break;
                }
            }

            start_Stop_cells.push_back(particles_Assigned);
            cells_Assigned++;
        }

    } while (particles_Assigned < num_Viral_particles);

    // cout << cells_Assigned - 1 << endl;
    srand(time(0));
    random_shuffle(parents_in_Tissue[0], parents_in_Tissue[0] + num_Viral_particles);

    if (gen_Phase == 1 || gen_Phase == 2)
    {
        int new_Parent_Count = -1;

        if (gen_Phase == 1)
        {
            cout << "Stationary phase\n";
            if (num_Viral_particles >= parents_Prev_generation)
            {
                normal_distribution<float> distribution(parents_Prev_generation, variable_1);
                new_Parent_Count = distribution(gen);
                // cout << "parents_Prev_generation: " << parents_Prev_generation << endl;
                // cout << new_Parent_Count << endl;
                if (new_Parent_Count < num_Viral_particles && new_Parent_Count >= 0)
                {
                    cout << "Parent population maintained at: " << new_Parent_Count << endl;
                }
                else
                {
                    new_Parent_Count = -1;
                }
            }
        }
        else if (gen_Phase == 2)
        {
            cout << "Depriciation phase\n";
            if (num_Viral_particles >= parents_Prev_generation)
            {
                new_Parent_Count = functions.beta_Distribution(variable_1, variable_2, gen) * parents_Prev_generation;
                new_Parent_Count = parents_Prev_generation - new_Parent_Count;
                cout << "Parent population reduced to: " << new_Parent_Count << endl;
            }
        }
        if (new_Parent_Count != -1)
        {
            int **temp = functions.create_INT_2D_arrays(2, new_Parent_Count);

            for (int parent = 0; parent < new_Parent_Count; parent++)
            {
                temp[0][parent] = parents_in_Tissue[0][parent];
                temp[1][parent] = parents_in_Tissue[1][parent];
            }

            functions.clear_Array_int_CPU(parents_in_Tissue, 2);

            parents_in_Tissue = functions.create_INT_2D_arrays(2, new_Parent_Count);

            for (int parent = 0; parent < new_Parent_Count; parent++)
            {
                parents_in_Tissue[0][parent] = temp[0][parent];
                parents_in_Tissue[1][parent] = temp[1][parent];
            }

            functions.clear_Array_int_CPU(temp, 2);

            cout << "Resizing parent cell array\n";

            vector<int> temp_Cells;
            if (new_Parent_Count > 0)
            {
                cout << "Configuring cell index\n";
                for (int cell = 0; cell < parents_in_Tissue[1][new_Parent_Count - 1] + 1; cell++)
                {
                    temp_Cells.push_back(start_Stop_cells[cell]);
                }

                if (temp_Cells[temp_Cells.size() - 1] > new_Parent_Count)
                {
                    temp_Cells[temp_Cells.size() - 1] = new_Parent_Count;
                }
            }
            else
            {
                temp_Cells.push_back(0);
            }

            start_Stop_cells.clear();
            start_Stop_cells = temp_Cells;

            parents_Prev_generation = new_Parent_Count;
        }
        else
        {
            parents_Prev_generation = num_Viral_particles;
        }
    }
    else
    {
        cout << "Neutral phase, all particles are viable\n";
        parents_Prev_generation = num_Viral_particles;
    }

    if ((start_Stop_cells.size() - 1) > cell_Limit[tissue])
    {
        vector<int> temp;
        for (int cell = 0; cell < cell_Limit[tissue] + 1; cell++)
        {
            temp.push_back(start_Stop_cells[cell]);
        }
        start_Stop_cells.clear();
        start_Stop_cells = temp;
        parents_Prev_generation = start_Stop_cells[start_Stop_cells.size() - 1];
    }

    return start_Stop_cells;
}

void node_within_host::intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions)
{
    current_Viral_load_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
    dead_Particle_count = (int *)malloc(sizeof(int) * num_Tissues);
    current_Generation = 0;

    set<int> init_removed_by_Transfer_Indexes;

    for (int tissue = 0; tissue < num_Tissues; tissue++)
    {
        vector<string> tissue_Sequence;
        functions.config_Folder(host_Folder + "/" + to_string(tissue), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(tissue + 1));

        tissue_Sequences.push_back(tissue_Sequence);
        current_Viral_load_per_Tissue[tissue] = 0;
        dead_Particle_count[tissue] = 0;

        removed_by_Transfer_Indexes.push_back(init_removed_by_Transfer_Indexes);
    }
}

int node_within_host::get_Load(int &num_tissues_Calc, int *tissue_array)
{
    int sum = 0;

    for (int tissue = 0; tissue < num_tissues_Calc; tissue++)
    {
        sum = sum + current_Viral_load_per_Tissue[tissue_array[tissue]];
    }

    return sum;
}

int node_within_host::infectious_status(int &num_tissues, int *tissue_array)
{
    if (get_Load(num_tissues, tissue_array) >= infectious_Load)
    {
        set_Infectious();
        return 1;
    }
    else
    {
        return 0;
    }
}
int node_within_host::terminal_status(int &num_tissues, int *tissue_array)
{
    if (get_Load(num_tissues, tissue_array) >= terminal_Load)
    {
        set_Dead();
        return 1;
    }
    else
    {
        return 0;
    }
}

string node_within_host::get_Name()
{
    return to_string(cave_ID) + "_" + to_string(host_ID);
}

string node_within_host::get_Status()
{
    return this->status;
}

int node_within_host::get_Profile()
{
    return profile_ID;
}

int node_within_host::get_host_Index()
{
    return host_Index;
}

int node_within_host::get_Generation()
{
    return current_Generation;
}

int *node_within_host::get_current_Viral_load_per_Tissue()
{
    return current_Viral_load_per_Tissue;
}

void node_within_host::set_Infected()
{
    this->status = "Infected";
}
void node_within_host::set_Infectious()
{
    this->status = "Infectious";
}
void node_within_host::set_Removed()
{
    this->status = "Removed";
}
void node_within_host::set_Dead()
{
    this->status = "Dead";
}