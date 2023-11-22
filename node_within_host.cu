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

void node_within_host::run_Generation(functions_library &functions,
                                      string source_sequence_Data_folder,
                                      vector<string> &tissue_Names,
                                      int terminal_tissues, int *terminal_array,
                                      int **cell_Distribution_Type, vector<pair<float, float>> &viral_distribution_per_Tissue_param,
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

                for (int tissue = 0; tissue < num_Tissues; tissue++)
                {
                    if (real_Particle_count_per_Tissue[tissue] > 0)
                    {
                        // real_Particle_count_per_Tissue[tissue] = 100;

                        cout << "\nSimulating " << real_Particle_count_per_Tissue[tissue] << " particle(s) for " << tissue_Names[tissue] << " tissue\n"
                             << endl;

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
                            cout << "\nIdentifying transfered viral indexe(s)\n";
                            for (auto it = removed_by_Transfer_Indexes[tissue].begin(); it != removed_by_Transfer_Indexes[tissue].end(); ++it)
                            {
                                int value = *it; // Dereference the iterator to get the value
                                check_to_Remove.insert(value);
                            }
                        }

                        int *parents_in_Tissue = (int *)malloc(sizeof(int) * real_Particle_count_per_Tissue[tissue]);
                        // int **parents_in_Tissue = functions.create_INT_2D_arrays(2, real_Particle_count_per_Tissue[tissue]);

                        // test
                        // check_to_Remove.insert(0);
                        // check_to_Remove.insert(1);
                        // check_to_Remove.insert(5);
                        // check_to_Remove.insert(99);

                        vector<int> start_Stop_cells = assign_Cells(parents_in_Tissue, real_Particle_count_per_Tissue[tissue], tissue,
                                                                    cell_Distribution_Type[profile_ID][tissue], viral_distribution_per_Tissue_param[tissue].first, viral_distribution_per_Tissue_param[tissue].second,
                                                                    check_to_Remove,
                                                                    gen);

                        check_to_Remove.clear();
                        removed_by_Transfer_Indexes[tissue].clear();
                        dead_Particle_count[tissue] = 0;

                        cout << "Number of cell(s) infected: " << start_Stop_cells.size() - 1 << endl;

                        if (start_Stop_cells.size() - 1 > 0)
                        {
                            for (int i = 0; i < start_Stop_cells.size() - 1; i++)
                            {
                                // cout << start_Stop_cells[i] << " : \t" << start_Stop_cells[i + 1] << endl;
                                for (int particle = start_Stop_cells[i]; particle < start_Stop_cells[i + 1]; particle++)
                                {
                                    // cout << parents_in_Tissue[0][particle] << " :\t" << parents_in_Tissue[1][particle] << endl;
                                    cout << parents_in_Tissue[particle] << endl;
                                }
                                cout << endl;
                            }
                        }

                        // functions.clear_Array_int_CPU(parents_in_Tissue, 2);
                        free(parents_in_Tissue);
                        exit(-1);
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

vector<int> node_within_host::assign_Cells(int *parents_in_Tissue, int num_Viral_particles, int &tissue,
                                           int distribution_Type, float &parameter_1, float &parameter_2,
                                           set<int> &check_to_Remove,
                                           mt19937 &gen)
{
    cout << "Assigning cell(s) their virulant particle(s)\n";

    vector<int> start_Stop_cells;

    // int cells_Assigned = 0;
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
                // parents_in_Tissue[1][particles_Assigned] = cells_Assigned;
                //// Account for dead
                if (index_Track_removed < removals.size())
                {
                    while (particle_ID == removals[index_Track_removed])
                    {
                        particle_ID++;
                        index_Track_removed++;
                    }
                }

                parents_in_Tissue[particles_Assigned] = particle_ID;

                particle_ID++;
                particles_Assigned++;

                if (particles_Assigned >= num_Viral_particles)
                {
                    break;
                }
            }

            start_Stop_cells.push_back(particles_Assigned);
            // cells_Assigned++;
        }

    } while (particles_Assigned < num_Viral_particles);

    // cout << cells_Assigned - 1 << endl;

    random_shuffle(&parents_in_Tissue[0], &parents_in_Tissue[num_Viral_particles]);

    if ((start_Stop_cells.size() - 1) > cell_Limit[tissue])
    {
        vector<int> temp;
        for (int cell = 0; cell < cell_Limit[tissue] + 1; cell++)
        {
            temp.push_back(start_Stop_cells[cell]);
        }
        start_Stop_cells.clear();
        start_Stop_cells = temp;
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