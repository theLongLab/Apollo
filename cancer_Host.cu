#include "cancer_Host.cuh"

cancer_Host::cancer_Host()
{
    cout << "\nSTEP 5: Cancer host intialization\n";
}

void cancer_Host::simulate_Generations(functions_library &functions,
                                       int &overall_Generations, float &date_Increment,
                                       int &stop_Type,
                                       int &stop_gen_Mode,
                                       int &stop_generations_Count, float &decimal_Date, float &stop_Date,
                                       vector<string> &tissue_Names,
                                       int &terminal_tissues, int *terminal_array,
                                       string source_sequence_Data_folder,
                                       string &enable_Folder_management, string &enable_Compression,
                                       int &terminal_Load,
                                       string &output_Node_location,
                                       vector<vector<float>> &time_Ratios_per_Tissue, vector<vector<string>> &phase_Type_per_tissue, vector<vector<pair<float, float>>> &phase_paramaters_per_Tissue,
                                       int &max_Cells_at_a_time,
                                       string &multi_Read, int &CPU_cores, int &num_Cuda_devices,
                                       int &genome_Length)
{
    cout << "\nSTEP 6: Conducting simulation\n";

    // this->terminal_Load = terminal_Load;

    random_device rd;
    mt19937 gen(rd());

    this->CPU_cores = CPU_cores;
    this->genome_Length = genome_Length;

    generational_Summary = output_Node_location + "/cancer_Host/node_generational_Summary.csv";
    functions.create_File(generational_Summary, "overall_Generation\tTissue\tnum_Parents\tnum_Progeny\tdead_Progeny");

    sequence_Profiles = output_Node_location + "/cancer_Host/sequence_Profiles.csv";
    // functions.create_File(sequence_Profiles, "Sequence_ID\tTissue");

    sequence_parent_Progeny_relationships = output_Node_location + "/cancer_Host/sequence_parent_Progeny_relationships.csv";
    // functions.create_File(sequence_parent_Progeny_relationships, "Source\tTarget\tType");

    // cells_of_parents_location = output_Node_location + "/cancer_Host/cells_of_Parents.csv";
    // functions.create_File(cells_of_parents_location, "Sequence_ID\tParent_Cell_ID");

    // cells_of_progeny_location = output_Node_location + "/cancer_Host/cells_of_Progeny.csv";
    // functions.create_File(cells_of_progeny_location, "Sequence_ID\tProgeny_Cell_ID");

    do
    {
        cout << "Calculating actual particles in each tissue: \n";
        ////clear array
        int *real_Particle_count_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
        int sum_Check = 0;
        for (int tissue = 0; tissue < num_Tissues; tissue++)
        {
            real_Particle_count_per_Tissue[tissue] = current_cell_load_per_Tissue[tissue] - removed_by_Transfer_Indexes[tissue].size() - dead_Particle_count[tissue];
            cout << tissue_Names[tissue] << " tissue: " << real_Particle_count_per_Tissue[tissue] << endl;
            sum_Check = sum_Check + real_Particle_count_per_Tissue[tissue];
        }

        if (sum_Check > 0)
        {
            if (terminal_status(terminal_tissues, terminal_array, source_sequence_Data_folder, enable_Folder_management, enable_Compression, terminal_Load) != 1)
            {
                vector<vector<pair<int, int>>> indexed_Source_Folders = functions.index_sequence_Folders(source_sequence_Data_folder, num_Tissues, overall_Generations, multi_Read);

                for (int tissue = 0; tissue < num_Tissues; tissue++)
                {
                    if (real_Particle_count_per_Tissue[tissue] > 0)
                    {
                        cout << "\nSimulating " << real_Particle_count_per_Tissue[tissue] << " particle(s) for " << tissue_Names[tissue] << " tissue\n"
                             << endl;

                        cout << "Identifying indexes to remove\n";
                        set<int> check_to_Remove;

                        if (dead_Particle_count[tissue] > 0)
                        {
                            cout << "\nIdentifying dead viral indexe(s)\n";

                            fstream dead_File;
                            dead_File.open(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) + "/dead_List.txt");
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
                                cout << "ERROR: UNABLE TO OPEN DEAD LIST FILE: " << source_sequence_Data_folder << "/" << tissue << "/generation_" << overall_Generations << "/dead_List.txt" << endl;
                                exit(-1);
                            }
                        }

                        if (removed_by_Transfer_Indexes[tissue].size() > 0)
                        {
                            cout << "Identifying transferred viral indexe(s)\n";
                            for (auto it = removed_by_Transfer_Indexes[tissue].begin(); it != removed_by_Transfer_Indexes[tissue].end(); ++it)
                            {
                                int value = *it; // Dereference the iterator to get the value
                                check_to_Remove.insert(value);
                            }
                            removed_by_Transfer_Indexes[tissue].clear();
                        }

                        cout << "\nFilling parent vector: ";
                        int *parents_in_Tissue = (int *)malloc(real_Particle_count_per_Tissue[tissue] * sizeof(int));

                        int fill_Count = 0;
                        int index = 0;
                        do
                        {
                            if (check_to_Remove.find(index) == check_to_Remove.end())
                            {
                                parents_in_Tissue[fill_Count] = index;
                                fill_Count++;
                            }
                            index++;
                        } while (fill_Count < real_Particle_count_per_Tissue[tissue]);

                        cout << "Completed\n";

                        cout << "Shuffling parent vector: ";
                        default_random_engine rng(time(nullptr)); // Seed the random number generator with current time
                        shuffle(parents_in_Tissue, parents_in_Tissue + real_Particle_count_per_Tissue[tissue], rng);
                        cout << "Complete\n";

                        float variable_1, variable_2;
                        string generation_Type = get_generation_Phase(overall_Generations,
                                                                      time_Ratios_per_Tissue[tissue],
                                                                      phase_Type_per_tissue[tissue],
                                                                      phase_paramaters_per_Tissue[tissue],
                                                                      variable_1, variable_2);

                        int parent_population_Count = real_Particle_count_per_Tissue[tissue];

                        if (overall_Generations != 0)
                        {
                            int new_Parent_Count = -1;
                            if (generation_Type == "STATIONARY")
                            {
                                cout << "\nStationary phase\n";
                                if (parent_population_Count >= parents_Prev_generation[tissue])
                                {
                                    normal_distribution<float> distribution(parents_Prev_generation[tissue], variable_1);
                                    new_Parent_Count = (int)(distribution(gen) + 0.5);

                                    if (new_Parent_Count < parent_population_Count && new_Parent_Count >= 0)
                                    {
                                        parent_population_Count = new_Parent_Count;
                                        cout << "Parent population maintained at: " << parent_population_Count << endl;
                                    }
                                }
                            }
                            else if (generation_Type == "DEPRICIATION")
                            {
                                cout << "\nDepriciation phase\n";
                                if (parent_population_Count >= parents_Prev_generation[tissue])
                                {
                                    new_Parent_Count = functions.beta_Distribution(variable_1, variable_2, gen) * parents_Prev_generation[tissue];
                                    new_Parent_Count = parents_Prev_generation[tissue] - new_Parent_Count;
                                    parent_population_Count = new_Parent_Count;
                                    cout << "Parent population reduced to: " << parent_population_Count << endl;
                                }
                            }
                        }

                        check_to_Remove.clear();
                        removed_by_Transfer_Indexes[tissue].clear();
                        dead_Particle_count[tissue] = 0;
                        current_cell_load_per_Tissue[tissue] = 0;

                        int last_index_Seq_Written = 0;

                        vector<pair<int, int>> cells_Rounds_start_stop = get_Rounds(parent_population_Count, max_Cells_at_a_time);

                        for (int cell_Round = 0; cell_Round < cells_Rounds_start_stop.size(); cell_Round++)
                        {
                            int num_of_Cells = cells_Rounds_start_stop[cell_Round].second - cells_Rounds_start_stop[cell_Round].first;
                            cout << "\nProcessing round " << cell_Round + 1 << " of " << cells_Rounds_start_stop.size() << ": " << num_of_Cells << " cell(s)" << endl;

                            simulate_cell_Round(functions, multi_Read, num_Cuda_devices,
                                                num_of_Cells, cells_Rounds_start_stop[cell_Round].first, cells_Rounds_start_stop[cell_Round].second,
                                                parents_in_Tissue, tissue, tissue_Names[tissue],
                                                indexed_Source_Folders[tissue],
                                                source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations),
                                                overall_Generations,
                                                last_index_Seq_Written,
                                                gen);
                        }

                        exit(-1);

                        free(parents_in_Tissue);
                    }
                }
            }
            else
            {
                stop_gen_Mode = 5;
            }
        }
        else
        {
            stop_gen_Mode = 4;
        }

        decimal_Date = decimal_Date + date_Increment;
        overall_Generations++;

        cout << "\nCompleted generation " << overall_Generations << " of time: " << decimal_Date << endl;

        if (stop_gen_Mode == 0)
        {
            if (overall_Generations >= stop_generations_Count)
            {
                stop_Type = 2;
            }
        }
        else
        {
            if (decimal_Date >= stop_Date)
            {
                stop_Type = 3;
            }
        }

        // ! STOP after testing
        stop_Type = 1;

    } while (stop_Type == 0);

    cout << "\nSimulation has concluded: ";
}

void cancer_Host::simulate_cell_Round(functions_library &functions, string &multi_Read, int &num_Cuda_devices,
                                      int &num_of_Cells, int &start, int &stop,
                                      int *parents_in_Tissue, int &tissue, string tissue_Name,
                                      vector<pair<int, int>> &indexed_Tissue_Folder,
                                      string source_sequence_Data_folder,
                                      int &overall_Generations,
                                      int &last_index_Seq_Written,
                                      mt19937 &gen)
{

    sort(parents_in_Tissue + start, parents_in_Tissue + stop);

    vector<int> parent_IDs;
    vector<string> collected_Sequences = find_Sequences_Master(start, tissue, tissue_Name, functions, source_sequence_Data_folder, parents_in_Tissue, num_of_Cells, indexed_Tissue_Folder, overall_Generations, parent_IDs, last_index_Seq_Written, gen);

    if (parent_IDs.size() == collected_Sequences.size())
    {
        int parent_Cells_Found = parent_IDs.size();
        cout << "Number of parents found: " << parent_Cells_Found << endl;
        cout << "Potential mitotic progeny produced: " << to_string(parent_Cells_Found * 2) << endl;

        cout << "\nMain Progeny generation\n";

        int **parent_Sequences;
        parent_Sequences = (int **)malloc(parent_Cells_Found * sizeof(int *));
        for (int row = 0; row < parent_Cells_Found; row++)
        {
            parent_Sequences[row] = (int *)malloc(genome_Length * sizeof(int));
        }

        string all_Sequences = "";

        for (int sequence = 0; sequence < parent_Cells_Found; sequence++)
        {
            all_Sequences.append(collected_Sequences[sequence]);
        }

        char *full_Char;
        full_Char = (char *)malloc((all_Sequences.size() + 1) * sizeof(char));
        strcpy(full_Char, all_Sequences.c_str());
        all_Sequences = "";

        cudaStream_t streams[num_Cuda_devices];
        cudaDeviceProp deviceProp;

        char *cuda_full_Char[num_Cuda_devices];

        cout << "\nGPU CPU parallel processing configuration\n";
        int cpus_to_GPU = 1;
        vector<pair<int, int>> perGPU_start_stop;

        if (num_Cuda_devices > 1 && multi_Read == "YES" && CPU_cores >= num_Cuda_devices)
        {
            cpus_to_GPU = num_Cuda_devices;
            int per_Block = parent_Cells_Found / cpus_to_GPU;
            int remainder_Last = parent_Cells_Found % cpus_to_GPU;

            for (int cpu_Index = 0; cpu_Index < cpus_to_GPU; cpu_Index++)
            {
                int start_From = (cpu_Index * per_Block) + start;
                int stop_To = start_From + per_Block;

                perGPU_start_stop.push_back(make_pair(start_From, stop_To));
            }
            perGPU_start_stop[perGPU_start_stop.size() - 1].second = perGPU_start_stop[perGPU_start_stop.size() - 1].second + remainder_Last;
        }
        else
        {
            perGPU_start_stop.push_back(make_pair(start, stop));
        }

        cout << "Using " << perGPU_start_stop.size() << " core(s) per GPU\n";

        vector<thread> threads_vec;

        for (int core_ID = 0; core_ID < cpus_to_GPU; core_ID++)
        {
            cout << "Initializing CPU core " << core_ID + 1 << endl;

            // threads_vec.push_back(thread{&cancer_Host::replication_Generation_thread, this, core_ID, parents_in_Tissue, ref(perGPU_start_stop[core_ID].first), ref(perGPU_start_stop[core_ID].second), perGPU_start_stop[core_ID].second - perGPU_start_stop[core_ID].first});
        }

        for (thread &t : threads_vec)
        {
            if (t.joinable())
            {
                t.join();
            }
        }

        threads_vec.clear();
    }
    else
    {
        cout << "CRITICAL ERROR: SEQUENCE IDs (" << parent_IDs.size() << ") AND SEQUENCES COLLECTED (" << collected_Sequences.size() << ") DO NOT MATCH\n";
        exit(-1);
    }

    exit(-1);
}

void cancer_Host::replication_Generation_thread(int core_ID,
                                                vector<int> &parent_IDs, vector<string> &collected_Sequences,
                                                int &start, int &stop, int cell_Count)
{
}

vector<string> cancer_Host::find_Sequences_Master(int &offset, int &tissue, string &tissue_Name, functions_library &functions, string &folder_Path, int *parents_in_Tissue, int &num_Sequences, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation, vector<int> &parent_IDs, int &last_index_Seq_Written, mt19937 &gen)
{

    cout << "Collecting " << num_Sequences << " sequence(s)\n";
    // string folder_Path = source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation);

    int num_per_Core = num_Sequences / this->CPU_cores;
    int remainder = num_Sequences % this->CPU_cores;

    vector<thread> threads_vec;

    for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    {
        int start_Cell = core_ID * num_per_Core;
        int stop_Cell = start_Cell + num_per_Core;

        threads_vec.push_back(thread{&cancer_Host::thread_find_Files, this, offset, start_Cell, stop_Cell, parents_in_Tissue, ref(indexed_Tissue_Folder)});
    }

    if (remainder != 0)
    {
        int start_Cell = num_Sequences - remainder;
        int stop_Cell = num_Sequences;

        threads_vec.push_back(thread{&cancer_Host::thread_find_Files, this, offset, start_Cell, stop_Cell, parents_in_Tissue, ref(indexed_Tissue_Folder)});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    vector<int> Tissue_files(found_Tissue_Folder_Indexes.begin(), found_Tissue_Folder_Indexes.end());
    found_Tissue_Folder_Indexes.clear();

    cout << Tissue_files.size() << " file(s) identified\n";

    // exit(-1);

    // vector<pair<int, int>> sequence_FileIndex_Position_list;
    vector<string> collected_Sequences;
    // for (int index = 0; index < num_Sequences; index++)
    // {
    //     sequence_FileIndex_Position_list.push_back(make_pair(parents_in_Tissue[index], index));
    //     collected_Sequences.push_back("");
    // }

    // sort(sequence_FileIndex_Position_list.begin(), sequence_FileIndex_Position_list.end());

    fstream nfasta;
    int index_Files = 0;
    int line_current = 0;

    cout << "Retrieving sequence(s)\n";

    // valid_Sequences = 0;

    nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].first) + "_" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].second) + ".nfasta", ios::in);

    uniform_real_distribution<float> check_Replicate_Dis(0.0, 1.0);

    fstream sequence_Profiles_File;
    sequence_Profiles_File.open(sequence_Profiles, ios::app);
    fstream sequence_parent_Progeny_relationships_File;
    sequence_parent_Progeny_relationships_File.open(sequence_parent_Progeny_relationships, ios::app);

    for (int find = 0; find < num_Sequences; find++)
    {
        // cout << "Looking for " << sequence_FileIndex_Position_list[find].first << "\n";

        while ((indexed_Tissue_Folder[Tissue_files[index_Files]].first <= parents_in_Tissue[find + offset] && indexed_Tissue_Folder[Tissue_files[index_Files]].second >= parents_in_Tissue[find + offset]) == 0)
        {
            nfasta.close();
            index_Files++;
            nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].first) + "_" + to_string(indexed_Tissue_Folder[Tissue_files[index_Files]].second) + ".nfasta", ios::in);
            line_current = 0;
        }

        if (nfasta.is_open())
        {
            int line_t0_check = (parents_in_Tissue[find + offset] - indexed_Tissue_Folder[Tissue_files[index_Files]].first) * 2;

            string line;
            string sequence = "";

            while (getline(nfasta, line))
            {
                if (line_t0_check == line_current)
                {
                    // cout << line << endl;
                    vector<string> line_Data;
                    functions.split(line_Data, line, '_');

                    int check_replicate = 0;

                    if (stof(line_Data[2]) < 1)
                    {
                        float val = check_Replicate_Dis(gen);
                        check_replicate = (val < stof(line_Data[2])) ? 0 : 1;
                    }

                    if (stoi(line_Data[0].substr(1)) == parents_in_Tissue[find + offset] && check_replicate == 0)
                    {
                        getline(nfasta, line);
                        collected_Sequences.push_back(line);
                        parent_IDs.push_back(parents_in_Tissue[find + offset]);
                        line_current++;
                    }
                    else if (stoi(line_Data[0].substr(1)) == parents_in_Tissue[find + offset] && check_replicate == 1)
                    {
                        // check if it survives to the next generation
                        int check_Survival = 0;
                        if (stof(line_Data[3]) < 1)
                        {
                            float val = check_Replicate_Dis(gen);
                            check_Survival = (val < stof(line_Data[3])) ? 0 : 1;
                        }
                        if (check_Survival == 0)
                        {
                            getline(nfasta, line);
                            to_write_Sequence_Store.push_back(make_pair(last_index_Seq_Written + "_A_" + line_Data[2] + "_" + line_Data[3], line));
                            sequence_Profiles_File << tissue_Name << "_" << to_string(current_Generation + 1) << "_" << last_index_Seq_Written << "\t" << tissue_Name << endl;
                            sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(current_Generation) << "_" << parents_in_Tissue[find + offset] << "\t"
                                                                       << tissue_Name << "_" << to_string(current_Generation + 1) << "_" << last_index_Seq_Written
                                                                       << "\tgeneration_Forward" << endl;

                            current_cell_load_per_Tissue[tissue] = current_cell_load_per_Tissue[tissue] + 1;
                            last_index_Seq_Written++;

                            line_current++;
                        }
                    }
                    else
                    {
                        cout << "ERROR: CORRECT SEQUENCE NOT FOUND AT INDEX\n";
                        cout << "Looking for: " << parents_in_Tissue[find + offset] << endl
                             << "Sequence ID at location: " << line << endl
                             << "File: " << folder_Path << "/" << indexed_Tissue_Folder[Tissue_files[index_Files]].first
                             << "_" << indexed_Tissue_Folder[Tissue_files[index_Files]].second << ".nfasta" << endl;
                        exit(-1);
                    }
                    line_current++;
                    break;
                }
                line_current++;
            }
        }
        else
        {
            cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << folder_Path << "/" << indexed_Tissue_Folder[Tissue_files[index_Files]].first << "_" << indexed_Tissue_Folder[Tissue_files[index_Files]].second << ".nfasta" << endl;
            exit(-1);
        }
    }

    nfasta.close();

    sequence_Profiles_File.close();
    sequence_parent_Progeny_relationships_File.close();

    // cout << valid_Sequences << " live sequence(s) collected\n";

    return collected_Sequences;
}

void cancer_Host::thread_find_Files(int offset, int start, int stop, int *parents_in_Tissue, vector<pair<int, int>> &indexed_Tissue_Folder)
{
    vector<int> caught_Indexes;
    for (int sequence = start; sequence < stop; sequence++)
    {
        int top = 0;
        int bottom = indexed_Tissue_Folder.size() - 1;
        int middle = top + ((bottom - top) / 2);

        while (top <= bottom)
        {
            if ((indexed_Tissue_Folder[middle].first <= parents_in_Tissue[sequence + offset]) && (indexed_Tissue_Folder[middle].second >= parents_in_Tissue[sequence + offset]))
            {
                caught_Indexes.push_back(middle);
                // found = 'Y';
                break;
            }
            else if (indexed_Tissue_Folder[middle].first < parents_in_Tissue[sequence + offset])
            {
                top = middle + 1;
            }
            else
            {
                bottom = middle - 1;
            }
            middle = top + ((bottom - top) / 2);
        }
    }

    // int index = 0;
    unique_lock<shared_mutex> ul(g_mutex);
    for (int index = 0; index < caught_Indexes.size(); index++)
    {
        found_Tissue_Folder_Indexes.insert(caught_Indexes[index]);
        //   index++;
    }
}

vector<pair<int, int>> cancer_Host::get_Rounds(int &total_Count, int &gpu_Max_Limit)
{
    vector<pair<int, int>> Rounds_start_stop;

    int full_Rounds = total_Count / gpu_Max_Limit;
    int partial_Rounds = total_Count % gpu_Max_Limit;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * gpu_Max_Limit;
        int stop = start + gpu_Max_Limit;
        Rounds_start_stop.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = total_Count - partial_Rounds;
        Rounds_start_stop.push_back(make_pair(start, total_Count));
    }

    return Rounds_start_stop;
}

string cancer_Host::get_generation_Phase(int &overall_Generations,
                                         vector<float> time_Generation,
                                         vector<string> phase_Type,
                                         vector<pair<float, float>> phase_paramaters,
                                         float &variable_1, float &variable_2)
{
    cout << "\nGetting current tissue generation phase\n";

    int index = -1;
    if (time_Generation[0] > overall_Generations)
    {
        index = 0;
    }
    else if (time_Generation[time_Generation.size() - 1] <= overall_Generations)
    {
        index = time_Generation.size() - 1;
    }
    else
    {
        for (int get_Index = 1; get_Index < time_Generation.size(); get_Index++)
        {
            if (overall_Generations >= time_Generation[get_Index - 1] && overall_Generations < time_Generation[get_Index])
            {
                index = get_Index;
                break;
            }
        }
    }

    if (index != -1)
    {

        cout << "Phase: " << phase_Type[index] << endl;

        variable_1 = phase_paramaters[index].first;
        variable_2 = phase_paramaters[index].second;

        return phase_Type[index];
    }
    else
    {
        cout << "\nERROR: Generation not found\n";
        exit(-1);
    }
}

int cancer_Host::terminal_status(int &num_tissues, int *tissue_array, string &source_sequence_Data_folder,
                                 string &enable_Folder_management, string &enable_Compression, int &terminal_Load)
{

    if (get_Load(num_tissues, tissue_array) >= terminal_Load)
    {
        if (enable_Folder_management == "YES")
        {
            compress_Folder(source_sequence_Data_folder, enable_Compression);
        }
        return 1;
    }
    else
    {
        return 0;
    }
}

void cancer_Host::compress_Folder(string path, string &enable_Compression)
{
    if (filesystem::exists(path))
    {
        cout << "\nCompressing folder: " << path << endl;

        string tar_Folder;
        string command_Tar;

        if (enable_Compression == "YES")
        {
            tar_Folder = path + ".tar.gz";
            command_Tar = "tar -czf " + tar_Folder + " " + path + " && rm -R " + path;
        }
        else
        {
            tar_Folder = path + ".tar";
            command_Tar = "tar -cf " + tar_Folder + " " + path + " && rm -R " + path;
        }

        int result = system(command_Tar.c_str());

        if (result == 0)
        {
            cout << "Compression successful" << endl;
        }
        else
        {
            cout << "Failed to compress the folder: " << path << endl;
            exit(-1);
        }
    }
    else
    {
        cout << "COMPRESSION ERROR: UNABLE TO FIND THE FOLDER: " << path << endl;
        exit(-1);
    }
}

int cancer_Host::get_Load(int &num_tissues_Calc, int *tissue_array)
{
    int sum = 0;

    for (int tissue = 0; tissue < num_tissues_Calc; tissue++)
    {
        sum = sum + (current_cell_load_per_Tissue[tissue_array[tissue]] - dead_Particle_count[tissue_array[tissue]] - removed_by_Transfer_Indexes[tissue_array[tissue]].size());
    }

    return sum;
}

void cancer_Host::initialize(functions_library &functions,
                             vector<string> &tissue_Names,
                             string &intermediary_Sequence_location, string &first_Infection,
                             int &current_Generation,
                             string &output_Node_location,
                             int &max_sequences_per_File)
{
    cout << "\nConfiguring cancer host\n";

    num_Tissues = tissue_Names.size();

    string host_Folder = intermediary_Sequence_location + "/cancer_Host";
    functions.config_Folder(host_Folder, "Cancer host sequence data");

    vector<vector<pair<string, string>>> tissue_Sequences;
    intialize_Tissues(host_Folder, tissue_Sequences, functions, current_Generation);

    string reference_Sequences = intermediary_Sequence_location + "/reference_Sequences";

    cout << "\nFirst infection mode: " << first_Infection << endl;

    if (first_Infection == "RANDOM")
    {
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

        vector<string> line_Data;

        cout << endl;
        // vector<char> seq_Status;
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
                        functions.split(line_Data, line, '_');
                        Sequence_IDs.push_back(line_Data[2] + "_" + line_Data[3]);
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
                nfasta.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << files[file] << endl;
                exit(-1);
            }
        }

        random_device rd; // Will be used to obtain a seed for the random number engine
        mt19937 gen(rd());
        uniform_int_distribution<int> entry_Tissue_select(0, tissue_Names.size() - 1);

        cout << endl;

        for (int sequence = 0; sequence < Sequences.size(); sequence++)
        {
            int tissue_Index = entry_Tissue_select(gen);
            cout << "Sequence " << sequence + 1 << " infects " << tissue_Names[tissue_Index] << " tissue of index: " << tissue_Index << endl;
            tissue_Sequences[tissue_Index].push_back(make_pair(Sequences[sequence], Sequence_IDs[sequence]));
        }

        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            if (tissue_Sequences[tissue].size() > 0)
            {
                current_cell_load_per_Tissue[tissue] = tissue_Sequences[tissue].size();
                functions.config_Folder(host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), "Tissue " + tissue_Names[tissue] + " Generation " + to_string(current_Generation));

                if (!filesystem::exists(output_Node_location + "/cancer_Host"))
                {
                    functions.config_Folder(output_Node_location + "/cancer_Host", "Cancer host node");
                    functions.create_File(output_Node_location + "/cancer_Host/sequence_Profiles.csv", "Sequence_ID\tTissue");
                    functions.create_File(output_Node_location + "/cancer_Host/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                }

                vector<pair<string, string>> sequence_Write_Store_All;
                int last_seq_Num = 0;
                sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[tissue],
                                            max_sequences_per_File, host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), last_seq_Num,
                                            output_Node_location + "/cancer_Host/sequence_Profiles.csv", tissue_Names[tissue], current_Generation);
                partial_Write_Check(sequence_Write_Store_All,
                                    host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), last_seq_Num,
                                    output_Node_location + "/cancer_Host/sequence_Profiles.csv", tissue_Names[tissue], current_Generation);
            }
        }
    }
    else
    {
        cout << "\nInfecting tissues\n";
        vector<string> line_Data;

        //// Write to sequence profiles file

        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            string reference_Sequences_tissue = reference_Sequences + "/" + tissue_Names[tissue];
            int last_seq_Num = 0;

            if (filesystem::exists(reference_Sequences_tissue) && filesystem::is_directory(reference_Sequences_tissue))
            {
                cout << "\nTissue: " << tissue_Names[tissue] << endl;
                for (const auto &entry : filesystem::directory_iterator(reference_Sequences_tissue))
                {
                    if (entry.is_regular_file() && entry.path().extension() == ".nfasta")
                    {
                        string file_Name = entry.path().stem();
                        functions.split(line_Data, file_Name, '_');

                        int num_Particles_Tissue = stoi(line_Data[1]) - stoi(line_Data[0]) + 1;

                        cout << "Sequences migrating: " << num_Particles_Tissue << endl;

                        current_cell_load_per_Tissue[tissue] = current_cell_load_per_Tissue[tissue] + num_Particles_Tissue;

                        functions.config_Folder(host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), "Tissue " + to_string(tissue) + " Generation " + to_string(current_Generation));

                        if (!filesystem::exists(output_Node_location + "/cancer_Host"))
                        {
                            functions.config_Folder(output_Node_location + "/cancer_Host", "Cancer host node");
                            functions.create_File(output_Node_location + "/cancer_Host" + "/sequence_Profiles.csv", "Sequence_ID\tTissue");
                            functions.create_File(output_Node_location + "/cancer_Host" + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                        }

                        filesystem::copy_file(entry.path().string(), host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation) + "/" + file_Name + ".nfasta");

                        fstream sequence_Profile;
                        sequence_Profile.open(output_Node_location + "/cancer_Host" + "/sequence_Profiles.csv", ios::app);

                        for (int sequence_Num = 0; sequence_Num < num_Particles_Tissue; sequence_Num++)
                        {
                            sequence_Profile << tissue_Names[tissue] << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue_Names[tissue] << endl;
                            last_seq_Num++;
                        }

                        sequence_Profile.close();
                    }
                }
            }
        }
    }
}

void cancer_Host::sequence_Write_Configurator(vector<pair<string, string>> &sequence_Write_Store_All, vector<pair<string, string>> &sequence_Write_Store,
                                              int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                              string sequence_Profiles_Location, string tissue, int current_Generation)
{
    for (int sequence_Collect = 0; sequence_Collect < sequence_Write_Store.size(); sequence_Collect++)
    {
        sequence_Write_Store_All.push_back(make_pair(sequence_Write_Store[sequence_Collect].first, sequence_Write_Store[sequence_Collect].second));
    }

    sequence_Write_Store.clear();

    if (sequence_Write_Store_All.size() >= max_sequences_per_File)
    {
        int full_Write_Count = sequence_Write_Store_All.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {

            string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);
            fstream sequence_Profile;
            //"Sequence_ID\tHost\tTissue"
            sequence_Profile.open(sequence_Profiles_Location, ios::app);

            if (fasta_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_seq_Num << "_A_";
                    // if (seq_Status.size() == 0)
                    // {
                    //     fasta_File << "_A";
                    // }
                    // else
                    // {
                    //     fasta_File << "_" << seq_Status[write_Seq];
                    // }
                    fasta_File << sequence_Write_Store_All[write_Seq].second;
                    fasta_File << endl;
                    fasta_File << sequence_Write_Store_All[write_Seq].first << endl;
                    sequence_Profile << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue << endl;
                    last_seq_Num++;
                }

                fasta_File.close();
                sequence_Profile.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                exit(-1);
            }
            // last_seq_Num = last_seq_Num + max_sequences_per_File;
        }

        // int parital_Write_Count = sequence_Write_Store_All.size() % max_sequences_per_File;
        vector<pair<string, string>> sequence_Write_Store_temp;
        vector<char> temp_seq_Status;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < sequence_Write_Store_All.size(); fill++)
        {
            sequence_Write_Store_temp.push_back(make_pair(sequence_Write_Store_All[fill].first, sequence_Write_Store_All[fill].second));
        }

        sequence_Write_Store_All.clear();
        sequence_Write_Store_All = sequence_Write_Store_temp;
    }
}

void cancer_Host::partial_Write_Check(vector<pair<string, string>> &sequence_Write_Store_All,
                                      const string &folder_Location, int &last_seq_Num,
                                      string sequence_Profiles_Location, string tissue, int current_Generation)
{
    if (sequence_Write_Store_All.size() > 0)
    {
        string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + sequence_Write_Store_All.size() - 1) + ".nfasta";
        fstream fasta_File;
        fasta_File.open(fasta_file_Location, ios::out);
        fstream sequence_Profile;
        sequence_Profile.open(sequence_Profiles_Location, ios::app);
        if (fasta_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < sequence_Write_Store_All.size(); write_Seq++)
            {
                fasta_File << ">" << last_seq_Num << "_A_";
                // if (seq_Status.size() == 0)
                // {
                //     fasta_File << "_A";
                // }
                // else
                // {
                //     fasta_File << "_" << seq_Status[write_Seq];
                // }
                fasta_File << sequence_Write_Store_All[write_Seq].second;
                fasta_File << endl;
                fasta_File << sequence_Write_Store_All[write_Seq].first << endl;
                sequence_Profile << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue << endl;
                last_seq_Num++;
            }

            fasta_File.close();
            sequence_Profile.close();
        }
        else
        {
            cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            exit(-1);
        }
        sequence_Write_Store_All.clear();
    }
}

void cancer_Host::intialize_Tissues(string &host_Folder, vector<vector<pair<string, string>>> &tissue_Sequences, functions_library &functions, int &current_Generation)
{
    current_cell_load_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
    dead_Particle_count = (int *)malloc(sizeof(int) * num_Tissues);
    parents_Prev_generation = (int *)malloc(sizeof(int) * num_Tissues);
    current_Generation = 0;

    set<int> init_removed_by_Transfer_Indexes;

    for (int tissue = 0; tissue < num_Tissues; tissue++)
    {
        vector<pair<string, string>> tissue_Sequence;
        functions.config_Folder(host_Folder + "/" + to_string(tissue), "Tissue " + to_string(tissue + 1));

        tissue_Sequences.push_back(tissue_Sequence);
        current_cell_load_per_Tissue[tissue] = 0;
        dead_Particle_count[tissue] = 0;

        parents_Prev_generation[tissue] = 0;

        removed_by_Transfer_Indexes.push_back(init_removed_by_Transfer_Indexes);
    }
}