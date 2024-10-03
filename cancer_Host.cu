#include "cancer_Host.cuh"

cancer_Host::cancer_Host()
{
    cout << "\nSTEP 5: Cancer host intialization\n";
}

void cancer_Host::cell_Migration_set(int &max_Limit, multiset<pair<float, int>> &migration_cell_List, pair<float, int> candidate_Cell)
{
    cout << "\nMigration set\n";
    migration_cell_List.insert(candidate_Cell);

    if (migration_cell_List.size() > max_Limit)
    {
        migration_cell_List.erase(prev(migration_cell_List.end()));
    }
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
                                       string &multi_Read, int &CPU_cores, int &num_Cuda_devices, int *CUDA_device_IDs,
                                       int &genome_Length,
                                       float *Reference_fitness_survivability_proof_reading, float *Reference_cancer_parameters,
                                       float **A_0_mutation,
                                       float **T_1_mutation,
                                       float **G_2_mutation,
                                       float **C_3_mutation,
                                       int &mutation_Hotspots,
                                       float **mutation_hotspot_parameters,
                                       int *num_effect_Segregating_sites,
                                       float **sequence_Survivability_changes,
                                       float **sequence_Proof_reading_changes,
                                       int *num_effect_Segregating_sites_Cancer,
                                       float **sequence_replication_factor_changes,
                                       float **sequence_mutation_rate_changes,
                                       float **sequence_generation_death_changes,
                                       float **sequence_replication_prob_changes,
                                       float **sequence_metastatic_prob_changes,
                                       int &max_sequences_per_File,
                                       string &viral_Migration, float **viral_Migration_Values, int *migration_start_Generation,
                                       int &count_tajima_Regions, int **tajima_regions_Start_Stop,
                                       string &reference_Genome_location,
                                       int *tissue_selection_Position_Count,
                                       int *Survivability_Positions,
                                       int *Proof_Positions,
                                       int *Replication_factor_Positions,
                                       int *Mutation_rate_factor_Positions,
                                       int *Generation_death_Positions,
                                       int *Replication_prob_Positions,
                                       int *Metastatic_Positions,
                                       float **tissues_ATGC_positions_Survivability,
                                       float **tissues_ATGC_positions_Proof,
                                       float **tissues_ATGC_positions_Replication_factor,
                                       float **tissues_ATGC_positions_Mutation_rate_factor,
                                       float **tissues_ATGC_positions_Generation_death,
                                       float **tissues_ATGC_positions_Replication_prob,
                                       float **tissues_ATGC_positions_Metastatic,
                                       int *profile_tissue_Limits)
{
    cout << "\nSTEP 6: Conducting simulation\n";

    // this->terminal_Load = terminal_Load;

    random_device rd;
    mt19937 gen(rd());

    this->CPU_cores = CPU_cores;
    this->genome_Length = genome_Length;

    generational_Summary = output_Node_location + "/cancer_Host/node_generational_Summary.csv";
    functions.create_File(generational_Summary, "Generation\tTissue\tPhase\tnum_Parents\tnum_Progeny\tdead_Progeny\trapid_Progeny");

    sequence_Profiles = output_Node_location + "/cancer_Host/sequence_Profiles.csv";
    // functions.create_File(sequence_Profiles, "Sequence_ID\tTissue");

    sequence_parent_Progeny_relationships = output_Node_location + "/cancer_Host/sequence_parent_Progeny_relationships.csv";

    string output_Tajima_File = output_Node_location + "/cancer_Host/tajimas_D_time_series.csv";
    if (count_tajima_Regions > 0)
    {
        // string columns_Tajima = "";
        // for (int region = 0; region < count_tajima_Regions; region++)
        // {
        //     columns_Tajima = columns_Tajima + "\tRegion_" + to_string(region + 1);
        // }
        functions.create_File(output_Tajima_File, "Type\tGeneration\tTissue\tN_cells\tseg_Sites\tpairwise_Diff\tCombinations\tPi\ttajima_D");
    }
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

        vector<int> tissue_Migration_Totals;
        vector<vector<pair<int, int>>> tissue_migration_Targets_amount;

        vector<vector<pair<int, int>>> indexed_Source_Folders = functions.index_sequence_Folders(source_sequence_Data_folder, num_Tissues, overall_Generations, multi_Read);

        for (int tissue = 0; tissue < num_Tissues; tissue++)
        {

            if (indexed_Source_Folders[tissue].size() > 0)
            {
                current_cell_load_per_Tissue[tissue] = indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1;
                cout << tissue << " Check: " << current_cell_load_per_Tissue[tissue] << endl;
            }

            real_Particle_count_per_Tissue[tissue] = current_cell_load_per_Tissue[tissue] - removed_by_Transfer_Indexes[tissue].size() - dead_Particle_count[tissue];
            cout << tissue_Names[tissue] << " tissue: " << real_Particle_count_per_Tissue[tissue] << endl;
            sum_Check = sum_Check + real_Particle_count_per_Tissue[tissue];

            tissue_Migration_Totals.push_back(0);
            vector<pair<int, int>> intialize_vec;
            tissue_migration_Targets_amount.push_back(intialize_vec);
        }

        if (viral_Migration == "YES")
        {
            cout << "\nPutative migrating particles per tissue: \n";

            for (int migration_Check = 0; migration_Check < (num_Tissues * (num_Tissues - 1)); migration_Check++)
            {
                if (viral_Migration_Values[migration_Check][0] != -1)
                {
                    if (overall_Generations >= migration_start_Generation[migration_Check])
                    {
                        int source = migration_Check / (num_Tissues - 1);
                        int destination = migration_Check % (num_Tissues - 1);

                        if (destination >= source)
                        {
                            destination = destination + 1;
                        }

                        binomial_distribution<int> num_Particles(viral_Migration_Values[migration_Check][0], viral_Migration_Values[migration_Check][1]);
                        int num_viruses_to_transfer = num_Particles(gen);

                        tissue_Migration_Totals[source] = tissue_Migration_Totals[source] + num_viruses_to_transfer;
                        tissue_migration_Targets_amount[source].push_back(make_pair(destination, num_viruses_to_transfer));
                    }
                }
            }

            for (int tissue = 0; tissue < num_Tissues; tissue++)
            {
                cout << "Tissue " << tissue_Names[tissue] << ": " << tissue_Migration_Totals[tissue] << endl;

                for (int path = 0; path < tissue_migration_Targets_amount[tissue].size(); path++)
                {
                    cout << tissue_migration_Targets_amount[tissue][path].second << " cells to " << tissue_Names[tissue_migration_Targets_amount[tissue][path].first] << endl;
                }
            }
        }

        // exit(-1);

        if (sum_Check > 0)
        {
            if (terminal_status(terminal_tissues, terminal_array, source_sequence_Data_folder, enable_Folder_management, enable_Compression, terminal_Load) != 1)
            {
                // vector<vector<pair<int, int>>> indexed_Source_Folders = functions.index_sequence_Folders(source_sequence_Data_folder, num_Tissues, overall_Generations, multi_Read);

                for (int tissue = 0; tissue < num_Tissues; tissue++)
                {
                    cout << "Tissue: " << tissue << endl;
                    if (real_Particle_count_per_Tissue[tissue] > 0)
                    {
                        // get last wrtten progeny number
                        cout << endl;
                        int last_Progeny_written_this_Gen = indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1;
                        string rapid_Progeny_Location = source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) + "/" + to_string(last_Progeny_written_this_Gen) + "_rapid_Progeny.nfasta";

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

                        // cout << "\nFilling parent vector: ";
                        // int *parents_in_Tissue = (int *)malloc(real_Particle_count_per_Tissue[tissue] * sizeof(int));

                        // int fill_Count = 0;
                        // int index = 0;
                        // do
                        // {
                        //     if (check_to_Remove.find(index) == check_to_Remove.end())
                        //     {
                        //         parents_in_Tissue[fill_Count] = index;
                        //         fill_Count++;
                        //     }
                        //     index++;
                        // } while (fill_Count < real_Particle_count_per_Tissue[tissue]);

                        // cout << "Completed\n";

                        // cout << "Shuffling parent vector: ";
                        // default_random_engine rng(time(nullptr)); // Seed the random number generator with current time
                        // shuffle(parents_in_Tissue, parents_in_Tissue + real_Particle_count_per_Tissue[tissue], rng);
                        // cout << "Complete\n";

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

                        if (profile_tissue_Limits[tissue] != -1 && parent_population_Count > profile_tissue_Limits[tissue])
                        {
                            cout << "\nMax cell count in tissue reached\n";
                            parent_population_Count = profile_tissue_Limits[tissue];
                        }

                        cout << "\nFilling parent set\n";

                        set<int> parents_to_get;
                        uniform_int_distribution<> distr(0, last_Progeny_written_this_Gen - 1);

                        while (parents_to_get.size() < parent_population_Count)
                        {
                            int random_parent = distr(gen);
                            if (check_to_Remove.find(random_parent) == check_to_Remove.end())
                            {
                                parents_to_get.insert(random_parent);
                            }
                        }
                        cout << "\nConverting parent vetor\n";
                        vector<int> parents_in_Tissue(parents_to_get.begin(), parents_to_get.end());
                        parents_to_get.clear();

                        check_to_Remove.clear();
                        removed_by_Transfer_Indexes[tissue].clear();
                        dead_Particle_count[tissue] = 0;
                        current_cell_load_per_Tissue[tissue] = 0;

                        parents_Prev_generation[tissue] = parent_population_Count;

                        int last_index_Seq_Written = 0;
                        // ! check if the new generation already exists and if so update;
                        string intermediary_Tissue_folder = source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations + 1);
                        string dead_List = intermediary_Tissue_folder + "/dead_List.txt";

                        // fstream this_Gen_progeny_parents;
                        // functions.create_File(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) + "/" + to_string(last_Progeny_written_this_Gen) + "_rapid_Progeny.nfasta");

                        if (filesystem::exists(intermediary_Tissue_folder) && filesystem::is_directory(intermediary_Tissue_folder))
                        {
                            cout << "Next generation already present\n";
                            vector<pair<int, int>> indexed_tissue_Folder = functions.index_sequence_Folder(intermediary_Tissue_folder);
                            if (indexed_tissue_Folder.size() > 0)
                            {
                                last_index_Seq_Written = indexed_tissue_Folder[indexed_tissue_Folder.size() - 1].second + 1;
                                current_cell_load_per_Tissue[tissue] = last_index_Seq_Written;

                                fstream dead_File;
                                dead_File.open(dead_List, ios::in);
                                if (dead_File.is_open())
                                {
                                    string line;
                                    int index = 0;
                                    while (getline(dead_File, line))
                                    {
                                        // check_to_Remove.insert(stoi(line));
                                        index++;
                                    }
                                    dead_Particle_count[tissue] = index;
                                    dead_File.close();
                                }
                                else
                                {
                                    cout << "ERROR: UNABLE TO OPEN DEAD LIST FILE: " << dead_List << endl;
                                    exit(-1);
                                }
                            }
                        }
                        else
                        {
                            functions.config_Folder(intermediary_Tissue_folder, to_string(overall_Generations + 1) + " generation Tissue " + tissue_Names[tissue] + " sequences");
                            functions.create_File(dead_List);
                        }

                        vector<pair<int, int>> cells_Rounds_start_stop = get_Rounds(parent_population_Count, max_Cells_at_a_time);

                        multiset<pair<float, int>> migration_cell_List;

                        for (int cell_Round = 0; cell_Round < cells_Rounds_start_stop.size(); cell_Round++)
                        {
                            int num_of_Cells = cells_Rounds_start_stop[cell_Round].second - cells_Rounds_start_stop[cell_Round].first;
                            cout << "\nProcessing round " << cell_Round + 1 << " of " << cells_Rounds_start_stop.size() << ": " << num_of_Cells << " cell(s)" << endl;

                            simulate_cell_Round(functions, multi_Read, num_Cuda_devices, CUDA_device_IDs,
                                                num_of_Cells, cells_Rounds_start_stop[cell_Round].first, cells_Rounds_start_stop[cell_Round].second,
                                                parents_in_Tissue, tissue, tissue_Names[tissue],
                                                indexed_Source_Folders[tissue],
                                                source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations),
                                                overall_Generations,
                                                last_index_Seq_Written,
                                                gen,
                                                Reference_fitness_survivability_proof_reading, Reference_cancer_parameters,
                                                A_0_mutation,
                                                T_1_mutation,
                                                G_2_mutation,
                                                C_3_mutation,
                                                mutation_Hotspots,
                                                mutation_hotspot_parameters,
                                                num_effect_Segregating_sites,
                                                sequence_Survivability_changes,
                                                sequence_Proof_reading_changes,
                                                num_effect_Segregating_sites_Cancer,
                                                sequence_replication_factor_changes,
                                                sequence_mutation_rate_changes,
                                                sequence_generation_death_changes,
                                                sequence_replication_prob_changes,
                                                sequence_metastatic_prob_changes,
                                                max_sequences_per_File, intermediary_Tissue_folder, source_sequence_Data_folder,
                                                last_Progeny_written_this_Gen, rapid_Progeny_Location,
                                                tissue_Migration_Totals[tissue], migration_cell_List,
                                                tissue_selection_Position_Count,
                                                Survivability_Positions,
                                                Proof_Positions,
                                                Replication_factor_Positions,
                                                Mutation_rate_factor_Positions,
                                                Generation_death_Positions,
                                                Replication_prob_Positions,
                                                Metastatic_Positions,
                                                tissues_ATGC_positions_Survivability,
                                                tissues_ATGC_positions_Proof,
                                                tissues_ATGC_positions_Replication_factor,
                                                tissues_ATGC_positions_Mutation_rate_factor,
                                                tissues_ATGC_positions_Generation_death,
                                                tissues_ATGC_positions_Replication_prob,
                                                tissues_ATGC_positions_Metastatic,
                                                viral_Migration);
                        }

                        // last_Progeny_written_this_Gen = indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1;

                        fstream gen_Summary;
                        gen_Summary.open(generational_Summary, ios::app);
                        if (gen_Summary.is_open())
                        {
                            //"Generation\tTissue\tPhase\tnum_Parents\tnum_Progeny\tdead_Progeny\trapid_Progeny"
                            gen_Summary << overall_Generations << "\t" << tissue_Names[tissue] << "\t" << generation_Type << "\t" << to_string(parent_population_Count)
                                        << "\t" << last_index_Seq_Written << "\t" << dead_Particle_count[tissue] << "\t";
                        }
                        else
                        {
                            cout << "ERROR: UNABLE TO OPEN GENERATIONAL SUMMARY FILE: " << generational_Summary << endl;
                            exit(-1);
                        }

                        string rename_rapid_Progeny = source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) +
                                                      "/" + to_string(indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1) + "_rapid_Progeny.nfasta";

                        if (filesystem::exists(rename_rapid_Progeny))
                        {
                            cout << "Renaming rapid progeny file: " << rename_rapid_Progeny << endl;
                            try
                            {
                                filesystem::rename(rename_rapid_Progeny, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) +
                                                                             "/" + to_string(indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1) + "_" + to_string(last_Progeny_written_this_Gen - 1) + ".nfasta");
                                gen_Summary << to_string(last_Progeny_written_this_Gen - (indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1));
                            }
                            catch (const filesystem::filesystem_error &e)
                            {
                                // Handle any errors that occur during renaming
                                std::cerr << "Error renaming file: " << e.what() << '\n';
                            }
                        }
                        else
                        {
                            gen_Summary << "0";
                        }

                        gen_Summary << endl;
                        gen_Summary.close();

                        remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
                        // remainder_Write_Sequences_NEXT_Generation(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

                        for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue].size(); forward++)
                        {
                            remainder_Write_Sequences_NEXT_Generation(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue][forward]);
                        }

                        // free(parents_in_Tissue);

                        // See which progeny qualify to migrate

                        cout << "\nBottleneck size for metastatsis: " << tissue_Migration_Totals[tissue] << endl;
                        cout << "Metastatic cells avaiable: " << migration_cell_List.size() << endl;

                        if (viral_Migration == "YES")
                        {
                            migration_of_Cells(source_sequence_Data_folder, tissue_Names,
                                               tissue, tissue_migration_Targets_amount[tissue], migration_cell_List,
                                               overall_Generations, functions);
                        }
                        // // ! Calculate Tajima's. Define parameters with gene regions for Tajima's

                        if (count_tajima_Regions > 0)
                        {
                            calculate_Tajima(functions,
                                             count_tajima_Regions, tajima_regions_Start_Stop,
                                             overall_Generations, tissue_Names[tissue], tissue,
                                             CUDA_device_IDs,
                                             source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations), max_Cells_at_a_time,
                                             output_Tajima_File,
                                             reference_Genome_location);
                        }

                        // exit(-1);
                        cout << "\nCompleted tissue: " << tissue_Names[tissue] << endl;
                    }
                    // overall_Generations++;
                }
            }
            else
            {
                stop_Type = 5;
            }
        }
        else
        {
            stop_Type = 4;
        }

        decimal_Date = decimal_Date + date_Increment;
        overall_Generations++;

        free(real_Particle_count_per_Tissue);

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
        // stop_Type = 1;

    } while (stop_Type == 0);

    cout << "\nSimulation has concluded: ";
}

__global__ void cuda_convert_Sequence(int genome_Length, char *sequence)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Length)
    {

        if (sequence[tid] == 'A' || sequence[tid] == 'a' || sequence[tid] == '0')
        {
            sequence[tid] = '0';
        }
        else if (sequence[tid] == 'T' || sequence[tid] == 't' || sequence[tid] == '1')
        {
            sequence[tid] = '1';
        }
        else if (sequence[tid] == 'G' || sequence[tid] == 'g' || sequence[tid] == '2')
        {
            sequence[tid] = '2';
        }
        else if (sequence[tid] == 'C' || sequence[tid] == 'c' || sequence[tid] == '3')
        {
            sequence[tid] = '3';
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void addToVariable(float *cuda_a_1, int N)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < (N + 1))
    {
        if (tid > 0)
        {
            atomicAdd(cuda_a_1, (float)1.0 / (float)(tid));
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void squared_addToVariable(float *cuda_a_2, int N)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < (N + 1))
    {
        if (tid > 0)
        {
            atomicAdd(cuda_a_2, (float)1.0 / (float)(tid * tid));
        }

        tid += blockDim.x * gridDim.x;
    }
}

void cancer_Host::calculate_Tajima(functions_library &functions,
                                   int &num_Regions, int **tajima_regions_Start_Stop,
                                   int &overall_Generations, string &tissue_Name, int &tissue_Index,
                                   int *CUDA_device_IDs,
                                   string sequence_Tissue_Folder, int &max_Cells_at_a_time,
                                   string output_Tajima_File,
                                   string reference_Genome_location)
{
    cout << "\nCalculating Tajima's for generation " << overall_Generations << " for tissue " << tissue_Name << endl;

    cout << "Reading reference genome:\n";
    fstream reference_Genome_file;
    reference_Genome_file.open(reference_Genome_location, ios::in);
    string reference_Genome = "";

    if (reference_Genome_file.is_open())
    {
        string line;
        while (getline(reference_Genome_file, line))
        {
            if (line.at(0) != '>')
            {
                reference_Genome.append(line);
            }
        }
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN REFERENCE GENOME FILE: " << reference_Genome_location << endl;
        exit(-1);
    }

    if (reference_Genome.size() == genome_Length)
    {
        cout << "Reference genome valid\n";
    }
    else
    {
        cout << "ERROR: Size of reference genome (" << reference_Genome.size() << ") is not equal to the genome length (" << genome_Length << ").";
        exit(-1);
    }

    cout << "Converting Reference genome to INT\n";
    // exit(-1);
    cudaSetDevice(CUDA_device_IDs[0]);

    char *full_Char;
    full_Char = (char *)malloc((reference_Genome.size() + 1) * sizeof(char));
    strcpy(full_Char, reference_Genome.c_str());

    char *cuda_full_Char;
    cudaMallocManaged(&cuda_full_Char, (reference_Genome.size() + 1) * sizeof(char));
    cudaMemcpy(cuda_full_Char, full_Char, (reference_Genome.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    free(full_Char);

    reference_Genome = "";
    cuda_convert_Sequence<<<functions.tot_Blocks_array[0], functions.tot_ThreadsperBlock_array[0]>>>(genome_Length, cuda_full_Char);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", 0, cudaGetErrorString(err));
        exit(-1);
    }

    cout << "Configuring common GPU memory\n";

    int **cuda_tajima_regions_Start_Stop;
    cudaMallocManaged(&cuda_tajima_regions_Start_Stop, num_Regions * sizeof(int *));
    for (int row = 0; row < num_Regions; row++)
    {
        cudaMalloc((void **)&(cuda_tajima_regions_Start_Stop[row]), 2 * sizeof(int));
        cudaMemcpy(cuda_tajima_regions_Start_Stop[row], tajima_regions_Start_Stop[row], 2 * sizeof(int), cudaMemcpyHostToDevice);
    }

    cout << "Reading folder: " << sequence_Tissue_Folder << endl;
    vector<pair<int, int>> indexed_Source_Folder = functions.index_sequence_Folder(sequence_Tissue_Folder);

    // cout << "\nCalculating pre-requisites: \n";

    int N_total = 0;
    int N_alive = 0;

    // exit(-1);

    string sequence_String = "";
    vector<string> line_Data;
    string dead_or_Alive = "";

    int count_Track = 0;

    int *per_Region = (int *)malloc(num_Regions * sizeof(int));
    int *per_Region_ALIVE = (int *)malloc(num_Regions * sizeof(int));

    for (int region = 0; region < num_Regions; region++)
    {
        per_Region[region] = 0;
        per_Region_ALIVE[region] = 0;
    }

    int *cuda_per_Region;
    cudaMallocManaged(&cuda_per_Region, num_Regions * sizeof(int));
    cudaMemcpy(cuda_per_Region, per_Region, num_Regions * sizeof(int), cudaMemcpyHostToDevice);

    int *cuda_per_Region_ALIVE;
    cudaMallocManaged(&cuda_per_Region_ALIVE, num_Regions * sizeof(int));
    cudaMemcpy(cuda_per_Region_ALIVE, per_Region_ALIVE, num_Regions * sizeof(int), cudaMemcpyHostToDevice);

    cout << "\nReading files\n";

    for (int file = 0; file < indexed_Source_Folder.size(); file++)
    {
        fstream nfasta;
        nfasta.open(sequence_Tissue_Folder + "/" + to_string(indexed_Source_Folder[file].first) + "_" + to_string(indexed_Source_Folder[file].second) + ".nfasta");
        if (nfasta.is_open())
        {
            string line;
            while (getline(nfasta, line))
            {
                if (line.at(0) != '>')
                {
                    sequence_String.append(line);
                    count_Track++;

                    if (count_Track == max_Cells_at_a_time)
                    {
                        process_Tajima_String(sequence_String, count_Track, num_Regions, cuda_tajima_regions_Start_Stop,
                                              cuda_full_Char, cuda_per_Region, functions, dead_or_Alive,
                                              cuda_per_Region_ALIVE);
                    }
                }
                else
                {
                    functions.split(line_Data, line, '_');
                    dead_or_Alive.append(line_Data[1]);
                    if (line_Data[1] == "A")
                    {
                        N_alive++;
                    }
                    N_total++;
                }
            }
            nfasta.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN SEQUENCE: "
                 << sequence_Tissue_Folder + "/" << indexed_Source_Folder[file].first + "_" << indexed_Source_Folder[file].second << ".nfasta\n";
            exit(-1);
        }
    }

    if (count_Track > 0)
    {
        process_Tajima_String(sequence_String, count_Track, num_Regions, cuda_tajima_regions_Start_Stop,
                              cuda_full_Char, cuda_per_Region, functions, dead_or_Alive, cuda_per_Region_ALIVE);
    }

    cout << "Completed reading files\n";

    cout << "Total number of cells (N): " << N_total << endl;

    cudaMemcpy(per_Region, cuda_per_Region, num_Regions * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(cuda_per_Region);

    cudaMemcpy(per_Region_ALIVE, cuda_per_Region_ALIVE, num_Regions * sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(cuda_per_Region_ALIVE);

    double b1;
    double b2;
    double c1;
    double c2;
    double e1;
    double e2;

    double out_a_1;

    double N_float = N_total;
    calc_pre_Requistes(b1,
                       b2,
                       c1,
                       c2,
                       e1,
                       e2,
                       out_a_1, N_total, functions);
    write_Tajima("ALL", per_Region, output_Tajima_File, overall_Generations, tissue_Name,
                 num_Regions, N_total, N_float, out_a_1, e1, e2);

    if (N_total != N_alive)
    {
        cout << "Processing ALIVE only\n";
        double N_float = N_alive;
        calc_pre_Requistes(b1,
                           b2,
                           c1,
                           c2,
                           e1,
                           e2,
                           out_a_1, N_alive, functions);
        write_Tajima("ALIVE", per_Region_ALIVE, output_Tajima_File, overall_Generations, tissue_Name,
                     num_Regions, N_alive, N_float, out_a_1, e1, e2);
    }

    // cout << "DONE\n";
    for (int row = 0; row < num_Regions; row++)
    {
        cudaFree(cuda_tajima_regions_Start_Stop[row]);
    }
    // cout << "DONE\n";
    cudaFree(cuda_tajima_regions_Start_Stop);
    cudaFree(cuda_full_Char);
    // cout << "DONE\n";
    free(per_Region);
    free(per_Region_ALIVE);
    // cout << "DONE\n";

    cout << "Completed Tajima calculation for generation: " << overall_Generations << endl;

    // exit(-1);
}

void cancer_Host::calc_pre_Requistes(double &b1,
                                     double &b2,
                                     double &c1,
                                     double &c2,
                                     double &e1,
                                     double &e2,
                                     double &out_a_1, int &N,
                                     functions_library &functions)
{
    cout << "\nCalculating pre-requistes\n";
    double N_float = N;

    if (N < 50000)
    {
        float a_1 = 0;
        float *cuda_a_1;

        // Allocate memory on the device
        cudaMalloc(&cuda_a_1, sizeof(float));

        // Copy the initial value from host to device
        cudaMemcpy(cuda_a_1, &a_1, sizeof(float), cudaMemcpyHostToDevice);

        addToVariable<<<functions.tot_Blocks_array[0], functions.tot_ThreadsperBlock_array[0]>>>(cuda_a_1, (N - 1));
        cudaDeviceSynchronize();

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", 0, cudaGetErrorString(err));
            exit(-1);
        }

        // Copy the result back to the host
        cudaMemcpy(&a_1, cuda_a_1, sizeof(float), cudaMemcpyDeviceToHost);
        // Free device memory
        cudaFree(cuda_a_1);

        cout << "a1: " << a_1 << endl;

        float a_2 = 0;
        float *cuda_a_2;

        cudaMalloc(&cuda_a_2, sizeof(float));
        cudaMemcpy(cuda_a_2, &a_2, sizeof(float), cudaMemcpyHostToDevice);

        squared_addToVariable<<<functions.tot_Blocks_array[0], functions.tot_ThreadsperBlock_array[0]>>>(cuda_a_2, (N - 1));
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", 0, cudaGetErrorString(err));
            exit(-1);
        }

        cudaMemcpy(&a_2, cuda_a_2, sizeof(float), cudaMemcpyDeviceToHost);
        cudaFree(cuda_a_2);

        cout << "a2: " << a_2 << endl;

        // double N_float = N;

        b1 = (N_float + 1.0) / (3.0 * (N_float - 1.0));
        b2 = (2.0 * ((N_float * N_float) + N_float + 3.0)) / (9.0 * N_float * (N_float - 1.0));
        c1 = b1 - (1.0 / a_1);
        c2 = b2 - ((N_float + 2.0) / (a_1 * N_float)) + (a_2 / (a_1 * a_1));
        e1 = c1 / a_1;
        e2 = c2 / ((a_1 * a_1) + a_2);

        out_a_1 = a_1;
    }
    else
    {
        double a_1 = 0;
        double a_2 = 0;

        for (double n = 1; n < N; n++)
        {
            a_1 = a_1 + (1.0 / n);
            a_2 = a_2 + (1.0 / (n * n));
        }
        cout << "a1: " << a_1 << endl;
        if (isinf(a_2))
        {
            // pi_squared /6
            a_2 = 1.64493128;
        }
        cout << "a2: " << a_2 << endl;

        b1 = (N_float + 1.0) / (3.0 * (N_float - 1.0));
        b2 = (2.0 * ((N_float * N_float) + N_float + 3.0)) / (9.0 * N_float * (N_float - 1.0));
        c1 = b1 - (1.0 / a_1);
        c2 = b2 - ((N_float + 2.0) / (a_1 * N_float)) + (a_2 / (a_1 * a_1));
        e1 = c1 / a_1;
        e2 = c2 / ((a_1 * a_1) + a_2);

        out_a_1 = a_1;
    }

    cout << "b1: " << b1 << endl;
    cout << "b2: " << b2 << endl;
    cout << "c1: " << c1 << endl;
    cout << "c2: " << c2 << endl;
    cout << "e1: " << e1 << endl;
    cout << "e2: " << e2 << endl;
}

void cancer_Host::write_Tajima(string type, int *per_Region, string &output_Tajima_File, int &overall_Generations, string &tissue_Name,
                               int &num_Regions, int &N, double &N_float, double &out_a_1, double &e1, double &e2)
{
    cout << "\nCalculating Tajima's D for ";
    cout << type << endl;
    // functions.create_File(output_Tajima_File, "Type\tGeneration\tTissue\tN_cells\tseg_Sites\tpairwise_Diff\tCombinations\tPi\ttajima_D");

    double tot_pairwise_Differences = 0;
    double seg_sites_Count = 0;
    for (int region = 0; region < num_Regions; region++)
    {
        double MAF = (double)per_Region[region] / N_float;
        if (MAF != 0 && MAF != 1)
        {
            if (MAF > 0.5)
            {
                MAF = 1.0 - MAF;
            }
            // cout << MAF << endl;
            tot_pairwise_Differences = tot_pairwise_Differences + (MAF * (1 - MAF) * pow(N_float, 2));
            seg_sites_Count++;
        }
    }

    if (seg_sites_Count > 0)
    {
        cout << "Total pairwise differences: " << tot_pairwise_Differences << endl;
        long int combinations = combos_N(N);
        cout << "Combinations: " << combinations << endl;
        double pi = (double)tot_pairwise_Differences / combinations;
        cout << "Pi: " << pi << endl;
        double D = (double)(pi - (seg_sites_Count / out_a_1)) / sqrt(((e1 * seg_sites_Count) + (e2 * seg_sites_Count * (seg_sites_Count - 1))));
        cout << "Tajima's D: " << D << endl;

        fstream tajima_Write;
        tajima_Write.open(output_Tajima_File, ios::app);
        if (tajima_Write.is_open())
        {
            cout << "Writing to file: " << output_Tajima_File << endl;
            tajima_Write << type << "\t" << to_string(overall_Generations) << "\t" << tissue_Name << "\t"
                         << to_string(N) << "\t" << seg_sites_Count << "\t" << to_string(tot_pairwise_Differences) << "\t" << to_string(combinations) << "\t" << to_string(pi) << "\t"
                         << to_string(D) << endl;
            tajima_Write.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN TAJIMA OUTPUT FILE: " << output_Tajima_File << endl;
            exit(-1);
        }
        // exit(-1);
    }
    else
    {
        cout << "No segregating sites due to fixation\n";

        fstream tajima_Write;
        tajima_Write.open(output_Tajima_File, ios::app);
        if (tajima_Write.is_open())
        {
            cout << "Writing to file: " << output_Tajima_File << endl;
            tajima_Write << type << "\t" << to_string(overall_Generations) << "\t" << tissue_Name << "\t"
                         << to_string(N) << "\t0\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t"
                         << "NA" << endl;
            tajima_Write.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN TAJIMA OUTPUT FILE: " << output_Tajima_File << endl;
            exit(-1);
        }
    }
}

long int cancer_Host::fact_half(int count)
{
    long int tot = 1;
    for (int i = count; i > count - 2; i--)
    {
        tot = tot * i;
    }
    return tot;
}

long int cancer_Host::combos_N(int count)
{
    long int combinations;

    combinations = fact_half(count) / 2;

    return combinations;
}

__global__ void cuda_tajima_calc(int num_Regions, int genome_Length, int **cuda_tajima_regions_Start_Stop,
                                 char *sites, char *cuda_reference_Genome, int *cuda_per_Region, int sequences_Total,
                                 char *cuda_char_dead_or_Alive, int *cuda_per_Region_ALIVE)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < sequences_Total)
    {
        int site_Start = tid * genome_Length;
        // int site_End = site_Start + genome_Length;

        for (int region = 0; region < num_Regions; region++)
        {
            int region_Count = 0;

            int region_Start = site_Start + cuda_tajima_regions_Start_Stop[region][0] - 1;
            int region_Stop = site_Start + cuda_tajima_regions_Start_Stop[region][1];

            int reference_Track = cuda_tajima_regions_Start_Stop[region][0] - 1;

            for (int site = region_Start; site < region_Stop; site++)
            {
                if (sites[site] != cuda_reference_Genome[reference_Track])
                {
                    region_Count++;
                    break;
                }
                reference_Track++;
            }

            if (region_Count != 0)
            {
                atomicAdd(&cuda_per_Region[region], 1);

                if (cuda_char_dead_or_Alive[tid] == 'A')
                {
                    atomicAdd(&cuda_per_Region_ALIVE[region], 1);
                }
                // cuda_per_Region[region] = 1;
            }
            // else
            // {
            //     cuda_per_Region[tid][region] = 0;
            // }
        }

        tid += blockDim.x * gridDim.x;
    }
}

void cancer_Host::process_Tajima_String(string &all_Sequences, int &count_Track, int &num_Regions, int **cuda_tajima_regions_Start_Stop,
                                        char *cuda_Reference_Genome, int *cuda_per_Region,
                                        functions_library &functions,
                                        string &dead_or_Alive, int *cuda_per_Region_ALIVE)
{
    cout << "Processing MAF for " << count_Track << " cells\n";

    // cout << "Dead alive size: " << dead_or_Alive.size() << endl;

    // exit(-1);

    char *full_Char;
    full_Char = (char *)malloc((all_Sequences.size() + 1) * sizeof(char));
    strcpy(full_Char, all_Sequences.c_str());

    char *cuda_full_Char;
    cudaMallocManaged(&cuda_full_Char, (all_Sequences.size() + 1) * sizeof(char));
    cudaMemcpy(cuda_full_Char, full_Char, (all_Sequences.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    free(full_Char);
    all_Sequences = "";

    char *char_dead_or_Alive;
    char_dead_or_Alive = (char *)malloc((dead_or_Alive.size() + 1) * sizeof(char));
    strcpy(char_dead_or_Alive, dead_or_Alive.c_str());

    char *cuda_char_dead_or_Alive;
    cudaMallocManaged(&cuda_char_dead_or_Alive, (dead_or_Alive.size() + 1) * sizeof(char));
    cudaMemcpy(cuda_char_dead_or_Alive, char_dead_or_Alive, (dead_or_Alive.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);
    free(char_dead_or_Alive);
    dead_or_Alive = "";

    cuda_tajima_calc<<<functions.tot_Blocks_array[0], functions.tot_ThreadsperBlock_array[0]>>>(num_Regions, genome_Length, cuda_tajima_regions_Start_Stop,
                                                                                                cuda_full_Char, cuda_Reference_Genome, cuda_per_Region, count_Track,
                                                                                                cuda_char_dead_or_Alive, cuda_per_Region_ALIVE);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", 0, cudaGetErrorString(err));
        exit(-1);
    }

    cudaFree(cuda_full_Char);
    cudaFree(cuda_char_dead_or_Alive);

    // exit(-1);

    count_Track = 0;
}

void cancer_Host::migration_of_Cells(string &source_sequence_Data_folder, vector<string> &tissue_Names,
                                     int source, vector<pair<int, int>> tissue_migration_Targets_amount, multiset<pair<float, int>> &migration_cell_List,
                                     int &overall_Generations,
                                     functions_library &functions)
{
    cout << "\nConfiguring metastatic cell migration from " << tissue_Names[source] << ": ";
    // check if migrating target tissue folder exists already
    int different_Tissues = tissue_migration_Targets_amount.size();
    cout << different_Tissues << " different tissue(s).\n";

    vector<vector<int>> tissue_cell_Index;
    for (int tissue = 0; tissue < different_Tissues; tissue++)
    {
        vector<int> indexes;
        tissue_cell_Index.push_back(indexes);
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, different_Tissues - 1);

    cout << "Assinging metastatic cells to new tissue\n";

    for (const auto &cell : migration_cell_List)
    {
        cout << "Cell: " << cell.second << " migrating to: ";
        int insert_Check = 0;
        do
        {
            int target_tissue = dis(gen);

            if (tissue_cell_Index[target_tissue].size() < tissue_migration_Targets_amount[target_tissue].second)
            {
                tissue_cell_Index[target_tissue].push_back(cell.second);
                insert_Check = 1;
                cout << tissue_Names[tissue_migration_Targets_amount[target_tissue].first] << endl;
                removed_by_Transfer_Indexes[source].insert(cell.second);
            }

        } while (insert_Check == 0);
    }

    cout << "\nAll cells assigned\nWriting migrations\n";

    string sequence_Tissue_Folder = source_sequence_Data_folder + "/" + to_string(source) + "/generation_" + to_string(overall_Generations + 1);
    vector<pair<int, int>> indexed_Source_Folder = functions.index_sequence_Folder(sequence_Tissue_Folder);

    for (int tissue = 0; tissue < different_Tissues; tissue++)
    {
        cout << "\nMoving sequences to tissue: " << tissue_Names[tissue_migration_Targets_amount[tissue].first] << endl;

        string sequence_target_Folder = source_sequence_Data_folder + "/" + to_string(tissue_migration_Targets_amount[tissue].first) + "/generation_" + to_string(overall_Generations + 1);
        int last_target_Index = 0;

        if (filesystem::exists(sequence_target_Folder))
        {
            vector<pair<int, int>> sequence_target_Index = functions.index_sequence_Folder(sequence_target_Folder);
            last_target_Index = sequence_target_Index[sequence_target_Index.size() - 1].second + 1;
        }
        else
        {
            string dead_List = sequence_target_Folder + "/dead_List.txt";
            functions.config_Folder(sequence_target_Folder, tissue_Names[tissue_migration_Targets_amount[tissue].first] + " Generation " + to_string(overall_Generations + 1));
            functions.create_File(dead_List);
        }

        int start_Index = last_target_Index;

        int num_Sequences = tissue_cell_Index[tissue].size();
        if (num_Sequences > 0)
        {
            sort(tissue_cell_Index[tissue].begin(), tissue_cell_Index[tissue].end());

            cout << "Collecting " << num_Sequences << " sequence(s)\n";
            // int num_per_Core = num_Sequences / this->CPU_cores;
            // int remainder = num_Sequences % this->CPU_cores;

            // vector<thread> threads_vec;

            // for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
            // {
            //     int start_Cell = core_ID * num_per_Core;
            //     int stop_Cell = start_Cell + num_per_Core;

            //     threads_vec.push_back(thread{&cancer_Host::thread_find_Files_2, this, start_Cell, stop_Cell, ref(tissue_cell_Index[tissue]), ref(indexed_Source_Folder)});
            // }

            // if (remainder != 0)
            // {
            //     int start_Cell = num_Sequences - remainder;
            //     int stop_Cell = num_Sequences;

            //     threads_vec.push_back(thread{&cancer_Host::thread_find_Files_2, this, start_Cell, stop_Cell, ref(tissue_cell_Index[tissue]), ref(indexed_Source_Folder)});
            // }

            // for (thread &t : threads_vec)
            // {
            //     if (t.joinable())
            //     {
            //         t.join();
            //     }
            // }

            // threads_vec.clear();

            // vector<int> Tissue_files(found_Tissue_Folder_Indexes.begin(), found_Tissue_Folder_Indexes.end());
            // found_Tissue_Folder_Indexes.clear();

            // sort(Tissue_files.begin(), Tissue_files.end());

            // cout << Tissue_files.size() << " file(s) identified\n";

            vector<pair<string, string>> all_Sequences;
            cout << "Retrieving sequence(s)\n";

            fstream nfasta;
            int index_Files = 0;
            int line_current = 0;
            nfasta.open(sequence_Tissue_Folder + "/" + to_string(indexed_Source_Folder[index_Files].first) + "_" + to_string(indexed_Source_Folder[index_Files].second) + ".nfasta", ios::in);
            cout << "File: " << sequence_Tissue_Folder << "/" << to_string(indexed_Source_Folder[index_Files].first) + "_" << to_string(indexed_Source_Folder[index_Files].second) << ".nfasta\n";

            fstream sequence_Profiles_File;
            sequence_Profiles_File.open(sequence_Profiles, ios::app);
            fstream sequence_parent_Progeny_relationships_File;
            sequence_parent_Progeny_relationships_File.open(sequence_parent_Progeny_relationships, ios::app);

            for (int find = 0; find < num_Sequences; find++)
            {
                cout << "Looking for " << find << " :" << tissue_cell_Index[tissue][find] << ": of " << num_Sequences << "\n";
                if (index_Files < indexed_Source_Folder.size())
                {
                    while ((indexed_Source_Folder[index_Files].first <= tissue_cell_Index[tissue][find] && indexed_Source_Folder[index_Files].second >= tissue_cell_Index[tissue][find]) == 0)
                    {
                        nfasta.close();
                        index_Files++;
                        nfasta.open(sequence_Tissue_Folder + "/" + to_string(indexed_Source_Folder[index_Files].first) + "_" + to_string(indexed_Source_Folder[index_Files].second) + ".nfasta", ios::in);

                        cout << "File: " << sequence_Tissue_Folder << "/" << to_string(indexed_Source_Folder[index_Files].first) + "_" << to_string(indexed_Source_Folder[index_Files].second) << ".nfasta\n";

                        line_current = 0;
                    }

                    if (nfasta.is_open())
                    {
                        int line_t0_check = (tissue_cell_Index[tissue][find] - indexed_Source_Folder[index_Files].first) * 2;

                        string line;
                        string sequence = "";

                        while (getline(nfasta, line))
                        {
                            if (line_t0_check == line_current)
                            {
                                cout << "Migration Sequence found\n";
                                vector<string> line_Data;
                                functions.split(line_Data, line, '_');
                                if (stoi(line_Data[0].substr(1)) == tissue_cell_Index[tissue][find])
                                {
                                    getline(nfasta, line);
                                    cout << "Write check\n";
                                    all_Sequences.push_back(make_pair(to_string(last_target_Index) + "_A_" + line_Data[2] + "_" + line_Data[3] + "_" + line_Data[4] + "_" + line_Data[5] + "_" + line_Data[6] + "_" + line_Data[7], line));
                                    cout << "Write check2\n";
                                    sequence_Profiles_File << tissue_Names[tissue_migration_Targets_amount[tissue].first] << "_" << to_string(overall_Generations + 1) << "_" << to_string(last_target_Index) << "\t" << tissue_Names[tissue_migration_Targets_amount[tissue].first] << "\t" << line_Data[5] << "\t" << line_Data[3] << "\t" << line_Data[2] << "\t" << line_Data[6] << "\t" << line_Data[7] << endl;
                                    cout << "Write check3\n";
                                    sequence_parent_Progeny_relationships_File << tissue_Names[source] << "_" << to_string(overall_Generations + 1) << "_" << tissue_cell_Index[tissue][find] << "\t"
                                                                               << tissue_Names[tissue_migration_Targets_amount[tissue].first] << "_" << to_string(overall_Generations + 1) << "_" << to_string(last_target_Index)
                                                                               << "\tMigration" << endl;
                                    last_target_Index++;
                                    line_current++;
                                }
                                else
                                {
                                    cout << "ERROR: CORRECT SEQUENCE NOT FOUND AT INDEX\n";
                                    cout << "Looking for: " << tissue_cell_Index[tissue][find] << endl
                                         << "Sequence ID at location: " << line << endl
                                         << "File: " << sequence_Tissue_Folder << "/" << indexed_Source_Folder[index_Files].first
                                         << "_" << indexed_Source_Folder[index_Files].second << ".nfasta" << endl;
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
                        cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << sequence_Tissue_Folder << "/" << indexed_Source_Folder[index_Files].first << "_" << indexed_Source_Folder[index_Files].second << ".nfasta" << endl;
                        exit(-1);
                    }
                }
                else
                {
                    cout << "ERROR: INDEX FILES OVERREACHED: " << index_Files << " SHOULD BE LESS THAN: " << indexed_Source_Folder.size();
                    exit(-1);
                }
            }
            nfasta.close();

            sequence_Profiles_File.close();
            sequence_parent_Progeny_relationships_File.close();

            cout << "Found " << all_Sequences.size() << " sequences\nWriting migrating sequences: ";
            fstream write_Migration;
            write_Migration.open(sequence_target_Folder + "/" + to_string(start_Index) + "_" + to_string(last_target_Index - 1) + ".nfasta", ios::out);

            if (write_Migration.is_open())
            {
                for (int sequence = 0; sequence < all_Sequences.size(); sequence++)
                {
                    write_Migration << ">" << all_Sequences[sequence].first << endl
                                    << all_Sequences[sequence].second << endl;
                }
                write_Migration.close();
                cout << "Complete\n";
            }
            else
            {
                cout << "ERROR: UNABLE TO CREATE METASTIC SEQUENCES FILE: " << sequence_target_Folder << "/" << to_string(start_Index) << "_" << to_string(last_target_Index - 1) << ".nfasta" << endl;
                exit(-1);
            }
        }
    }
}

void cancer_Host::thread_find_Files_2(int start, int stop, vector<int> &cell_Indexes, vector<pair<int, int>> &indexed_Tissue_Folder)
{
    vector<int> caught_Indexes;

    for (int sequence = start; sequence < stop; sequence++)
    {
        int top = 0;
        int bottom = indexed_Tissue_Folder.size() - 1;
        int middle = top + ((bottom - top) / 2);

        while (top <= bottom)
        {
            if ((indexed_Tissue_Folder[middle].first <= cell_Indexes[sequence]) && (indexed_Tissue_Folder[middle].second >= cell_Indexes[sequence]))
            {
                caught_Indexes.push_back(middle);
                // found = 'Y';
                break;
            }
            else if (indexed_Tissue_Folder[middle].first < cell_Indexes[sequence])
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

    unique_lock<shared_mutex> ul(g_mutex);
    for (int index = 0; index < caught_Indexes.size(); index++)
    {
        found_Tissue_Folder_Indexes.insert(caught_Indexes[index]);
    }
}

__global__ void cuda_replicate_Progeny_Main(int num_Parents,
                                            int **progeny_sequences_INT, int genome_Length, char *sites,
                                            float *cuda_Reference_fitness_survivability_proof_reading, float *cuda_Reference_cancer_parameters,
                                            float **cuda_A_0_mutation,
                                            float **cuda_T_1_mutation,
                                            float **cuda_G_2_mutation,
                                            float **cuda_C_3_mutation,
                                            int mutation_Hotspots,
                                            float **cuda_mutation_hotspot_parameters,
                                            int *cuda_num_effect_Segregating_sites,
                                            float **cuda_sequence_Survivability_changes,
                                            float **cuda_sequence_Proof_reading_changes,
                                            int *cuda_num_effect_Segregating_sites_Cancer,
                                            float **cuda_sequence_replication_factor_changes,
                                            float **cuda_sequence_mutation_rate_changes,
                                            float **cuda_sequence_generation_death_changes,
                                            float **cuda_sequence_replication_prob_changes,
                                            float **cuda_sequence_metastatic_prob_changes,
                                            float **progeny_Configuration_Cancer,
                                            float *cuda_parents_Elapsed, float *cuda_progeny_Elapsed,
                                            int *cuda_tissue_selection_Position_Count,
                                            float **cuda_tissues_ATGC_positions_Survivability,
                                            float **cuda_tissues_ATGC_positions_Proof,
                                            float **cuda_tissues_ATGC_positions_Replication_factor,
                                            float **cuda_tissues_ATGC_positions_Mutation_rate_factor,
                                            float **cuda_tissues_ATGC_positions_Generation_death,
                                            float **cuda_tissues_ATGC_positions_Replication_prob,
                                            float **cuda_tissues_ATGC_positions_Metastatic,
                                            int *cuda_Survivability_Positions,
                                            int *cuda_Proof_Positions,
                                            int *cuda_Replication_factor_Positions,
                                            int *cuda_Mutation_rate_factor_Positions,
                                            int *cuda_Generation_death_Positions,
                                            int *cuda_Replication_prob_Positions,
                                            int *cuda_Metastatic_Positions,
                                            int tissue)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Parents)
    {
        int site_Start = tid * genome_Length;
        int site_End = site_Start + genome_Length;

        int bp_Pos = 0;

        int progeny_1 = tid * 2;
        // int progeny_2 = progeny_1 + 1;

        for (int site = site_Start; site < site_End; site++)
        {
            if (sites[site] == '0')
            {
                progeny_sequences_INT[progeny_1][bp_Pos] = 0;
                progeny_sequences_INT[progeny_1 + 1][bp_Pos] = 0;
            }
            else if (sites[site] == '1')
            {
                progeny_sequences_INT[progeny_1][bp_Pos] = 1;
                progeny_sequences_INT[progeny_1 + 1][bp_Pos] = 1;
            }
            else if (sites[site] == '2')
            {
                progeny_sequences_INT[progeny_1][bp_Pos] = 2;
                progeny_sequences_INT[progeny_1 + 1][bp_Pos] = 2;
            }
            else if (sites[site] == '3')
            {
                progeny_sequences_INT[progeny_1][bp_Pos] = 3;
                progeny_sequences_INT[progeny_1 + 1][bp_Pos] = 3;
            }

            bp_Pos++;
        }

        curandState localState;
        curand_init(clock64(), tid, 0, &localState);

        for (int progeny = 0; progeny < 2; progeny++)
        {
            int progeny_Index = progeny_1 + progeny;

            float replication_Factor = cuda_Reference_cancer_parameters[0];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[0]; pos++)
            {
                replication_Factor = replication_Factor * cuda_sequence_replication_factor_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_factor_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[2]; pos++)
            {
                replication_Factor = replication_Factor * cuda_tissues_ATGC_positions_Replication_factor[progeny_sequences_INT[progeny_Index][cuda_Replication_factor_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            cuda_progeny_Elapsed[progeny_Index] = cuda_parents_Elapsed[tid] + replication_Factor;

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

                        float mutation_Factor = cuda_Reference_cancer_parameters[1];

                        for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[1]; pos++)
                        {
                            mutation_Factor = mutation_Factor * cuda_sequence_mutation_rate_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_mutation_rate_changes[pos][0] - 1] + 1];
                        }

                        for (int pos = 0; pos < cuda_tissue_selection_Position_Count[3]; pos++)
                        {
                            mutation_Factor = mutation_Factor * cuda_tissues_ATGC_positions_Mutation_rate_factor[progeny_sequences_INT[progeny_Index][cuda_Mutation_rate_factor_Positions[pos] - 1] + (tissue * 4)][pos];
                        }

                        num_Mutations = (int)((num_Mutations * mutation_Factor) + 0.5);
                    }

                    if (num_Mutations > 0 && cuda_Reference_fitness_survivability_proof_reading[2] != -1)
                    {
                        float proof_Reading = cuda_Reference_fitness_survivability_proof_reading[2];

                        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[2]; pos++)
                        {
                            proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] + 1];
                        }

                        for (int pos = 0; pos < cuda_tissue_selection_Position_Count[1]; pos++)
                        {
                            proof_Reading = proof_Reading + cuda_tissues_ATGC_positions_Proof[progeny_sequences_INT[progeny_Index][cuda_Proof_Positions[pos] - 1] + (tissue * 4)][pos];
                        }

                        proof_Reading = (proof_Reading > 1) ? 1 : (proof_Reading < 0) ? 0
                                                                                      : proof_Reading;

                        int count = 0;

                        for (int trial = 0; trial < num_Mutations; trial++)
                        {
                            if (curand_uniform(&localState) < proof_Reading)
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
                            int position = (int)(curand_uniform(&localState) * (((int)cuda_mutation_hotspot_parameters[hotspot][1] - 1) - ((int)cuda_mutation_hotspot_parameters[hotspot][0] - 1) + 1)) + ((int)cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                            float rand_num = curand_uniform(&localState);
                            float cumulative_prob = 0.0f;

                            int original_BASE = progeny_sequences_INT[progeny_Index][position];
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

                            progeny_sequences_INT[progeny_Index][position] = new_Base;
                        }
                    }
                }
            }
            // Progeny configuration
            replication_Factor = cuda_Reference_cancer_parameters[0];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[0]; pos++)
            {
                replication_Factor = replication_Factor * cuda_sequence_replication_factor_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_factor_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[2]; pos++)
            {
                replication_Factor = replication_Factor * cuda_tissues_ATGC_positions_Replication_factor[progeny_sequences_INT[progeny_Index][cuda_Replication_factor_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            progeny_Configuration_Cancer[progeny_Index][0] = replication_Factor;

            float gen_Death_prob = cuda_Reference_cancer_parameters[2];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[2]; pos++)
            {
                gen_Death_prob = gen_Death_prob + cuda_sequence_generation_death_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_generation_death_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[4]; pos++)
            {
                gen_Death_prob = gen_Death_prob + cuda_tissues_ATGC_positions_Generation_death[progeny_sequences_INT[progeny_Index][cuda_Generation_death_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            gen_Death_prob = (gen_Death_prob > 1) ? 1 : (gen_Death_prob < 0) ? 0
                                                                             : gen_Death_prob;
            progeny_Configuration_Cancer[progeny_Index][1] = gen_Death_prob;

            float replication_prob = cuda_Reference_cancer_parameters[3];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[3]; pos++)
            {
                replication_prob = replication_prob + cuda_sequence_replication_prob_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_prob_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[5]; pos++)
            {
                replication_prob = replication_prob + cuda_tissues_ATGC_positions_Replication_prob[progeny_sequences_INT[progeny_Index][cuda_Replication_prob_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            replication_prob = (replication_prob > 1) ? 1 : (replication_prob < 0) ? 0
                                                                                   : replication_prob;
            progeny_Configuration_Cancer[progeny_Index][2] = replication_prob;

            float metatstatic_prob = cuda_Reference_cancer_parameters[4];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[4]; pos++)
            {
                metatstatic_prob = metatstatic_prob + cuda_sequence_metastatic_prob_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_metastatic_prob_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[6]; pos++)
            {
                metatstatic_prob = metatstatic_prob + cuda_tissues_ATGC_positions_Metastatic[progeny_sequences_INT[progeny_Index][cuda_Metastatic_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            metatstatic_prob = (metatstatic_prob > 1) ? 1 : (metatstatic_prob < 0) ? 0
                                                                                   : metatstatic_prob;
            progeny_Configuration_Cancer[progeny_Index][3] = metatstatic_prob;

            float survivability = cuda_Reference_fitness_survivability_proof_reading[1];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites[1]; pos++)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Survivability_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[0]; pos++)
            {
                survivability = survivability + cuda_tissues_ATGC_positions_Survivability[progeny_sequences_INT[progeny_Index][cuda_Survivability_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            survivability = (survivability > 1) ? 1 : (survivability < 0) ? 0
                                                                          : survivability;
            progeny_Configuration_Cancer[progeny_Index][4] = survivability;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void cancer_Host::remainder_Write_Sequences_NEXT_Generation(string next_Generation_location,
                                                            functions_library &functions, vector<pair<string, string>> &to_write_Sequence_Store)
{
    if (to_write_Sequence_Store.size() > 0)
    {
        string start_String = to_write_Sequence_Store[0].first;
        string stop_String = to_write_Sequence_Store[to_write_Sequence_Store.size() - 1].first;

        vector<string> line_Data;

        functions.split(line_Data, start_String, '_');
        start_String = line_Data[0];

        functions.split(line_Data, stop_String, '_');
        stop_String = line_Data[0];

        fstream nFASTA_File;
        nFASTA_File.open(next_Generation_location + "/" + start_String + "_" + stop_String + ".nfasta", ios::out);

        cout << "Writing sequences to file: " << next_Generation_location << "/" << start_String << "_" + stop_String << ".nfasta" << endl;

        if (nFASTA_File.is_open())
        {
            for (int sequence_Index = 0; sequence_Index < to_write_Sequence_Store.size(); sequence_Index++)
            {
                nFASTA_File << ">" << to_write_Sequence_Store[sequence_Index].first << endl
                            << to_write_Sequence_Store[sequence_Index].second << endl;
            }
            nFASTA_File.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO CREATE NFASTA SEQUENCE FILE: " << next_Generation_location << "/" << start_String << "_" << stop_String << ".nfasta";
        }

        to_write_Sequence_Store.clear();
    }
}

void cancer_Host::full_Write_Sequences_NEXT_Generation(int &max_sequences_per_File, string next_Generation_location,
                                                       functions_library &functions, vector<pair<string, string>> &to_write_Sequence_Store)
{
    if (to_write_Sequence_Store.size() > max_sequences_per_File)
    {
        int full_Write_Count = to_write_Sequence_Store.size() / max_sequences_per_File;

        for (int batch = 0; batch < full_Write_Count; batch++)
        {
            string start_String = to_write_Sequence_Store[batch * max_sequences_per_File].first;
            string stop_String = to_write_Sequence_Store[(batch * max_sequences_per_File) + max_sequences_per_File - 1].first;

            vector<string> line_Data;

            functions.split(line_Data, start_String, '_');
            start_String = line_Data[0];

            functions.split(line_Data, stop_String, '_');
            stop_String = line_Data[0];

            fstream nFASTA_File;
            nFASTA_File.open(next_Generation_location + "/" + start_String + "_" + stop_String + ".nfasta", ios::out);

            cout << "Writing sequences to file: " << next_Generation_location << "/" << start_String << "_" + stop_String << ".nfasta" << endl;

            if (nFASTA_File.is_open())
            {
                for (int sequence_Index = (batch * max_sequences_per_File); sequence_Index < ((batch * max_sequences_per_File) + max_sequences_per_File); sequence_Index++)
                {
                    nFASTA_File << ">" << to_write_Sequence_Store[sequence_Index].first << endl
                                << to_write_Sequence_Store[sequence_Index].second << endl;
                }
                nFASTA_File.close();
            }
            else
            {
                cout << "ERROR: UNABLE TO CREATE NFASTA SEQUENCE FILE: " << next_Generation_location << "/" << start_String << "_" << stop_String << ".nfasta";
            }
        }

        if (to_write_Sequence_Store.size() % max_sequences_per_File != 0)
        {
            vector<pair<string, string>> to_write_Sequence_Store_NEXT_Gen_TEMP;

            for (int temp = (to_write_Sequence_Store.size() - (to_write_Sequence_Store.size() % max_sequences_per_File)); temp < to_write_Sequence_Store.size(); temp++)
            {
                to_write_Sequence_Store_NEXT_Gen_TEMP.push_back(make_pair(to_write_Sequence_Store[temp].first, to_write_Sequence_Store[temp].second));
            }
            to_write_Sequence_Store.clear();

            to_write_Sequence_Store = to_write_Sequence_Store_NEXT_Gen_TEMP;
        }
    }
}

__global__ void cuda_replicate_Progeny_reRun(int num_Parents,
                                             int **cuda_parent_sequences_INT, int **cuda_progeny_sequences_INT, int genome_Length,
                                             float *cuda_parents_Elapsed, float *cuda_progeny_Elapsed,
                                             float *cuda_Reference_fitness_survivability_proof_reading, float *cuda_Reference_cancer_parameters,
                                             float **cuda_sequence_replication_factor_changes,
                                             int mutation_Hotspots,
                                             float **cuda_mutation_hotspot_parameters,
                                             float **cuda_A_0_mutation,
                                             float **cuda_T_1_mutation,
                                             float **cuda_G_2_mutation,
                                             float **cuda_C_3_mutation,
                                             int *cuda_num_effect_Segregating_sites,
                                             int *cuda_num_effect_Segregating_sites_Cancer,
                                             float **cuda_sequence_Survivability_changes,
                                             float **cuda_sequence_Proof_reading_changes,
                                             float **cuda_sequence_mutation_rate_changes,
                                             float **cuda_sequence_generation_death_changes,
                                             float **cuda_sequence_replication_prob_changes,
                                             float **cuda_sequence_metastatic_prob_changes,
                                             float **progeny_Configuration_Cancer,
                                             int *cuda_tissue_selection_Position_Count,
                                             float **cuda_tissues_ATGC_positions_Survivability,
                                             float **cuda_tissues_ATGC_positions_Proof,
                                             float **cuda_tissues_ATGC_positions_Replication_factor,
                                             float **cuda_tissues_ATGC_positions_Mutation_rate_factor,
                                             float **cuda_tissues_ATGC_positions_Generation_death,
                                             float **cuda_tissues_ATGC_positions_Replication_prob,
                                             float **cuda_tissues_ATGC_positions_Metastatic,
                                             int *cuda_Survivability_Positions,
                                             int *cuda_Proof_Positions,
                                             int *cuda_Replication_factor_Positions,
                                             int *cuda_Mutation_rate_factor_Positions,
                                             int *cuda_Generation_death_Positions,
                                             int *cuda_Replication_prob_Positions,
                                             int *cuda_Metastatic_Positions,
                                             int tissue)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Parents)
    {
        // int index_Parent = cuda_rerun_Progeny_Indexes[tid];

        int progeny_1 = tid * 2;

        for (int base = 0; base < genome_Length; base++)
        {
            cuda_progeny_sequences_INT[progeny_1][base] = cuda_parent_sequences_INT[tid][base];
            cuda_progeny_sequences_INT[progeny_1 + 1][base] = cuda_parent_sequences_INT[tid][base];
        }

        curandState localState;
        curand_init(clock64(), tid, 0, &localState);

        for (int progeny = 0; progeny < 2; progeny++)
        {
            int progeny_Index = progeny_1 + progeny;

            float replication_Factor = cuda_Reference_cancer_parameters[0];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[0]; pos++)
            {
                replication_Factor = replication_Factor * cuda_sequence_replication_factor_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_factor_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[2]; pos++)
            {
                replication_Factor = replication_Factor * cuda_tissues_ATGC_positions_Replication_factor[cuda_progeny_sequences_INT[progeny_Index][cuda_Replication_factor_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            cuda_progeny_Elapsed[progeny_Index] = cuda_parents_Elapsed[tid] + replication_Factor;

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

                        float mutation_Factor = cuda_Reference_cancer_parameters[1];

                        for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[1]; pos++)
                        {
                            mutation_Factor = mutation_Factor * cuda_sequence_mutation_rate_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_mutation_rate_changes[pos][0] - 1] + 1];
                        }

                        for (int pos = 0; pos < cuda_tissue_selection_Position_Count[3]; pos++)
                        {
                            mutation_Factor = mutation_Factor * cuda_tissues_ATGC_positions_Mutation_rate_factor[cuda_progeny_sequences_INT[progeny_Index][cuda_Mutation_rate_factor_Positions[pos] - 1] + (tissue * 4)][pos];
                        }

                        num_Mutations = (int)((num_Mutations * mutation_Factor) + 0.5);
                    }

                    if (num_Mutations > 0 && cuda_Reference_fitness_survivability_proof_reading[2] != -1)
                    {
                        float proof_Reading = cuda_Reference_fitness_survivability_proof_reading[2];

                        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[2]; pos++)
                        {
                            proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] + 1];
                        }

                        for (int pos = 0; pos < cuda_tissue_selection_Position_Count[1]; pos++)
                        {
                            proof_Reading = proof_Reading + cuda_tissues_ATGC_positions_Proof[cuda_progeny_sequences_INT[progeny_Index][cuda_Proof_Positions[pos] - 1] + (tissue * 4)][pos];
                        }

                        proof_Reading = (proof_Reading > 1) ? 1 : (proof_Reading < 0) ? 0
                                                                                      : proof_Reading;

                        int count = 0;

                        for (int trial = 0; trial < num_Mutations; trial++)
                        {
                            if (curand_uniform(&localState) < proof_Reading)
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
                            int position = (int)(curand_uniform(&localState) * (((int)cuda_mutation_hotspot_parameters[hotspot][1] - 1) - ((int)cuda_mutation_hotspot_parameters[hotspot][0] - 1) + 1)) + ((int)cuda_mutation_hotspot_parameters[hotspot][0] - 1);

                            float rand_num = curand_uniform(&localState);
                            float cumulative_prob = 0.0f;

                            int original_BASE = cuda_progeny_sequences_INT[progeny_Index][position];
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

                            cuda_progeny_sequences_INT[progeny_Index][position] = new_Base;
                        }
                    }
                }
            }
            // Progeny configuration
            replication_Factor = cuda_Reference_cancer_parameters[0];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[0]; pos++)
            {
                replication_Factor = replication_Factor * cuda_sequence_replication_factor_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_factor_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[2]; pos++)
            {
                replication_Factor = replication_Factor * cuda_tissues_ATGC_positions_Replication_factor[cuda_progeny_sequences_INT[progeny_Index][cuda_Replication_factor_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            progeny_Configuration_Cancer[progeny_Index][0] = replication_Factor;

            float gen_Death_prob = cuda_Reference_cancer_parameters[2];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[2]; pos++)
            {
                gen_Death_prob = gen_Death_prob + cuda_sequence_generation_death_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_generation_death_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[4]; pos++)
            {
                gen_Death_prob = gen_Death_prob + cuda_tissues_ATGC_positions_Generation_death[cuda_progeny_sequences_INT[progeny_Index][cuda_Generation_death_Positions[pos] - 1] + (tissue * 4)][pos];
            }
            gen_Death_prob = (gen_Death_prob > 1) ? 1 : (gen_Death_prob < 0) ? 0
                                                                             : gen_Death_prob;
            progeny_Configuration_Cancer[progeny_Index][1] = gen_Death_prob;

            float replication_prob = cuda_Reference_cancer_parameters[3];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[3]; pos++)
            {
                replication_prob = replication_prob + cuda_sequence_replication_prob_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_prob_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[5]; pos++)
            {
                replication_prob = replication_prob + cuda_tissues_ATGC_positions_Replication_prob[cuda_progeny_sequences_INT[progeny_Index][cuda_Replication_prob_Positions[pos] - 1] + (tissue * 4)][pos];
            }
            replication_prob = (replication_prob > 1) ? 1 : (replication_prob < 0) ? 0
                                                                                   : replication_prob;
            progeny_Configuration_Cancer[progeny_Index][2] = replication_prob;

            float metatstatic_prob = cuda_Reference_cancer_parameters[4];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[4]; pos++)
            {
                metatstatic_prob = metatstatic_prob + cuda_sequence_metastatic_prob_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_metastatic_prob_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[6]; pos++)
            {
                metatstatic_prob = metatstatic_prob + cuda_tissues_ATGC_positions_Metastatic[cuda_progeny_sequences_INT[progeny_Index][cuda_Metastatic_Positions[pos] - 1] + (tissue * 4)][pos];
            }
            metatstatic_prob = (metatstatic_prob > 1) ? 1 : (metatstatic_prob < 0) ? 0
                                                                                   : metatstatic_prob;
            progeny_Configuration_Cancer[progeny_Index][3] = metatstatic_prob;

            float survivability = cuda_Reference_fitness_survivability_proof_reading[1];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites[1]; pos++)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Survivability_changes[pos][0] - 1] + 1];
            }

            for (int pos = 0; pos < cuda_tissue_selection_Position_Count[0]; pos++)
            {
                survivability = survivability + cuda_tissues_ATGC_positions_Survivability[cuda_progeny_sequences_INT[progeny_Index][cuda_Survivability_Positions[pos] - 1] + (tissue * 4)][pos];
            }

            survivability = (survivability > 1) ? 1 : (survivability < 0) ? 0
                                                                          : survivability;
            progeny_Configuration_Cancer[progeny_Index][4] = survivability;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void cancer_Host::simulate_cell_Round(functions_library &functions, string &multi_Read, int &num_Cuda_devices, int *CUDA_device_IDs,
                                      int &num_of_Cells, int &start, int &stop,
                                      vector<int> &parents_in_Tissue, int &tissue, string tissue_Name,
                                      vector<pair<int, int>> &indexed_Tissue_Folder,
                                      string this_Gen_intermediary_Sequences,
                                      int &overall_Generations,
                                      int &last_index_Seq_Written,
                                      mt19937 &gen,
                                      float *Reference_fitness_survivability_proof_reading, float *Reference_cancer_parameters,
                                      float **A_0_mutation,
                                      float **T_1_mutation,
                                      float **G_2_mutation,
                                      float **C_3_mutation,
                                      int &mutation_Hotspots,
                                      float **mutation_hotspot_parameters,
                                      int *num_effect_Segregating_sites,
                                      float **sequence_Survivability_changes,
                                      float **sequence_Proof_reading_changes,
                                      int *num_effect_Segregating_sites_Cancer,
                                      float **sequence_replication_factor_changes,
                                      float **sequence_mutation_rate_changes,
                                      float **sequence_generation_death_changes,
                                      float **sequence_replication_prob_changes,
                                      float **sequence_metastatic_prob_changes,
                                      int &max_sequences_per_File, string &intermediary_Tissue_folder, string &source_sequence_Data_folder,
                                      int &last_Progeny_written_this_Gen, string rapid_Progeny_Location,
                                      int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List,
                                      int *tissue_selection_Position_Count,
                                      int *Survivability_Positions,
                                      int *Proof_Positions,
                                      int *Replication_factor_Positions,
                                      int *Mutation_rate_factor_Positions,
                                      int *Generation_death_Positions,
                                      int *Replication_prob_Positions,
                                      int *Metastatic_Positions,
                                      float **tissues_ATGC_positions_Survivability,
                                      float **tissues_ATGC_positions_Proof,
                                      float **tissues_ATGC_positions_Replication_factor,
                                      float **tissues_ATGC_positions_Mutation_rate_factor,
                                      float **tissues_ATGC_positions_Generation_death,
                                      float **tissues_ATGC_positions_Replication_prob,
                                      float **tissues_ATGC_positions_Metastatic,
                                      string &viral_Migration)
{

    // sort(parents_in_Tissue + start, parents_in_Tissue + stop);

    vector<int> parent_IDs;
    // vector<float> parents_Elapsed;

    float *parents_Elapsed = (float *)malloc(sizeof(float) * num_of_Cells);

    string all_Sequences = find_Sequences_Master(start, tissue, tissue_Name, functions, this_Gen_intermediary_Sequences, parents_in_Tissue, num_of_Cells, indexed_Tissue_Folder, overall_Generations, parent_IDs, parents_Elapsed, last_index_Seq_Written, gen, tissue_Migration_Total, migration_cell_List, viral_Migration);
    full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
    // remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions);

    // for (int test = 0; test < parent_IDs.size(); test++)
    // {
    //     cout << parents_Elapsed[test] << endl;
    // }

    // exit(-1);

    parents_in_Tissue.clear();

    if (parent_IDs.size() != 0)
    {
        if (all_Sequences != "")
        {
            int parent_Cells_Found = parent_IDs.size();
            cout << "\nNumber of parents found: " << parent_Cells_Found << endl;
            cout << "Potential mitotic progeny produced: " << to_string(parent_Cells_Found * 2) << endl;

            cout << "\nMain Progeny generation\n";

            int standard_num_per_GPU = parent_Cells_Found / num_Cuda_devices;
            int remainder = parent_Cells_Found % num_Cuda_devices;

            vector<pair<int, int>> start_stop_Per_GPU;

            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                int start_GPU = gpu * standard_num_per_GPU;
                int stop_GPU = start_GPU + standard_num_per_GPU;

                start_stop_Per_GPU.push_back(make_pair(start_GPU, stop_GPU));
            }

            start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

            // string all_Sequences = "";

            // for (int sequence = 0; sequence < parent_Cells_Found; sequence++)
            // {
            //     all_Sequences.append(collected_Sequences[sequence]);
            // }

            char *full_Char;
            full_Char = (char *)malloc((all_Sequences.size() + 1) * sizeof(char));
            strcpy(full_Char, all_Sequences.c_str());
            all_Sequences = "";

            cudaStream_t streams[num_Cuda_devices];

            char *cuda_full_Char[num_Cuda_devices];
            int **cuda_progeny_Sequences[num_Cuda_devices];

            float *cuda_Reference_fitness_survivability_proof_reading[num_Cuda_devices];
            float *cuda_Reference_cancer_parameters[num_Cuda_devices];

            float **cuda_A_0_mutation[num_Cuda_devices];
            float **cuda_T_1_mutation[num_Cuda_devices];
            float **cuda_G_2_mutation[num_Cuda_devices];
            float **cuda_C_3_mutation[num_Cuda_devices];

            float **cuda_mutation_hotspot_parameters[num_Cuda_devices];

            int *cuda_num_effect_Segregating_sites[num_Cuda_devices];
            float **cuda_sequence_Survivability_changes[num_Cuda_devices];
            float **cuda_sequence_Proof_reading_changes[num_Cuda_devices];

            int *cuda_num_effect_Segregating_sites_Cancer[num_Cuda_devices];
            float **cuda_sequence_replication_factor_changes[num_Cuda_devices];
            float **cuda_sequence_mutation_rate_changes[num_Cuda_devices];
            float **cuda_sequence_generation_death_changes[num_Cuda_devices];
            float **cuda_sequence_replication_prob_changes[num_Cuda_devices];
            float **cuda_sequence_metastatic_prob_changes[num_Cuda_devices];

            float **cuda_progeny_Configuration_Cancer[num_Cuda_devices];

            float *cuda_parents_Elapsed[num_Cuda_devices];

            float *cuda_progeny_Elapsed[num_Cuda_devices];

            // Tisse specific selection arrays
            int *cuda_tissue_selection_Position_Count[num_Cuda_devices];

            float **cuda_tissues_ATGC_positions_Survivability[num_Cuda_devices];
            float **cuda_tissues_ATGC_positions_Proof[num_Cuda_devices];
            float **cuda_tissues_ATGC_positions_Replication_factor[num_Cuda_devices];
            float **cuda_tissues_ATGC_positions_Mutation_rate_factor[num_Cuda_devices];
            float **cuda_tissues_ATGC_positions_Generation_death[num_Cuda_devices];
            float **cuda_tissues_ATGC_positions_Replication_prob[num_Cuda_devices];
            float **cuda_tissues_ATGC_positions_Metastatic[num_Cuda_devices];

            int *cuda_Survivability_Positions[num_Cuda_devices];
            int *cuda_Proof_Positions[num_Cuda_devices];
            int *cuda_Replication_factor_Positions[num_Cuda_devices];
            int *cuda_Mutation_rate_factor_Positions[num_Cuda_devices];
            int *cuda_Generation_death_Positions[num_Cuda_devices];
            int *cuda_Replication_prob_Positions[num_Cuda_devices];
            int *cuda_Metastatic_Positions[num_Cuda_devices];

            // int tissue;

            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                cudaDeviceProp deviceProp;
                cudaSetDevice(CUDA_device_IDs[gpu]);
                cudaGetDeviceProperties(&deviceProp, gpu);

                cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

                int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

                cudaMalloc(&cuda_full_Char[gpu], cell_Count * genome_Length * sizeof(char));
                cudaMemcpy(cuda_full_Char[gpu], full_Char + (start_stop_Per_GPU[gpu].first * genome_Length), cell_Count * genome_Length * sizeof(char), cudaMemcpyHostToDevice);

                cudaMallocManaged(&cuda_progeny_Sequences[gpu], (cell_Count * 2) * sizeof(int *));
                for (int row = 0; row < (cell_Count * 2); row++)
                {
                    cudaMalloc((void **)&cuda_progeny_Sequences[gpu][row], genome_Length * sizeof(int));
                }

                cudaMallocManaged(&cuda_Reference_fitness_survivability_proof_reading[gpu], 3 * sizeof(float));
                cudaMemcpy(cuda_Reference_fitness_survivability_proof_reading[gpu], Reference_fitness_survivability_proof_reading, 3 * sizeof(float), cudaMemcpyHostToDevice);

                cudaMallocManaged(&cuda_Reference_cancer_parameters[gpu], 5 * sizeof(float));
                cudaMemcpy(cuda_Reference_cancer_parameters[gpu], Reference_cancer_parameters, 5 * sizeof(float), cudaMemcpyHostToDevice);

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

                cudaMallocManaged(&cuda_num_effect_Segregating_sites[gpu], 3 * sizeof(int));
                cudaMemcpy(cuda_num_effect_Segregating_sites[gpu], num_effect_Segregating_sites, 3 * sizeof(int), cudaMemcpyHostToDevice);

                cudaMallocManaged(&cuda_sequence_Survivability_changes[gpu], num_effect_Segregating_sites[1] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_Survivability_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_Survivability_changes[gpu][row], sequence_Survivability_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_sequence_Proof_reading_changes[gpu], num_effect_Segregating_sites[2] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_Proof_reading_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_Proof_reading_changes[gpu][row], sequence_Proof_reading_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_num_effect_Segregating_sites_Cancer[gpu], 5 * sizeof(int));
                cudaMemcpy(cuda_num_effect_Segregating_sites_Cancer[gpu], num_effect_Segregating_sites_Cancer, 5 * sizeof(int), cudaMemcpyHostToDevice);

                cudaMallocManaged(&cuda_sequence_replication_factor_changes[gpu], num_effect_Segregating_sites_Cancer[0] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_replication_factor_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_replication_factor_changes[gpu][row], sequence_replication_factor_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_sequence_mutation_rate_changes[gpu], num_effect_Segregating_sites_Cancer[1] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites_Cancer[1]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_mutation_rate_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_mutation_rate_changes[gpu][row], sequence_mutation_rate_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_sequence_generation_death_changes[gpu], num_effect_Segregating_sites_Cancer[2] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites_Cancer[2]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_generation_death_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_generation_death_changes[gpu][row], sequence_generation_death_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_sequence_replication_prob_changes[gpu], num_effect_Segregating_sites_Cancer[3] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites_Cancer[3]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_replication_prob_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_replication_prob_changes[gpu][row], sequence_replication_prob_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_sequence_metastatic_prob_changes[gpu], num_effect_Segregating_sites_Cancer[4] * sizeof(float *));
                for (int row = 0; row < num_effect_Segregating_sites_Cancer[4]; row++)
                {
                    cudaMalloc((void **)&cuda_sequence_metastatic_prob_changes[gpu][row], 5 * sizeof(float));
                    cudaMemcpy(cuda_sequence_metastatic_prob_changes[gpu][row], sequence_metastatic_prob_changes[row], 5 * sizeof(float), cudaMemcpyHostToDevice);
                }

                cudaMallocManaged(&cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float));

                cudaMalloc(&cuda_parents_Elapsed[gpu], cell_Count * sizeof(float));
                cudaMemcpy(cuda_parents_Elapsed[gpu], parents_Elapsed + start_stop_Per_GPU[gpu].first, cell_Count * sizeof(float), cudaMemcpyHostToDevice);

                cudaMallocManaged(&cuda_progeny_Configuration_Cancer[gpu], (cell_Count * 2) * sizeof(float *));
                for (int row = 0; row < (cell_Count * 2); row++)
                {
                    cudaMalloc((void **)&cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float));
                }
                // tissue selection
                // cout << "Test 1\n";
                cudaMallocManaged(&cuda_tissue_selection_Position_Count[gpu], 7 * sizeof(int));
                cudaMemcpy(cuda_tissue_selection_Position_Count[gpu], tissue_selection_Position_Count, 7 * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 2\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Survivability[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[0] > 0)
                {
                    // cout << "Test 2.1\n";
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Survivability[gpu][row], tissue_selection_Position_Count[0] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Survivability[gpu][row], tissues_ATGC_positions_Survivability[row], tissue_selection_Position_Count[0] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                    // cout << "Test 2.2\n";
                }
                // cout << "Test 3\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Proof[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[1] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Proof[gpu][row], tissue_selection_Position_Count[1] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Proof[gpu][row], tissues_ATGC_positions_Proof[row], tissue_selection_Position_Count[1] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                }
                // cout << "Test 4\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Replication_factor[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[2] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Replication_factor[gpu][row], tissue_selection_Position_Count[2] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Replication_factor[gpu][row], tissues_ATGC_positions_Replication_factor[row], tissue_selection_Position_Count[2] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                }
                // cout << "Test 5\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[3] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu][row], tissue_selection_Position_Count[3] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu][row], tissues_ATGC_positions_Mutation_rate_factor[row], tissue_selection_Position_Count[3] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                }
                // cout << "Test 6\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Generation_death[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[4] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Generation_death[gpu][row], tissue_selection_Position_Count[4] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Generation_death[gpu][row], tissues_ATGC_positions_Generation_death[row], tissue_selection_Position_Count[4] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                }
                // cout << "Test 7\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Replication_prob[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[5] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Replication_prob[gpu][row], tissue_selection_Position_Count[5] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Replication_prob[gpu][row], tissues_ATGC_positions_Replication_prob[row], tissue_selection_Position_Count[5] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                }
                // cout << "Test 8\n";
                cudaMallocManaged(&cuda_tissues_ATGC_positions_Metastatic[gpu], (4 * num_Tissues) * sizeof(float *));
                if (tissue_selection_Position_Count[6] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaMalloc((void **)&cuda_tissues_ATGC_positions_Metastatic[gpu][row], tissue_selection_Position_Count[6] * sizeof(float));
                        cudaMemcpy(cuda_tissues_ATGC_positions_Metastatic[gpu][row], tissues_ATGC_positions_Metastatic[row], tissue_selection_Position_Count[6] * sizeof(float), cudaMemcpyHostToDevice);
                    }
                }
                // cout << "Test 9\n";
                cudaMallocManaged(&cuda_Survivability_Positions[gpu], tissue_selection_Position_Count[0] * sizeof(int));
                cudaMemcpy(cuda_Survivability_Positions[gpu], Survivability_Positions, tissue_selection_Position_Count[0] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 10\n";
                cudaMallocManaged(&cuda_Proof_Positions[gpu], tissue_selection_Position_Count[1] * sizeof(int));
                cudaMemcpy(cuda_Proof_Positions[gpu], Proof_Positions, tissue_selection_Position_Count[1] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 11\n";
                cudaMallocManaged(&cuda_Replication_factor_Positions[gpu], tissue_selection_Position_Count[2] * sizeof(int));
                cudaMemcpy(cuda_Replication_factor_Positions[gpu], Replication_factor_Positions, tissue_selection_Position_Count[2] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 12\n";
                cudaMallocManaged(&cuda_Mutation_rate_factor_Positions[gpu], tissue_selection_Position_Count[3] * sizeof(int));
                cudaMemcpy(cuda_Mutation_rate_factor_Positions[gpu], Mutation_rate_factor_Positions, tissue_selection_Position_Count[3] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 13\n";
                cudaMallocManaged(&cuda_Generation_death_Positions[gpu], tissue_selection_Position_Count[4] * sizeof(int));
                cudaMemcpy(cuda_Generation_death_Positions[gpu], Generation_death_Positions, tissue_selection_Position_Count[4] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 14\n";
                cudaMallocManaged(&cuda_Replication_prob_Positions[gpu], tissue_selection_Position_Count[5] * sizeof(int));
                cudaMemcpy(cuda_Replication_prob_Positions[gpu], Replication_prob_Positions, tissue_selection_Position_Count[5] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 15\n";
                cudaMallocManaged(&cuda_Metastatic_Positions[gpu], tissue_selection_Position_Count[6] * sizeof(int));
                cudaMemcpy(cuda_Metastatic_Positions[gpu], Metastatic_Positions, tissue_selection_Position_Count[6] * sizeof(int), cudaMemcpyHostToDevice);
                // cout << "Test 16\n";
                cudaStreamCreate(&streams[gpu]);
            }

            free(full_Char);
            free(parents_Elapsed);

            cout << "GPUs intialized\n";

            cout << "Loaded " << parent_Cells_Found << " parent sequence(s) and all pre-requisites to the GPU(s)\nInitiating GPU(s) execution\n";

            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                cudaSetDevice(CUDA_device_IDs[gpu]);
                cuda_replicate_Progeny_Main<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first,
                                                                                                                                            cuda_progeny_Sequences[gpu], genome_Length, cuda_full_Char[gpu],
                                                                                                                                            cuda_Reference_fitness_survivability_proof_reading[gpu], cuda_Reference_cancer_parameters[gpu],
                                                                                                                                            cuda_A_0_mutation[gpu],
                                                                                                                                            cuda_T_1_mutation[gpu],
                                                                                                                                            cuda_G_2_mutation[gpu],
                                                                                                                                            cuda_C_3_mutation[gpu],
                                                                                                                                            mutation_Hotspots,
                                                                                                                                            cuda_mutation_hotspot_parameters[gpu],
                                                                                                                                            cuda_num_effect_Segregating_sites[gpu],
                                                                                                                                            cuda_sequence_Survivability_changes[gpu],
                                                                                                                                            cuda_sequence_Proof_reading_changes[gpu],
                                                                                                                                            cuda_num_effect_Segregating_sites_Cancer[gpu],
                                                                                                                                            cuda_sequence_replication_factor_changes[gpu],
                                                                                                                                            cuda_sequence_mutation_rate_changes[gpu],
                                                                                                                                            cuda_sequence_generation_death_changes[gpu],
                                                                                                                                            cuda_sequence_replication_prob_changes[gpu],
                                                                                                                                            cuda_sequence_metastatic_prob_changes[gpu],
                                                                                                                                            cuda_progeny_Configuration_Cancer[gpu],
                                                                                                                                            cuda_parents_Elapsed[gpu], cuda_progeny_Elapsed[gpu],
                                                                                                                                            cuda_tissue_selection_Position_Count[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Survivability[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Proof[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Replication_factor[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Generation_death[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Replication_prob[gpu],
                                                                                                                                            cuda_tissues_ATGC_positions_Metastatic[gpu],
                                                                                                                                            cuda_Survivability_Positions[gpu],
                                                                                                                                            cuda_Proof_Positions[gpu],
                                                                                                                                            cuda_Replication_factor_Positions[gpu],
                                                                                                                                            cuda_Mutation_rate_factor_Positions[gpu],
                                                                                                                                            cuda_Generation_death_Positions[gpu],
                                                                                                                                            cuda_Replication_prob_Positions[gpu],
                                                                                                                                            cuda_Metastatic_Positions[gpu],
                                                                                                                                            tissue);
            }

            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                cudaSetDevice(CUDA_device_IDs[gpu]);
                cudaStreamSynchronize(streams[gpu]);

                cudaError_t err = cudaGetLastError();
                if (err != cudaSuccess)
                {
                    fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", gpu, cudaGetErrorString(err));
                    exit(-1);
                }
            }

            cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

            //// CLEAR ARRAYS
            int **progeny_Sequences;
            float **progeny_Configuration_Cancer;
            float *progeny_Elapsed = (float *)malloc(sizeof(float) * parent_Cells_Found * 2);

            progeny_Sequences = (int **)malloc(parent_Cells_Found * 2 * sizeof(int *));
            progeny_Configuration_Cancer = (float **)malloc(parent_Cells_Found * 2 * sizeof(float *));
            for (int row = 0; row < (parent_Cells_Found * 2); row++)
            {
                progeny_Sequences[row] = (int *)malloc(genome_Length * sizeof(int));
                progeny_Configuration_Cancer[row] = (float *)malloc(5 * sizeof(float));
                // converted_Sequences.push_back("");
            }

            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                cudaSetDevice(CUDA_device_IDs[gpu]);
                int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

                for (int row = 0; row < (cell_Count * 2); row++)
                {
                    cudaMemcpy(progeny_Sequences[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Sequences[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
                    cudaMemcpy(progeny_Configuration_Cancer[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float), cudaMemcpyDeviceToHost);
                }

                cudaMemcpy(progeny_Elapsed + (start_stop_Per_GPU[gpu].first * 2), cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float), cudaMemcpyDeviceToHost);
            }

            cout << "Data received by host\n";

            cout << "Partial termination of GPU streams: ";
            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                cudaSetDevice(CUDA_device_IDs[gpu]);
                cudaFree(cuda_full_Char[gpu]);

                cudaFree(cuda_parents_Elapsed[gpu]);

                int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

                for (int row = 0; row < (cell_Count * 2); row++)
                {
                    cudaFree(cuda_progeny_Sequences[gpu][row]);
                }
                cudaFree(cuda_progeny_Sequences[gpu]);

                // cuda_progeny_Configuration_Cancer;
                for (int row = 0; row < (cell_Count * 2); row++)
                {
                    cudaFree(cuda_progeny_Configuration_Cancer[gpu][row]);
                }
                cudaFree(cuda_progeny_Configuration_Cancer[gpu]);

                cudaFree(cuda_progeny_Elapsed[gpu]);

                cudaStreamDestroy(streams[gpu]);
            }

            cout << " GPU(s) released\n";

            for (int test = 0; test < parent_Cells_Found; test++)
            {
                cout << test << ": " << endl;
                for (int progeny = (test * 2); progeny < ((test * 2) + 2); progeny++)
                {
                    cout << progeny_Sequences[progeny][0];
                    cout << " ";
                    for (int col = 0; col < 5; col++)
                    {
                        cout << progeny_Configuration_Cancer[progeny][col] << " ";
                    }
                    cout << "| " << progeny_Elapsed[progeny];
                    cout << endl;
                }
            }

            // exit(-1);

            // cout << "\nConverting sequences to String\n";

            // int num_Progeny_being_Processed = parent_Cells_Found * 2;

            // int num_per_Core = num_Progeny_being_Processed / CPU_cores;
            // int remainder_Core = num_Progeny_being_Processed % CPU_cores;

            // vector<thread> threads_vec;

            // for (int core_ID = 0; core_ID < CPU_cores; core_ID++)
            // {
            //     int start_Node = core_ID * num_per_Core;
            //     int stop_Node = start_Node + num_per_Core;

            //     threads_vec.push_back(thread{&cancer_Host::thread_Sequence_to_String_Cancer, this, start_Node, stop_Node, progeny_Sequences});
            // }

            // if (remainder_Core != 0)
            // {
            //     int start_Node = num_Progeny_being_Processed - remainder_Core;
            //     int stop_Node = num_Progeny_being_Processed;

            //     threads_vec.push_back(thread{&cancer_Host::thread_Sequence_to_String_Cancer, this, start_Node, stop_Node, progeny_Sequences});
            // }

            // for (thread &t : threads_vec)
            // {
            //     if (t.joinable())
            //     {
            //         t.join();
            //     }
            // }

            // threads_vec.clear();

            // for (int test = 0; test < parent_Cells_Found; test++)
            // {
            //     for (int progeny = (test * 2); progeny < ((test * 2) + 2); progeny++)
            //     {
            //         string sequence = "";
            //         for (int base = 0; base < genome_Length; base++)
            //         {
            //             sequence.append(to_string(progeny_Sequences[progeny][base]));
            //             // cout << progeny_Sequences[progeny][base];
            //         }
            //         converted_Sequences.push_back(sequence);
            //     }
            // }

            // exit(-1);

            // cout << "Converted sequences to String\n";

            vector<pair<int, int>> rerun_Progeny = compile_Progeny(functions,
                                                                   intermediary_Tissue_folder, rapid_Progeny_Location,
                                                                   parent_Cells_Found,
                                                                   progeny_Elapsed, progeny_Configuration_Cancer,
                                                                   last_index_Seq_Written, overall_Generations, tissue_Name,
                                                                   progeny_Sequences, tissue,
                                                                   gen, parent_IDs,
                                                                   source_sequence_Data_folder,
                                                                   last_Progeny_written_this_Gen,
                                                                   tissue_Migration_Total, migration_cell_List, viral_Migration);

            full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
            //  full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

            for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue].size(); forward++)
            {
                full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue][forward]);
            }

            cout << "\nClearing arrays\n";
            for (int row = 0; row < parent_Cells_Found * 2; row++)
            {
                free(progeny_Configuration_Cancer[row]);
            }
            free(progeny_Configuration_Cancer);

            if (rerun_Progeny.size() > 0)
            {
                cout << "\nRerun progeny: " << rerun_Progeny.size() << endl;

                cout << "\nProcessing reRun progeny\n";
                int rounds_reRun = 0;

                int tot_Parents = parent_Cells_Found;

                do
                {
                    cout << "reRun Round: " << (rounds_reRun + 1) << endl;
                    parent_Cells_Found = rerun_Progeny.size();

                    start_stop_Per_GPU.clear();

                    standard_num_per_GPU = parent_Cells_Found / num_Cuda_devices;
                    remainder = parent_Cells_Found % num_Cuda_devices;

                    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
                    {
                        int start = gpu * standard_num_per_GPU;
                        int stop_GPU = start + standard_num_per_GPU;

                        start_stop_Per_GPU.push_back(make_pair(start, stop_GPU));
                    }

                    start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

                    cout << "Configuring datapoints: ";
                    // int *rerun_Progeny_Indexes = (int *)malloc(sizeof(int) * parent_Cells_Found);
                    parents_Elapsed = (float *)malloc(sizeof(float) * parent_Cells_Found);
                    parent_IDs.clear();
                    // vector<int> parent_IDs;
                    for (int parent = 0; parent < parent_Cells_Found; parent++)
                    {
                        // rerun_Progeny_Indexes[parent] = rerun_Progeny[parent];
                        parents_Elapsed[parent] = progeny_Elapsed[rerun_Progeny[parent].first];
                        parent_IDs.push_back(rerun_Progeny[parent].second);
                    }
                    rerun_Progeny.clear();
                    free(progeny_Elapsed);

                    cout << "Done\n";

                    int **cuda_progeny_Sequences_INT[num_Cuda_devices];
                    // float *cuda_progeny_Elapsed[num_Cuda_devices];
                    // float **cuda_progeny_Configuration_Cancer[num_Cuda_devices];

                    int **cuda_parent_sequences_INT[num_Cuda_devices];
                    // cuda_parents_Elapsed[num_Cuda_devices];

                    cudaStream_t streams[num_Cuda_devices];

                    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
                    {
                        cudaDeviceProp deviceProp;
                        cudaSetDevice(CUDA_device_IDs[gpu]);
                        cudaGetDeviceProperties(&deviceProp, gpu);

                        cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

                        int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

                        cudaMallocManaged(&cuda_parent_sequences_INT[gpu], cell_Count * sizeof(int *));
                        for (int row = 0; row < cell_Count; row++)
                        {
                            cudaMalloc((void **)&cuda_parent_sequences_INT[gpu][row], genome_Length * sizeof(int));
                            cudaMemcpy(cuda_parent_sequences_INT[gpu][row], progeny_Sequences[rerun_Progeny[row + start_stop_Per_GPU[gpu].first].first], genome_Length * sizeof(int), cudaMemcpyHostToDevice);
                        }

                        cudaMalloc(&cuda_parents_Elapsed[gpu], cell_Count * sizeof(float));
                        cudaMemcpy(cuda_parents_Elapsed[gpu], parents_Elapsed + start_stop_Per_GPU[gpu].first, cell_Count * sizeof(float), cudaMemcpyHostToDevice);

                        // cudaMalloc(&cuda_rerun_Progeny_Indexes[gpu], cell_Count * sizeof(int));
                        // cudaMemcpy(cuda_rerun_Progeny_Indexes[gpu], rerun_Progeny_Indexes + start_stop_Per_GPU[gpu].first, cell_Count * sizeof(int), cudaMemcpyHostToDevice);

                        cudaMallocManaged(&cuda_progeny_Sequences_INT[gpu], (cell_Count * 2) * sizeof(int *));
                        for (int row = 0; row < (cell_Count * 2); row++)
                        {
                            cudaMalloc((void **)&cuda_progeny_Sequences_INT[gpu][row], genome_Length * sizeof(int));
                        }

                        cudaMallocManaged(&cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float));

                        cudaMallocManaged(&cuda_progeny_Configuration_Cancer[gpu], (cell_Count * 2) * sizeof(float *));
                        for (int row = 0; row < (cell_Count * 2); row++)
                        {
                            cudaMalloc((void **)&cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float));
                        }

                        cudaStreamCreate(&streams[gpu]);
                    }

                    free(parents_Elapsed);

                    for (int row = 0; row < tot_Parents * 2; row++)
                    {
                        free(progeny_Sequences[row]);
                    }
                    free(progeny_Sequences);

                    tot_Parents = parent_Cells_Found;

                    cout << "GPUs intialized\n";
                    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
                    {
                        cudaSetDevice(CUDA_device_IDs[gpu]);

                        cuda_replicate_Progeny_reRun<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first,
                                                                                                                                                     cuda_parent_sequences_INT[gpu], cuda_progeny_Sequences_INT[gpu], genome_Length,
                                                                                                                                                     cuda_parents_Elapsed[gpu], cuda_progeny_Elapsed[gpu],
                                                                                                                                                     cuda_Reference_fitness_survivability_proof_reading[gpu], cuda_Reference_cancer_parameters[gpu],
                                                                                                                                                     cuda_sequence_replication_factor_changes[gpu],
                                                                                                                                                     mutation_Hotspots,
                                                                                                                                                     cuda_mutation_hotspot_parameters[gpu],
                                                                                                                                                     cuda_A_0_mutation[gpu],
                                                                                                                                                     cuda_T_1_mutation[gpu],
                                                                                                                                                     cuda_G_2_mutation[gpu],
                                                                                                                                                     cuda_C_3_mutation[gpu],
                                                                                                                                                     cuda_num_effect_Segregating_sites[gpu],
                                                                                                                                                     cuda_num_effect_Segregating_sites_Cancer[gpu],
                                                                                                                                                     cuda_sequence_Survivability_changes[gpu],
                                                                                                                                                     cuda_sequence_Proof_reading_changes[gpu],
                                                                                                                                                     cuda_sequence_mutation_rate_changes[gpu],
                                                                                                                                                     cuda_sequence_generation_death_changes[gpu],
                                                                                                                                                     cuda_sequence_replication_prob_changes[gpu],
                                                                                                                                                     cuda_sequence_metastatic_prob_changes[gpu],
                                                                                                                                                     cuda_progeny_Configuration_Cancer[gpu],
                                                                                                                                                     cuda_tissue_selection_Position_Count[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Survivability[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Proof[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Replication_factor[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Generation_death[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Replication_prob[gpu],
                                                                                                                                                     cuda_tissues_ATGC_positions_Metastatic[gpu],
                                                                                                                                                     cuda_Survivability_Positions[gpu],
                                                                                                                                                     cuda_Proof_Positions[gpu],
                                                                                                                                                     cuda_Replication_factor_Positions[gpu],
                                                                                                                                                     cuda_Mutation_rate_factor_Positions[gpu],
                                                                                                                                                     cuda_Generation_death_Positions[gpu],
                                                                                                                                                     cuda_Replication_prob_Positions[gpu],
                                                                                                                                                     cuda_Metastatic_Positions[gpu],
                                                                                                                                                     tissue);
                    }

                    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
                    {
                        cudaSetDevice(CUDA_device_IDs[gpu]);
                        cudaStreamSynchronize(streams[gpu]);

                        cudaError_t err = cudaGetLastError();
                        if (err != cudaSuccess)
                        {
                            fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", gpu, cudaGetErrorString(err));
                            exit(-1);
                        }
                    }

                    cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

                    // float **progeny_Configuration_Cancer;
                    progeny_Elapsed = (float *)malloc(sizeof(float) * parent_Cells_Found * 2);
                    progeny_Sequences = (int **)malloc(parent_Cells_Found * 2 * sizeof(int *));

                    progeny_Configuration_Cancer = (float **)malloc(parent_Cells_Found * 2 * sizeof(float *));
                    for (int row = 0; row < (parent_Cells_Found * 2); row++)
                    {
                        progeny_Sequences[row] = (int *)malloc(genome_Length * sizeof(int));
                        progeny_Configuration_Cancer[row] = (float *)malloc(5 * sizeof(float));
                        // converted_Sequences.push_back("");
                    }

                    cout << "Receiving data from GPU: ";
                    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
                    {
                        cudaSetDevice(CUDA_device_IDs[gpu]);
                        int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

                        for (int row = 0; row < (cell_Count * 2); row++)
                        {
                            cudaMemcpy(progeny_Sequences[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Sequences_INT[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
                            cudaMemcpy(progeny_Configuration_Cancer[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float), cudaMemcpyDeviceToHost);
                        }

                        cudaMemcpy(progeny_Elapsed + (start_stop_Per_GPU[gpu].first * 2), cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float), cudaMemcpyDeviceToHost);
                    }

                    cout << "Data received by host\n";

                    for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
                    {
                        cudaSetDevice(CUDA_device_IDs[gpu]);

                        cudaFree(cuda_parents_Elapsed[gpu]);
                        cudaFree(cuda_progeny_Elapsed[gpu]);

                        int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

                        for (int row = 0; row < (cell_Count * 2); row++)
                        {
                            cudaFree(cuda_progeny_Configuration_Cancer[gpu][row]);
                            cudaFree(cuda_progeny_Sequences_INT[gpu][row]);
                        }
                        cudaFree(cuda_progeny_Configuration_Cancer[gpu]);
                        cudaFree(cuda_progeny_Sequences_INT[gpu]);

                        for (int row = 0; row < cell_Count; row++)
                        {
                            cudaFree(cuda_parent_sequences_INT[gpu][row]);
                        }
                        cudaFree(cuda_parent_sequences_INT[gpu]);

                        cudaStreamDestroy(streams[gpu]);
                    }

                    cout << " GPU(s) released\n";

                    for (int test = 0; test < parent_Cells_Found; test++)
                    {
                        cout << test << ": " << endl;
                        for (int progeny = (test * 2); progeny < ((test * 2) + 2); progeny++)
                        {
                            cout << progeny_Sequences[progeny][0];
                            cout << " ";
                            for (int col = 0; col < 5; col++)
                            {
                                cout << progeny_Configuration_Cancer[progeny][col] << " ";
                            }
                            cout << "| " << progeny_Elapsed[progeny];
                            cout << endl;
                        }
                    }

                    rerun_Progeny = compile_Progeny(functions,
                                                    intermediary_Tissue_folder, rapid_Progeny_Location,
                                                    parent_Cells_Found,
                                                    progeny_Elapsed, progeny_Configuration_Cancer,
                                                    last_index_Seq_Written, overall_Generations, tissue_Name,
                                                    progeny_Sequences, tissue,
                                                    gen, parent_IDs,
                                                    source_sequence_Data_folder,
                                                    last_Progeny_written_this_Gen,
                                                    tissue_Migration_Total, migration_cell_List, viral_Migration);

                    for (int row = 0; row < parent_Cells_Found * 2; row++)
                    {
                        free(progeny_Configuration_Cancer[row]);
                    }
                    free(progeny_Configuration_Cancer);

                    rounds_reRun++;

                    full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
                    // full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

                    for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue].size(); forward++)
                    {
                        full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue][forward]);
                    }

                } while (rerun_Progeny.size() > 0);

                for (int row = 0; row < parent_Cells_Found * 2; row++)
                {
                    free(progeny_Sequences[row]);
                }
                free(progeny_Sequences);
                free(progeny_Elapsed);

                // remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);

                cout << "Reruns completed after: " << rounds_reRun << " rounds\n";

                // exit(-1);
                // rerun_Progeny_THIS_gen(functions, rerun_Progeny, start_stop_Per_GPU, num_Cuda_devices, parent_Cells_Found,
                //                        progeny_Sequences, progeny_Elapsed,
                //                        cuda_Reference_fitness_survivability_proof_reading,
                //                        cuda_Reference_cancer_parameters,
                //                        cuda_sequence_replication_factor_changes,
                //                        mutation_Hotspots, cuda_mutation_hotspot_parameters,
                //                        cuda_A_0_mutation,
                //                        cuda_T_1_mutation,
                //                        cuda_G_2_mutation,
                //                        cuda_C_3_mutation,
                //                        cuda_num_effect_Segregating_sites,
                //                        cuda_num_effect_Segregating_sites_Cancer,
                //                        cuda_sequence_Survivability_changes,
                //                        cuda_sequence_Proof_reading_changes,
                //                        cuda_sequence_mutation_rate_changes,
                //                        cuda_sequence_generation_death_changes,
                //                        cuda_sequence_replication_prob_changes,
                //                        cuda_sequence_metastatic_prob_changes,
                //                        CUDA_device_IDs,
                //                        intermediary_Tissue_folder, rapid_Progeny_Location,
                //                        last_index_Seq_Written, overall_Generations, tissue_Name,
                //                        tissue, gen, source_sequence_Data_folder, last_Progeny_written_this_Gen,
                //                        max_sequences_per_File,
                //                        tissue_Migration_Total, migration_cell_List);
            }
            else
            {
                cout << "Clearing arrays\n";
                // remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions);
                // cout << "Check 1\n";
                // functions.clear_Array_int_CPU(progeny_Sequences, parent_Cells_Found * 2);
                for (int row = 0; row < parent_Cells_Found * 2; row++)
                {
                    free(progeny_Sequences[row]);
                }
                free(progeny_Sequences);
                // cout << "Check 2\n";
                free(progeny_Elapsed);
                // cout << "Check 3\n";
            }

            full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
            // full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

            for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue].size(); forward++)
            {
                full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue][forward]);
            }

            cout << "Complete termination of GPU streams: ";
            for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
            {
                cudaSetDevice(CUDA_device_IDs[gpu]);

                cudaFree(cuda_Reference_fitness_survivability_proof_reading[gpu]);
                cudaFree(cuda_Reference_cancer_parameters[gpu]);

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

                cudaFree(cuda_num_effect_Segregating_sites[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
                {
                    cudaFree(cuda_sequence_Survivability_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_Survivability_changes[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
                {
                    cudaFree(cuda_sequence_Proof_reading_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_Proof_reading_changes[gpu]);

                cudaFree(cuda_num_effect_Segregating_sites_Cancer[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
                {
                    cudaFree(cuda_sequence_replication_factor_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_replication_factor_changes[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites_Cancer[1]; row++)
                {
                    cudaFree(cuda_sequence_mutation_rate_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_mutation_rate_changes[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites_Cancer[2]; row++)
                {
                    cudaFree(cuda_sequence_generation_death_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_generation_death_changes[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites_Cancer[3]; row++)
                {
                    cudaFree(cuda_sequence_replication_prob_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_replication_prob_changes[gpu]);

                for (int row = 0; row < num_effect_Segregating_sites_Cancer[4]; row++)
                {
                    cudaFree(cuda_sequence_metastatic_prob_changes[gpu][row]);
                }
                cudaFree(cuda_sequence_metastatic_prob_changes[gpu]);

                cudaFree(cuda_tissue_selection_Position_Count[gpu]);

                cudaFree(cuda_Survivability_Positions[gpu]);
                cudaFree(cuda_Proof_Positions[gpu]);
                cudaFree(cuda_Replication_factor_Positions[gpu]);
                cudaFree(cuda_Mutation_rate_factor_Positions[gpu]);
                cudaFree(cuda_Generation_death_Positions[gpu]);
                cudaFree(cuda_Replication_prob_Positions[gpu]);
                cudaFree(cuda_Metastatic_Positions[gpu]);

                if (tissue_selection_Position_Count[0] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Survivability[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Survivability[gpu]);

                if (tissue_selection_Position_Count[1] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Proof[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Proof[gpu]);

                if (tissue_selection_Position_Count[2] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Replication_factor[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Replication_factor[gpu]);

                if (tissue_selection_Position_Count[3] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Mutation_rate_factor[gpu]);

                if (tissue_selection_Position_Count[4] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Generation_death[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Generation_death[gpu]);

                if (tissue_selection_Position_Count[5] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Replication_prob[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Replication_prob[gpu]);

                if (tissue_selection_Position_Count[6] > 0)
                {
                    for (int row = 0; row < (4 * num_Tissues); row++)
                    {
                        cudaFree(cuda_tissues_ATGC_positions_Metastatic[gpu][row]);
                    }
                }
                cudaFree(cuda_tissues_ATGC_positions_Metastatic[gpu]);
            }

            cout << "Round Complete\n";

            //  exit(-1);
        }
        else
        {
            cout << "CRITICAL ERROR: SEQUENCES ARE BLANK\n";
            // exit(-1);
        }
    }
    else
    {
        cout << "\nNo parents undergoing mitosis\n";
    }

    // exit(-1);
}

vector<pair<int, int>> cancer_Host::compile_Progeny(functions_library &functions,
                                                    string &intermediary_Tissue_folder, string &rapid_Progeny_Location,
                                                    int &parent_Cells_Found,
                                                    float *progeny_Elapsed, float **progeny_Configuration_Cancer,
                                                    int &last_index_Seq_Written, int &overall_Generations, string &tissue_Name,
                                                    int **progeny_Sequences, int &tissue_Index,
                                                    mt19937 &gen, vector<int> &parent_IDs,
                                                    string &source_sequence_Data_folder,
                                                    int &last_Progeny_written_this_Gen,
                                                    int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List, string &viral_Migration)
{
    vector<pair<int, int>> rerun_Progeny;

    uniform_real_distribution<float> check_Survival_Dis(0.0, 1.0);
    int check_Survival = 0;

    fstream sequence_Profiles_File;
    sequence_Profiles_File.open(sequence_Profiles, ios::app);
    fstream sequence_parent_Progeny_relationships_File;
    sequence_parent_Progeny_relationships_File.open(sequence_parent_Progeny_relationships, ios::app);

    // string intermediary_Tissue_folder = source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations + 1);
    string dead_List = intermediary_Tissue_folder + "/dead_List.txt";

    fstream dead_List_File;
    dead_List_File.open(dead_List, ios::app);

    fstream rapid_Progeny_File;
    if (filesystem::exists(rapid_Progeny_Location))
    {
        rapid_Progeny_File.open(rapid_Progeny_Location, ios::app);
    }
    cout << "\nCompiling progeny: ";

    for (int parent = 0; parent < parent_Cells_Found; parent++)
    {
        for (int progeny_Index = (parent * 2); progeny_Index < ((parent * 2) + 2); progeny_Index++)
        {
            if (progeny_Elapsed[progeny_Index] >= 1)
            {
                int gen_Displacement = int(progeny_Elapsed[progeny_Index]);
                // get new_Elapsed time
                progeny_Elapsed[progeny_Index] = progeny_Elapsed[progeny_Index] - gen_Displacement;

                if (gen_Displacement == 1)
                {
                    check_Survival = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][4]) ? 0 : 1;
                    string survival_Status = "_A_";
                    if (check_Survival == 1)
                    {
                        survival_Status = "_D_";
                        dead_List_File << last_index_Seq_Written << endl;
                        dead_Particle_count[tissue_Index] = dead_Particle_count[tissue_Index] + 1;
                    }
                    sequence_Profiles_File << tissue_Name << "_" << to_string(overall_Generations + 1) << "_" << to_string(last_index_Seq_Written) << "\t" << tissue_Name;
                    for (int col = 0; col < 5; col++)
                    {
                        sequence_Profiles_File << "\t" << to_string(progeny_Configuration_Cancer[progeny_Index][col]);
                    }
                    sequence_Profiles_File << endl;
                    sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(overall_Generations) << "_" << parent_IDs[parent] << "\t"
                                                               << tissue_Name << "_" << to_string(overall_Generations + 1) << "_" << to_string(last_index_Seq_Written)
                                                               << "\tprimary_Parent" << endl;
                    string sequence = "";
                    for (int base = 0; base < genome_Length; base++)
                    {
                        sequence.append(to_string(progeny_Sequences[progeny_Index][base]));
                    }
                    to_write_Sequence_Store_NEXT_Gen.push_back(make_pair(to_string(last_index_Seq_Written) + survival_Status + to_string(progeny_Configuration_Cancer[progeny_Index][2]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][1]) + "_" + to_string(progeny_Elapsed[progeny_Index]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][0]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][3]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][4]), sequence));

                    if (viral_Migration == "YES")
                    {
                        int migration_Check = 0;
                        migration_Check = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][3]) ? 0 : 1;
                        if (migration_Check == 0)
                        {
                            cell_Migration_set(tissue_Migration_Total, migration_cell_List, make_pair(progeny_Configuration_Cancer[progeny_Index][3], last_index_Seq_Written));
                        }
                    }

                    last_index_Seq_Written++;
                    current_cell_load_per_Tissue[tissue_Index] = current_cell_load_per_Tissue[tissue_Index] + 1;
                }
                else
                {
                    int new_Genenerations = gen_Displacement + overall_Generations;

                    string OTHER_intermediary_Tissue_folder = source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(new_Genenerations);
                    string OTHER_dead_List = OTHER_intermediary_Tissue_folder + "/dead_List.txt";

                    if (to_write_Sequence_Store_OTHER_Gens[tissue_Index].size() < (new_Genenerations + 1))
                    {
                        int create_Stores = (new_Genenerations + 1) - to_write_Sequence_Store_OTHER_Gens[tissue_Index].size();
                        for (int gen = 0; gen < create_Stores; gen++)
                        {
                            vector<pair<string, string>> store_Initialize;
                            to_write_Sequence_Store_OTHER_Gens[tissue_Index].push_back(store_Initialize);
                            last_index_Seq_Written_OTHERs[tissue_Index].push_back(0);
                        }
                    }

                    if (!filesystem::exists(OTHER_intermediary_Tissue_folder))
                    {
                        functions.config_Folder(OTHER_intermediary_Tissue_folder, to_string(new_Genenerations) + " generation Tissue " + tissue_Name + " sequences");
                        functions.create_File(OTHER_dead_List);
                    }

                    check_Survival = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][4]) ? 0 : 1;
                    string survival_Status = "_A_";
                    if (check_Survival == 1)
                    {
                        fstream dead_List_File_Others;
                        dead_List_File_Others.open(OTHER_dead_List, ios::app);

                        survival_Status = "_D_";
                        dead_List_File_Others << last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations] << endl;

                        dead_List_File_Others.close();
                    }

                    sequence_Profiles_File << tissue_Name << "_" << to_string(new_Genenerations) << "_" << to_string(last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations]) << "\t" << tissue_Name;
                    for (int col = 0; col < 5; col++)
                    {
                        sequence_Profiles_File << "\t" << to_string(progeny_Configuration_Cancer[progeny_Index][col]);
                    }
                    sequence_Profiles_File << endl;
                    sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(overall_Generations) << "_" << parent_IDs[parent] << "\t"
                                                               << tissue_Name << "_" << to_string(new_Genenerations) << "_" << to_string(last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations])
                                                               << "\tprimary_Parent" << endl;
                    string sequence = "";
                    for (int base = 0; base < genome_Length; base++)
                    {
                        sequence.append(to_string(progeny_Sequences[progeny_Index][base]));
                    }
                    to_write_Sequence_Store_OTHER_Gens[tissue_Index][new_Genenerations].push_back(make_pair(to_string(last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations]) + survival_Status + to_string(progeny_Configuration_Cancer[progeny_Index][2]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][1]) + "_" + to_string(progeny_Elapsed[progeny_Index]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][0]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][3]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][4]), sequence));

                    last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations]++;
                }
            }
            else if (progeny_Elapsed[progeny_Index] < 1)
            {
                float check_Rep_Time = progeny_Elapsed[progeny_Index] + progeny_Configuration_Cancer[progeny_Index][0];
                if (check_Rep_Time > 1)
                {
                    int gen_Displacement = (int)(progeny_Elapsed[progeny_Index] + progeny_Configuration_Cancer[progeny_Index][0]);
                    if (gen_Displacement == 1)
                    {
                        check_Survival = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][4]) ? 0 : 1;
                        string survival_Status = "_A_";
                        if (check_Survival == 1)
                        {
                            survival_Status = "_D_";
                            dead_List_File << last_index_Seq_Written << endl;
                            dead_Particle_count[tissue_Index] = dead_Particle_count[tissue_Index] + 1;
                        }
                        sequence_Profiles_File << tissue_Name << "_" << to_string(overall_Generations + 1) << "_" << to_string(last_index_Seq_Written) << "\t" << tissue_Name;
                        for (int col = 0; col < 5; col++)
                        {
                            sequence_Profiles_File << "\t" << to_string(progeny_Configuration_Cancer[progeny_Index][col]);
                        }
                        sequence_Profiles_File << endl;
                        sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(overall_Generations) << "_" << parent_IDs[parent] << "\t"
                                                                   << tissue_Name << "_" << to_string(overall_Generations + 1) << "_" << to_string(last_index_Seq_Written)
                                                                   << "\tprimary_Parent" << endl;
                        string sequence = "";
                        for (int base = 0; base < genome_Length; base++)
                        {
                            sequence.append(to_string(progeny_Sequences[progeny_Index][base]));
                        }
                        to_write_Sequence_Store_NEXT_Gen.push_back(make_pair(to_string(last_index_Seq_Written) + survival_Status + to_string(progeny_Configuration_Cancer[progeny_Index][2]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][1]) + "_" + to_string(progeny_Elapsed[progeny_Index] - 1), sequence));

                        if (viral_Migration == "YES")
                        {
                            int migration_Check = 0;
                            migration_Check = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][3]) ? 0 : 1;
                            if (migration_Check == 0)
                            {
                                cell_Migration_set(tissue_Migration_Total, migration_cell_List, make_pair(progeny_Configuration_Cancer[progeny_Index][3], last_index_Seq_Written));
                            }
                        }

                        last_index_Seq_Written++;
                        current_cell_load_per_Tissue[tissue_Index] = current_cell_load_per_Tissue[tissue_Index] + 1;
                    }
                    else
                    {
                        int new_Genenerations = gen_Displacement + overall_Generations;

                        string OTHER_intermediary_Tissue_folder = source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(new_Genenerations);
                        string OTHER_dead_List = OTHER_intermediary_Tissue_folder + "/dead_List.txt";

                        if (to_write_Sequence_Store_OTHER_Gens[tissue_Index].size() < (new_Genenerations + 1))
                        {
                            int create_Stores = (new_Genenerations + 1) - to_write_Sequence_Store_OTHER_Gens[tissue_Index].size();
                            for (int gen = 0; gen < create_Stores; gen++)
                            {
                                vector<pair<string, string>> store_Initialize;
                                to_write_Sequence_Store_OTHER_Gens[tissue_Index].push_back(store_Initialize);
                                last_index_Seq_Written_OTHERs[tissue_Index].push_back(0);
                            }
                        }

                        if (!filesystem::exists(OTHER_intermediary_Tissue_folder))
                        {
                            functions.config_Folder(OTHER_intermediary_Tissue_folder, to_string(new_Genenerations) + " generation Tissue " + tissue_Name + " sequences");
                            functions.create_File(OTHER_dead_List);
                        }

                        check_Survival = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][4]) ? 0 : 1;
                        string survival_Status = "_A_";
                        if (check_Survival == 1)
                        {
                            fstream dead_List_File_Others;
                            dead_List_File_Others.open(OTHER_dead_List, ios::app);

                            survival_Status = "_D_";
                            dead_List_File_Others << last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations] << endl;

                            dead_List_File_Others.close();
                        }

                        sequence_Profiles_File << tissue_Name << "_" << to_string(new_Genenerations) << "_" << to_string(last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations]) << "\t" << tissue_Name;
                        for (int col = 0; col < 5; col++)
                        {
                            sequence_Profiles_File << "\t" << to_string(progeny_Configuration_Cancer[progeny_Index][col]);
                        }
                        sequence_Profiles_File << endl;
                        sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(overall_Generations) << "_" << parent_IDs[parent] << "\t"
                                                                   << tissue_Name << "_" << to_string(new_Genenerations) << "_" << to_string(last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations])
                                                                   << "\tprimary_Parent" << endl;
                        string sequence = "";
                        for (int base = 0; base < genome_Length; base++)
                        {
                            sequence.append(to_string(progeny_Sequences[progeny_Index][base]));
                        }
                        to_write_Sequence_Store_OTHER_Gens[tissue_Index][new_Genenerations].push_back(make_pair(to_string(last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations]) + survival_Status + to_string(progeny_Configuration_Cancer[progeny_Index][2]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][1]) + "_" + to_string(progeny_Elapsed[progeny_Index] - 1) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][0]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][3]) + "_" + to_string(progeny_Configuration_Cancer[progeny_Index][4]), sequence));

                        last_index_Seq_Written_OTHERs[tissue_Index][new_Genenerations]++;
                    }
                }
                else
                {
                    sequence_Profiles_File << tissue_Name << "_" << to_string(overall_Generations) << "_" << to_string(last_Progeny_written_this_Gen) << "\t" << tissue_Name;
                    for (int col = 0; col < 5; col++)
                    {
                        sequence_Profiles_File << "\t" << to_string(progeny_Configuration_Cancer[progeny_Index][col]);
                    }
                    sequence_Profiles_File << endl;
                    sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(overall_Generations) << "_" << parent_IDs[parent] << "\t"
                                                               << tissue_Name << "_" << to_string(overall_Generations) << "_" << to_string(last_Progeny_written_this_Gen)
                                                               << "\tprimary_Parent" << endl;
                    string survival_Status = "_A_";
                    check_Survival = (check_Survival_Dis(gen) < progeny_Configuration_Cancer[progeny_Index][4]) ? 0 : 1;
                    if (check_Survival == 1)
                    {
                        fstream dead_List_File_this_Gen;
                        survival_Status = "_D_";
                        dead_List_File_this_Gen.open(source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(overall_Generations) + "/dead_List.txt", ios::app);
                        dead_List_File_this_Gen << last_Progeny_written_this_Gen << endl;
                        dead_List_File_this_Gen.close();
                    }
                    else
                    {
                        rerun_Progeny.push_back(make_pair(progeny_Index, last_Progeny_written_this_Gen));
                    }

                    if (!filesystem::exists(rapid_Progeny_Location))
                    {
                        functions.create_File(rapid_Progeny_Location);
                        rapid_Progeny_File.open(rapid_Progeny_Location, ios::app);
                    }

                    rapid_Progeny_File << ">" << to_string(last_Progeny_written_this_Gen) << survival_Status << to_string(progeny_Configuration_Cancer[progeny_Index][2]) << "_" << to_string(progeny_Configuration_Cancer[progeny_Index][1]) << "_" << to_string(progeny_Elapsed[progeny_Index]) << "_" << to_string(progeny_Configuration_Cancer[progeny_Index][0]) << "_" << to_string(progeny_Configuration_Cancer[progeny_Index][3]) << "_" << to_string(progeny_Configuration_Cancer[progeny_Index][4]) << endl;
                    string sequence = "";
                    for (int base = 0; base < genome_Length; base++)
                    {
                        sequence.append(to_string(progeny_Sequences[progeny_Index][base]));
                    }
                    rapid_Progeny_File << sequence << endl;
                    last_Progeny_written_this_Gen++;
                }
            }
        }
    }

    if (filesystem::exists(rapid_Progeny_Location))
    {
        rapid_Progeny_File.close();
    }
    sequence_Profiles_File.close();
    sequence_parent_Progeny_relationships_File.close();
    dead_List_File.close();

    cout << "Completed\n";

    return rerun_Progeny;
}

// void cancer_Host::rerun_Progeny_THIS_gen(functions_library &functions, vector<pair<int, int>> &rerun_Progeny, vector<pair<int, int>> &start_stop_Per_GPU, int &num_Cuda_devices, int &tot_Parents,
//                                          int **parent_sequences_INT, float *parents_Elapsed_full,
//                                          float *cuda_Reference_fitness_survivability_proof_reading[],
//                                          float *cuda_Reference_cancer_parameters[],
//                                          float **cuda_sequence_replication_factor_changes[],
//                                          int mutation_Hotspots, float **cuda_mutation_hotspot_parameters[],
//                                          float **cuda_A_0_mutation[],
//                                          float **cuda_T_1_mutation[],
//                                          float **cuda_G_2_mutation[],
//                                          float **cuda_C_3_mutation[],
//                                          int *cuda_num_effect_Segregating_sites[],
//                                          int *cuda_num_effect_Segregating_sites_Cancer[],
//                                          float **cuda_sequence_Survivability_changes[],
//                                          float **cuda_sequence_Proof_reading_changes[],
//                                          float **cuda_sequence_mutation_rate_changes[],
//                                          float **cuda_sequence_generation_death_changes[],
//                                          float **cuda_sequence_replication_prob_changes[],
//                                          float **cuda_sequence_metastatic_prob_changes[],
//                                          int *CUDA_device_IDs,
//                                          string &intermediary_Tissue_folder, string &rapid_Progeny_Location,
//                                          int &last_index_Seq_Written, int &overall_Generations, string &tissue_Name,
//                                          int &tissue_Index, mt19937 &gen, string &source_sequence_Data_folder, int &last_Progeny_written_this_Gen,
//                                          int &max_sequences_per_File,
//                                          int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List)
// {
//     cout << "\nProcessing reRun progeny\n";
//     int rounds_reRun = 0;
//     do
//     {
//         // CONVERT THE TOT cuda_parent_sequences_INT PARENTS TO ONE SO NOT SHARED BETWEEN GPUs
//         cout << "reRun Round: " << (rounds_reRun + 1) << endl;
//         int parent_Cells_Found = rerun_Progeny.size();

//         start_stop_Per_GPU.clear();

//         int standard_num_per_GPU = parent_Cells_Found / num_Cuda_devices;
//         int remainder = parent_Cells_Found % num_Cuda_devices;

//         for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
//         {
//             int start = gpu * standard_num_per_GPU;
//             int stop_GPU = start + standard_num_per_GPU;

//             start_stop_Per_GPU.push_back(make_pair(start, stop_GPU));
//         }

//         start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

//         cout << "Configuring datapoints: ";
//         // int *rerun_Progeny_Indexes = (int *)malloc(sizeof(int) * parent_Cells_Found);
//         float *parents_Elapsed = (float *)malloc(sizeof(float) * parent_Cells_Found);
//         vector<int> parent_IDs;
//         for (int parent = 0; parent < parent_Cells_Found; parent++)
//         {
//             // rerun_Progeny_Indexes[parent] = rerun_Progeny[parent];
//             parents_Elapsed[parent] = parents_Elapsed_full[rerun_Progeny[parent].first];
//             parent_IDs.push_back(rerun_Progeny[parent].second);
//         }
//         rerun_Progeny.clear();
//         free(parents_Elapsed_full);

//         cout << "Done\n";

//         // exit(-1);

//         // int *cuda_rerun_Progeny_Indexes[num_Cuda_devices];
//         int **cuda_progeny_Sequences_INT[num_Cuda_devices];
//         float *cuda_progeny_Elapsed[num_Cuda_devices];
//         float **cuda_progeny_Configuration_Cancer[num_Cuda_devices];

//         int **cuda_parent_sequences_INT[num_Cuda_devices];
//         float *cuda_parents_Elapsed[num_Cuda_devices];

//         cudaStream_t streams[num_Cuda_devices];

//         for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
//         {
//             cudaDeviceProp deviceProp;
//             cudaSetDevice(CUDA_device_IDs[gpu]);
//             cudaGetDeviceProperties(&deviceProp, gpu);

//             cout << "Intializing GPU " << CUDA_device_IDs[gpu] << "'s stream: " << deviceProp.name << endl;

//             int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

//             cudaMallocManaged(&cuda_parent_sequences_INT[gpu], cell_Count * sizeof(int *));
//             for (int row = 0; row < cell_Count; row++)
//             {
//                 cudaMalloc((void **)&cuda_parent_sequences_INT[gpu][row], genome_Length * sizeof(int));
//                 cudaMemcpy(cuda_parent_sequences_INT[gpu][row], parent_sequences_INT[rerun_Progeny[row + start_stop_Per_GPU[gpu].first].first], genome_Length * sizeof(int), cudaMemcpyHostToDevice);
//             }

//             cudaMalloc(&cuda_parents_Elapsed[gpu], cell_Count * sizeof(float));
//             cudaMemcpy(cuda_parents_Elapsed[gpu], parents_Elapsed + start_stop_Per_GPU[gpu].first, cell_Count * sizeof(float), cudaMemcpyHostToDevice);

//             // cudaMalloc(&cuda_rerun_Progeny_Indexes[gpu], cell_Count * sizeof(int));
//             // cudaMemcpy(cuda_rerun_Progeny_Indexes[gpu], rerun_Progeny_Indexes + start_stop_Per_GPU[gpu].first, cell_Count * sizeof(int), cudaMemcpyHostToDevice);

//             cudaMallocManaged(&cuda_progeny_Sequences_INT[gpu], (cell_Count * 2) * sizeof(int *));
//             for (int row = 0; row < (cell_Count * 2); row++)
//             {
//                 cudaMalloc((void **)&cuda_progeny_Sequences_INT[gpu][row], genome_Length * sizeof(int));
//             }

//             cudaMallocManaged(&cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float));

//             cudaMallocManaged(&cuda_progeny_Configuration_Cancer[gpu], (cell_Count * 2) * sizeof(float *));
//             for (int row = 0; row < (cell_Count * 2); row++)
//             {
//                 cudaMalloc((void **)&cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float));
//             }

//             cudaStreamCreate(&streams[gpu]);
//         }

//         // functions.clear_Array_int_CPU(parent_sequences_INT, tot_Parents * 2);
//         for (int row = 0; row < tot_Parents * 2; row++)
//         {
//             free(parent_sequences_INT[row]);
//         }
//         free(parent_sequences_INT);

//         tot_Parents = parent_Cells_Found;

//         cout << "GPUs intialized\n";

//         for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
//         {
//             cudaSetDevice(CUDA_device_IDs[gpu]);

//             cuda_replicate_Progeny_reRun<<<functions.tot_Blocks_array[gpu], functions.tot_ThreadsperBlock_array[gpu], 0, streams[gpu]>>>(start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first,
//                                                                                                                                          cuda_parent_sequences_INT[gpu], cuda_progeny_Sequences_INT[gpu], genome_Length,
//                                                                                                                                          cuda_parents_Elapsed[gpu], cuda_progeny_Elapsed[gpu],
//                                                                                                                                          cuda_Reference_fitness_survivability_proof_reading[gpu], cuda_Reference_cancer_parameters[gpu],
//                                                                                                                                          cuda_sequence_replication_factor_changes[gpu],
//                                                                                                                                          mutation_Hotspots,
//                                                                                                                                          cuda_mutation_hotspot_parameters[gpu],
//                                                                                                                                          cuda_A_0_mutation[gpu],
//                                                                                                                                          cuda_T_1_mutation[gpu],
//                                                                                                                                          cuda_G_2_mutation[gpu],
//                                                                                                                                          cuda_C_3_mutation[gpu],
//                                                                                                                                          cuda_num_effect_Segregating_sites[gpu],
//                                                                                                                                          cuda_num_effect_Segregating_sites_Cancer[gpu],
//                                                                                                                                          cuda_sequence_Survivability_changes[gpu],
//                                                                                                                                          cuda_sequence_Proof_reading_changes[gpu],
//                                                                                                                                          cuda_sequence_mutation_rate_changes[gpu],
//                                                                                                                                          cuda_sequence_generation_death_changes[gpu],
//                                                                                                                                          cuda_sequence_replication_prob_changes[gpu],
//                                                                                                                                          cuda_sequence_metastatic_prob_changes[gpu],
//                                                                                                                                          cuda_progeny_Configuration_Cancer[gpu],);
//         }

//         for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
//         {
//             cudaSetDevice(CUDA_device_IDs[gpu]);
//             cudaStreamSynchronize(streams[gpu]);

//             cudaError_t err = cudaGetLastError();
//             if (err != cudaSuccess)
//             {
//                 fprintf(stderr, "ERROR: CUDA error after synchronizing stream on GPU %d: %s\n", gpu, cudaGetErrorString(err));
//                 exit(-1);
//             }
//         }

//         cout << "GPU(s) streams completed and synchronized\nCopying data from GPU to Host memory\n";

//         float **progeny_Configuration_Cancer;
//         parents_Elapsed_full = (float *)malloc(sizeof(float) * parent_Cells_Found * 2);
//         parent_sequences_INT = (int **)malloc(parent_Cells_Found * 2 * sizeof(int *));

//         progeny_Configuration_Cancer = (float **)malloc(parent_Cells_Found * 2 * sizeof(float *));
//         for (int row = 0; row < (parent_Cells_Found * 2); row++)
//         {
//             parent_sequences_INT[row] = (int *)malloc(genome_Length * sizeof(int));
//             progeny_Configuration_Cancer[row] = (float *)malloc(5 * sizeof(float));
//             // converted_Sequences.push_back("");
//         }

//         cout << "Receiving data from GPU: ";
//         for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
//         {
//             cudaSetDevice(CUDA_device_IDs[gpu]);
//             int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

//             for (int row = 0; row < (cell_Count * 2); row++)
//             {
//                 cudaMemcpy(parent_sequences_INT[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Sequences_INT[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
//                 cudaMemcpy(progeny_Configuration_Cancer[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float), cudaMemcpyDeviceToHost);
//             }

//             cudaMemcpy(parents_Elapsed_full + (start_stop_Per_GPU[gpu].first * 2), cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float), cudaMemcpyDeviceToHost);
//         }

//         cout << "Data received by host\n";

//         cout << "Partial termination of GPU streams: ";
//         for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
//         {
//             cudaSetDevice(CUDA_device_IDs[gpu]);

//             cudaFree(cuda_parents_Elapsed[gpu]);
//             cudaFree(cuda_progeny_Elapsed[gpu]);

//             int cell_Count = start_stop_Per_GPU[gpu].second - start_stop_Per_GPU[gpu].first;

//             for (int row = 0; row < (cell_Count * 2); row++)
//             {
//                 cudaFree(cuda_progeny_Configuration_Cancer[gpu][row]);
//                 cudaFree(cuda_progeny_Sequences_INT[gpu][row]);
//             }
//             cudaFree(cuda_progeny_Configuration_Cancer[gpu]);
//             cudaFree(cuda_progeny_Sequences_INT[gpu]);

//             for (int row = 0; row < cell_Count; row++)
//             {
//                 cudaFree(cuda_parent_sequences_INT[gpu][row]);
//             }
//             cudaFree(cuda_parent_sequences_INT[gpu]);

//             cudaStreamDestroy(streams[gpu]);
//         }

//         cout << " GPU(s) released\n";

//         for (int test = 0; test < parent_Cells_Found; test++)
//         {
//             cout << test << ": " << endl;
//             for (int progeny = (test * 2); progeny < ((test * 2) + 2); progeny++)
//             {
//                 cout << parent_sequences_INT[progeny][0];
//                 cout << " ";
//                 for (int col = 0; col < 5; col++)
//                 {
//                     cout << progeny_Configuration_Cancer[progeny][col] << " ";
//                 }
//                 cout << "| " << parents_Elapsed_full[progeny];
//                 cout << endl;
//             }
//         }

//         rerun_Progeny = compile_Progeny(functions,
//                                         intermediary_Tissue_folder, rapid_Progeny_Location,
//                                         parent_Cells_Found,
//                                         parents_Elapsed_full, progeny_Configuration_Cancer,
//                                         last_index_Seq_Written, overall_Generations, tissue_Name,
//                                         parent_sequences_INT, tissue_Index,
//                                         gen, parent_IDs,
//                                         source_sequence_Data_folder,
//                                         last_Progeny_written_this_Gen,
//                                         tissue_Migration_Total, migration_cell_List);

//         functions.clear_Array_float_CPU(progeny_Configuration_Cancer, parent_Cells_Found * 2);

//         rounds_reRun++;

//         full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
//         // full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

//         for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue_Index].size(); forward++)
//         {
//             full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue_Index][forward]);
//         }

//     } while (rerun_Progeny.size() > 0);

//     // functions.clear_Array_int_CPU(parent_sequences_INT, tot_Parents * 2);
//     for (int row = 0; row < tot_Parents * 2; row++)
//     {
//         free(parent_sequences_INT[row]);
//     }
//     free(parent_sequences_INT);
//     free(parents_Elapsed_full);

//     cout << "Reruns completed after: " << rounds_reRun << " rounds\n";
// }

// void cancer_Host::thread_Sequence_to_String_Cancer(int start, int stop, int **progeny_Sequences)
// {
//     vector<string> converted_Sequences_Store;

//     for (int progeny = start; progeny < stop; progeny++)
//     {
//         string sequence = "";
//         for (int base = 0; base < genome_Length; base++)
//         {
//             sequence.append(to_string(progeny_Sequences[progeny][base]));
//         }
//         converted_Sequences_Store.push_back(sequence);
//     }

//     unique_lock<shared_mutex> ul(g_mutex);
//     int index = 0;
//     for (int progeny = start; progeny < stop; progeny++)
//     {
//         converted_Sequences[progeny] = converted_Sequences_Store[index];
//         index++;
//     }
// }

string cancer_Host::find_Sequences_Master(int &offset, int &tissue, string &tissue_Name, functions_library &functions, string &folder_Path, vector<int> &parents_in_Tissue, int &num_Sequences, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation, vector<int> &parent_IDs, float *parents_Elapsed, int &last_index_Seq_Written, mt19937 &gen,
                                          int &tissue_Migration_Total, multiset<pair<float, int>> &migration_cell_List,
                                          string &viral_Migration)
{

    cout << "Master collecting " << num_Sequences << " sequence(s)\n";
    // string folder_Path = source_Target_file_Location + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation);

    // int num_per_Core = num_Sequences / this->CPU_cores;
    // int remainder = num_Sequences % this->CPU_cores;

    // vector<thread> threads_vec;

    // for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
    // {
    //     int start_Cell = core_ID * num_per_Core;
    //     int stop_Cell = start_Cell + num_per_Core;

    //     threads_vec.push_back(thread{&cancer_Host::thread_find_Files, this, offset, start_Cell, stop_Cell, ref(parents_in_Tissue), ref(indexed_Tissue_Folder)});
    // }

    // if (remainder != 0)
    // {
    //     int start_Cell = num_Sequences - remainder;
    //     int stop_Cell = num_Sequences;

    //     threads_vec.push_back(thread{&cancer_Host::thread_find_Files, this, offset, start_Cell, stop_Cell, ref(parents_in_Tissue), ref(indexed_Tissue_Folder)});
    // }

    // for (thread &t : threads_vec)
    // {
    //     if (t.joinable())
    //     {
    //         t.join();
    //     }
    // }

    // threads_vec.clear();

    // vector<int> Tissue_files(found_Tissue_Folder_Indexes.begin(), found_Tissue_Folder_Indexes.end());
    // found_Tissue_Folder_Indexes.clear();

    // sort(Tissue_files.begin(), Tissue_files.end());

    // cout << Tissue_files.size() << " file(s) identified\n";

    // exit(-1);

    // vector<pair<int, int>> sequence_FileIndex_Position_list;
    string all_Sequences = "";
    // vector<string> collected_Sequences;
    //  for (int index = 0; index < num_Sequences; index++)
    //  {
    //      sequence_FileIndex_Position_list.push_back(make_pair(parents_in_Tissue[index], index));
    //      collected_Sequences.push_back("");
    //  }

    // sort(sequence_FileIndex_Position_list.begin(), sequence_FileIndex_Position_list.end());

    fstream nfasta;
    int index_Files = 0;
    int line_current = 0;

    cout << "Retrieving sequence(s)\n";

    // valid_Sequences = 0;

    nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[index_Files].first) + "_" + to_string(indexed_Tissue_Folder[index_Files].second) + ".nfasta", ios::in);

    cout << "File: " << folder_Path << "/" << to_string(indexed_Tissue_Folder[index_Files].first) + "_" << to_string(indexed_Tissue_Folder[index_Files].second) << ".nfasta\n";

    uniform_real_distribution<float> check_Replicate_Dis(0.0, 1.0);

    fstream sequence_Profiles_File;
    sequence_Profiles_File.open(sequence_Profiles, ios::app);
    fstream sequence_parent_Progeny_relationships_File;
    sequence_parent_Progeny_relationships_File.open(sequence_parent_Progeny_relationships, ios::app);

    uniform_real_distribution<float> check_Survival_Dis(0.0, 1.0);

    for (int find = 0; find < num_Sequences; find++)
    {
        cout << "Looking for " << find << " :" << parents_in_Tissue[find + offset] << ": of " << num_Sequences << "\n";

        if (index_Files < indexed_Tissue_Folder.size())
        {

            while ((indexed_Tissue_Folder[index_Files].first <= parents_in_Tissue[find + offset] && indexed_Tissue_Folder[index_Files].second >= parents_in_Tissue[find + offset]) == 0)
            {
                nfasta.close();
                index_Files++;
                nfasta.open(folder_Path + "/" + to_string(indexed_Tissue_Folder[index_Files].first) + "_" + to_string(indexed_Tissue_Folder[index_Files].second) + ".nfasta", ios::in);

                cout << "File: " << folder_Path << "/" << to_string(indexed_Tissue_Folder[index_Files].first) + "_" << to_string(indexed_Tissue_Folder[index_Files].second) << ".nfast\n";

                line_current = 0;
            }

            if (nfasta.is_open())
            {
                int line_t0_check = (parents_in_Tissue[find + offset] - indexed_Tissue_Folder[index_Files].first) * 2;

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
                            all_Sequences.append(line);
                            // collected_Sequences.push_back(line);
                            parent_IDs.push_back(parents_in_Tissue[find + offset]);
                            cout << "Parent idetified\n";
                            parents_Elapsed[parent_IDs.size() - 1] = stof(line_Data[4]);
                            line_current++;
                        }
                        else if (stoi(line_Data[0].substr(1)) == parents_in_Tissue[find + offset] && check_replicate == 1)
                        {
                            cout << "Not replicating\n";
                            // check if it survives to the next generation
                            int check_Survival = 0;
                            if (stof(line_Data[3]) < 1)
                            {
                                float val = check_Replicate_Dis(gen);
                                check_Survival = (val < stof(line_Data[3])) ? 0 : 1;
                            }
                            if (check_Survival == 1)
                            {
                                cout << "Moving\n";
                                getline(nfasta, line);
                                cout << "Line_got\n";
                                to_write_Sequence_Store_NEXT_Gen.push_back(make_pair(to_string(last_index_Seq_Written) + "_A_" + line_Data[2] + "_" + line_Data[3] + "_0_" + line_Data[5] + "_" + line_Data[6] + "_" + line_Data[7], line));
                                cout << "check\n";
                                sequence_Profiles_File << tissue_Name << "_" << to_string(current_Generation + 1) << "_" << to_string(last_index_Seq_Written) << "\t" << tissue_Name << "\t" << line_Data[5] << "\t" << line_Data[3] << "\t" << line_Data[2] << "\t" << line_Data[6] << "\t" << line_Data[7] << endl;
                                cout << "check\n";
                                sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(current_Generation) << "_" << parents_in_Tissue[find + offset] << "\t"
                                                                           << tissue_Name << "_" << to_string(current_Generation + 1) << "_" << to_string(last_index_Seq_Written)
                                                                           << "\tgeneration_Forward" << endl;
                                cout << "MOVED\n";
                                if (viral_Migration == "YES")
                                {
                                    int migration_Check = 0;
                                    migration_Check = (check_Survival_Dis(gen) < stof(line_Data[6])) ? 0 : 1;
                                    if (migration_Check == 0)
                                    {
                                        cell_Migration_set(tissue_Migration_Total, migration_cell_List, make_pair(stof(line_Data[6]), last_index_Seq_Written));
                                    }
                                }
                                cout << "check\n";
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
                                 << "File: " << folder_Path << "/" << indexed_Tissue_Folder[index_Files].first
                                 << "_" << indexed_Tissue_Folder[index_Files].second << ".nfasta" << endl;
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
                cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << folder_Path << "/" << indexed_Tissue_Folder[index_Files].first << "_" << indexed_Tissue_Folder[index_Files].second << ".nfasta" << endl;
                exit(-1);
            }
        }
        else
        {
            cout << "ERROR: INDEX FILES OVERREACHED: " << index_Files << " SHOULD BE LESS THAN: " << indexed_Tissue_Folder.size();
            exit(-1);
        }
    }

    nfasta.close();

    sequence_Profiles_File.close();
    sequence_parent_Progeny_relationships_File.close();

    // cout << valid_Sequences << " live sequence(s) collected\n";

    return all_Sequences;
}

void cancer_Host::thread_find_Files(int offset, int start, int stop, vector<int> &parents_in_Tissue, vector<pair<int, int>> &indexed_Tissue_Folder)
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
    vector<vector<string>> profile_Lines_Tissues;

    intialize_Tissues(host_Folder, tissue_Sequences, profile_Lines_Tissues, functions, current_Generation);

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

        vector<string> profile_Lines;

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
                        Sequence_IDs.push_back(line_Data[2] + "_" + line_Data[3] + "_" + line_Data[4] + "_" + line_Data[5] + "_" + line_Data[6] + "_" + line_Data[7]);
                        profile_Lines.push_back(line_Data[5] + "\t" + line_Data[3] + "\t" + line_Data[2] + "\t" + line_Data[6] + "\t" + line_Data[7]);
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
            profile_Lines_Tissues[tissue_Index].push_back(profile_Lines[sequence]);
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
                    functions.create_File(output_Node_location + "/cancer_Host/sequence_Profiles.csv", "Sequence_ID\tTissue\treplication_Factor\tgenerational_death_Prob\treplication_Prob\tmetastatic_Prob\tSurvivability");
                    functions.create_File(output_Node_location + "/cancer_Host/sequence_parent_Progeny_relationships.csv", "Source\tTarget\tType");
                }

                vector<pair<string, string>> sequence_Write_Store_All;
                int last_seq_Num = 0;
                sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[tissue], profile_Lines_Tissues[tissue],
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

void cancer_Host::sequence_Write_Configurator(vector<pair<string, string>> &sequence_Write_Store_All, vector<pair<string, string>> &sequence_Write_Store, vector<string> &profile_Lines,
                                              int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                              string sequence_Profiles_Location, string tissue, int current_Generation)
{
    fstream sequence_Profile;
    sequence_Profile.open(sequence_Profiles_Location, ios::app);
    int temp_Count = last_seq_Num;
    for (int sequence_Collect = 0; sequence_Collect < sequence_Write_Store.size(); sequence_Collect++)
    {
        sequence_Write_Store_All.push_back(make_pair(sequence_Write_Store[sequence_Collect].first, sequence_Write_Store[sequence_Collect].second));
        sequence_Profile << tissue << "_" << current_Generation << "_" << temp_Count << "\t" << tissue << "\t" << profile_Lines[sequence_Collect] << endl;
        temp_Count++;
    }
    sequence_Profile.close();
    sequence_Write_Store.clear();

    if (sequence_Write_Store_All.size() >= max_sequences_per_File)
    {
        int full_Write_Count = sequence_Write_Store_All.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {

            string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);

            //"Sequence_ID\tHost\tTissue"

            if (fasta_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_seq_Num << "_A_";
                    fasta_File << sequence_Write_Store_All[write_Seq].second;
                    fasta_File << endl;
                    fasta_File << sequence_Write_Store_All[write_Seq].first << endl;

                    last_seq_Num++;
                }

                fasta_File.close();
                // sequence_Profile.close();
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
        // fstream sequence_Profile;
        // sequence_Profile.open(sequence_Profiles_Location, ios::app);
        if (fasta_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < sequence_Write_Store_All.size(); write_Seq++)
            {
                fasta_File << ">" << last_seq_Num << "_A_";
                fasta_File << sequence_Write_Store_All[write_Seq].second;
                fasta_File << endl;
                fasta_File << sequence_Write_Store_All[write_Seq].first << endl;
                //      sequence_Profile << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue << endl;
                last_seq_Num++;
            }

            fasta_File.close();
            // sequence_Profile.close();
        }
        else
        {
            cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            exit(-1);
        }
        sequence_Write_Store_All.clear();
    }
}

void cancer_Host::intialize_Tissues(string &host_Folder, vector<vector<pair<string, string>>> &tissue_Sequences, vector<vector<string>> &profile_Lines_Tissues, functions_library &functions, int &current_Generation)
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

        vector<string> profile_Lines;
        profile_Lines_Tissues.push_back(profile_Lines);

        tissue_Sequences.push_back(tissue_Sequence);
        current_cell_load_per_Tissue[tissue] = 0;
        dead_Particle_count[tissue] = 0;

        parents_Prev_generation[tissue] = 0;

        vector<int> tissue_initalize;
        last_index_Seq_Written_OTHERs.push_back(tissue_initalize);

        vector<vector<pair<string, string>>> others_Initalize;
        to_write_Sequence_Store_OTHER_Gens.push_back(others_Initalize);

        removed_by_Transfer_Indexes.push_back(init_removed_by_Transfer_Indexes);
    }
}