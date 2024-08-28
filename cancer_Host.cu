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
                                       int &max_sequences_per_File)
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
                        // get last wrtten progeny number
                        int last_Progeny_written_this_Gen = indexed_Source_Folders[tissue][indexed_Source_Folders[tissue].size() - 1].second + 1;
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
                        // ! check if the new generation already exists and if so update;
                        string intermediary_Tissue_folder = source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations + 1);
                        string dead_List = intermediary_Tissue_folder + "/dead_List.txt";

                        // fstream this_Gen_progeny_parents;
                        functions.create_File(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) + "/" + to_string(last_Progeny_written_this_Gen) + "_rapid_Progeny.nFASTA");

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
                                                last_Progeny_written_this_Gen, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations) + "/" + to_string(last_Progeny_written_this_Gen) + "_rapid_Progeny.nFASTA");
                        }

                        remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
                        // remainder_Write_Sequences_NEXT_Generation(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

                        for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue].size(); forward++)
                        {
                            remainder_Write_Sequences_NEXT_Generation(source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue][forward]);
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
                                            float *cuda_parents_Elapsed, float *cuda_progeny_Elapsed)
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

                        num_Mutations = (int)((num_Mutations * mutation_Factor) + 0.5);
                    }

                    if (num_Mutations > 0 && cuda_Reference_fitness_survivability_proof_reading[2] != -1)
                    {
                        float proof_Reading = cuda_Reference_fitness_survivability_proof_reading[2];

                        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[2]; pos++)
                        {
                            proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] + 1];
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
            progeny_Configuration_Cancer[progeny_Index][0] = replication_Factor;

            float gen_Death_prob = cuda_Reference_cancer_parameters[2];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[2]; pos++)
            {
                gen_Death_prob = gen_Death_prob + cuda_sequence_generation_death_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_generation_death_changes[pos][0] - 1] + 1];
            }
            gen_Death_prob = (gen_Death_prob > 1) ? 1 : (gen_Death_prob < 0) ? 0
                                                                             : gen_Death_prob;
            progeny_Configuration_Cancer[progeny_Index][1] = gen_Death_prob;

            float replication_prob = cuda_Reference_cancer_parameters[3];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[3]; pos++)
            {
                replication_prob = replication_prob + cuda_sequence_replication_prob_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_prob_changes[pos][0] - 1] + 1];
            }
            replication_prob = (replication_prob > 1) ? 1 : (replication_prob < 0) ? 0
                                                                                   : replication_prob;
            progeny_Configuration_Cancer[progeny_Index][2] = replication_prob;

            float metatstatic_prob = cuda_Reference_cancer_parameters[4];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[4]; pos++)
            {
                metatstatic_prob = metatstatic_prob + cuda_sequence_metastatic_prob_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_metastatic_prob_changes[pos][0] - 1] + 1];
            }
            metatstatic_prob = (metatstatic_prob > 1) ? 1 : (metatstatic_prob < 0) ? 0
                                                                                   : metatstatic_prob;
            progeny_Configuration_Cancer[progeny_Index][3] = metatstatic_prob;

            float survivability = cuda_Reference_fitness_survivability_proof_reading[1];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites[1]; pos++)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Survivability_changes[pos][0] - 1] + 1];
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

void cancer_Host::simulate_cell_Round(functions_library &functions, string &multi_Read, int &num_Cuda_devices, int *CUDA_device_IDs,
                                      int &num_of_Cells, int &start, int &stop,
                                      int *parents_in_Tissue, int &tissue, string tissue_Name,
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
                                      int &last_Progeny_written_this_Gen, string rapid_Progeny_Location)
{

    sort(parents_in_Tissue + start, parents_in_Tissue + stop);

    vector<int> parent_IDs;
    // vector<float> parents_Elapsed;

    float *parents_Elapsed = (float *)malloc(sizeof(float) * num_Tissues);

    string all_Sequences = find_Sequences_Master(start, tissue, tissue_Name, functions, this_Gen_intermediary_Sequences, parents_in_Tissue, num_of_Cells, indexed_Tissue_Folder, overall_Generations, parent_IDs, parents_Elapsed, last_index_Seq_Written, gen);
    full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
    // remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions);

    // for (int test = 0; test < parent_IDs.size(); test++)
    // {
    //     cout << parents_Elapsed[test] << endl;
    // }

    // exit(-1);

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

            cudaStreamCreate(&streams[gpu]);
        }

        free(full_Char);

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
                                                                                                                                        cuda_parents_Elapsed[gpu], cuda_progeny_Elapsed[gpu]);
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
                                                               last_Progeny_written_this_Gen);

        full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
        //  full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

        for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue].size(); forward++)
        {
            full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue][forward]);
        }

        cout << "\nClearing arrays\n";
        functions.clear_Array_float_CPU(progeny_Configuration_Cancer, parent_Cells_Found * 2);

        if (rerun_Progeny.size() > 0)
        {
            cout << "\nRerun progeny: " << rerun_Progeny.size() << endl;
            rerun_Progeny_THIS_gen(functions, rerun_Progeny, start_stop_Per_GPU, num_Cuda_devices, parent_Cells_Found,
                                   progeny_Sequences, progeny_Elapsed,
                                   cuda_Reference_fitness_survivability_proof_reading,
                                   cuda_Reference_cancer_parameters,
                                   cuda_sequence_replication_factor_changes,
                                   mutation_Hotspots, cuda_mutation_hotspot_parameters,
                                   cuda_A_0_mutation,
                                   cuda_T_1_mutation,
                                   cuda_G_2_mutation,
                                   cuda_C_3_mutation,
                                   cuda_num_effect_Segregating_sites,
                                   cuda_num_effect_Segregating_sites_Cancer,
                                   cuda_sequence_Survivability_changes,
                                   cuda_sequence_Proof_reading_changes,
                                   cuda_sequence_mutation_rate_changes,
                                   cuda_sequence_generation_death_changes,
                                   cuda_sequence_replication_prob_changes,
                                   cuda_sequence_metastatic_prob_changes,
                                   CUDA_device_IDs,
                                   intermediary_Tissue_folder, rapid_Progeny_Location,
                                   last_index_Seq_Written, overall_Generations, tissue_Name,
                                   tissue, gen, source_sequence_Data_folder, last_Progeny_written_this_Gen,
                                   max_sequences_per_File);
        }
        else
        {
            cout << "Clearing arrays\n";
            // remainder_Write_Sequences_NEXT_Generation(intermediary_Tissue_folder, functions);
            // cout << "Check 1\n";
            functions.clear_Array_int_CPU(progeny_Sequences, parent_Cells_Found * 2);
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
        }

        cout << "Round Complete\n";

        //  exit(-1);
    }
    else
    {
        cout << "CRITICAL ERROR: SEQUENCES ARE BLANK\n";
        exit(-1);
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
                                                    int &last_Progeny_written_this_Gen)
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
    rapid_Progeny_File.open(rapid_Progeny_Location, ios::app);

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
                    last_index_Seq_Written++;
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
                        last_index_Seq_Written++;
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

    rapid_Progeny_File.close();
    sequence_Profiles_File.close();
    sequence_parent_Progeny_relationships_File.close();
    dead_List_File.close();

    cout << "Completed\n";

    return rerun_Progeny;
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
                                             float **progeny_Configuration_Cancer)
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

                        num_Mutations = (int)((num_Mutations * mutation_Factor) + 0.5);
                    }

                    if (num_Mutations > 0 && cuda_Reference_fitness_survivability_proof_reading[2] != -1)
                    {
                        float proof_Reading = cuda_Reference_fitness_survivability_proof_reading[2];

                        for (int pos = 0; pos < cuda_num_effect_Segregating_sites[2]; pos++)
                        {
                            proof_Reading = proof_Reading + cuda_sequence_Proof_reading_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Proof_reading_changes[pos][0] - 1] + 1];
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
            progeny_Configuration_Cancer[progeny_Index][0] = replication_Factor;

            float gen_Death_prob = cuda_Reference_cancer_parameters[2];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[2]; pos++)
            {
                gen_Death_prob = gen_Death_prob + cuda_sequence_generation_death_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_generation_death_changes[pos][0] - 1] + 1];
            }
            gen_Death_prob = (gen_Death_prob > 1) ? 1 : (gen_Death_prob < 0) ? 0
                                                                             : gen_Death_prob;
            progeny_Configuration_Cancer[progeny_Index][1] = gen_Death_prob;

            float replication_prob = cuda_Reference_cancer_parameters[3];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[3]; pos++)
            {
                replication_prob = replication_prob + cuda_sequence_replication_prob_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_replication_prob_changes[pos][0] - 1] + 1];
            }
            replication_prob = (replication_prob > 1) ? 1 : (replication_prob < 0) ? 0
                                                                                   : replication_prob;
            progeny_Configuration_Cancer[progeny_Index][2] = replication_prob;

            float metatstatic_prob = cuda_Reference_cancer_parameters[4];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites_Cancer[4]; pos++)
            {
                metatstatic_prob = metatstatic_prob + cuda_sequence_metastatic_prob_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_metastatic_prob_changes[pos][0] - 1] + 1];
            }
            metatstatic_prob = (metatstatic_prob > 1) ? 1 : (metatstatic_prob < 0) ? 0
                                                                                   : metatstatic_prob;
            progeny_Configuration_Cancer[progeny_Index][3] = metatstatic_prob;

            float survivability = cuda_Reference_fitness_survivability_proof_reading[1];
            for (int pos = 0; pos < cuda_num_effect_Segregating_sites[1]; pos++)
            {
                survivability = survivability + cuda_sequence_Survivability_changes[pos][cuda_progeny_sequences_INT[progeny_Index][(int)cuda_sequence_Survivability_changes[pos][0] - 1] + 1];
            }
            survivability = (survivability > 1) ? 1 : (survivability < 0) ? 0
                                                                          : survivability;
            progeny_Configuration_Cancer[progeny_Index][4] = survivability;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void cancer_Host::rerun_Progeny_THIS_gen(functions_library &functions, vector<pair<int, int>> &rerun_Progeny, vector<pair<int, int>> &start_stop_Per_GPU, int &num_Cuda_devices, int &tot_Parents,
                                         int **parent_sequences_INT, float *parents_Elapsed_full,
                                         float *cuda_Reference_fitness_survivability_proof_reading[],
                                         float *cuda_Reference_cancer_parameters[],
                                         float **cuda_sequence_replication_factor_changes[],
                                         int mutation_Hotspots, float **cuda_mutation_hotspot_parameters[],
                                         float **cuda_A_0_mutation[],
                                         float **cuda_T_1_mutation[],
                                         float **cuda_G_2_mutation[],
                                         float **cuda_C_3_mutation[],
                                         int *cuda_num_effect_Segregating_sites[],
                                         int *cuda_num_effect_Segregating_sites_Cancer[],
                                         float **cuda_sequence_Survivability_changes[],
                                         float **cuda_sequence_Proof_reading_changes[],
                                         float **cuda_sequence_mutation_rate_changes[],
                                         float **cuda_sequence_generation_death_changes[],
                                         float **cuda_sequence_replication_prob_changes[],
                                         float **cuda_sequence_metastatic_prob_changes[],
                                         int *CUDA_device_IDs,
                                         string &intermediary_Tissue_folder, string &rapid_Progeny_Location,
                                         int &last_index_Seq_Written, int &overall_Generations, string &tissue_Name,
                                         int &tissue_Index, mt19937 &gen, string &source_sequence_Data_folder, int &last_Progeny_written_this_Gen,
                                         int &max_sequences_per_File)
{
    cout << "\nProcessing reRun progeny\n";
    int rounds_reRun = 0;
    do
    {
        // CONVERT THE TOT cuda_parent_sequences_INT PARENTS TO ONE SO NOT SHARED BETWEEN GPUs
        cout << "reRun Round: " << (rounds_reRun + 1) << endl;
        int parent_Cells_Found = rerun_Progeny.size();

        start_stop_Per_GPU.clear();

        int standard_num_per_GPU = parent_Cells_Found / num_Cuda_devices;
        int remainder = parent_Cells_Found % num_Cuda_devices;

        for (int gpu = 0; gpu < num_Cuda_devices; gpu++)
        {
            int start = gpu * standard_num_per_GPU;
            int stop_GPU = start + standard_num_per_GPU;

            start_stop_Per_GPU.push_back(make_pair(start, stop_GPU));
        }

        start_stop_Per_GPU[num_Cuda_devices - 1].second = start_stop_Per_GPU[num_Cuda_devices - 1].second + remainder;

        cout << "Configuring datapoints: ";
        // int *rerun_Progeny_Indexes = (int *)malloc(sizeof(int) * parent_Cells_Found);
        float *parents_Elapsed = (float *)malloc(sizeof(float) * parent_Cells_Found);
        vector<int> parent_IDs;
        for (int parent = 0; parent < parent_Cells_Found; parent++)
        {
            // rerun_Progeny_Indexes[parent] = rerun_Progeny[parent];
            parents_Elapsed[parent] = parents_Elapsed_full[rerun_Progeny[parent].first];
            parent_IDs.push_back(rerun_Progeny[parent].second);
        }
        rerun_Progeny.clear();
        free(parents_Elapsed_full);

        cout << "Done\n";

        // exit(-1);

        // int *cuda_rerun_Progeny_Indexes[num_Cuda_devices];
        int **cuda_progeny_Sequences_INT[num_Cuda_devices];
        float *cuda_progeny_Elapsed[num_Cuda_devices];
        float **cuda_progeny_Configuration_Cancer[num_Cuda_devices];

        int **cuda_parent_sequences_INT[num_Cuda_devices];
        float *cuda_parents_Elapsed[num_Cuda_devices];

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
                cudaMemcpy(cuda_parent_sequences_INT[gpu][row], parent_sequences_INT[rerun_Progeny[row + start_stop_Per_GPU[gpu].first].first], genome_Length * sizeof(int), cudaMemcpyHostToDevice);
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

        functions.clear_Array_int_CPU(parent_sequences_INT, tot_Parents * 2);
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
                                                                                                                                         cuda_progeny_Configuration_Cancer[gpu]);
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

        float **progeny_Configuration_Cancer;
        parents_Elapsed_full = (float *)malloc(sizeof(float) * parent_Cells_Found * 2);
        parent_sequences_INT = (int **)malloc(parent_Cells_Found * 2 * sizeof(int *));

        progeny_Configuration_Cancer = (float **)malloc(parent_Cells_Found * 2 * sizeof(float *));
        for (int row = 0; row < (parent_Cells_Found * 2); row++)
        {
            parent_sequences_INT[row] = (int *)malloc(genome_Length * sizeof(int));
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
                cudaMemcpy(parent_sequences_INT[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Sequences_INT[gpu][row], genome_Length * sizeof(int), cudaMemcpyDeviceToHost);
                cudaMemcpy(progeny_Configuration_Cancer[(start_stop_Per_GPU[gpu].first * 2) + row], cuda_progeny_Configuration_Cancer[gpu][row], 5 * sizeof(float), cudaMemcpyDeviceToHost);
            }

            cudaMemcpy(parents_Elapsed_full + (start_stop_Per_GPU[gpu].first * 2), cuda_progeny_Elapsed[gpu], (cell_Count * 2) * sizeof(float), cudaMemcpyDeviceToHost);
        }

        cout << "Data received by host\n";

        cout << "Partial termination of GPU streams: ";
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
                cout << parent_sequences_INT[progeny][0];
                cout << " ";
                for (int col = 0; col < 5; col++)
                {
                    cout << progeny_Configuration_Cancer[progeny][col] << " ";
                }
                cout << "| " << parents_Elapsed_full[progeny];
                cout << endl;
            }
        }

        rerun_Progeny = compile_Progeny(functions,
                                        intermediary_Tissue_folder, rapid_Progeny_Location,
                                        parent_Cells_Found,
                                        parents_Elapsed_full, progeny_Configuration_Cancer,
                                        last_index_Seq_Written, overall_Generations, tissue_Name,
                                        parent_sequences_INT, tissue_Index,
                                        gen, parent_IDs,
                                        source_sequence_Data_folder,
                                        last_Progeny_written_this_Gen);

        functions.clear_Array_float_CPU(progeny_Configuration_Cancer, parent_Cells_Found * 2);

        rounds_reRun++;

        full_Write_Sequences_NEXT_Generation(max_sequences_per_File, intermediary_Tissue_folder, functions, to_write_Sequence_Store_NEXT_Gen);
        // full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(overall_Generations), functions, to_write_Sequence_Store_THIS_Gen);

        for (int forward = 0; forward < to_write_Sequence_Store_OTHER_Gens[tissue_Index].size(); forward++)
        {
            full_Write_Sequences_NEXT_Generation(max_sequences_per_File, source_sequence_Data_folder + "/" + to_string(tissue_Index) + "/generation_" + to_string(forward), functions, to_write_Sequence_Store_OTHER_Gens[tissue_Index][forward]);
        }

    } while (rerun_Progeny.size() > 0);

    functions.clear_Array_int_CPU(parent_sequences_INT, tot_Parents * 2);
    free(parents_Elapsed_full);

    cout << "Reruns completed after: " << rounds_reRun << " rounds\n";
}

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

string cancer_Host::find_Sequences_Master(int &offset, int &tissue, string &tissue_Name, functions_library &functions, string &folder_Path, int *parents_in_Tissue, int &num_Sequences, vector<pair<int, int>> &indexed_Tissue_Folder, int &current_Generation, vector<int> &parent_IDs, float *parents_Elapsed, int &last_index_Seq_Written, mt19937 &gen)
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
                        all_Sequences.append(line);
                        // collected_Sequences.push_back(line);
                        parent_IDs.push_back(parents_in_Tissue[find + offset]);
                        parents_Elapsed[parent_IDs.size() - 1] = stof(line_Data[4]);
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
                        if (check_Survival == 1)
                        {
                            getline(nfasta, line);
                            to_write_Sequence_Store_NEXT_Gen.push_back(make_pair(to_string(last_index_Seq_Written) + "_A_" + line_Data[2] + "_" + line_Data[3] + "_0_" + line_Data[5] + "_" + line_Data[6] + "_" + line_Data[7], line));
                            sequence_Profiles_File << tissue_Name << "_" << to_string(current_Generation + 1) << "_" << to_string(last_index_Seq_Written) << "\t" << tissue_Name << "\t" << line_Data[5] << "\t" << line_Data[3] << "\t" << line_Data[2] << "\t" << line_Data[6] << "\t" << line_Data[7] << endl;
                            sequence_parent_Progeny_relationships_File << tissue_Name << "_" << to_string(current_Generation) << "_" << parents_in_Tissue[find + offset] << "\t"
                                                                       << tissue_Name << "_" << to_string(current_Generation + 1) << "_" << to_string(last_index_Seq_Written)
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

    return all_Sequences;
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