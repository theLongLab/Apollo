#include "cancer.cuh"

cancer::cancer(string parameter_Master_Location)
{
    cout << "Initiating cancer simulation\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"CUDA Device IDs\"",
        "\"CPU cores\"",
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Multi read\"",
        "\"Nodes master profile\"",
        "\"Sequence master profile\"",
        "\"Intermediate Sequences per file\"",
        "\"Process cell rate\"",
        "\"Start date\"",
        "\"Enable folder management\"",
        "\"First infection\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    cout << "\nConfiguring folders:\n";

    intermediate_Folder_location = Parameters.get_STRING(found_Parameters[2]);
    output_Folder_location = Parameters.get_STRING(found_Parameters[3]);

    function.config_Folder(intermediate_Folder_location, "Intermediate");
    function.config_Folder(output_Folder_location, "Output");

    max_sequences_per_File = Parameters.get_INT(found_Parameters[7]);
    max_Cells_at_a_time = Parameters.get_INT(found_Parameters[8]);
    start_Date = Parameters.get_STRING(found_Parameters[9]);
    first_Infection = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[11]));

    if (max_sequences_per_File <= 0)
    {
        cout << "ERROR: Intermediate Sequences per file PARAMETER MUST BE GREATER THAN ZERO.\n";
        exit(-1);
    }

    if (max_Cells_at_a_time <= 0)
    {
        cout << "ERROR: Process cell rate PARAMETER MUST BE GREATER THAN ZERO.\n";
        exit(-1);
    }

    cout << "\nFolder management: ";

    if (function.to_Upper_Case(Parameters.get_STRING(found_Parameters[10])) == "YES")
    {
        enable_Folder_management = "YES";

        vector<string> folder_Management = {"\"Compress folders\""};
        vector<string> folder_management_Parameters = Parameters.get_parameters(parameter_Master_Location, folder_Management);

        if (function.to_Upper_Case(Parameters.get_STRING(folder_Management[0])) == "YES")
        {
            enable_Compression = "YES";
        }
    }
    cout << enable_Folder_management << endl;
    cout << "Folder compression: " << enable_Compression << endl;

    cout << "\nConfiguring node master profiles:\n";
    node_Master_location = Parameters.get_STRING(found_Parameters[5]);

    cout << "\nConfiguring sequence master profiles:\n";
    sequence_Master_location = Parameters.get_STRING(found_Parameters[6]);

    // stop_after_generations = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[10]));

    cout << "\nConfiguring simulation termination by overall generations run\n";

    vector<string> parameters_List_stop_Gen_Mode = {"\"Mode to stop\""};
    vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

    if (function.to_Upper_Case(Parameters.get_STRING(stop_Gen_Parameters[0])) == "GENERATIONS")
    {
        stop_gen_Mode = 0;
        cout << "Stop after given number of generations: ";
        parameters_List_stop_Gen_Mode.clear();
        stop_Gen_Parameters.clear();

        vector<string> parameters_List_stop_Gen_Mode = {"\"Number of generations\""};
        vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

        stop_generations_Count = Parameters.get_INT(stop_Gen_Parameters[0]);
        cout << stop_generations_Count << endl;
    }
    else
    {
        cout << "Stop after given date: ";
        stop_gen_Mode = 1;
        parameters_List_stop_Gen_Mode.clear();
        stop_Gen_Parameters.clear();

        vector<string> parameters_List_stop_Gen_Mode = {"\"End date\""};
        vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

        string stop_Date_String = Parameters.get_STRING(stop_Gen_Parameters[0]);

        cout << stop_Date_String << endl;
        vector<string> split_Date;
        function.split(split_Date, stop_Date_String, '-');

        stop_Date = function.date_to_Decimal(stoi(split_Date[0]), stoi(split_Date[1]), stoi(split_Date[2]));
        cout << "Decimal date: " << stop_Date << endl;
    }

    cout << "\nConfiguring hardware resources:\n\n";
    this->CPU_cores = Parameters.get_INT(found_Parameters[1]);
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    cout << "Maximum cells processed at a time: " << this->max_Cells_at_a_time << endl
         << endl;

    cout << "Maximum sequences per file (Intermediary): " << this->max_sequences_per_File << endl
         << endl;
    this->multi_Read = Parameters.get_STRING(found_Parameters[4]);
    transform(multi_Read.begin(), multi_Read.end(), multi_Read.begin(), ::toupper);
    cout << "Multiple read and write: " << this->multi_Read << endl
         << endl;

    string cuda_IDs_String = Parameters.get_STRING(found_Parameters[0]);
    vector<string> cuda_IDs;
    function.split(cuda_IDs, cuda_IDs_String, ',');
    if (cuda_IDs.size() > 0)
    {
        this->num_Cuda_devices = cuda_IDs.size();
        CUDA_device_IDs = (int *)malloc(sizeof(int) * num_Cuda_devices);
        tot_Blocks = (int *)malloc(sizeof(int) * num_Cuda_devices);
        tot_ThreadsperBlock = (int *)malloc(sizeof(int) * num_Cuda_devices);
        function.print_Cuda_devices(cuda_IDs, this->CUDA_device_IDs, num_Cuda_devices, this->tot_Blocks, this->tot_ThreadsperBlock);
    }
    else
    {
        cout << "ERROR: THERE HAS TO BE AT LEAST ONE CUDA DEVICE SELECTED\n";
        exit(-1);
    }
}

void cancer::ingress()
{
    cout << "\nConfiguring parameters\n";

    functions_library functions = functions_library(tot_Blocks, tot_ThreadsperBlock, CUDA_device_IDs, num_Cuda_devices, gpu_Limit, CPU_cores);

    cout << "STEP 1: Configuring node profile and within host mechanics\n\n";
    node_Master_Manager(functions);

    // exit(-1);

    cout << "STEP 2: Configuring sequence profiles\n\n";
    sequence_Master_Manager(functions);
}

void cancer::node_Master_Manager(functions_library &functions)
{
    parameter_load Parameters = parameter_load();
    cout << "Loading nodes master profile: " << this->node_Master_location << endl;

    vector<string> parameters_List = {
        "\"Shape replication time\"",
        "\"Scale replication time\""};

    vector<string> found_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
    parameters_List.clear();

    random_device rd;
    mt19937 gen(rd());

    gamma_distribution<float> generation_Time_dis(Parameters.get_FLOAT(found_Parameters[0]), Parameters.get_FLOAT(found_Parameters[1]));
    float generation_Time = generation_Time_dis(gen);

    if (generation_Time > 0)
    {
        cout << "\n";
        cout << "Time taken for cell replication (days): " << generation_Time << endl;
    }
    else
    {
        cout << "\nERROR: GENERATION SHOULE BE A NON ZERO POSITIVE VALUE. PLEASE CHECK THE REPLICATION PARAMETERS\n";
        exit(-1);
    }

    vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");

    // for (size_t i = 0; i < Tissue_profiles_block_Data.size(); i++)
    // {
    //     cout << Tissue_profiles_block_Data[i].first << " : " << Tissue_profiles_block_Data[i].second << endl;
    // }

    cout << "\nConfiguring tissues:\n";
    num_tissues_per_Node = Parameters.get_INT(Tissue_profiles_block_Data, "Number of tissues");

    if (num_tissues_per_Node > 0)
    {
        cout << "\nNumber of tissues in a node: " << num_tissues_per_Node << endl;

        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
        {
            string check_Tissue = "Tissue " + to_string(tissue + 1) + " Name";
            tissue_Names.push_back(Parameters.get_STRING(Tissue_profiles_block_Data, check_Tissue));
            cout << check_Tissue << ": " << tissue_Names[tissue] << endl;
        }
    }
    else
    {
        cout << "ERROR: TISSUE NUMBER HAS TO BE GREATER THAN ZERO.\n\n";
    }

    cout << endl;

    string viral_Terminal = Parameters.get_STRING(Tissue_profiles_block_Data, "Terminal load tissues");

    vector<string> viral_tissue_Split;

    functions.split(viral_tissue_Split, viral_Terminal, ',');
    this->terminal_tissues = viral_tissue_Split.size();
    this->terminal_array = (int *)malloc(terminal_tissues * sizeof(int));
    cout << this->terminal_tissues << " tissue(s) determine node termination (death): ";
    for (size_t i = 0; i < terminal_tissues; i++)
    {
        terminal_array[i] = stoi(viral_tissue_Split[i]) - 1;
        cout << tissue_Names[terminal_array[i]];

        if ((i + 1) != terminal_tissues)
        {
            cout << ", ";
        }
        else
        {
            cout << endl;
        }
    }

    viral_Migration = Parameters.get_STRING(Tissue_profiles_block_Data, "Metastatic migration");
    transform(viral_Migration.begin(), viral_Migration.end(), viral_Migration.begin(), ::toupper);

    if (viral_Migration == "YES")
    {
        cout << "Metastatic migration: Activated\n";

        viral_Migration_Values = functions.create_Fill_2D_array_FLOAT(num_tissues_per_Node * (num_tissues_per_Node - 1), 2, -1);
        migration_start_Generation = (int *)malloc(num_tissues_per_Node * (num_tissues_per_Node - 1) * sizeof(int));

        for (int fill_mig = 0; fill_mig < num_tissues_per_Node * (num_tissues_per_Node - 1); fill_mig++)
        {
            migration_start_Generation[fill_mig] = -1;
        }

        vector<pair<string, string>> Viral_migration_block_Data = Parameters.get_block_from_block(Tissue_profiles_block_Data, "Cell migration");

        cout << "Configuring tissue to tissue migrations:\n\n";

        for (int migration_Check = 0; migration_Check < (num_tissues_per_Node * (num_tissues_per_Node - 1)); migration_Check++)
        {
            int source = migration_Check / (num_tissues_per_Node - 1);
            int destination = migration_Check % (num_tissues_per_Node - 1);

            if (destination >= source)
            {
                destination = destination + 1;
            }

            string check_source_destination = to_string(source + 1) + "_" + to_string(destination + 1);

            vector<pair<string, string>> block_Migration = Parameters.check_block_from_block(Viral_migration_block_Data, check_source_destination);

            if (block_Migration.size() > 0)
            {
                cout << "From " << tissue_Names[source] << " to " << tissue_Names[destination] << endl;
                for (int i = 0; i < block_Migration.size(); i++)
                {
                    if (Parameters.get_STRING(block_Migration[i].first) == "Cell migration Binomial trials")
                    {
                        viral_Migration_Values[migration_Check][0] = Parameters.get_INT(block_Migration[i].second);
                    }
                    else if (Parameters.get_STRING(block_Migration[i].first) == "Cell migration Binomial probability")
                    {
                        viral_Migration_Values[migration_Check][1] = Parameters.get_FLOAT(block_Migration[i].second);
                    }
                    else if (Parameters.get_STRING(block_Migration[i].first) == "Start generation")
                    {
                        migration_start_Generation[migration_Check] = Parameters.get_INT(block_Migration[i].second);
                    }
                    else
                    {
                        cout << "ERROR INVALID ENTRY AT " << check_source_destination << endl;
                        exit(-1);
                    }
                }
                cout << "Migration start generation: " << migration_start_Generation[migration_Check] << endl;
                cout << "Cell migration Binomial trials: " << viral_Migration_Values[migration_Check][0] << endl;
                cout << "Cell migration Binomial probability: " << viral_Migration_Values[migration_Check][1] << endl;
                cout << endl;
                // exit(-1);
            }
        }
    }
    else
    {
        cout << "Viral migration: Does not occur\n";
    }

    // EXTRACT NODE PROFILE INFORMATION
    cout << "\nExtracting host profile information\n";
    parameters_List = {"\"Location of node profile\""};

    found_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
    parameters_List.clear();

    string node_Profile_Location = Parameters.get_STRING(found_Parameters[0]);
    cout << "Extracting profiles from: " << node_Profile_Location << endl
         << endl;

    if (filesystem::exists(node_Profile_Location))
    {
        cout << "Configuring profile:\n";

        parameters_List = {"\"Profile name\""};

        vector<string> profile_Parameters = Parameters.get_parameters(node_Profile_Location, parameters_List);

        profile_Name = Parameters.get_STRING(profile_Parameters[0]);
        cout << "Profile name: " << profile_Name << endl;

        if (terminal_tissues > 0)
        {
            parameters_List = {"\"Terminal load distribution type\""};

            terminal_load_Profiles_param = functions.create_Fill_2D_array_FLOAT(1, 3, -1);

            vector<string> terminal_Tissue_distribution = Parameters.get_parameters(node_Profile_Location, parameters_List);
            transform(terminal_Tissue_distribution[0].begin(), terminal_Tissue_distribution[0].end(), terminal_Tissue_distribution[0].begin(), ::toupper);

            cout << "\nTerminal viral load distribution: " << terminal_Tissue_distribution[0] << endl;

            if (terminal_Tissue_distribution[0] == "\"BINOMIAL\"")
            {
                parameters_List = {"\"Terminal load Binomial trials\"",
                                   "\"Terminal load Binomial probability\""};
                terminal_Tissue_distribution = Parameters.get_parameters(node_Profile_Location, parameters_List);

                terminal_load_Profiles_param[0][0] = 0;
                terminal_load_Profiles_param[0][1] = Parameters.get_INT(terminal_Tissue_distribution[0]);
                terminal_load_Profiles_param[0][2] = Parameters.get_FLOAT(terminal_Tissue_distribution[1]);

                cout << "Terminal load Binomial trials: " << terminal_load_Profiles_param[0][1] << endl;
                cout << "Terminal load Binomial probability: " << terminal_load_Profiles_param[0][2] << endl;
            }
            else if (terminal_Tissue_distribution[0] == "\"FIXED\"")
            {
                parameters_List = {"\"Terminal load Fixed\""};
                terminal_Tissue_distribution = Parameters.get_parameters(node_Profile_Location, parameters_List);

                terminal_load_Profiles_param[0][0] = -1;
                terminal_load_Profiles_param[0][1] = Parameters.get_INT(terminal_Tissue_distribution[0]);

                cout << "Terminal load Fixed: " << terminal_load_Profiles_param[0][1] << endl;
            }
            else
            {
                cout << "PROFILE ERROR TERMINAL TISSUE DISTRIBUTION SHOULD BE FIXED OR BINOMIAL.\n";
                exit(-1);
            }
        }

        cout << "\nCollecting tissue data\n";
        vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Profile_Location, "Tissue profiles");

        profile_tissue_Limits = (int *)malloc(sizeof(int) * num_tissues_per_Node);

        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
        {
            vector<float> time_Ratios;
            vector<string> phase_Type;
            vector<pair<float, float>> phase_paramaters;

            string get_Tissue = "Tissue " + to_string(tissue + 1);
            cout << "Processing: " << get_Tissue << endl;

            vector<pair<string, string>> current_tissue_Profile_block_Data = Parameters.get_block_from_block(Tissue_profiles_block_Data, get_Tissue);

            string cell_Limit = Parameters.get_STRING(current_tissue_Profile_block_Data, "Cell limit");
            transform(cell_Limit.begin(), cell_Limit.end(), cell_Limit.begin(), ::toupper);

            cout << tissue_Names[tissue] << " tissue cell limit: ";

            if (cell_Limit == "YES")
            {
                cout << "YES\n";
                // profile_tissue_Limits[tissue] = 1;
                // cout << functions.to_Upper_Case(Parameters.get_STRING(current_tissue_Profile_block_Data, "Cell limit")) << endl;

                if (functions.to_Upper_Case(Parameters.get_STRING(current_tissue_Profile_block_Data, "Cell limit type")) == "BINOMIAL")
                {
                    cout << "Distribution type: BINOMIAL\n";
                    int trials_Cell_count = Parameters.get_INT(current_tissue_Profile_block_Data, "Cell limit Binomial trials");
                    float prob_Cell_count = Parameters.get_FLOAT(current_tissue_Profile_block_Data, "Cell limit Binomial probability");

                    cout << "Cell limit Binomial trials: " << trials_Cell_count << endl
                         << "Cell limit Binomial probability: " << prob_Cell_count << endl;

                    binomial_distribution<int> distribution(trials_Cell_count, prob_Cell_count);
                    profile_tissue_Limits[tissue] = distribution(gen);
                }
                else
                {
                    cout << "Distribution type: FIXED\n";
                    profile_tissue_Limits[tissue] = Parameters.get_INT(current_tissue_Profile_block_Data, "Cell limit Fixed");
                }

                cout << "Max cell count: " << profile_tissue_Limits[tissue] << endl;
            }
            else
            {
                cout << " NO\n";
                profile_tissue_Limits[tissue] = -1;
            }

            cout << endl;

            cout << "Configuring Tissue " << tissue + 1 << " replication phases: \n";
            vector<pair<string, string>> replication_Phases_Block = Parameters.get_block_from_block(current_tissue_Profile_block_Data, "Replication phases");
            int num_Phases = Parameters.get_INT(replication_Phases_Block, "Number of phases");

            if (num_Phases > 0)
            {
                replication_phases_tissues.push_back(num_Phases);
                cout << "Number of phases: " << replication_phases_tissues[replication_phases_tissues.size() - 1] << endl;

                float time_Check = 0;

                for (int rep_Phase = 0; rep_Phase < num_Phases; rep_Phase++)
                {
                    string phase_keyword = "Phase " + to_string(rep_Phase + 1);
                    string phase_Mode = phase_keyword + " Mode";
                    string phase_Time_ratio = phase_keyword + " Time ratio";

                    phase_Type.push_back(Parameters.get_STRING(replication_Phases_Block, phase_Mode));
                    time_Ratios.push_back(Parameters.get_FLOAT(replication_Phases_Block, phase_Time_ratio));

                    time_Check = time_Check + time_Ratios[time_Ratios.size() - 1];

                    transform(phase_Type[phase_Type.size() - 1].begin(), phase_Type[phase_Type.size() - 1].end(), phase_Type[phase_Type.size() - 1].begin(), ::toupper);

                    if (phase_Type[phase_Type.size() - 1] == "NEUTRAL")
                    {
                        phase_paramaters.push_back(make_pair(-1, -1));
                    }
                    else if (phase_Type[phase_Type.size() - 1] == "STATIONARY")
                    {
                        phase_paramaters.push_back(make_pair(Parameters.get_FLOAT(replication_Phases_Block, phase_keyword + " Variance"), -1));
                    }
                    else if (phase_Type[phase_Type.size() - 1] == "DEPRICIATION")
                    {
                        phase_paramaters.push_back(make_pair(Parameters.get_FLOAT(replication_Phases_Block, phase_keyword + " Alpha"), Parameters.get_FLOAT(replication_Phases_Block, phase_keyword + " Beta")));
                    }
                    else
                    {
                        cout << "ERROR  PROFILE TISSUE: " << tissue + 1 << "TISSUE REPLICATION MODE HAS TO BE ONE OF NEUTRAL, STATIONARY OT DEPRICIATION.\n";
                        exit(-1);
                    }

                    cout << "\nPhase " << rep_Phase + 1 << ": \n";
                    cout << "Mode: " << phase_Type[phase_Type.size() - 1] << endl;
                    // CHECK TIME RATIO ADDITIONS
                    cout << "Time ratio: " << time_Ratios[time_Ratios.size() - 1] << endl;

                    if (phase_paramaters[phase_paramaters.size() - 1].first != -1 && phase_paramaters[phase_paramaters.size() - 1].second == -1)
                    {
                        cout << "Variance of Stationary: " << phase_paramaters[phase_paramaters.size() - 1].first << endl;
                    }
                    else if (phase_paramaters[phase_paramaters.size() - 1].first != -1 && phase_paramaters[phase_paramaters.size() - 1].second != -1)
                    {
                        cout << "Alpha of Depriciation: " << phase_paramaters[phase_paramaters.size() - 1].first;
                        cout << " ; Beta of Depriciation: " << phase_paramaters[phase_paramaters.size() - 1].second << endl;
                    }
                }

                cout << "\nPerforming time ratio check:";
                if (time_Check != 1)
                {
                    cout << " FAILED\nERROR: TIME RATIOS MUST SUM UPTO 1, NOW THEY SUM UPTO: " << time_Check << endl;
                    exit(-1);
                }
                else
                {
                    phase_Type_per_tissue.push_back(phase_Type);
                    phase_paramaters_per_Tissue.push_back(phase_paramaters);
                    time_Ratios_per_Tissue.push_back(time_Ratios);
                    cout << " PASSED\n";
                }
            }
            else
            {
                cout << "ERROR: PROFILE TISSUE: " << tissue + 1 << " CANNOT HAVE LESS THAN ONE PHASES." << endl;
                exit(-1);
            }
        }
    }
    else
    {
        cout << "ERROR: PROFILE NOT FOUND, PLEASE CHECK LOCATION: " << node_Profile_Location << endl;
        exit(-1);
    }

    cout << "\nCompleted host profile configuration\n";
    cout << endl;

    // exit(-1);
}

void cancer::sequence_Master_Manager(functions_library &functions)
{
    parameter_load Parameters = parameter_load();
    cout << "Loading sequence master profile: " << this->sequence_Master_location << endl;

    vector<string> parameters_List = {
        "\"Parent sequences folder\"",
        "\"Mutation availability\"",
        "\"Reference Fitness\"",
        "\"Reference Survivability\""};

    vector<string> found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

    parent_Sequence_Folder = Parameters.get_STRING(found_Parameters[0]);
    cout << "\nParent sequences folder: " << parent_Sequence_Folder << endl;

    cout << "\nConfiguring reference genome parameters:\n";
    Reference_fitness_survivability_proof_reading = (float *)malloc(sizeof(float) * 3);

    cout << "Reference Fitness: ";
    Reference_fitness_survivability_proof_reading[0] = Parameters.get_FLOAT(found_Parameters[2]);
    cout << Reference_fitness_survivability_proof_reading[0] << endl;

    cout << "Reference Survivability: ";
    Reference_fitness_survivability_proof_reading[1] = Parameters.get_FLOAT(found_Parameters[3]);
    cout << Reference_fitness_survivability_proof_reading[1] << endl;

    mutation_proof_Reading_availability = (int *)malloc(sizeof(int) * 2);
    cout << "\nSequence mechanisms: \n";

    string status = Parameters.get_STRING(found_Parameters[1]);
    transform(status.begin(), status.end(), status.begin(), ::toupper);

    cout << "Mutations: ";

    if (status == "YES")
    {
        mutation_proof_Reading_availability[0] = 1;
        cout << "Active\n";

        vector<string> parameter_Proof_Reading = {"\"Proof reading availability\""};
        vector<string> found_Proof_Reading = Parameters.get_parameters(sequence_Master_location, parameter_Proof_Reading);

        transform(found_Proof_Reading[0].begin(), found_Proof_Reading[0].end(), found_Proof_Reading[0].begin(), ::toupper);

        cout << "Proof reading: ";
        if (Parameters.get_STRING(found_Proof_Reading[0]) == "YES")
        {
            mutation_proof_Reading_availability[1] = 1;
            cout << "Active\n";

            parameter_Proof_Reading = {"\"Reference Proof Reading\""};
            found_Proof_Reading = Parameters.get_parameters(sequence_Master_location, parameter_Proof_Reading);
            cout << "Reference Proof reading: ";
            Reference_fitness_survivability_proof_reading[2] = Parameters.get_FLOAT(found_Proof_Reading[0]);
            cout << Reference_fitness_survivability_proof_reading[2] << endl;
        }
        else
        {
            mutation_proof_Reading_availability[1] = 0;
            Reference_fitness_survivability_proof_reading[2] = -1;
            cout << "Not active\n";
        }
    }
    else
    {
        mutation_proof_Reading_availability[0] = 0;
        mutation_proof_Reading_availability[1] = 0;

        Reference_fitness_survivability_proof_reading[2] = -1;

        cout << "Not active\n";
    }

    if (mutation_proof_Reading_availability[0] == 1)
    {
        vector<pair<string, string>> mutations_Block = Parameters.get_block_from_File(sequence_Master_location, "Mutations");

        mutation_Hotspots = Parameters.get_INT(mutations_Block, "Number of hotspots");

        if (mutation_Hotspots > 0)
        {
            cout << "\nProcessing " << this->mutation_Hotspots << " mutation hotspots: \n";

            A_0_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);
            T_1_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);
            G_2_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);
            C_3_mutation = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 4, -1);

            mutation_hotspot_parameters = functions.create_Fill_2D_array_FLOAT(mutation_Hotspots, 5, -1);

            for (int mutation_hotspot = 0; mutation_hotspot < mutation_Hotspots; mutation_hotspot++)
            {
                string hotspot_ID = "Hotspot " + to_string(mutation_hotspot + 1);
                cout << "\nProcessing: " << hotspot_ID << endl;
                vector<pair<string, string>> mutations_hotspot_Block = Parameters.get_block_from_block(mutations_Block, hotspot_ID);

                string region = Parameters.get_STRING(mutations_hotspot_Block, "Region");
                string clock_Model = Parameters.get_STRING(mutations_hotspot_Block, "Clock model");
                transform(clock_Model.begin(), clock_Model.end(), clock_Model.begin(), ::toupper);

                vector<string> split_Region;
                functions.split(split_Region, region, '_');

                mutation_hotspot_parameters[mutation_hotspot][0] = stof(split_Region[0]);
                mutation_hotspot_parameters[mutation_hotspot][1] = stof(split_Region[1]);
                cout << "Region: From " << mutation_hotspot_parameters[mutation_hotspot][0] << " to " << mutation_hotspot_parameters[mutation_hotspot][1] << endl;

                cout << "Clock model: " << clock_Model << endl;

                if (clock_Model == "POISSON")
                {
                    mutation_hotspot_parameters[mutation_hotspot][2] = 0;
                    mutation_hotspot_parameters[mutation_hotspot][3] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Poisson mean");

                    cout << "Poisson mean: " << mutation_hotspot_parameters[mutation_hotspot][3] << endl;
                }
                else if (clock_Model == "NEGATIVE BINOMIAL")
                {
                    mutation_hotspot_parameters[mutation_hotspot][2] = 1;
                    mutation_hotspot_parameters[mutation_hotspot][3] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Negative Binomial sucesses");
                    mutation_hotspot_parameters[mutation_hotspot][4] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Negative Binomial probability");

                    cout << "Negative Binomial sucesses: " << mutation_hotspot_parameters[mutation_hotspot][3] << endl;
                    cout << "Negative Binomial probability: " << mutation_hotspot_parameters[mutation_hotspot][4] << endl;
                }
                else if (clock_Model == "FIXED")
                {
                    mutation_hotspot_parameters[mutation_hotspot][2] = 2;
                    mutation_hotspot_parameters[mutation_hotspot][3] = Parameters.get_FLOAT(mutations_hotspot_Block, "Clock model Fixed probability");

                    cout << "Fixed probability: " << mutation_hotspot_parameters[mutation_hotspot][3] << endl;
                }
                else
                {
                    cout << "HOTSPOT " << mutation_hotspot + 1 << "'S CLOCK MODEL HAS TO BE POISSON, NEGATIVE BINOMIAL OR FIXED.\n";
                    exit(-1);
                }

                vector<string> bases = {"A", "T", "G", "C"};

                cout << "\nConfiguring mutation probabilities:\n";

                for (int reference_Base = 0; reference_Base < bases.size(); reference_Base++)
                {
                    for (int mutation_Base = 0; mutation_Base < bases.size(); mutation_Base++)
                    {
                        string search_Query = bases[reference_Base] + bases[mutation_Base];
                        cout << search_Query << ": ";

                        if (reference_Base == 0)
                        {
                            A_0_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                            cout << A_0_mutation[mutation_hotspot][mutation_Base];
                        }
                        else if (reference_Base == 1)
                        {
                            T_1_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                            cout << T_1_mutation[mutation_hotspot][mutation_Base];
                        }
                        else if (reference_Base == 2)
                        {
                            G_2_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                            cout << G_2_mutation[mutation_hotspot][mutation_Base];
                        }
                        else
                        {
                            C_3_mutation[mutation_hotspot][mutation_Base] = Parameters.get_FLOAT(mutations_hotspot_Block, search_Query);
                            cout << C_3_mutation[mutation_hotspot][mutation_Base];
                        }

                        if ((mutation_Base + 1) != bases.size())
                        {
                            cout << " | ";
                        }
                    }
                    cout << endl;
                }
            }
        }
        else
        {
            cout << "No mutational hotspots present to configure.\n";
        }
    }

    parameters_List = {
        "\"Fitness profile file\"",
        "\"Survivability profile file\""};

    found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

    string fitness_Profile_Location = Parameters.get_STRING(found_Parameters[0]);
    string survivability_Profile_Location = Parameters.get_STRING(found_Parameters[1]);

    num_effect_Segregating_sites = (int *)malloc(sizeof(int) * 3);

    if (fitness_Profile_Location != "NA")
    {
        sequence_Fitness_changes = Parameters.get_Profile_Array(fitness_Profile_Location, num_effect_Segregating_sites[0], functions);

        cout << "Printing Fitness matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites[0]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_Fitness_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No fitness profile\n\n";
        num_effect_Segregating_sites[0] = 0;
    }

    if (survivability_Profile_Location != "NA")
    {
        sequence_Survivability_changes = Parameters.get_Profile_Array(survivability_Profile_Location, num_effect_Segregating_sites[1], functions);

        cout << "Printing Survivability matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites[1]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_Survivability_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No survivability profile\n\n";
        num_effect_Segregating_sites[1] = 0;
    }

    if (mutation_proof_Reading_availability[1] == 1)
    {
        parameters_List = {
            "\"Proof reading profile file\""};

        found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);
        string proof_Reading_Profile_Location = Parameters.get_STRING(found_Parameters[0]);

        if (proof_Reading_Profile_Location != "NA")
        {
            sequence_Proof_reading_changes = Parameters.get_Profile_Array(proof_Reading_Profile_Location, num_effect_Segregating_sites[2], functions);

            cout << "Printing Proof reading matrix:\n";

            for (int row = 0; row < num_effect_Segregating_sites[2]; row++)
            {
                for (int col = 0; col < 5; col++)
                {
                    cout << sequence_Proof_reading_changes[row][col] << "\t";
                }
                cout << endl;
            }
        }
        else
        {
            cout << "No proof reading profile\n\n";
            num_effect_Segregating_sites[2] = 0;
        }
    }
    else
    {
        num_effect_Segregating_sites[2] = 0;
    }

    cout << "\nCancer parameters\n\n";

    parameters_List = {
        "\"Reference_replication_factor\"",
        "\"Reference_mutation_rate_factor\"",
        "\"Reference_generation_death_prob\"",
        "\"Reference_replication_prob\"",
        "\"Reference_metastatic_prob\""};

    found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

    Reference_cancer_parameters = (float *)malloc(sizeof(float) * 5);

    Reference_cancer_parameters[0] = Parameters.get_FLOAT(found_Parameters[0]);
    cout << "Reference replication factor: " << Reference_cancer_parameters[0] << endl;

    Reference_cancer_parameters[1] = Parameters.get_FLOAT(found_Parameters[1]);
    cout << "Reference mutation rate factor: " << Reference_cancer_parameters[1] << endl;

    Reference_cancer_parameters[2] = Parameters.get_FLOAT(found_Parameters[2]);
    cout << "Reference generational death: " << Reference_cancer_parameters[2] << endl;

    Reference_cancer_parameters[3] = Parameters.get_FLOAT(found_Parameters[3]);
    cout << "Reference replication probability: " << Reference_cancer_parameters[3] << endl;

    Reference_cancer_parameters[4] = Parameters.get_FLOAT(found_Parameters[4]);
    cout << "Reference metastatic probability: " << Reference_cancer_parameters[4] << endl;

    parameters_List = {
        "\"Replication_factor profile file\"",
        "\"Mutation_rate_factor profile file\"",
        "\"Generation_death profile file\"",
        "\"Replication_prob profile file\"",
        "\"Reference_metastatic profile file\""};

    found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

    string Replication_factor_Profile_Location = Parameters.get_STRING(found_Parameters[0]);
    string Mutation_rate_factor_Profile_Location = Parameters.get_STRING(found_Parameters[1]);
    string Generation_death_Profile_Location = Parameters.get_STRING(found_Parameters[2]);
    string Replication_prob_Profile_Location = Parameters.get_STRING(found_Parameters[3]);
    string Reference_metastatic_Profile_Location = Parameters.get_STRING(found_Parameters[4]);

    num_effect_Segregating_sites_Cancer = (int *)malloc(sizeof(int) * 5);

    if (Replication_factor_Profile_Location != "NA")
    {
        sequence_replication_factor_changes = Parameters.get_Profile_Array(Replication_factor_Profile_Location, num_effect_Segregating_sites_Cancer[0], functions);

        cout << "\nPrinting replication factor matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_replication_factor_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No replication factor profile\n\n";
        num_effect_Segregating_sites_Cancer[0] = 0;
    }

    if (Mutation_rate_factor_Profile_Location != "NA")
    {
        sequence_mutation_rate_changes = Parameters.get_Profile_Array(Mutation_rate_factor_Profile_Location, num_effect_Segregating_sites_Cancer[1], functions);

        cout << "Printing mutation rateF factor matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_mutation_rate_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No mutation rate factor profile\n\n";
        num_effect_Segregating_sites_Cancer[0] = 0;
    }

    if (Generation_death_Profile_Location != "NA")
    {
        sequence_generation_death_changes = Parameters.get_Profile_Array(Generation_death_Profile_Location, num_effect_Segregating_sites_Cancer[2], functions);

        cout << "Printing geerational death matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_generation_death_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No generational death profile\n\n";
        num_effect_Segregating_sites_Cancer[0] = 0;
    }

    if (Replication_prob_Profile_Location != "NA")
    {
        sequence_replication_prob_changes = Parameters.get_Profile_Array(Replication_prob_Profile_Location, num_effect_Segregating_sites_Cancer[3], functions);

        cout << "Printing replication probability matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_replication_prob_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No replication probability profile\n\n";
        num_effect_Segregating_sites_Cancer[0] = 0;
    }

    if (Reference_metastatic_Profile_Location != "NA")
    {
        sequence_metastatic_prob_changes = Parameters.get_Profile_Array(Reference_metastatic_Profile_Location, num_effect_Segregating_sites_Cancer[4], functions);

        cout << "Printing metastasis probability matrix:\n";

        for (int row = 0; row < num_effect_Segregating_sites_Cancer[0]; row++)
        {
            for (int col = 0; col < 5; col++)
            {
                cout << sequence_metastatic_prob_changes[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    else
    {
        cout << "No metastasis probability profile\n\n";
        num_effect_Segregating_sites_Cancer[0] = 0;
    }

    cout << "\nCompleted sequence configuration\n";
}