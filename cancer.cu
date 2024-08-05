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

    cout << "STEP 1: Configuring sequence profiles\n\n";
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

    cout << endl;

    exit(-1);
}

void cancer::sequence_Master_Manager(functions_library &functions)
{
    parameter_load Parameters = parameter_load();
    cout << "Loading sequence master profile: " << this->sequence_Master_location << endl;

    vector<string> parameters_List = {
        "\"Parent sequences folder\"",
        "\"Mutation availability\"",
        "\"Recombination availability\"",
        "\"Reference Fitness\"",
        "\"Reference Survivability\""};

    vector<string> found_Parameters = Parameters.get_parameters(sequence_Master_location, parameters_List);

    parent_Sequence_Folder = Parameters.get_STRING(found_Parameters[0]);
    cout << "\nParent sequences folder: " << parent_Sequence_Folder << endl;

    cout << "\nConfiguring reference genome parameters:\n";
    Reference_fitness_survivability_proof_reading = (float *)malloc(sizeof(float) * 3);

    cout << "Reference Fitness: ";
    Reference_fitness_survivability_proof_reading[0] = Parameters.get_FLOAT(found_Parameters[3]);
    cout << Reference_fitness_survivability_proof_reading[0] << endl;

    cout << "Reference Survivability: ";
    Reference_fitness_survivability_proof_reading[1] = Parameters.get_FLOAT(found_Parameters[4]);
    cout << Reference_fitness_survivability_proof_reading[1] << endl;

    mutation_proof_Reading_availability = (int *)malloc(sizeof(int) * 3);
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
}