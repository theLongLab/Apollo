#include "simulator_Master.cuh"

simulator_Master::simulator_Master(string parameter_Master_Location)
{
    cout << "Intializing Simulator (based on CATE's architecture)\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"CUDA Device ID\"",
        "\"CPU cores\"",
        "\"GPU max units\"",
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Multi read\"",
        "\"Network profile\"",
        "\"Nodes master profile\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    cout << "Configuring folders:\n";

    output_Folder_location = Parameters.get_STRING(found_Parameters[4]);
    intermediate_Folder_location = Parameters.get_STRING(found_Parameters[3]);

    function.config_Folder(intermediate_Folder_location, "Intermediate");
    function.config_Folder(output_Folder_location, "Output");

    output_Network_location = this->output_Folder_location + "/network_Data";
    function.config_Folder(output_Network_location, "Network");
    network_File_location = output_Network_location + "/node_node_Relationships.csv";
    function.create_File(network_File_location, "Source\tTarget");

    cout << "\nConfiguring parameter profiles:\n";
    node_Master_location = Parameters.get_STRING(found_Parameters[7]);

    cout << "\nConfiguring hardware resources:\n\n";
    this->CPU_cores = Parameters.get_INT(found_Parameters[1]);
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    this->multi_Read = Parameters.get_STRING(found_Parameters[5]);
    cout << "Multiple read and write: " << this->multi_Read << endl
         << endl;

    this->CUDA_device_number = Parameters.get_INT(found_Parameters[0]);
    function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

    this->gpu_Limit = Parameters.get_INT(found_Parameters[2]);

    cout << "Per round GPU max unit: " << this->gpu_Limit << endl
         << endl;

    configure_Network_Profile(Parameters.get_STRING(found_Parameters[6]), Parameters);
    cout << "\n";
}

void simulator_Master::configure_Network_Profile(string network_Profile_File, parameter_load &Parameters)
{
    cout << "Configuring network profile: " << network_Profile_File << endl;

    vector<string> parameters_List = {"\"Network type\""};
    vector<string> found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

    transform(found_Parameters[0].begin(), found_Parameters[0].end(), found_Parameters[0].begin(), ::toupper);

    parameters_List.clear();
    found_Parameters.clear();

    if (Parameters.get_STRING(found_Parameters[0]) == "BA MODEL")
    {
        cout << "\nBarabsi Albert model selected: \n";
        network_Model = "BA";

        parameters_List = {"\"BA model number of nodes\"",
                           "\"BA model standard new connections\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        Total_number_of_Nodes = Parameters.get_INT(found_Parameters[0]);
        cout << "Number of nodes: " << Total_number_of_Nodes << endl;
        connection_Model = Parameters.get_STRING(found_Parameters[1]);
        transform(connection_Model.begin(), connection_Model.end(), connection_Model.begin(), ::toupper);

        cout << "Node connection type: " << connection_Model << endl;

        if (connection_Model == "FIXED")
        {
            parameters_List = {"\"BA model fixed new connections\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            BA_FIXED = Parameters.get_INT(found_Parameters[0]);
            cout << "Fixed new connections: " << BA_FIXED << endl;
        }
        else if (connection_Model == "NEGATIVE BINOMIAL")
        {
            parameters_List = {"\"BA model Negative binomial sucesses\"",
                               "\"BA model Negative binomial probability\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            BA_NB_sucesses = Parameters.get_INT(found_Parameters[0]);
            BA_NB_probability = Parameters.get_FLOAT(found_Parameters[1]);

            cout << "Negative Binomial sucesses: " << BA_NB_sucesses << endl;
            cout << "Negative Binomial probability: " << BA_NB_probability << endl;
        }
        else if (connection_Model == "POISSON")
        {
            parameters_List = {"\"BA model Poisson mean\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            BA_Poisson_mean = Parameters.get_FLOAT(found_Parameters[0]);

            cout << "Poisson mean: " << BA_Poisson_mean << endl;
        }
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "SC MODEL")
    {
        cout << "\nStandard Caveman model selected: \n";
        network_Model = "SCM";

        parameters_List = {"\"SC_model number of caves\"",
                           "\"SC_model number of nodes per caves\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        SCM_number_of_caves = Parameters.get_INT(found_Parameters[0]);
        SCM_number_of_nodes_per_cave = Parameters.get_INT(found_Parameters[1]);

        cout << "Number of caves: " << SCM_number_of_caves;
        cout << "\nNumber of nodes per cave: " << SCM_number_of_nodes_per_cave;
        Total_number_of_Nodes = SCM_number_of_caves * SCM_number_of_nodes_per_cave;
        cout << "\nTotal nodes in network: " << Total_number_of_Nodes << endl;
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "DC MODEL")
    {
        cout << "\nDynamic Caveman model selected: \n";
        network_Model = "DCM";

        parameters_List = {"\"DC_model number of caves\"",
                           "\"DC_model node distribution\"",
                           "\"DC_model neighbouring nodes percentage\"",
                           "\"DC_model global nodes percentage\""};

        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        DCM_number_of_caves = Parameters.get_INT(found_Parameters[0]);
        cout << "Number of caves: " << DCM_number_of_caves << endl;

        connection_Model = Parameters.get_STRING(found_Parameters[1]);
        transform(connection_Model.begin(), connection_Model.end(), connection_Model.begin(), ::toupper);

        DC_percent_Neighbouring = Parameters.get_FLOAT(found_Parameters[2]);
        DC_percent_Global_freedom = Parameters.get_FLOAT(found_Parameters[3]);

        cout << "Node connection type: " << connection_Model << endl;

        if (connection_Model == "NEGATIVE BINOMIAL")
        {
            parameters_List = {"\"DC_model node Negative binomial sucesses\"",
                               "\"DC_model node Negative binomial probability\""};

            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            DC_ND_sucesses = Parameters.get_INT(found_Parameters[0]);
            DC_ND_probability = Parameters.get_FLOAT(found_Parameters[1]);

            cout << "Negative Binomial sucesses: " << DC_ND_sucesses << endl;
            cout << "Negative Binomial probability: " << DC_ND_probability << endl;
        }
        else if (connection_Model == "POISSON")
        {
            parameters_List = {"\"DC_model node Poisson mean\""};
            found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

            DC_Poisson_mean = Parameters.get_FLOAT(found_Parameters[0]);

            cout << "Poisson mean: " << DC_Poisson_mean << endl;
        }
    }
    else if (Parameters.get_STRING(found_Parameters[0]) == "RANDOM MODEL")
    {
        cout << "\nRandom model selected: \n";
        network_Model = "RANDOM";

        parameters_List = {"\"Random model number of nodes\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        Total_number_of_Nodes = Parameters.get_INT(found_Parameters[0]);
    }
    else
    {
        cout << "ERROR Incorrect network selected. Please check \"Network type\" in the network parameter file: \"" << network_Profile_File << "\"";
        exit(-1);
    }
}

void simulator_Master::ingress()
{
    functions_library functions = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    cout << "STEP 1: Configuring population network\n\n";

    // ! Compatible for both BA and Caveman
    // INT, INT = Cave_ID and Node, for BA CaveID is 0 for all.
    vector<vector<pair<int, int>>> each_Nodes_Connection;
    // ! Only for CAVEMAN models, to keep track of each nodes Caves

    network_Manager(each_Nodes_Connection, functions);

    cout << "STEP 2: Configuring within host mechanics\n\n";
    node_Master_Manager(functions);

    cout << "STEP 3: Configuring of Population profiles\n\n";
}

void simulator_Master::node_Master_Manager(functions_library &functions)
{
    parameter_load Parameters = parameter_load();
    cout << "Loading nodes master profile: " << this->node_Master_location << endl;

    vector<string> parameters_List = {
        "\"Shape replication time\"",
        "\"Scale replication time\"",
        "\"Shape days in host\"",
        "\"Scale days in host\"",
        "\"Sampling present\"",
        "\"Progeny distribution type\"",
        "\"Number of node profiles\"",
        "\"Location of node profiles\""};

    vector<string> found_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
    parameters_List.clear();

    random_device rd;
    mt19937 gen(rd());

    gamma_distribution<float> days_in_Host(Parameters.get_FLOAT(found_Parameters[0]), Parameters.get_FLOAT(found_Parameters[1]));
    this->generation_Time = days_in_Host(gen);

    if (this->generation_Time > 0)
    {
        cout << "\n";
        cout << "Time taken for viral replication (from cell attachment to release) (days): " << generation_Time << endl;

        this->shape_days_in_Host = Parameters.get_FLOAT(found_Parameters[2]);
        this->scale_days_in_Host = Parameters.get_FLOAT(found_Parameters[3]);
        cout << "Days in host gamma distribution: Shape: " << this->shape_days_in_Host << " Scale: " << this->scale_days_in_Host << endl;

        cout << "\nSampling status: ";
        if (functions.to_Upper_Case(Parameters.get_STRING(found_Parameters[4])) == "YES")
        {
            cout << "Active\n";
            parameters_List = {"\"Sampling effect\"",
                               "\"Sampling rate Binomial trials\"",
                               "\"Sampling rate Binomial probability\""};

            vector<string> sampling_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);

            sampling_effect = Parameters.get_STRING(sampling_Parameters[0]);
            transform(sampling_effect.begin(), sampling_effect.end(), sampling_effect.begin(), ::toupper);

            cout << "Effect of sampling: " << sampling_effect << endl;

            sampling_trials = Parameters.get_INT(sampling_Parameters[1]);
            sampling_probability = Parameters.get_FLOAT(sampling_Parameters[2]);

            cout << "Sampling rate distribution trials: " << sampling_trials << endl;
            cout << "Sampling rate distribution probability: " << sampling_probability << endl;
        }
        else
        {
            cout << "Inactive\n";
        }

        progeny_distribution_Model = Parameters.get_STRING(found_Parameters[5]);
        transform(progeny_distribution_Model.begin(), progeny_distribution_Model.end(), progeny_distribution_Model.begin(), ::toupper);

        vector<string> progeny_distribution_Parameters;

        if (progeny_distribution_Model == "NEGATIVE BINOMIAL")
        {
            parameters_List = {"\"Progeny Negative Binomial sucesses\"",
                               "\"Progeny Negative Binomial probability\""};

            progeny_distribution_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);

            this->progeny_NB_sucesses = Parameters.get_INT(progeny_distribution_Parameters[0]);
            this->progeny_NB_probability = Parameters.get_FLOAT(progeny_distribution_Parameters[1]);

            cout << "\nProgeny generation via a Negative Binomial disribution\n";
            cout << "Negative Binomial sucesses: " << progeny_NB_sucesses
                 << "\nNegative Binomial probability: " << progeny_NB_probability << endl;
        }
        else if (progeny_distribution_Model == "GAMMA")
        {
            parameters_List = {"\"Progeny Gamma shape\"",
                               "\"Progeny Gamma scale\""};

            progeny_distribution_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
            this->progeny_Gamma_shape = Parameters.get_FLOAT(progeny_distribution_Parameters[0]);
            this->progeny_Gamma_scale = Parameters.get_FLOAT(progeny_distribution_Parameters[1]);

            cout << "\nProgeny generation via a Gamma disribution\n";
            cout << "Gamma shape: " << progeny_Gamma_shape
                 << "\nGamma scale: " << progeny_Gamma_scale << endl;
        }
        else if (progeny_distribution_Model == "POISSON")
        {
            parameters_List = {"\"Progeny Poisson mean\""};

            progeny_distribution_Parameters = Parameters.get_parameters(node_Master_location, parameters_List);
            this->progeny_Poisson_mean = Parameters.get_FLOAT(progeny_distribution_Parameters[0]);

            cout << "\nProgeny generation via a Poisson disribution\n";
            cout << "Poisson mean: " << progeny_Poisson_mean << endl;
        }
        else
        {
            cout << "\nERROR: PROGENY GENERATION DISTRIBUTION MODEL TYPE\n";
            exit(-1);
        }
        progeny_distribution_Parameters.clear();

        // Configure tissue profile read
        vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");

        for (size_t i = 0; i < Tissue_profiles_block_Data.size(); i++)
        {
            cout << Tissue_profiles_block_Data[i].first << " : " << Tissue_profiles_block_Data[i].second << endl;
        }

        cout << "******\n\n";

        vector<pair<string, string>> Viral_migration_block_Data = Parameters.get_block_from_block(Tissue_profiles_block_Data, "Viral migration");

        for (size_t i = 0; i < Viral_migration_block_Data.size(); i++)
        {
            cout << Viral_migration_block_Data[i].first << " : " << Viral_migration_block_Data[i].second << endl;
        }
    }
    else
    {
        cout << "\nERROR: GENERATION SHOULE BE A NON ZERO POSITIVE VALUE. PLEASE CHECK THE REPLICATION PARAMETERS\n";
        exit(-1);
    }
    cout << endl;
}

void simulator_Master::network_Manager(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions)
{
    for (int initialize = 0; initialize < Total_number_of_Nodes; initialize++)
    {
        vector<pair<int, int>> intialize_Vector;
        each_Nodes_Connection.push_back(intialize_Vector);
    }

    if (network_Model == "BA")
    {
        BA_Model_Engine(each_Nodes_Connection, functions);

        // TEST node connections = DONE
        // for (int test = 0; test < each_Nodes_Connection.size(); test++)
        // {
        //     for (size_t i = 0; i < each_Nodes_Connection[test].size(); i++)
        //     {
        //         cout << each_Nodes_Connection[test][i].second << " ";
        //     }
        //     cout << endl;
        // }

        cout << "Completed Barbasi Albert model network engine\n\n";
    }
    else if (network_Model == "SCM")
    {
        SCM_Model_Engine(each_Nodes_Connection, functions);
    }
    else if (network_Model == "DCM")
    {
        DCM_Model_Engine(each_Nodes_Connection, functions);
    }
    else if (network_Model == "RANDOM")
    {
        RANDOM_Model_Engine(each_Nodes_Connection, functions);
    }
    else
    {
        cout << "ERROR Incorrect network selected. Please check \"Network type\" in the network parameter file.\n";
        exit(-1);
    }
}

void simulator_Master::RANDOM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions)
{
    cout << "Intializing Random model network engine\n";

    random_device rd;
    mt19937 gen(rd());

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    cout << "Forming node to node relationships\n";

    cout << "Configuring node 1 of " << Total_number_of_Nodes << " node(s)" << endl;

    for (int node = 1; node < Total_number_of_Nodes; node++)
    {
        cout << "Configuring node " << (node + 1) << " of " << Total_number_of_Nodes << " node(s)" << endl;

        uniform_int_distribution<int> distribution_Neighbour(0, node - 1);

        int attach_Node = -1;
        do
        {
            attach_Node = distribution_Neighbour(gen);
        } while (attach_Node == node);

        each_Nodes_Connection[node].push_back(make_pair(0, attach_Node));
        each_Nodes_Connection[attach_Node].push_back(make_pair(0, node));

        network_File << "0"
                     << "_" << node << "\t"
                     << "0"
                     << "_" << attach_Node << "\n";
    }

    network_File.close();
    cout << endl;
}

void simulator_Master::DCM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions)
{
    cout << "Intializing Dynamic Caveman model network engine\n";

    random_device rd;
    mt19937 gen(rd());

    cout << "Determining per cave node counts\n";

    per_cave_Stride = (int *)malloc((DCM_number_of_caves + 1) * sizeof(int));
    per_cave_Stride[0] = 0;

    int **network_Array = functions.create_INT_2D_arrays(DCM_number_of_caves, 3);
    vector<vector<int>> neighbour_Nodes;
    vector<vector<int>> global_Nodes;

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        int number_of_Nodes;

        if (connection_Model == "NEGATIVE BINOMIAL")
        {
            negative_binomial_distribution<int> distribution(DC_ND_sucesses, DC_ND_probability);
            number_of_Nodes = distribution(gen);
        }
        else if (connection_Model == "POISSON")
        {
            poisson_distribution<int> distribution(BA_Poisson_mean);
            number_of_Nodes = distribution(gen);
        }

        number_of_Nodes = (number_of_Nodes < 1) ? 1 : number_of_Nodes;

        // if (number_of_Nodes < 1)
        // {
        //     number_of_Nodes = 1;
        // }

        network_Array[cave_ID][0] = number_of_Nodes;
        per_cave_Stride[cave_ID + 1] = per_cave_Stride[cave_ID] + network_Array[cave_ID][0];

        vector<int> intialize;
        neighbour_Nodes.push_back(intialize);
        global_Nodes.push_back(intialize);
    }

    Total_number_of_Nodes = per_cave_Stride[DCM_number_of_caves];

    cout << "Total nodes in network: " << Total_number_of_Nodes << endl;

    cout << "Determining relationships for neighour and global nodes\n";

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        binomial_distribution<int> binomialDist_Neighbour(network_Array[cave_ID][0], DC_percent_Neighbouring);
        network_Array[cave_ID][1] = binomialDist_Neighbour(gen);

        uniform_int_distribution<int> distribution_Neighbour(0, network_Array[cave_ID][0] - 1);

        if (network_Array[cave_ID][1] < 1)
        {
            network_Array[cave_ID][1] = 1;
        }

        vector<int> no_repetition;

        while (neighbour_Nodes[cave_ID].size() < network_Array[cave_ID][1])
        {
            int check_indicator = 0;
            int get_Node = distribution_Neighbour(gen);
            for (int check = 0; check < no_repetition.size(); check++)
            {
                if (get_Node == no_repetition[check])
                {
                    check_indicator = 1;
                    break;
                }
            }
            if (check_indicator == 0)
            {
                neighbour_Nodes[cave_ID].push_back(get_Node);
                no_repetition.push_back(get_Node);
            }
        }

        binomial_distribution<int> binomialDist_Global(network_Array[cave_ID][1], DC_percent_Global_freedom);
        network_Array[cave_ID][2] = binomialDist_Global(gen);

        uniform_int_distribution<int> distribution_Global(0, network_Array[cave_ID][1] - 1);

        no_repetition.clear();

        while (global_Nodes[cave_ID].size() < network_Array[cave_ID][2])
        {
            int check_indicator = 0;
            int get_Index = distribution_Global(gen);

            for (int check = 0; check < no_repetition.size(); check++)
            {
                if (get_Index == no_repetition[check])
                {
                    check_indicator = 1;
                    break;
                }
            }

            if (check_indicator == 0)
            {
                global_Nodes[cave_ID].push_back(neighbour_Nodes[cave_ID][get_Index]);
                no_repetition.push_back(get_Index);
            }
        }
    }

    vector<vector<pair<int, int>>> each_Nodes_Connections;

    for (size_t i = 0; i < Total_number_of_Nodes; i++)
    {
        vector<pair<int, int>> nodes_Intialize;
        each_Nodes_Connections.push_back(nodes_Intialize);
    }

    cout << "Configuring cave cohort relationships\n";

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        int nodes_Total = network_Array[cave_ID][0];

        // cout << cave_ID << endl;

        for (int node_parent = 0; node_parent < nodes_Total; node_parent++)
        {

            int node_Count = per_cave_Stride[cave_ID] + node_parent;
            // cout << node_Count << endl;

            // node_Summary << cave_ID << "_" << node_parent << "\t" << cave_ID << "\n";
            for (int node_child = node_parent + 1; node_child < nodes_Total; node_child++)
            {
                int node_child_Main = node_Count + (node_child - node_parent);
                network_File << cave_ID << "_" << node_parent << "\t" << cave_ID << "_" << node_child << "\n";
                each_Nodes_Connections[node_Count].push_back(make_pair(cave_ID, node_child));
                each_Nodes_Connections[node_child_Main].push_back(make_pair(cave_ID, node_parent));
                // cout << node_Count << "\t" << node_child_Main << endl;
            }
        }
        // cout << endl;
    }

    cout << "Configuring neigbouring relationships\n";

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        uniform_int_distribution<int> distribution_Neighbour_parent(0, neighbour_Nodes[cave_ID].size() - 1);

        int parent_Node_Index = neighbour_Nodes[cave_ID][distribution_Neighbour_parent(gen)];

        int attach_Cave = cave_ID + 1;
        int attach_Node_Index;
        if (attach_Cave != DCM_number_of_caves)
        {
            uniform_int_distribution<int> distribution_Neighbour_attach(0, neighbour_Nodes[attach_Cave].size() - 1);
            attach_Node_Index = neighbour_Nodes[attach_Cave][distribution_Neighbour_attach(gen)];
        }
        else
        {
            attach_Cave = 0;
            uniform_int_distribution<int> distribution_Neighbour_attach(0, neighbour_Nodes[attach_Cave].size() - 1);
            attach_Node_Index = neighbour_Nodes[attach_Cave][distribution_Neighbour_attach(gen)];
        }

        network_File << cave_ID << "_" << parent_Node_Index << "\t" << attach_Cave << "_" << attach_Node_Index << "\n";

        int parent_Node_all_Index = per_cave_Stride[cave_ID] + parent_Node_Index;
        int attach_Node_all_Index = per_cave_Stride[attach_Cave] + attach_Node_Index;

        each_Nodes_Connections[parent_Node_all_Index].push_back(make_pair(attach_Cave, attach_Node_Index));
        each_Nodes_Connections[attach_Node_all_Index].push_back(make_pair(cave_ID, parent_Node_Index));
    }

    cout << "Configuring Global networks attachments\n";

    vector<string> repeat_Check;

    for (int cave_ID = 0; cave_ID < DCM_number_of_caves; cave_ID++)
    {
        // cout << cave_ID;
        // if (cave_ID + 1 != DCM_number_of_caves)
        // {
        //     cout << ", ";
        // }
        for (size_t i = 0; i < global_Nodes[cave_ID].size(); i++)
        {
            int cave_Global;
            do
            {
                uniform_int_distribution<int> distribution_global_CAVEs(0, DCM_number_of_caves - 1);
                cave_Global = distribution_global_CAVEs(gen);
                // cout << cave_ID << ":\t" << cave_Global << endl;
            } while (cave_Global == cave_ID);

            if (global_Nodes[cave_Global].size() > 0)
            {
                // cout << "Nodes global: " << global_Nodes[cave_Global].size() << endl;
                uniform_int_distribution<int> distribution_Parent_attach(0, global_Nodes[cave_ID].size() - 1);
                uniform_int_distribution<int> distribution_Global_attach(0, global_Nodes[cave_Global].size() - 1);

                int parent_Cave_Node = global_Nodes[cave_ID][distribution_Parent_attach(gen)];
                int Global_Cave_Node = global_Nodes[cave_Global][distribution_Global_attach(gen)];

                string current = to_string(cave_ID) + "_" + to_string(parent_Cave_Node) + "\t" + to_string(cave_Global) + "_" + to_string(Global_Cave_Node);
                string repeat = to_string(cave_Global) + "_" + to_string(Global_Cave_Node) + "\t" + to_string(cave_ID) + "_" + to_string(parent_Cave_Node);

                int repeat_Index = -1;

                for (int check = 0; check < repeat_Check.size(); check++)
                {
                    if (repeat_Check[check] == repeat || repeat_Check[check] == current)
                    {
                        repeat_Index = 0;
                        break;
                    }
                }

                if (repeat_Index == -1)
                {
                    repeat_Check.push_back(current);
                    repeat_Check.push_back(repeat);
                    network_File << current << "\n";

                    int parent_Node_all_Index = per_cave_Stride[cave_ID] + parent_Cave_Node;
                    int attach_Node_all_Index = per_cave_Stride[cave_Global] + Global_Cave_Node;

                    each_Nodes_Connections[parent_Node_all_Index].push_back(make_pair(cave_Global, Global_Cave_Node));
                    each_Nodes_Connections[attach_Node_all_Index].push_back(make_pair(cave_ID, parent_Cave_Node));
                }
            }
        }
    }

    functions.clear_Array_int_CPU(network_Array, DCM_number_of_caves);
    network_File.close();

    cout << endl;
}

void simulator_Master::SCM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions)
{
    cout << "Intializing Standard Caveman model network engine\n";

    random_device rd;
    mt19937 gen(rd());

    cout << "Configuring cave cohort relationships\n";

    vector<int> neighbour_Node;

    uniform_int_distribution<int> distribution(0, SCM_number_of_nodes_per_cave - 1);

    cout << "Configuring cave neighbour nodes" << endl;
    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        neighbour_Node.push_back(distribution(gen));
        // cout << cave_ID << ": " << neighbour_Node[cave_ID] << endl;
    }

    cout << "Configuring cave cohort nodes" << endl;
    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;

        for (int node_Index_Parent = 0; node_Index_Parent < SCM_number_of_nodes_per_cave; node_Index_Parent++)
        {
            for (int node_Index_child = node_Index_Parent + 1; node_Index_child < SCM_number_of_nodes_per_cave; node_Index_child++)
            {
                each_Nodes_Connection[node_start_Index + node_Index_Parent].push_back(make_pair(cave_ID, node_Index_child));
                each_Nodes_Connection[node_start_Index + node_Index_child].push_back(make_pair(cave_ID, node_Index_Parent));
            }
        }
    }

    cout << "Re-routing cave cohort nodes to neighbours" << endl;

    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;
        int rerouting_Node = neighbour_Node[cave_ID] + node_start_Index;

        uniform_int_distribution<int> rerouting_Connection_Draw(0, each_Nodes_Connection[rerouting_Node].size() - 1);

        int rerouting_Connection = rerouting_Connection_Draw(gen);
        int rerouting_Connection_Node = each_Nodes_Connection[rerouting_Node][rerouting_Connection].second + node_start_Index;

        int index = -1;
        for (int find = 0; find < each_Nodes_Connection[rerouting_Connection_Node].size(); find++)
        {
            if (each_Nodes_Connection[rerouting_Connection_Node][find].second == neighbour_Node[cave_ID])
            {
                index = find;
                break;
            }
        }

        if (index == -1)
        {
            cout << "ERROR in node rerouting removal\n";
        }
        else
        {
            each_Nodes_Connection[rerouting_Connection_Node].erase(each_Nodes_Connection[rerouting_Connection_Node].begin() + index);
        }

        if ((cave_ID + 1) != SCM_number_of_caves)
        {
            each_Nodes_Connection[rerouting_Node][rerouting_Connection].first = cave_ID + 1;
            each_Nodes_Connection[rerouting_Node][rerouting_Connection].second = neighbour_Node[cave_ID + 1];
        }
        else
        {
            each_Nodes_Connection[rerouting_Node][rerouting_Connection].first = 0;
            each_Nodes_Connection[rerouting_Node][rerouting_Connection].second = neighbour_Node[0];
        }
    }

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    {
        cout << "Finalizing cave " << cave_ID + 1 << " of " << SCM_number_of_caves << endl;
        int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;

        for (int node = 0; node < SCM_number_of_nodes_per_cave; node++)
        {
            for (int connections = 0; connections < each_Nodes_Connection[node_start_Index + node].size(); connections++)
            {
                if (each_Nodes_Connection[node_start_Index + node][connections].first == cave_ID)
                {
                    if (each_Nodes_Connection[node_start_Index + node][connections].second > node)
                    {
                        network_File << cave_ID << "_" << to_string(node) << "\t"
                                     << cave_ID << "_" << to_string(each_Nodes_Connection[node_start_Index + node][connections].second) << "\n";
                    }
                }
                else
                {
                    network_File << to_string(cave_ID) << "_" << to_string(node) << "\t"
                                 << to_string(each_Nodes_Connection[node_start_Index + node][connections].first) << "_" << to_string(each_Nodes_Connection[node_start_Index + node][connections].second) << "\n";
                }
            }
        }
    }

    network_File.close();

    // for (int cave_ID = 0; cave_ID < SCM_number_of_caves; cave_ID++)
    // {
    //     int node_start_Index = cave_ID * SCM_number_of_nodes_per_cave;
    //     for (size_t i = node_start_Index; i < (node_start_Index + SCM_number_of_nodes_per_cave); i++)
    //     {
    //         for (int x = 0; x < each_Nodes_Connection[i].size(); x++)
    //         {
    //             cout << each_Nodes_Connection[i][x].first << "_" << each_Nodes_Connection[i][x].second << ", ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    cout << endl;
}

void simulator_Master::BA_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions)
{
    cout << "Intializing Barbasi Albert model network engine\n";

    // string output_Network_location = this->output_Folder_location + "/network_Data";
    // functions.config_Folder(output_Network_location, "Network");
    // network_File_location = output_Network_location + "/node_node_Relationships.csv";
    // functions.create_File(network_File_location, "Source\tTarget");

    // for (int initialize = 0; initialize < Total_number_of_Nodes; initialize++)
    // {
    //     vector<pair<int, int>> intialize_Vector;
    //     each_Nodes_Connection.push_back(intialize_Vector);
    // }

    cout << "Forming node to node relationships\n";

    vector<int> connections_per_Node;
    connections_per_Node.push_back(1);
    cout << "Configuring node 1 of " << Total_number_of_Nodes << " node(s)" << endl;

    int tot_connections = 2;

    // int iterative_Attach = BA_FIXED;

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    random_device rd;
    mt19937 gen(rd());

    for (int node = 1; node < Total_number_of_Nodes; node++)
    {
        cout << "Configuring node " << (node + 1) << " of " << Total_number_of_Nodes << " node(s)" << endl;
        connections_per_Node.push_back(0);

        int iterative_Attach;
        if (connection_Model == "FIXED")
        {
            iterative_Attach = BA_FIXED;
        }
        else if (connection_Model == "NEGATIVE BINOMIAL")
        {
            negative_binomial_distribution<int> distribution(BA_NB_sucesses, BA_NB_probability);
            iterative_Attach = distribution(gen);
        }
        else if (connection_Model == "POISSON")
        {
            poisson_distribution<int> distribution(BA_Poisson_mean);
            iterative_Attach = distribution(gen);
        }
        else
        {
            cout << "ERROR in network parameter \"BA model fixed new connections\". CHECK PARAMETER\n";
            exit(-1);
        }

        iterative_Attach = (iterative_Attach < 1) ? 1 : iterative_Attach;

        for (int interative = 0; interative < iterative_Attach; interative++)
        {
            tot_connections++;
            int attach_Node = node;

            while (attach_Node == node)
            {
                float randomNum = static_cast<float>(std::rand()) / RAND_MAX;
                float cum_Prob = 0;

                for (int check_Node = 0; check_Node < connections_per_Node.size(); check_Node++)
                {
                    cum_Prob += ((float)(connections_per_Node[check_Node] + 1) / (float)tot_connections);
                    if (randomNum < cum_Prob)
                    {
                        attach_Node = check_Node;
                        break;
                    }
                }
            }

            if (attach_Node != -1)
            {
                int check_Present = 0;

                for (int check = 0; check < each_Nodes_Connection[node].size(); check++)
                {
                    if (attach_Node == each_Nodes_Connection[node][check].second)
                    {
                        check_Present = 1;
                        break;
                    }
                }

                if (check_Present == 0)
                {
                    // cout << "Node " << node + 1 << " attached to " << attach_Node + 1 << endl;
                    network_File << "0_" << to_string(attach_Node + 1) << "\t"
                                 << "0_" << to_string(node + 1) << "\n";

                    connections_per_Node[attach_Node] = connections_per_Node[attach_Node] + 1;

                    each_Nodes_Connection[attach_Node].push_back(make_pair(0, node));
                    each_Nodes_Connection[node].push_back(make_pair(0, attach_Node));

                    tot_connections++;
                }
            }
            else
            {
                cout << "ERROR in node configuration\n";
                exit(-1);
            }
        }
    }

    network_File.close();
}