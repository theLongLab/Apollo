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
        "\"Network profile\""};

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
    vector<int> node_cave_IDs;

    network_Manager(each_Nodes_Connection, functions, node_cave_IDs);

    cout << "STEP 2: Infection of Population\n\n";
}

void simulator_Master::network_Manager(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions, vector<int> &node_cave_IDs)
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
        SCM_Model_Engine(each_Nodes_Connection, functions, node_cave_IDs);
    }
    else
    {
        cout << "ERROR Incorrect network selected. Please check \"Network type\" in the network parameter file.\n";
        exit(-1);
    }
}

void simulator_Master::SCM_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions, vector<int> &node_cave_IDs)
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
