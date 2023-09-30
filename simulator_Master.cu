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

    if (Parameters.get_STRING(found_Parameters[0]) == "BA MODEL")
    {
        cout << "\nBarabsi Albert model selected: \n";
        network_Model = "BA";

        parameters_List = {"\"BA model number of nodes\"",
                           "\"BA model standard new connections\""};
        found_Parameters = Parameters.get_parameters(network_Profile_File, parameters_List);

        number_of_Nodes_BA = Parameters.get_INT(found_Parameters[0]);
        cout << "Number of nodes: " << number_of_Nodes_BA << endl;
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

    string network_Summary_Location, network_node_Location;

    // ! Compatible for both BA and Caveman
    // INT, INT = Cave_ID and Node, for BA CaveID is 0 for all.
    vector<vector<pair<int, int>>> each_Nodes_Connection;

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
    else
    {
        cout << "ERROR Incorrect network selected. Please check \"Network type\" in the network parameter file.\n";
        exit(-1);
    }
}

void simulator_Master::BA_Model_Engine(vector<vector<pair<int, int>>> &each_Nodes_Connection, functions_library &functions)
{
    cout << "Intializing Barbasi Albert model network engine\n";

    string output_Network_location = this->output_Folder_location + "/network_Data";
    functions.config_Folder(output_Network_location, "Network");
    network_File_location = output_Network_location + "/node_node_Relationships.csv";
    functions.create_File(network_File_location, "Source\tTarget");

    for (int initialize = 0; initialize < number_of_Nodes_BA; initialize++)
    {
        vector<pair<int, int>> intialize_Vector;
        each_Nodes_Connection.push_back(intialize_Vector);
    }

    cout << "Forming node to node relationships\n";

    vector<int> connections_per_Node;
    connections_per_Node.push_back(1);
    cout << "Configuring node 1 of " << number_of_Nodes_BA << " node(s)" << endl;

    int tot_connections = 2;

    // int iterative_Attach = BA_FIXED;

    fstream network_File;
    network_File.open(network_File_location, ios::app);

    random_device rd;
    mt19937 gen(rd());

    for (int node = 1; node < number_of_Nodes_BA; node++)
    {
        cout << "Configuring node " << (node + 1) << " of " << number_of_Nodes_BA << " node(s)" << endl;
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
                    network_File << to_string(attach_Node + 1) << "\t" << to_string(node + 1) << "\n";

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
