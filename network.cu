#include "network.cuh"
#include "functions_library.cuh"
#include "parameter_load.h"

network::network(int CUDA_device_number, int CPU_cores, int gpu_Limit, string multi_READ)
{
    cout << "\nNetwork Sim\n";

    functions_library function = functions_library();

    this->CPU_cores = CPU_cores;
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    this->multi_READ = multi_READ;
    cout << "Multiple read and write: " << this->multi_READ << endl
         << endl;

    this->CUDA_device_number = CUDA_device_number;
    function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

    this->gpu_Limit = gpu_Limit;

    cout << "Per round GPU max unit: " << this->gpu_Limit << endl
         << endl;
}

void network::ingress()
{
    functions_library function = functions_library();

    int nodes = 100;

    cout << "Barabasi Albert Model\n\nNodes: " << nodes << endl;

    vector<int> connections_per_Node;

    int tot_connections = 0;

    string generation_population_Network = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/generation_Summary_File.csv";
    function.config_File_delete_create(generation_population_Network, "Source\tTarget");
    fstream network_Write;
    network_Write.open(generation_population_Network, ios::app);

    for (int node = 0; node < nodes; node++)
    {
        connections_per_Node.push_back(0);
        tot_connections++;

        // randomNum = 0.99;

        // cout << randomNum << endl;

        int attach_Node = -1;

        do
        {
            float randomNum = static_cast<float>(std::rand()) / RAND_MAX;
            float cum_Prob = 0;
            for (int check_Node = 0; check_Node < connections_per_Node.size(); check_Node++)
            {
                cum_Prob += ((float)(connections_per_Node[check_Node] + 1) / (float)tot_connections);
                // cout << check_Node + 1 << ": " << connections_per_Node[check_Node] + 1 << "/" << tot_connections << "=" << ((float)(connections_per_Node[check_Node] + 1) / (float)tot_connections) << endl;
                // cout << "Tot_prob: " << cum_Prob << " : " << randomNum << endl;
                if (randomNum < cum_Prob)
                {
                    attach_Node = check_Node;
                    break;
                }
            }
            // cum_Prob = 0;
        } while (attach_Node == node && node != 0);

        if (attach_Node != -1)
        {
            if (node != 0)
            {
                cout << "Node " << node + 1 << " attached to " << attach_Node + 1 << endl;
                network_Write << to_string(attach_Node + 1) << "\t" << to_string(node + 1) << "\n";
            }
            connections_per_Node[attach_Node] = connections_per_Node[attach_Node] + 1;
            tot_connections++;
        }
        else
        {
            cout << "ERROR\n";
            exit(-1);
        }
    }

    network_Write.close();

    cout << endl;

    string node_Summary_Path = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/node_Summary_File.csv";
    function.config_File_delete_create(node_Summary_Path, "ID\tNumber_of_links");

    fstream node_Summary;
    node_Summary.open(node_Summary_Path, ios::app);

    node_Summary << "1\t" << to_string(connections_per_Node[0] - 1) << "\n";

    for (size_t i = 1; i < connections_per_Node.size(); i++)
    {
        node_Summary << to_string(i + 1) << "\t" << to_string(connections_per_Node[i] + 1) << "\n";
    }

    node_Summary.close();
}