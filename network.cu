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
    int nodes = 10;

    cout << "Barabasi Albert Model\n\nNodes: " << nodes << endl;

    vector<int> connections_per_Node;

    int tot_connections = 0;

    for (int node = 0; node < nodes; node++)
    {
        connections_per_Node.push_back(0);
        tot_connections++;

        float randomNum = static_cast<float>(std::rand()) / RAND_MAX;

        // cout << randomNum << endl;

        int attach_Node = -1;
        float cum_Prob = 0;

        for (int check_Node = 0; check_Node < connections_per_Node.size(); check_Node++)
        {
            cum_Prob += cum_Prob + ((float)(connections_per_Node[check_Node] + 1) / (float)tot_connections);
            // cout << tot_connections << endl;
            if (randomNum < cum_Prob)
            {
                attach_Node = check_Node;
                break;
            }
        }

        if (attach_Node != -1)
        {
            cout << "Node " << node + 1 << " attached to " << attach_Node + 1 << endl;
            connections_per_Node[attach_Node] + connections_per_Node[attach_Node] + 1;
            tot_connections++;
        }
    }
}