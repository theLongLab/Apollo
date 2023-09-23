#include "network.cuh"
#include "functions_library.cuh"
#include "parameter_load.h"

network::network(int CUDA_device_number, int CPU_cores, int gpu_Limit, string multi_READ, int activate)
{
    cout << "\nNetwork Sim\n";

    functions_library function = functions_library();

    this->CPU_cores = CPU_cores;
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    this->multi_READ = multi_READ;
    cout << "Multiple read and write: " << this->multi_READ << endl
         << endl;

    if (activate == 1)
    {
        this->CUDA_device_number = CUDA_device_number;
        function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

        this->gpu_Limit = gpu_Limit;

        cout << "Per round GPU max unit: " << this->gpu_Limit << endl
             << endl;
    }
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

void network::ingress_2()
{
    functions_library function = functions_library();

    int caves = 10;
    int avg_size_of_cave = 10;
    // int outside_world_nodes = 1;

    cout << "Connected Caveman models\n\nNodes: " << caves << endl;

    // parent_child
    vector<vector<pair<int, int>>> connections_All;

    // for (size_t i = 0; i < 100; i++)
    // {
    //     cout << dist(gen) << endl;
    // }

    // exit(-1);

    // create cave connections
    string node_Summary_Path = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/node_Summary_File.csv";
    function.config_File_delete_create(node_Summary_Path, "ID\tCave_ID");
    fstream node_Summary;
    node_Summary.open(node_Summary_Path, ios::app);

    for (int cave = 0; cave < caves; cave++)
    {

        vector<pair<int, int>> connections_of_cave;

        for (int node_parent = 0; node_parent < avg_size_of_cave; node_parent++)
        {
            node_Summary << cave << "_" << node_parent << "\t" << cave << "\n";
            for (int node_child = node_parent + 1; node_child < avg_size_of_cave; node_child++)
            {
                connections_of_cave.push_back(make_pair(node_parent, node_child));
            }
        }

        connections_All.push_back(connections_of_cave);
    }

    node_Summary.close();

    // Define the distribution
    // Seed the random number generator
    random_device rd;
    mt19937 gen(rd());

    vector<int> parent_outside;

    vector<pair<int, int>> connections_of_cave;
    connections_of_cave = connections_All[0];

    uniform_int_distribution<int> dist(0, connections_of_cave.size() - 1);
    int select_Connection = dist(gen);

    parent_outside.push_back(connections_of_cave[select_Connection].first);
    connections_of_cave.erase(connections_of_cave.begin() + select_Connection);
    connections_All[0] = connections_of_cave;

    vector<pair<int, int>> redirect_outside_Connections;

    for (int redirect_Cave = 1; redirect_Cave < connections_All.size(); redirect_Cave++)
    {
        connections_of_cave = connections_All[redirect_Cave];

        uniform_int_distribution<int> dist(0, connections_of_cave.size() - 1);
        select_Connection = dist(gen);

        parent_outside.push_back(connections_of_cave[select_Connection].first);
        connections_of_cave.erase(connections_of_cave.begin() + select_Connection);
        connections_All[redirect_Cave] = connections_of_cave;

        redirect_outside_Connections.push_back(make_pair(parent_outside[redirect_Cave - 1], parent_outside[redirect_Cave]));
    }

    redirect_outside_Connections.push_back(make_pair(parent_outside[parent_outside.size() - 1], parent_outside[0]));

    string generation_population_Network = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/generation_Summary_File.csv";
    function.config_File_delete_create(generation_population_Network, "Source\tTarget");

    fstream network_Write;
    network_Write.open(generation_population_Network, ios::app);

    for (size_t i = 0; i < connections_All.size(); i++)
    {
        vector<pair<int, int>> connections_of_cave;
        connections_of_cave = connections_All[i];

        for (size_t x = 0; x < connections_of_cave.size(); x++)
        {
            network_Write << i << "_" << connections_of_cave[x].first << "\t" << i << "_" << connections_of_cave[x].second << "\n";
            cout << i << "_" << connections_of_cave[x].first << "\t" << i << "_" << connections_of_cave[x].second << endl;
        }
        cout << endl;

        int outside = i + 1;

        if (outside >= connections_All.size())
        {
            outside = 0;
        }

        cout << "outside: ";
        cout << i << "_" << redirect_outside_Connections[i].first << "\t" << outside << "_" << redirect_outside_Connections[i].second << endl;
        network_Write << i << "_" << redirect_outside_Connections[i].first << "\t" << outside << "_" << redirect_outside_Connections[i].second << "\n";

        cout << endl;
    }

    network_Write.close();
}

void network::ingress_flexible_caveman()
{

    functions_library function = functions_library();

    int cave_number = 10;

    int *per_cave_Stride = (int *)malloc((cave_number + 1) * sizeof(int));
    int **network_Array = function.create_INT_2D_arrays(cave_number, 3);

    vector<vector<int>> neighbour_Nodes;
    vector<vector<int>> global_Nodes;

    per_cave_Stride[0] = 0;

    cout << "Dynamic Caveman models\n\nNodes: " << cave_number << endl
         << endl;

    random_device rd;
    mt19937 gen(rd());

    float shape_Nodes_per_cave = 10;
    float scale_Nodes_per_cave = 2;

    gamma_distribution<float> gamma(shape_Nodes_per_cave, scale_Nodes_per_cave);

    for (size_t i = 0; i < cave_number; i++)
    {
        int number_of_Nodes = (int)round(gamma(gen));
        // cout << number_of_Nodes << endl;
        network_Array[i][0] = number_of_Nodes;
        per_cave_Stride[i + 1] = per_cave_Stride[i] + network_Array[i][0];

        vector<int> intialize;
        neighbour_Nodes.push_back(intialize);
        global_Nodes.push_back(intialize);
    }

    // min 1 needed
    float percent_Neighbouring = 0.25;
    float percent_Global_freedom = 0.35;

    cout << "Determining node counts and neighour and global nodes\n\n";

    for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    {
        binomial_distribution<int> binomialDist_Neighbour(network_Array[cave_ID][0], percent_Neighbouring);
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

        // for (int test_Get = 0; test_Get < network_Array[cave_ID][1]; test_Get++)
        // {
        //     neighbour_Nodes[cave_ID].push_back(distribution_Neighbour(gen));
        // }

        binomial_distribution<int> binomialDist_Global(network_Array[cave_ID][1], percent_Global_freedom);
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

        // for (int test_Get = 0; test_Get < network_Array[cave_ID][2]; test_Get++)
        // {
        //     int get_Index = distribution_Global(gen);
        //     global_Nodes.push_back(neighbour_Nodes[cave_ID][get_Index]);
        // }
    }

    // for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    // {
    //     for (int column = 0; column < 3; column++)
    //     {
    //         cout << network_Array[cave_ID][column] << "\t";
    //     }
    //     cout << "\n";
    // }

    // cout << "\n";

    // for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    // {
    //     cout << "Local nodes: ";
    //     for (int local = 0; local < neighbour_Nodes[cave_ID].size(); local++)
    //     {
    //         cout << cave_ID << "_" << neighbour_Nodes[cave_ID][local] << " ";
    //     }
    //     cout << endl;

    //     cout << "Global nodes: ";
    //     for (int global = 0; global < global_Nodes[cave_ID].size(); global++)
    //     {
    //         cout << cave_ID << "_" << global_Nodes[cave_ID][global] << " ";
    //     }
    //     cout << endl;
    //     cout << endl;
    // }

    // exit(-1);

    vector<vector<pair<int, int>>> each_Nodes_Connections;

    int total_Nodes = per_cave_Stride[cave_number];
    cout << "Total nodes in the network: " << total_Nodes << endl
         << endl;

    for (size_t i = 0; i < total_Nodes; i++)
    {
        vector<pair<int, int>> nodes_Intialize;
        each_Nodes_Connections.push_back(nodes_Intialize);
    }

    // cout << each_Nodes_Connections.size() << endl;

    // exit(-1);

    string generation_population_Network = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/generation_Summary_File_flex.csv";
    function.config_File_delete_create(generation_population_Network, "Source\tTarget");
    fstream network_Write;
    network_Write.open(generation_population_Network, ios::app);

    string node_Summary_Path = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/node_Summary_File_flex.csv";
    function.config_File_delete_create(node_Summary_Path, "ID\tCave_ID");
    fstream node_Summary;
    node_Summary.open(node_Summary_Path, ios::app);

    // 1: 1 to 2
    // 2: 2 to 1

    // interconnect the caves nodes
    // Assign nodes their relationships

    cout << "Configuring local networks\n";

    for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    {
        int nodes_Total = network_Array[cave_ID][0];

        // cout << cave_ID << endl;

        for (int node_parent = 0; node_parent < nodes_Total; node_parent++)
        {

            int node_Count = per_cave_Stride[cave_ID] + node_parent;
            // cout << node_Count << endl;

            node_Summary << cave_ID << "_" << node_parent << "\t" << cave_ID << "\n";
            for (int node_child = node_parent + 1; node_child < nodes_Total; node_child++)
            {
                int node_child_Main = node_Count + (node_child - node_parent);
                network_Write << cave_ID << "_" << node_parent << "\t" << cave_ID << "_" << node_child << "\n";
                each_Nodes_Connections[node_Count].push_back(make_pair(cave_ID, node_child));
                each_Nodes_Connections[node_child_Main].push_back(make_pair(cave_ID, node_parent));
                // cout << node_Count << "\t" << node_child_Main << endl;
            }
        }
        // cout << endl;
    }

    node_Summary.close();

    cout << "Configuring Neighbouring networks attachments\n";

    for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    {
        uniform_int_distribution<int> distribution_Neighbour_parent(0, neighbour_Nodes[cave_ID].size() - 1);

        int parent_Node_Index = neighbour_Nodes[cave_ID][distribution_Neighbour_parent(gen)];

        int attach_Cave = cave_ID + 1;
        int attach_Node_Index;
        if (attach_Cave != cave_number)
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

        network_Write << cave_ID << "_" << parent_Node_Index << "\t" << attach_Cave << "_" << attach_Node_Index << "\n";

        int parent_Node_all_Index = per_cave_Stride[cave_ID] + parent_Node_Index;
        int attach_Node_all_Index = per_cave_Stride[attach_Cave] + attach_Node_Index;

        each_Nodes_Connections[parent_Node_all_Index].push_back(make_pair(attach_Cave, attach_Node_Index));
        each_Nodes_Connections[attach_Node_all_Index].push_back(make_pair(cave_ID, parent_Node_Index));
    }

    cout << "Configuring Global networks attachments\n";
    cout << "Global processing caves: ";
    for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    {
        cout << cave_ID;
        if (cave_ID + 1 != cave_number)
        {
            cout << ", ";
        }
        for (size_t i = 0; i < global_Nodes[cave_ID].size(); i++)
        {

            int cave_Global;
            do
            {
                uniform_int_distribution<int> distribution_global_CAVEs(0, cave_number - 1);
                cave_Global = distribution_global_CAVEs(gen);
                // cout << cave_Global << endl;
            } while (cave_Global == cave_ID);

            if (global_Nodes[cave_Global].size() > 0)
            {
                uniform_int_distribution<int> distribution_Parent_attach(0, global_Nodes[cave_ID].size() - 1);
                uniform_int_distribution<int> distribution_Global_attach(0, global_Nodes[cave_Global].size() - 1);

                int parent_Cave_Node = global_Nodes[cave_ID][distribution_Parent_attach(gen)];
                int Global_Cave_Node = global_Nodes[cave_Global][distribution_Parent_attach(gen)];

                network_Write << cave_ID << "_" << parent_Cave_Node << "\t" << cave_Global << "_" << Global_Cave_Node << "\n";

                int parent_Node_all_Index = per_cave_Stride[cave_ID] + parent_Cave_Node;
                int attach_Node_all_Index = per_cave_Stride[cave_Global] + Global_Cave_Node;

                each_Nodes_Connections[parent_Node_all_Index].push_back(make_pair(cave_Global, Global_Cave_Node));
                each_Nodes_Connections[attach_Node_all_Index].push_back(make_pair(cave_ID, parent_Cave_Node));
            }
        }
    }

    cout << endl;

    network_Write.close();

    // cout << "hello\n";

    // exit(-1);

    // for (int cave_ID = 0; cave_ID < cave_number; cave_ID++)
    // {
    //     cout << "CAVE ID: " << cave_ID << endl
    //          << endl;
    //     for (size_t i = per_cave_Stride[cave_ID]; i < per_cave_Stride[cave_ID + 1]; i++)
    //     {
    //         for (int count = 0; count < each_Nodes_Connections[i].size(); count++)
    //         {
    //             cout << each_Nodes_Connections[i][count].first << "_" << each_Nodes_Connections[i][count].second << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
}

void network::sim_cLD()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    random_device rd;
    mt19937 generator(rd());

    // default_random_engine generator;

    cout << "cLD simulator\n"
         << endl;

    int generations = 50;

    int eff_Population = 10;

    // PROGENY BINOMIAL
    int n = 4;
    float prob = 0.75;
    binomial_distribution<int> binomialDist(n, prob);

    // cout << x << endl;

    float mutation_Rate = 0.01;
    poisson_distribution<int> dist_Poisson(mutation_Rate);

    int mutation_points = 10;

    float recombination_Prob = 0.01;
    int interactions = 0;

    int *parent_IDs = (int *)malloc((eff_Population) * sizeof(int));

    float **pop_GeneA_GeneB_Parent = function.create_Fill_2D_array_FLOAT(eff_Population, 3, 0);

    vector<vector<pair<int, int>>> mutation_points_All_Parents;

    cout << "Configure ideal population\n";
    for (size_t i = 0; i < eff_Population / 2; i++)
    {
        pop_GeneA_GeneB_Parent[i][0] = 0;
        pop_GeneA_GeneB_Parent[i][1] = 0;

        vector<pair<int, int>> mutation_points_store;

        for (int point = 0; point < mutation_points; point++)
        {
            mutation_points_store.push_back(make_pair(0, 0));
        }

        mutation_points_All_Parents.push_back(mutation_points_store);
    }

    uniform_int_distribution<int> mutation_point_Draw(0, mutation_points - 1);

    for (size_t i = eff_Population / 2; i < eff_Population; i++)
    {
        pop_GeneA_GeneB_Parent[i][0] = 1;
        pop_GeneA_GeneB_Parent[i][1] = 1;

        vector<pair<int, int>> mutation_points_store;

        for (int point = 0; point < mutation_points; point++)
        {
            mutation_points_store.push_back(make_pair(0, 0));
        }

        for (int location_Draw = 0; location_Draw < 1; location_Draw++)
        {
            int location = mutation_point_Draw(generator);
            mutation_points_store[location].first = mutation_points_store[location].first + 1;
        }

        for (int location_Draw = 0; location_Draw < 1; location_Draw++)
        {
            int location = mutation_point_Draw(generator);
            mutation_points_store[location].second = mutation_points_store[location].second + 1;
        }

        mutation_points_All_Parents.push_back(mutation_points_store);
    }

    float Pa = 0;
    float Pb = 0;
    float Pab = 0;

    for (size_t i = 0; i < eff_Population; i++)
    {
        pop_GeneA_GeneB_Parent[i][2] = 0.80;
        parent_IDs[i] = i;

        if (pop_GeneA_GeneB_Parent[i][0] != 0)
        {
            Pa++;
        }
        if (pop_GeneA_GeneB_Parent[i][1] != 0)
        {
            Pb++;
        }
        if (pop_GeneA_GeneB_Parent[i][0] != 0 && pop_GeneA_GeneB_Parent[i][1] != 0)
        {
            Pab++;
        }
    }

    Pa = Pa / eff_Population;
    Pb = Pb / eff_Population;
    Pab = Pab / eff_Population;

    // Pa = Pa / sum_Progeny;
    // Pb = Pb / sum_Progeny;
    // Pab = Pab / sum_Progeny;

    string cLD_write = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/cLD_check_0.010_No.csv";
    function.config_File_delete_create(cLD_write, "Generation\tPa\tPb\tPab\tcLD");

    fstream cLD_writer;
    cLD_writer.open(cLD_write, ios::app);
    cLD_writer.precision(9);

    cout << "\nPa: " << Pa << "\tPb: " << Pb << "\tPab: " << Pab << endl;

    float cLD = pow((Pab - (Pa * Pb)), 2) / (Pa * (1 - Pa) * Pb * (1 - Pb));

    cout << "cLD: " << cLD << endl
         << endl;

    cLD_writer << to_string(0) << "\t" << to_string(Pa) << "\t" << to_string(Pb) << "\t" << to_string(Pab) << "\t" << to_string(cLD) << "\n";

    for (int gen = 0; gen < generations; gen++)
    {
        cout << "Processing generation " << gen + 1 << " of " << generations << " generations\n"
             << endl;

        random_shuffle(&parent_IDs[0], &parent_IDs[eff_Population]);

        int **parent_Pairs = function.create_INT_2D_arrays(eff_Population / 2, 3);
        int sum_Progeny = 0;

        cout << "Assigning parent pairs and children numbers\n";
        int increment = eff_Population / 2;
        for (size_t i = 0; i < eff_Population / 2; i++)
        {
            parent_Pairs[i][0] = parent_IDs[i];
            parent_Pairs[i][1] = parent_IDs[i + increment];
            int x = binomialDist(generator);
            sum_Progeny = sum_Progeny + x;
            parent_Pairs[i][2] = x;
        }

        cout << "Number of children to be simulated: " << sum_Progeny << endl;
        float **pop_GeneA_GeneB_Progeny = function.create_Fill_2D_array_FLOAT(sum_Progeny, 3, 0);
        // cout << "TEST 1" << endl;
        int progeny_Fill_count = 0;

        vector<int> progeny_surviving_ID;

        Pa = 0;
        Pb = 0;
        Pab = 0;

        vector<vector<pair<int, int>>> mutation_points_All_Progeny;

        // cout << "TEST 2" << endl;

        for (int parent_pair = 0; parent_pair < eff_Population / 2; parent_pair++)
        {
            // cout << "\nProcessing parent pair: " << parent_pair << " of " << eff_Population / 2 << endl;

            for (int progeny = 0; progeny < parent_Pairs[parent_pair][2]; progeny++)
            {
                vector<pair<int, int>> mutation_points_store_Progeny;

                // cout << "\nProcessing progeny: " << progeny << " of " << parent_Pairs[parent_pair][2] / 2;
                uniform_real_distribution<> dis(0.0, 1.0);

                int mom_0_dad_1 = 0;
                if (dis(generator) < 0.5)
                {
                    // cout << dis(generator) << endl;
                    mom_0_dad_1 = 1;
                }

                // pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][0];
                // pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][1];

                vector<pair<int, int>> mutation_points_store_Parent = mutation_points_All_Parents[parent_Pairs[parent_pair][mom_0_dad_1]];

                // cout << parent_Pairs[parent_pair][mom_0_dad_1] << endl;
                // cout << mutation_points_All_Parents.size() << endl;
                //  cout << mutation_points_store_Parent.size() << endl;

                // cout << "TEST 3" << endl;
                for (int mutation_Fill = 0; mutation_Fill < mutation_points; mutation_Fill++)
                {
                    // cout << mutation_points << endl;
                    int A_value = mutation_points_store_Parent[mutation_Fill].first;
                    int B_value = mutation_points_store_Parent[mutation_Fill].second;
                    // cout << A_value << "\t" << B_value << endl;
                    mutation_points_store_Progeny.push_back(make_pair(A_value, B_value));
                    // cout << mutation_Fill << endl;
                }

                mutation_points_store_Parent.clear();
                // cout << "TEST 4" << endl;
                //  RECOMBINATION A

                for (int mutation_Fill = 0; mutation_Fill < mutation_points; mutation_Fill++)
                {
                    if (dis(generator) < recombination_Prob)
                    {
                        mom_0_dad_1 = 0;
                        if (dis(generator) < 0.5)
                        {
                            // cout << "dad recombinant A" << endl;
                            mom_0_dad_1 = 1;
                        }

                        mutation_points_store_Parent = mutation_points_All_Parents[parent_Pairs[parent_pair][mom_0_dad_1]];
                        int A_value = mutation_points_store_Parent[mutation_Fill].first;
                        mutation_points_store_Progeny[mutation_Fill].first = A_value;
                        mutation_points_store_Parent.clear();

                        // pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][0];
                    }

                    // RECOMBINATION B
                    if (dis(generator) < recombination_Prob)
                    {
                        mom_0_dad_1 = 0;
                        if (dis(generator) < 0.5)
                        {
                            // cout << "dad recombinant B" << endl;
                            mom_0_dad_1 = 1;
                        }

                        mutation_points_store_Parent = mutation_points_All_Parents[parent_Pairs[parent_pair][mom_0_dad_1]];
                        int B_value = mutation_points_store_Parent[mutation_Fill].second;
                        mutation_points_store_Progeny[mutation_Fill].second = B_value;
                        mutation_points_store_Parent.clear();

                        // pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][1];
                    }
                }

                // cout << "TEST 5" << endl;

                pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = 0.80;

                if (interactions != 1)
                {
                    for (int mutation_Fill = 0; mutation_Fill < mutation_points; mutation_Fill++)
                    {
                        if (mutation_points_store_Progeny[mutation_Fill].first > 0 && mutation_points_store_Progeny[mutation_Fill].second > 0)
                        {
                            if (mutation_points_store_Progeny[mutation_Fill].first == mutation_points_store_Progeny[mutation_Fill].second)
                            {
                                // if (interactions == 1)
                                // {
                                //     pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 0.05;
                                // }
                                // else
                                // {
                                pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.1;
                                //}
                            }
                        }
                        else if (mutation_points_store_Progeny[mutation_Fill].first > 0 && mutation_points_store_Progeny[mutation_Fill].second == 0)
                        {
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.05;
                        }
                        else if (mutation_points_store_Progeny[mutation_Fill].first == 0 && mutation_points_store_Progeny[mutation_Fill].second > 0)
                        {
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.05;
                        }
                    }

                    // if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] > 0.9)
                    // {
                    //     pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = 0.9;
                    // }
                    // else
                    if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] < 0)
                    {
                        pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = 0;
                    }
                }

                // if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] == pop_GeneA_GeneB_Progeny[progeny_Fill_count][1])
                // {
                //     if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] + pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] != 0)
                //     {
                //         if (interactions == 1)
                //         {
                //             pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 0.10;
                //         }
                //         else
                //         {
                //             pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                //         }
                //     }
                // }
                // else
                // {
                //     pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                // }

                // mutate
                // reduce those with mismatc mutations survivability

                vector<int> A_mutations_NEW;
                vector<int> B_mutations_NEW;

                for (size_t mutation_Fill = 0; mutation_Fill < mutation_points; mutation_Fill++)
                {
                    A_mutations_NEW.push_back(0);
                    B_mutations_NEW.push_back(0);
                }

                int mutations_A = dist_Poisson(generator);
                int mutations_B = dist_Poisson(generator);

                for (int mutation_Fill = 0; mutation_Fill < mutations_A; mutation_Fill++)
                {
                    int mut_Loc = mutation_point_Draw(generator);
                    mutation_points_store_Progeny[mut_Loc].first = mutation_points_store_Progeny[mut_Loc].first + 1;
                    A_mutations_NEW[mut_Loc] = A_mutations_NEW[mut_Loc] + 1;
                    // mutation_points_store_Progeny[mutation_Fill].second = mutation_points_store_Progeny[mutation_Fill].second + mutations_B;
                }

                for (int mutation_Fill = 0; mutation_Fill < mutations_B; mutation_Fill++)
                {
                    int mut_Loc = mutation_point_Draw(generator);
                    mutation_points_store_Progeny[mut_Loc].second = mutation_points_store_Progeny[mut_Loc].second + 1;
                    B_mutations_NEW[mut_Loc] = B_mutations_NEW[mut_Loc] + 1;
                    // mutation_points_store_Progeny[mutation_Fill].second = mutation_points_store_Progeny[mutation_Fill].second + mutations_B;
                }

                // pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] + mutations_A;
                // pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] + mutations_B;
                //  int check_A_B = 0;

                // if (mutations_A != 0)
                // {
                //     Pa++;
                //     // pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] + mutations_A;
                //     //  check_A_B++;
                // }

                // if (mutations_B != 0)
                // {
                //     Pb++;
                //     // pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] + mutations_B;
                //     //  check_A_B++;
                // }

                // if (mutations_A != 0 && mutations_B != 0)
                // {
                //     Pab++;
                // }

                // if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] != 0)
                // {
                //     Pa++;
                // }

                // if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] != 0)
                // {
                //     Pb++;
                // }

                // if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] != 0 && pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] != 0)
                // {
                //     Pab++;
                // }

                // ADD regions for mutations and use that for interactions

                // if (interactions == 1)
                // {
                //     if (mutations_A == mutations_B)
                //     {
                //         if ((mutations_A + mutations_B) != 0)
                //         {
                //             pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 0.10;
                //         }
                //     }
                //     else
                //     {
                //         pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                //     }
                // }
                // else
                // {
                //     if ((mutations_A + mutations_B) != 0)
                //     {
                //         pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                //     }
                // }

                for (int mutation_Fill = 0; mutation_Fill < mutation_points; mutation_Fill++)
                {
                    if (A_mutations_NEW[mutation_Fill] > 0 && B_mutations_NEW[mutation_Fill] > 0)
                    {
                        if (interactions == 1)
                        {
                            // cout << "HIT 1" << endl;
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 20;
                        }
                        else
                        {
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.1;
                        }
                    }
                    else if (A_mutations_NEW[mutation_Fill] > 0 && B_mutations_NEW[mutation_Fill] == 0)
                    {
                        // cout << "HIT 2" << endl;
                        pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.05;
                    }
                    else if (A_mutations_NEW[mutation_Fill] == 0 && B_mutations_NEW[mutation_Fill] > 0)
                    {
                        // cout << "HIT 3" << endl;
                        pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.05;
                    }
                }

                if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] > 1)
                {
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = 1;
                }
                else if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] < 0)
                {
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = 0.00;
                }

                pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = 0;
                pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = 0;

                for (int mutation_Fill = 0; mutation_Fill < mutation_points; mutation_Fill++)
                {
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] + mutation_points_store_Progeny[mutation_Fill].first;
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] + mutation_points_store_Progeny[mutation_Fill].second;
                }

                // if (mutations_A != 0 || mutations_B != 0)
                // {
                //     if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] == pop_GeneA_GeneB_Progeny[progeny_Fill_count][1])
                //     {
                //         if (interactions == 1)
                //         {
                //             pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 0.25;
                //         }
                //         else
                //         {
                //             pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.25;
                //         }
                //         // cout << mutations_A << "\t" << mutations_B << "\n";
                //     }
                //     else
                //     {
                //         pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.25;
                //     }
                // }

                // cout << pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] << "\t";
                bernoulli_distribution distribution(pop_GeneA_GeneB_Progeny[progeny_Fill_count][2]);

                // Flip the coin and output the result
                bool result = distribution(generator);
                if (result)
                {
                    progeny_surviving_ID.push_back(progeny_Fill_count);
                    // cout << "Yes\n";
                    //  cout << result << " Heads" << endl;
                    if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] != 0)
                    {
                        Pa++;
                    }
                    if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] != 0)
                    {
                        Pb++;
                    }
                    if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] != 0 && pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] != 0)
                    {
                        Pab++;
                    }
                }
                // else
                // {
                //     // survive = 0;
                //     //cout << "No\n";
                //     // cout << result << " Tails" << endl;
                // }

                progeny_Fill_count++;
                mutation_points_All_Progeny.push_back(mutation_points_store_Progeny);
            }
        }

        Pa = Pa / progeny_surviving_ID.size();
        Pb = Pb / progeny_surviving_ID.size();
        Pab = Pab / progeny_surviving_ID.size();

        // Pa = Pa / sum_Progeny;
        // Pb = Pb / sum_Progeny;
        // Pab = Pab / sum_Progeny;

        cout << "\nPa: " << Pa << "\tPb: " << Pb << "\tPab: " << Pab << endl;

        cLD = pow((Pab - (Pa * Pb)), 2) / (Pa * (1 - Pa) * Pb * (1 - Pb));

        cout << "cLD: " << cLD << endl
             << endl;

        cLD_writer << to_string(gen + 1) << "\t" << (Pa) << "\t" << (Pb) << "\t" << (Pab) << "\t" << (cLD) << "\n";

        cLD_writer.flush();

        cout << "Progeny moving to next generation: " << progeny_surviving_ID.size();

        cout << "\nConfiguring next generation\n\n";
        function.clear_Array_float_CPU(pop_GeneA_GeneB_Parent, eff_Population);
        pop_GeneA_GeneB_Parent = function.create_Fill_2D_array_FLOAT(progeny_surviving_ID.size(), 3, 0);

        free(parent_IDs);
        parent_IDs = (int *)malloc((progeny_surviving_ID.size()) * sizeof(int));

        mutation_points_All_Parents.clear();

        for (size_t i = 0; i < progeny_surviving_ID.size(); i++)
        {
            pop_GeneA_GeneB_Parent[i][0] = pop_GeneA_GeneB_Progeny[progeny_surviving_ID[i]][0];
            pop_GeneA_GeneB_Parent[i][1] = pop_GeneA_GeneB_Progeny[progeny_surviving_ID[i]][1];
            pop_GeneA_GeneB_Parent[i][2] = pop_GeneA_GeneB_Progeny[progeny_surviving_ID[i]][2];
            parent_IDs[i] = i;
            mutation_points_All_Parents.push_back(mutation_points_All_Progeny[progeny_surviving_ID[i]]);
        }

        mutation_points_All_Progeny.clear();

        function.clear_Array_float_CPU(pop_GeneA_GeneB_Progeny, progeny_surviving_ID.size());
        function.clear_Array_int_CPU(parent_Pairs, eff_Population / 2);

        eff_Population = progeny_surviving_ID.size();

        // for (size_t i = 0; i < sum_Progeny; i++)
        // {
        //     cout << pop_GeneA_GeneB_Progeny[i][0] << "\t";
        //     cout << pop_GeneA_GeneB_Progeny[i][1] << "\t";
        //     cout << pop_GeneA_GeneB_Progeny[i][2] << "\n";
        // }

        // REMOVE

        // clear 2D arrays
        // if (gen == 1)
        // {
        //     break;
        // }
    }
    cLD_writer.close();
}

void network::ncbi_find_conserved()
{
    functions_library function = functions_library();

    cout << "Indetifying conserved regions\n\n";

    string parent_Folder = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/gene_Sequences/";
    string summary_File = parent_Folder + "/" + "summary_excel.csv";

    fstream summary_File_read;
    summary_File_read.open(summary_File, ios::in);

    if (summary_File_read.is_open())
    {
        string line;
        getline(summary_File_read, line);

        vector<string> line_Data;

        while (getline(summary_File_read, line))
        {
            function.split(line_Data, line, ',');
            int sequences_Present = stoi(line_Data[1]);
            if (sequences_Present >= 15)
            {
                string gene_Name = parent_Folder + "/alignments/" + line_Data[0] + ".fasta.align";
                cout << "Processing: " << gene_Name << endl;
                cout << "Sequences present: " << sequences_Present << endl;

                fstream align_File;
                align_File.open(gene_Name, ios::in);

                if (align_File.is_open())
                {
                    string align_Line;
                    vector<string> line_Data_Align;
                    getline(align_File, align_Line);

                    vector<string> sequence_Reconstruction;

                    for (int i = 0; i < (sequences_Present + 1); i++)
                    {
                        sequence_Reconstruction.push_back("");
                    }

                    int count = 0;
                    int char_Length = 0;
                    while (getline(align_File, align_Line))
                    {
                        if (align_Line.length() > 0)
                        {
                            // cout << align_Line << endl;
                            if (count < sequences_Present)
                            {
                                function.split(line_Data_Align, align_Line, ' ');
                                // cout << line_Data_Align[1] << endl;
                                sequence_Reconstruction[count] = sequence_Reconstruction[count].append(line_Data_Align[1]);
                                if (char_Length == 0)
                                {
                                    char_Length = line_Data_Align[0].size() + 1;
                                }
                                count++;
                            }
                            else
                            {
                                // cout << align_Line.substr(char_Length, align_Line.size()) << endl;
                                sequence_Reconstruction[count] = sequence_Reconstruction[count].append(align_Line.substr(char_Length, align_Line.size()));
                                count = 0;
                            }
                        }
                    }
                    align_File.close();

                    for (int i = 0; i < (sequences_Present + 1); i++)
                    {
                        cout << sequence_Reconstruction[i] << endl;
                    }
                }

                cout << endl;

                // REMOVE later
                exit(-1);
            }
        }
        summary_File_read.close();
    }
}

void network::ncbi_Read()
{
    functions_library function = functions_library();

    cout << "NCBI files\n\n";

    string primary_Location = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Dr_Koul/New test";

    string output_Folder = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/gene_Sequences/";
    function.config_Folder("/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/gene_Sequences/", "Sequence folder");

    vector<string> folders;

    cout << "Looking for folders in: " << primary_Location << endl
         << endl;

    for (const auto &entry : filesystem::directory_iterator(primary_Location))
    {
        if (filesystem::is_directory(entry.status()))
        {
            string folder_Name = entry.path().filename().string();
            folders.push_back(folder_Name);
        }
    }

    cout << folders.size() << " folders found\n\n";

    vector<pair<string, int>> gene_Count;
    vector<string> invalid_Genes;

    vector<vector<pair<string, int>>> gene_Count_Folder;

    for (int folder_ID = 0; folder_ID < folders.size(); folder_ID++)
    {
        cout << "Reading folder: " << folders[folder_ID] << endl
             << endl;

        string nest_Location = primary_Location + "/" + folders[folder_ID] + "/ncbi_dataset/data";
        vector<string> sequence_GFF_folders;

        for (const auto &entry : filesystem::directory_iterator(nest_Location))
        {
            if (filesystem::is_directory(entry.status()))
            {
                string folder_Name = entry.path().filename().string();
                sequence_GFF_folders.push_back(folder_Name);
                // cout << folder_Name << endl;
            }
        }

        cout << sequence_GFF_folders.size() << " sequence folders found\n\n";

        vector<pair<string, int>> gene_Count_Folder_ONLY;

        for (int seq_Folder_ID = 0; seq_Folder_ID < sequence_GFF_folders.size(); seq_Folder_ID++)
        {
            cout << "Processing: " << sequence_GFF_folders[seq_Folder_ID] << endl
                 << endl;

            string sequence_Folder_location = nest_Location + "/" + sequence_GFF_folders[seq_Folder_ID];

            string fna_File = "";
            string gff_File = "";

            for (const auto &entry : filesystem::directory_iterator(sequence_Folder_location))
            {
                if (filesystem::is_regular_file(entry.status()))
                {
                    string file_Query = entry.path().filename().string();
                    // cout << file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) << endl;
                    if (file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".fna" || file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".fasta")
                    {
                        fna_File = file_Query;
                    }
                    else if (file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".gff")
                    {
                        gff_File = file_Query;
                    }
                }
            }
            cout << "GFF file\t: " << gff_File << "\n"
                 << "Sequence file: " << fna_File << endl
                 << endl;

            if (gff_File != "" && fna_File != "")
            {

                fstream gff_Read;
                gff_Read.open(sequence_Folder_location + "/" + gff_File, ios::in);

                if (gff_Read.is_open())
                {
                    cout << "Processing GFF file\n";

                    string line;

                    while (getline(gff_Read, line))
                    {
                        if (line.at(0) != '#')
                        {
                            break;
                        }
                    }

                    do
                    {
                        // cout << line << endl;
                        vector<string> line_Split;
                        function.split(line_Split, line, '\t');

                        if (line_Split[1] == "RefSeq")
                        {
                            if (line_Split[2] == "gene")
                            {
                                string seq_ID = line_Split[0];
                                int start = stoi(line_Split[3]) - 1;
                                int stop = stoi(line_Split[4]);

                                string gene_Name;

                                // cout << seq_ID << "\t";
                                // cout << start << "\t";
                                // cout << stop << endl;

                                vector<string> description_Split;
                                function.split(description_Split, line_Split[line_Split.size() - 1], ';');
                                for (string sub_Split : description_Split)
                                {
                                    vector<string> sub_Split_vector;
                                    function.split(sub_Split_vector, sub_Split, '=');
                                    if (sub_Split_vector[0] == "Name")
                                    {
                                        gene_Name = sub_Split_vector[1];
                                        break;
                                    }
                                }

                                // cout << "Check\n";
                                int already_Present = 0;
                                for (int check_Name = 0; check_Name < gene_Count_Folder_ONLY.size(); check_Name++)
                                {
                                    if (gene_Count_Folder_ONLY[check_Name].first == gene_Name)
                                    {
                                        gene_Count_Folder_ONLY[check_Name].second = gene_Count_Folder_ONLY[check_Name].second + 1;
                                        already_Present = 1;
                                        break;
                                    }
                                }

                                if (already_Present == 0)
                                {
                                    gene_Count_Folder_ONLY.push_back(make_pair(gene_Name, 1));
                                }

                                // cout << "Check\n";
                                // cout << "Processing gene: " << gene_Name << endl;
                            }
                        }

                        getline(gff_Read, line);
                    } while (line.at(0) != '#');

                    gff_Read.close();
                }

                // exit(-1);

                // LOAD Sequence file to memory
                // vector<pair<string, string>> ID_sequences;

                // fstream sequence_Read;
                // sequence_Read.open(sequence_Folder_location + "/" + fna_File, ios::in);

                // if (sequence_Read.is_open())
                // {
                //     cout << "Processing sequence file\n";
                //     string line;
                //     int sequence_Catch = 0;

                //     string sequence_Full = "";
                //     string sequence_Name = "";

                //     while (getline(sequence_Read, line))
                //     {
                //         if (line.at(0) == '>')
                //         {
                //             // cout << line << endl;

                //             if (sequence_Catch == 1)
                //             {
                //                 ID_sequences.push_back(make_pair(sequence_Name, sequence_Full));
                //                 sequence_Name = line.substr(1, line.length());
                //                 // sequence_Catch = 0;
                //                 sequence_Full = "";
                //             }
                //             else
                //             {
                //                 sequence_Catch = 1;
                //                 sequence_Name = line.substr(1, line.length());
                //             }
                //         }
                //         else
                //         {
                //             sequence_Full.append(line);
                //         }
                //     }

                //     ID_sequences.push_back(make_pair(sequence_Name, sequence_Full));
                //     sequence_Read.close();
                // }

                // cout << ID_sequences.size() << " sequence(s) found\n\n";

                // // cout << ID_sequences[0].first << "\n"
                // //      << ID_sequences[0].second.size() << endl;

                // // GET genes from GFF
                // fstream gff_Read;
                // gff_Read.open(sequence_Folder_location + "/" + gff_File, ios::in);

                // if (gff_Read.is_open())
                // {
                //     cout << "Processing GFF file\n";

                //     string line;

                //     while (getline(gff_Read, line))
                //     {
                //         if (line.at(0) != '#')
                //         {
                //             break;
                //         }
                //     }

                //     do
                //     {
                //         // cout << line << endl;
                //         vector<string> line_Split;
                //         function.split(line_Split, line, '\t');

                //         if (line_Split[1] == "RefSeq")
                //         {
                //             if (line_Split[2] == "gene")
                //             {
                //                 string seq_ID = line_Split[0];
                //                 int start = stoi(line_Split[3]) - 1;
                //                 int stop = stoi(line_Split[4]);

                //                 string gene_Name;

                //                 // cout << seq_ID << "\t";
                //                 // cout << start << "\t";
                //                 // cout << stop << endl;

                //                 vector<string> description_Split;
                //                 function.split(description_Split, line_Split[line_Split.size() - 1], ';');
                //                 for (string sub_Split : description_Split)
                //                 {
                //                     vector<string> sub_Split_vector;
                //                     function.split(sub_Split_vector, sub_Split, '=');
                //                     if (sub_Split_vector[0] == "Name")
                //                     {
                //                         gene_Name = sub_Split_vector[1];
                //                         break;
                //                     }
                //                 }

                //                 // cout << "Processing gene: " << gene_Name << endl;

                //                 for (int f_Seq = 0; f_Seq < ID_sequences.size(); f_Seq++)
                //                 {
                //                     // cout << ID_sequences[f_Seq].first << endl;
                //                     // cout << "Sequence: " << seq_ID << endl;
                //                     if (ID_sequences[f_Seq].first.find(seq_ID) != string::npos)
                //                     {
                //                         cout << ID_sequences[f_Seq].first << endl;
                //                         cout << "Found sequence: " << seq_ID << endl;
                //                         string sequence = ID_sequences[f_Seq].second.substr(start, (stop - start));

                //                         // Write seq to file
                //                         string file_Write_Location = output_Folder + gene_Name + ".fasta";

                //                         if (!filesystem::exists(file_Write_Location))
                //                         {
                //                             function.create_File(file_Write_Location);
                //                             gene_Count.push_back(make_pair(gene_Name, 1));
                //                         }
                //                         else
                //                         {
                //                             for (int count_Genes = 0; count_Genes < gene_Count.size(); count_Genes++)
                //                             {
                //                                 if (gene_Count[count_Genes].first == gene_Name)
                //                                 {
                //                                     gene_Count[count_Genes].second = gene_Count[count_Genes].second + 1;
                //                                     break;
                //                                 }
                //                             }
                //                         }

                //                         fstream gene_Write;
                //                         gene_Write.open(file_Write_Location, ios::app);

                //                         gene_Write << ">" << folders[folder_ID] << " " << seq_ID << " " << gene_Name << " " << to_string(start + 1) << "_" << to_string(stop) << "\n";
                //                         gene_Write << sequence << "\n";

                //                         gene_Write.close();
                //                         break;
                //                     }
                //                 }
                //             }
                //         }

                //         getline(gff_Read, line);
                //     } while (line.at(0) != '#');

                //     gff_Read.close();
                // }
            }

            // REMOVE
            // break;
        }

        // cout << "Check\n";
        // for (int add_Check = 0; add_Check < gene_Count_Folder.size(); add_Check++)
        // {
        //     if (gene_Count_Folder[add_Check].second == 1)
        //     {
        //         int check_Caught = 0;
        //         for (int invalid_Check = 0; invalid_Check < invalid_Genes.size(); invalid_Check++)
        //         {
        //             if (invalid_Genes[invalid_Check] == gene_Count_Folder[add_Check].first)
        //             {
        //                 check_Caught = 1;
        //                 break;
        //             }
        //         }
        //         if (check_Caught == 0)
        //         {
        //             int check_Present = 0;
        //             for (int check_Present_int = 0; check_Present_int < gene_Count.size(); check_Present_int++)
        //             {
        //                 if (gene_Count[check_Present_int].first == gene_Count_Folder[add_Check].first)
        //                 {
        //                     gene_Count[check_Present_int].second = gene_Count[check_Present_int].second + 1;
        //                     check_Present = 1;
        //                     break;
        //                 }
        //             }
        //             if (check_Present == 0)
        //             {
        //                 gene_Count.push_back(make_pair(gene_Count_Folder[add_Check].first, gene_Count_Folder[add_Check].second));
        //             }
        //         }
        //     }
        //     else
        //     {
        //         invalid_Genes.push_back(gene_Count_Folder[add_Check].first);
        //     }
        // }
        // REMOVE
        // break;

        gene_Count_Folder.push_back(gene_Count_Folder_ONLY);
    }

    for (size_t i = 0; i < gene_Count_Folder.size(); i++)
    {

        vector<pair<string, int>> gene_Count_Folder_ONLY;
        gene_Count_Folder_ONLY = gene_Count_Folder[i];

        for (int add_Check = 0; add_Check < gene_Count_Folder_ONLY.size(); add_Check++)
        {
            if (gene_Count_Folder_ONLY[add_Check].second != 1)
            {
                int check_Caught = 0;
                for (int invalid_Check = 0; invalid_Check < invalid_Genes.size(); invalid_Check++)
                {
                    if (invalid_Genes[invalid_Check] == gene_Count_Folder_ONLY[add_Check].first)
                    {
                        check_Caught = 1;
                        break;
                    }
                }
                if (check_Caught == 0)
                {
                    invalid_Genes.push_back(gene_Count_Folder_ONLY[add_Check].first);
                }
                // if (check_Caught == 0)
                // {
                //     int check_Present = 0;
                //     for (int check_Present_int = 0; check_Present_int < gene_Count.size(); check_Present_int++)
                //     {
                //         if (gene_Count[check_Present_int].first == gene_Count_Folder[add_Check].first)
                //         {
                //             gene_Count[check_Present_int].second = gene_Count[check_Present_int].second + 1;
                //             check_Present = 1;
                //             break;
                //         }
                //     }
                //     if (check_Present == 0)
                //     {
                //         gene_Count.push_back(make_pair(gene_Count_Folder[add_Check].first, gene_Count_Folder[add_Check].second));
                //     }
                // }
            }
            // else
            // {
            //     invalid_Genes.push_back(gene_Count_Folder[add_Check].first);
            // }
        }
    }

    for (size_t i = 0; i < gene_Count_Folder.size(); i++)
    {

        vector<pair<string, int>> gene_Count_Folder_ONLY;
        gene_Count_Folder_ONLY = gene_Count_Folder[i];

        for (int add_Check = 0; add_Check < gene_Count_Folder_ONLY.size(); add_Check++)
        {
            if (gene_Count_Folder_ONLY[add_Check].second == 1)
            {
                int check_Caught = 0;
                for (int invalid_Check = 0; invalid_Check < invalid_Genes.size(); invalid_Check++)
                {
                    if (invalid_Genes[invalid_Check] == gene_Count_Folder_ONLY[add_Check].first)
                    {
                        check_Caught = 1;
                        break;
                    }
                }
                if (check_Caught == 0)
                {
                    int check_Present = 0;
                    for (int check_Present_int = 0; check_Present_int < gene_Count.size(); check_Present_int++)
                    {
                        if (gene_Count[check_Present_int].first == gene_Count_Folder_ONLY[add_Check].first)
                        {
                            gene_Count[check_Present_int].second = gene_Count[check_Present_int].second + 1;
                            check_Present = 1;
                            break;
                        }
                    }
                    if (check_Present == 0)
                    {
                        gene_Count.push_back(make_pair(gene_Count_Folder_ONLY[add_Check].first, gene_Count_Folder_ONLY[add_Check].second));
                    }
                }
            }
            // else
            // {
            //     invalid_Genes.push_back(gene_Count_Folder[add_Check].first);
            // }
        }
    }
    // cout << "Done\n\n";

    // for (int print = 0; print < gene_Count.size(); print++)
    // {
    //     if (gene_Count[print].second == 5)
    //     {
    //         cout << gene_Count[print].first << "\t" << gene_Count[print].second << "\n";
    //     }
    // }

    cout << "Getting genes out\n\n";

    for (int folder_ID = 0; folder_ID < folders.size(); folder_ID++)
    {
        cout << "Reading folder: " << folders[folder_ID] << endl
             << endl;

        string nest_Location = primary_Location + "/" + folders[folder_ID] + "/ncbi_dataset/data";
        vector<string> sequence_GFF_folders;

        for (const auto &entry : filesystem::directory_iterator(nest_Location))
        {
            if (filesystem::is_directory(entry.status()))
            {
                string folder_Name = entry.path().filename().string();
                sequence_GFF_folders.push_back(folder_Name);
                // cout << folder_Name << endl;
            }
        }

        cout << sequence_GFF_folders.size() << " sequence folders found\n\n";

        vector<pair<string, int>> gene_Count_Folder;

        for (int seq_Folder_ID = 0; seq_Folder_ID < sequence_GFF_folders.size(); seq_Folder_ID++)
        {
            cout << "Processing: " << sequence_GFF_folders[seq_Folder_ID] << endl
                 << endl;

            string sequence_Folder_location = nest_Location + "/" + sequence_GFF_folders[seq_Folder_ID];

            string fna_File = "";
            string gff_File = "";

            for (const auto &entry : filesystem::directory_iterator(sequence_Folder_location))
            {
                if (filesystem::is_regular_file(entry.status()))
                {
                    string file_Query = entry.path().filename().string();
                    // cout << file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) << endl;
                    if (file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".fna" || file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".fasta")
                    {
                        fna_File = file_Query;
                    }
                    else if (file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".gff")
                    {
                        gff_File = file_Query;
                    }
                }
            }

            cout << "GFF file\t: " << gff_File << "\n"
                 << "Sequence file: " << fna_File << endl
                 << endl;

            if (gff_File != "" && fna_File != "")
            {
                vector<pair<string, string>> ID_sequences;

                fstream sequence_Read;
                sequence_Read.open(sequence_Folder_location + "/" + fna_File, ios::in);

                if (sequence_Read.is_open())
                {
                    cout << "Processing sequence file\n";
                    string line;
                    int sequence_Catch = 0;

                    string sequence_Full = "";
                    string sequence_Name = "";

                    while (getline(sequence_Read, line))
                    {
                        if (line.at(0) == '>')
                        {
                            // cout << line << endl;

                            if (sequence_Catch == 1)
                            {
                                ID_sequences.push_back(make_pair(sequence_Name, sequence_Full));
                                sequence_Name = line.substr(1, line.length());
                                // sequence_Catch = 0;
                                sequence_Full = "";
                            }
                            else
                            {
                                sequence_Catch = 1;
                                sequence_Name = line.substr(1, line.length());
                            }
                        }
                        else
                        {
                            sequence_Full.append(line);
                        }
                    }

                    ID_sequences.push_back(make_pair(sequence_Name, sequence_Full));
                    sequence_Read.close();
                }

                cout << ID_sequences.size() << " sequence(s) found\n\n";

                fstream gff_Read;
                gff_Read.open(sequence_Folder_location + "/" + gff_File, ios::in);

                if (gff_Read.is_open())
                {
                    cout << "Processing GFF file\n";

                    string line;

                    while (getline(gff_Read, line))
                    {
                        if (line.at(0) != '#')
                        {
                            break;
                        }
                    }
                    do
                    {
                        // cout << line << endl;
                        vector<string> line_Split;
                        function.split(line_Split, line, '\t');

                        if (line_Split[1] == "RefSeq")
                        {
                            if (line_Split[2] == "gene")
                            {
                                string seq_ID = line_Split[0];
                                int start = stoi(line_Split[3]) - 1;
                                int stop = stoi(line_Split[4]);

                                string gene_Name;

                                // cout << seq_ID << "\t";
                                // cout << start << "\t";
                                // cout << stop << endl;

                                vector<string> description_Split;
                                function.split(description_Split, line_Split[line_Split.size() - 1], ';');
                                for (string sub_Split : description_Split)
                                {
                                    vector<string> sub_Split_vector;
                                    function.split(sub_Split_vector, sub_Split, '=');
                                    if (sub_Split_vector[0] == "Name")
                                    {
                                        gene_Name = sub_Split_vector[1];
                                        break;
                                    }
                                }

                                int check_Exists = 0;
                                for (int check = 0; check < gene_Count.size(); check++)
                                {
                                    if (gene_Count[check].first == gene_Name)
                                    {
                                        check_Exists = 1;
                                        break;
                                    }
                                }

                                if (check_Exists == 1)
                                {
                                    for (int f_Seq = 0; f_Seq < ID_sequences.size(); f_Seq++)
                                    {
                                        if (ID_sequences[f_Seq].first.find(seq_ID) != string::npos)
                                        {
                                            cout << ID_sequences[f_Seq].first << endl;
                                            cout << "Found sequence: " << seq_ID << endl;
                                            string sequence = ID_sequences[f_Seq].second.substr(start, (stop - start));

                                            // Write seq to file
                                            string file_Write_Location = output_Folder + gene_Name + ".fasta";

                                            if (!filesystem::exists(file_Write_Location))
                                            {
                                                function.create_File(file_Write_Location);
                                                // gene_Count.push_back(make_pair(gene_Name, 1));
                                            }
                                            // else
                                            // {
                                            //     for (int count_Genes = 0; count_Genes < gene_Count.size(); count_Genes++)
                                            //     {
                                            //         if (gene_Count[count_Genes].first == gene_Name)
                                            //         {
                                            //             gene_Count[count_Genes].second = gene_Count[count_Genes].second + 1;
                                            //             break;
                                            //         }
                                            //     }
                                            // }
                                            fstream gene_Write;
                                            gene_Write.open(file_Write_Location, ios::app);

                                            gene_Write << ">" << folders[folder_ID] << "_" << seq_ID << "_" << gene_Name << "_" << to_string(start + 1) << "-" << to_string(stop) << "\n";
                                            gene_Write << sequence << "\n";

                                            gene_Write.close();
                                            break;
                                        }
                                    }
                                }

                                // cout << "Check\n";

                                // cout << "Check\n";
                                // cout << "Processing gene: " << gene_Name << endl;
                            }
                        }

                        getline(gff_Read, line);
                    } while (line.at(0) != '#');

                    gff_Read.close();
                }
            }
        }
    }

    cout << "\n\nWriing summary\n\n";

    string final_Summary = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/gene_Sequences/summary.csv";

    function.config_File_delete_create(final_Summary, "Gene_name\tCount");

    fstream summary_Write;
    summary_Write.open(final_Summary, ios::app);

    for (int i = 0; i < gene_Count.size(); i++)
    {
        summary_Write << gene_Count[i].first << "\t" << to_string(gene_Count[i].second);

        if (gene_Count[i].second < folders.size())
        {
            string gene_File = output_Folder + gene_Count[i].first + ".fasta";

            fstream gene_Read;
            gene_Read.open(gene_File, ios::in);

            if (gene_Read.is_open())
            {
                vector<string> ID_names;
                vector<string> not_Found;

                cout << "Processing sequence file: " << gene_File << "\n";

                string line;

                while (getline(gene_Read, line))
                {
                    if (line.at(0) == '>')
                    {
                        vector<string> line_Data;
                        function.split(line_Data, line, '_');
                        ID_names.push_back(line_Data[0].substr(1, line_Data[0].length()));
                        // cout << line_Data[0].substr(1, line_Data[0].length()) << endl;
                    }
                }

                for (int folder = 0; folder < folders.size(); folder++)
                {
                    int found = -1;
                    for (int IDs = 0; IDs < ID_names.size(); IDs++)
                    {
                        if (ID_names[IDs] == folders[folder])
                        {
                            found = IDs;
                            break;
                        }
                    }
                    if (found == -1)
                    {
                        summary_Write << "\t" << folders[folder];
                    }
                }

                gene_Read.close();
            }
        }

        // for (int folder_ID = 0; folder_ID < folders.size(); folder_ID++)
        // {
        //     int check_Exists = 0;
        //     cout << "Reading folder: " << folders[folder_ID] << endl
        //          << endl;

        //     string nest_Location = primary_Location + "/" + folders[folder_ID] + "/ncbi_dataset/data";
        //     vector<string> sequence_GFF_folders;

        //     for (const auto &entry : filesystem::directory_iterator(nest_Location))
        //     {
        //         if (filesystem::is_directory(entry.status()))
        //         {
        //             string folder_Name = entry.path().filename().string();
        //             sequence_GFF_folders.push_back(folder_Name);
        //             // cout << folder_Name << endl;
        //         }
        //     }

        //     cout << sequence_GFF_folders.size() << " sequence folders found\n\n";

        //     vector<pair<string, int>> gene_Count_Folder;

        //     for (int seq_Folder_ID = 0; seq_Folder_ID < sequence_GFF_folders.size(); seq_Folder_ID++)
        //     {
        //         cout << "Processing: " << sequence_GFF_folders[seq_Folder_ID] << endl
        //              << endl;

        //         string sequence_Folder_location = nest_Location + "/" + sequence_GFF_folders[seq_Folder_ID];

        //         string fna_File = "";
        //         string gff_File = "";

        //         for (const auto &entry : filesystem::directory_iterator(sequence_Folder_location))
        //         {
        //             if (filesystem::is_regular_file(entry.status()))
        //             {
        //                 string file_Query = entry.path().filename().string();
        //                 // cout << file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) << endl;
        //                 if (file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".fna" || file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".fasta")
        //                 {
        //                     fna_File = file_Query;
        //                 }
        //                 else if (file_Query.substr(file_Query.find_last_of('.'), file_Query.length()) == ".gff")
        //                 {
        //                     gff_File = file_Query;
        //                 }
        //             }
        //         }

        //         cout << "GFF file\t: " << gff_File << "\n"
        //              << "Sequence file: " << fna_File << endl
        //              << endl;

        //         if (gff_File != "" && fna_File != "")
        //         {
        //             fstream gff_Read;
        //             gff_Read.open(sequence_Folder_location + "/" + gff_File, ios::in);

        //             if (gff_Read.is_open())
        //             {
        //                 cout << "Processing GFF file\n";

        //                 string line;

        //                 while (getline(gff_Read, line))
        //                 {
        //                     if (line.at(0) != '#')
        //                     {
        //                         break;
        //                     }
        //                 }
        //                 do
        //                 {
        //                     // cout << line << endl;
        //                     vector<string> line_Split;
        //                     function.split(line_Split, line, '\t');

        //                     if (line_Split[1] == "RefSeq")
        //                     {
        //                         if (line_Split[2] == "gene")
        //                         {
        //                             string seq_ID = line_Split[0];
        //                             int start = stoi(line_Split[3]) - 1;
        //                             int stop = stoi(line_Split[4]);

        //                             string gene_Name;

        //                             // cout << seq_ID << "\t";
        //                             // cout << start << "\t";
        //                             // cout << stop << endl;

        //                             vector<string> description_Split;
        //                             function.split(description_Split, line_Split[line_Split.size() - 1], ';');
        //                             for (string sub_Split : description_Split)
        //                             {
        //                                 vector<string> sub_Split_vector;
        //                                 function.split(sub_Split_vector, sub_Split, '=');
        //                                 if (sub_Split_vector[0] == "Name")
        //                                 {
        //                                     gene_Name = sub_Split_vector[1];
        //                                     break;
        //                                 }
        //                             }

        //                             if (gene_Count[i].first == gene_Name)
        //                             {
        //                                 check_Exists = 1;
        //                                 break;
        //                             }

        //                             // cout << "Check\n";

        //                             // cout << "Check\n";
        //                             // cout << "Processing gene: " << gene_Name << endl;
        //                         }
        //                     }

        //                     getline(gff_Read, line);
        //                 } while (line.at(0) != '#');

        //                 gff_Read.close();

        //                 if (check_Exists == 0)
        //                 {
        //                     summary_Write << "\t" << folders[folder_ID];
        //                 }
        //             }
        //         }
        //     }
        // }
        summary_Write << "\n";
    }

    summary_Write.close();

    cout << "Done\n\n";
}