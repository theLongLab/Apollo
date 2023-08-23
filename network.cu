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
    per_cave_Stride[0] = 0;

    cout << "Dynamic Caveman models\n\nNodes: " << cave_number << endl
         << endl;

    random_device rd;
    mt19937 gen(rd());

    float shape_Nodes_per_cave = 10;
    float scale_Nodes_per_cave = 1;

    gamma_distribution<float> gamma(shape_Nodes_per_cave, scale_Nodes_per_cave);

    for (size_t i = 0; i < cave_number; i++)
    {
        int number_of_Nodes = (int)round(gamma(gen));
        cout << number_of_Nodes << endl;
        per_cave_Stride[i + 1] = per_cave_Stride[i] + number_of_Nodes;
    }

    // float percent_Outside = 0.2;
    int **network_Array = function.create_INT_2D_arrays(cave_number, 3);

    vector<vector<int>> global_Nodes_per_Cave;
}

void network::sim_cLD()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    // random_device rd;
    // mt19937 gen(rd());

    default_random_engine generator;

    cout << "cLD simulator\n"
         << endl;

    int generations = 10;

    int eff_Population = 1000;

    // PROGENY BINOMIAL
    int n = 10;
    float prob = 0.70;
    binomial_distribution<int> binomialDist(n, prob);

    // cout << x << endl;

    float mutation_Rate = 0.8;
    poisson_distribution<int> dist_Poisson(mutation_Rate);

    // int mutation_points = 10;
    float recombination_Prob = 0.005;
    int interactions = 1;

    int *parent_IDs = (int *)malloc((eff_Population) * sizeof(int));

    float **pop_GeneA_GeneB_Parent = function.create_Fill_2D_array_FLOAT(eff_Population, 3, 0);

    cout << "Configure ideal population\n";
    for (size_t i = 0; i < eff_Population / 2; i++)
    {
        pop_GeneA_GeneB_Parent[i][0] = 0;
        pop_GeneA_GeneB_Parent[i][1] = 0;
    }

    for (size_t i = eff_Population / 2; i < eff_Population; i++)
    {
        pop_GeneA_GeneB_Parent[i][0] = 1;
        pop_GeneA_GeneB_Parent[i][1] = 1;
    }

    float Pa = 0;
    float Pb = 0;
    float Pab = 0;

    for (size_t i = 0; i < eff_Population; i++)
    {
        pop_GeneA_GeneB_Parent[i][2] = 0.70;
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

    string cLD_write = "/mnt/d/Deshan/Books/University of Calgary/Experiments/Simulator_Linux/results_of_Simulation/cLD.csv";
    function.config_File_delete_create(cLD_write, "Generation\tPa\tPb\tPab\tcLD");

    fstream cLD_writer;
    cLD_writer.open(cLD_write, ios::app);

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

        int progeny_Fill_count = 0;

        vector<int> progeny_surviving_ID;

        Pa = 0;
        Pb = 0;
        Pab = 0;

        for (int parent_pair = 0; parent_pair < eff_Population / 2; parent_pair++)
        {
            // cout << "\nProcessing parent pair: " << parent_pair << " of " << eff_Population / 2 << endl;

            for (int progeny = 0; progeny < parent_Pairs[parent_pair][2]; progeny++)
            {
                // cout << "\nProcessing progeny: " << progeny << " of " << parent_Pairs[parent_pair][2] / 2;
                uniform_real_distribution<> dis(0.0, 1.0);

                int mom_0_dad_1 = 0;
                if (dis(generator) < 0.5)
                {
                    // cout << dis(generator) << endl;
                    mom_0_dad_1 = 1;
                }

                pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][0];
                pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][1];

                // RECOMBINATION A

                if (dis(generator) < recombination_Prob)
                {
                    mom_0_dad_1 = 0;
                    if (dis(generator) < 0.5)
                    {
                        // cout << "dad recombinant A" << endl;
                        mom_0_dad_1 = 1;
                    }
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][0];
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
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Parent[parent_Pairs[parent_pair][mom_0_dad_1]][1];
                }

                pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = 0.70;

                if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] == pop_GeneA_GeneB_Progeny[progeny_Fill_count][1])
                {
                    if (pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] + pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] != 0)
                    {
                        if (interactions == 1)
                        {
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 0.10;
                        }
                        else
                        {
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                        }
                    }
                }
                else
                {
                    pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                }

                // mutate
                // reduce those with mismatc mutations survivability

                int mutations_A = dist_Poisson(generator);
                int mutations_B = dist_Poisson(generator);

                pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][0] + mutations_A;
                pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][1] + mutations_B;
                // int check_A_B = 0;

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

                if (interactions == 1)
                {
                    if (mutations_A == mutations_B)
                    {
                        if ((mutations_A + mutations_B) != 0)
                        {
                            pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] + 0.10;
                        }
                    }
                    else
                    {
                        pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                    }
                }
                else
                {
                    if ((mutations_A + mutations_B) != 0)
                    {
                        pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] = pop_GeneA_GeneB_Progeny[progeny_Fill_count][2] - 0.15;
                    }
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

        cLD_writer << to_string(gen + 1) << "\t" << to_string(Pa) << "\t" << to_string(Pb) << "\t" << to_string(Pab) << "\t" << to_string(cLD) << "\n";

        cout << "Progeny moving to next generation: " << progeny_surviving_ID.size();

        cout << "\nConfiguring next generation\n\n";
        function.clear_Array_float_CPU(pop_GeneA_GeneB_Parent, eff_Population);
        pop_GeneA_GeneB_Parent = function.create_Fill_2D_array_FLOAT(progeny_surviving_ID.size(), 3, 0);

        free(parent_IDs);
        parent_IDs = (int *)malloc((progeny_surviving_ID.size()) * sizeof(int));

        for (size_t i = 0; i < progeny_surviving_ID.size(); i++)
        {
            pop_GeneA_GeneB_Parent[i][0] = pop_GeneA_GeneB_Progeny[progeny_surviving_ID[i]][0];
            pop_GeneA_GeneB_Parent[i][1] = pop_GeneA_GeneB_Progeny[progeny_surviving_ID[i]][1];
            pop_GeneA_GeneB_Parent[i][2] = pop_GeneA_GeneB_Progeny[progeny_surviving_ID[i]][2];
            parent_IDs[i] = i;
        }

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