#include "node_within_host.cuh"

node_within_host::node_within_host()
{
    cout << "Intializing host: ";
}

void node_within_host::setHost(int host_Index, int cave_ID, int host_ID, int profile_ID)
{
    this->host_Index = host_Index;
    this->cave_ID = cave_ID;
    this->host_ID = host_ID;
    this->profile_ID = profile_ID;
    cout << this->cave_ID << "_" << host_ID << endl;
}

void node_within_host::setNum_Generation(int num_Generation)
{
    this->num_Generation = num_Generation;
}

void node_within_host::setInfectious_Load(int infectious_Load)
{
    this->infectious_Load = infectious_Load;
}

void node_within_host::setTerminal_Load(int terminal_Load)
{
    this->terminal_Load = terminal_Load;
}

void node_within_host::setSampling_Effect(float sampling_Effect)
{
    this->sampling_Effect = sampling_Effect;
}

void node_within_host::setCell_Limit(vector<int> cell_Limit_vec)
{
    num_Tissues = cell_Limit_vec.size();
    this->cell_Limit = (int *)malloc(sizeof(int) * num_Tissues);

    for (size_t i = 0; i < num_Tissues; i++)
    {
        cell_Limit[i] = cell_Limit_vec[i];
    }
}

void node_within_host::print_All()
{
    cout << host_Index << "\t"
         << cave_ID << "_" << host_ID << "\t"
         << profile_ID << "\t"
         << num_Generation << "\t"
         << infectious_Load << "\t"
         << terminal_Load << "\t"
         << sampling_Effect;

    for (size_t i = 0; i < num_Tissues; i++)
    {
        cout << "\t" << cell_Limit[i];
    }
    cout << endl;
}