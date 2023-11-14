#include "node_within_host.cuh"

node_within_host::node_within_host()
{
    cout << "Intializing host: ";
}

void node_within_host::setHost(int host_Index, int cave_ID, int host_ID, int profile_ID, int num_Tissues)
{
    this->host_Index = host_Index;
    this->cave_ID = cave_ID;
    this->host_ID = host_ID;
    this->profile_ID = profile_ID;
    this->num_Tissues = num_Tissues;
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

void node_within_host::begin_Infection(functions_library &functions, string &intermediary_Sequence_location,
                                       int entry_tissues, int *entry_array, int &max_sequences_per_File)
{
    // FIRST NODE OF INFECTION IN THE HOST
    string host_Folder = intermediary_Sequence_location + "/" + to_string(host_Index);
    functions.config_Folder(host_Folder, to_string(cave_ID) + "_" + to_string(host_ID));

    vector<vector<string>> tissue_Sequences;
    intialize_Tissues(host_Folder, tissue_Sequences, functions);

    string reference_Sequences = intermediary_Sequence_location + "/reference_Sequences";

    vector<string> files;

    for (const auto &entry : filesystem::directory_iterator(reference_Sequences))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".nfasta")
        {
            files.push_back(entry.path().string());
        }
    }

    vector<string> Sequences;
    vector<string> Sequence_IDs;

    cout << endl;
    for (int file = 0; file < files.size(); file++)
    {
        cout << "Reading file: " << files[file] << endl;
        fstream nfasta;
        nfasta.open(files[file], ios::in);

        if (nfasta.is_open())
        {
            string line;
            string sequence = "";

            while (getline(nfasta, line))
            {
                if (line.at(0) != '>')
                {
                    sequence.append(line);
                }
                else
                {
                    Sequence_IDs.push_back(line);
                    if (sequence != "")
                    {
                        Sequences.push_back(sequence);
                        sequence = "";
                    }
                }
            }

            if (sequence != "")
            {
                Sequences.push_back(sequence);
                sequence = "";
            }

            random_device rd; // Will be used to obtain a seed for the random number engine
            mt19937 gen(rd());
            uniform_int_distribution<int> entry_Tissue_select(0, entry_tissues - 1);

            cout << endl;

            for (int sequence = 0; sequence < Sequences.size(); sequence++)
            {
                int tissue_Index = entry_array[entry_Tissue_select(gen)];
                cout << "Sequence " << sequence + 1 << " infects tissue: " << tissue_Index << endl;
                tissue_Sequences[tissue_Index].push_back(Sequences[sequence]);
            }

            for (int tissue = 0; tissue < entry_tissues; tissue++)
            {
                if (tissue_Sequences[entry_array[tissue]].size() > 0)
                {
                    current_Viral_load_per_Tissue[entry_array[tissue]] = tissue_Sequences[entry_array[tissue]].size();
                    functions.config_Folder(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(entry_array[tissue]) + " Generation 0");

                    vector<string> sequence_Write_Store_All;
                    int last_seq_Num = 0;
                    functions.sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[entry_array[tissue]],
                                                          max_sequences_per_File, host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num);
                    functions.partial_Write_Check(sequence_Write_Store_All,
                                                  host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num);
                }
            }
        }
        else
        {
            cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << files[file] << endl;
            exit(-1);
        }
    }
    status = "Infected";
    // for (int tissue = 0; tissue < num_Tissues; tissue++)
    // {
    //     cout << current_Viral_load_per_Tissue[tissue] << endl;
    // }
    // exit(-1);
}

void node_within_host::begin_Infection(functions_library &functions, string &intermediary_Sequence_location, string &source_Target_file_Location,
                                       int &source_Index, int &source_Generation, string &source_Name, int *source_current_Viral_load_per_Tissue,
                                       int num_viruses_to_transfer,
                                       int &entry_tissues, int *entry_array, int &exit_Load, int &exit_tissues, int *exit_array,
                                       int &max_sequences_per_File,
                                       mt19937 &gen)
{
    cout << "Node " << this->cave_ID << "_" << this->host_ID << " is being infected by " << source_Name;

    if (num_viruses_to_transfer > exit_Load)
    {
        num_viruses_to_transfer = exit_Load;
    }

    cout << "Transferring " << num_viruses_to_transfer << " viral particle(s)\n";

    uniform_int_distribution<> distribution_exit_Tissue(0, exit_tissues - 1);

    for (int particle = 0; particle < num_viruses_to_transfer; particle++)
    {
        int exit_tissue_Index = distribution_exit_Tissue(gen);

        uniform_int_distribution<> distribution_particle(0, exit_array[exit_tissue_Index] - 1);

        int viral_Index = distribution_particle(gen);
    }
}

void node_within_host::intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions)
{
    current_Viral_load_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
    for (int tissue = 0; tissue < num_Tissues; tissue++)
    {
        vector<string> tissue_Sequence;
        functions.config_Folder(host_Folder + "/" + to_string(tissue), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(tissue + 1));

        tissue_Sequences.push_back(tissue_Sequence);
        current_Viral_load_per_Tissue[tissue] = 0;
    }
}

int node_within_host::get_Load(int &num_tissues_Calc, int *tissue_array)
{
    int sum = 0;

    for (int tissue = 0; tissue < num_tissues_Calc; tissue++)
    {
        sum = sum + current_Viral_load_per_Tissue[tissue_array[tissue]];
    }

    return sum;
}

int node_within_host::infectious_status(int &num_tissues, int *tissue_array)
{
    if (get_Load(num_tissues, tissue_array) >= infectious_Load)
    {
        status = "Infectious";
        return 1;
    }
    else
    {
        return 0;
    }
}
int node_within_host::terminal_status(int &num_tissues, int *tissue_array)
{
    if (get_Load(num_tissues, tissue_array) >= terminal_Load)
    {
        status = "Dead";
        return 1;
    }
    else
    {
        return 0;
    }
}

string node_within_host::get_Name()
{
    return to_string(cave_ID) + "_" + to_string(host_ID);
}

string node_within_host::get_Status()
{
    return this->status;
}

int node_within_host::get_Profile()
{
    return profile_ID;
}