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
                    functions.config_Folder(host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), to_string(cave_ID) + "_" + to_string(host_ID) + " Tissue " + to_string(entry_array[tissue]) + " Generation 0");

                    vector<string> sequence_Write_Store_All;
                    int last_seq_Num = 0;
                    functions.sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[entry_array[tissue]],
                                                          max_sequences_per_File, host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num);
                    functions.partial_Write_Check(sequence_Write_Store_All,
                                                  host_Folder + "/" + to_string(entry_array[tissue]) + "/generation_" + to_string(current_Generation), last_seq_Num);

                    current_Viral_load_per_Tissue[tissue] = tissue_Sequences[entry_array[tissue]].size();
                }
            }
        }
        else
        {
            cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << files[file] << endl;
            exit(-1);
        }
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