#include "cancer_Host.cuh"

cancer_Host::cancer_Host()
{
    cout << "\nCancer host intialization\n";
}

void cancer_Host::initialize(functions_library &functions,
                             vector<string> &tissue_Names,
                             string &intermediary_Sequence_location, string &first_Infection,
                             int &current_Generation,
                             string &output_Node_location,
                             int &max_sequences_per_File)
{
    cout << "\nConfiguring cancer host\n";

    num_Tissues = tissue_Names.size();

    string host_Folder = intermediary_Sequence_location + "/cancer_Host";
    functions.config_Folder(host_Folder, "Cancer host sequence data");

    vector<vector<string>> tissue_Sequences;
    intialize_Tissues(host_Folder, tissue_Sequences, functions, current_Generation);

    string reference_Sequences = intermediary_Sequence_location + "/reference_Sequences";

    cout << "\nFirst infection mode: " << first_Infection << endl;

    if (first_Infection == "RANDOM")
    {
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
        vector<char> seq_Status;
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
                nfasta.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN NFATSA FILE: " << files[file] << endl;
                exit(-1);
            }
        }

        random_device rd; // Will be used to obtain a seed for the random number engine
        mt19937 gen(rd());
        uniform_int_distribution<int> entry_Tissue_select(0, tissue_Names.size() - 1);

        cout << endl;

        for (int sequence = 0; sequence < Sequences.size(); sequence++)
        {
            int tissue_Index = entry_Tissue_select(gen);
            cout << "Sequence " << sequence + 1 << " infects tissue: " << tissue_Index << endl;
            tissue_Sequences[tissue_Index].push_back(Sequences[sequence]);
        }

        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            if (tissue_Sequences[tissue].size() > 0)
            {
                current_cell_load_per_Tissue[tissue] = tissue_Sequences[tissue].size();
                functions.config_Folder(host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), "Tissue " + tissue_Names[tissue] + " Generation " + to_string(current_Generation));

                if (!filesystem::exists(output_Node_location + "/cancer_Host"))
                {
                    functions.config_Folder(output_Node_location + "/cancer_Host", "Cancer host node");
                    functions.create_File(output_Node_location + "/cancer_Host/sequence_Profiles.csv", "Sequence_ID\tTissue");
                    functions.create_File(output_Node_location + "/cancer_Host/sequence_parent_Progeny_relationships.csv", "Source\tTarget");
                }

                vector<string> sequence_Write_Store_All;
                int last_seq_Num = 0;
                sequence_Write_Configurator(sequence_Write_Store_All, tissue_Sequences[tissue],
                                            max_sequences_per_File, host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), last_seq_Num, seq_Status,
                                            output_Node_location + "/cancer_Host/sequence_Profiles.csv", tissue_Names[tissue], current_Generation);
                partial_Write_Check(sequence_Write_Store_All,
                                    host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), last_seq_Num, seq_Status,
                                    output_Node_location + "/cancer_Host/sequence_Profiles.csv", tissue_Names[tissue], current_Generation);
            }
        }
    }
    else
    {
        cout << "\nInfecting tissues\n";
        vector<string> line_Data;

        //// Write to sequence profiles file

        for (int tissue = 0; tissue < tissue_Names.size(); tissue++)
        {
            string reference_Sequences_tissue = reference_Sequences + "/" + tissue_Names[tissue];
            int last_seq_Num = 0;

            if (filesystem::exists(reference_Sequences_tissue) && filesystem::is_directory(reference_Sequences_tissue))
            {
                cout << "\nTissue: " << tissue_Names[tissue] << endl;
                for (const auto &entry : filesystem::directory_iterator(reference_Sequences_tissue))
                {
                    if (entry.is_regular_file() && entry.path().extension() == ".nfasta")
                    {
                        string file_Name = entry.path().stem();
                        functions.split(line_Data, file_Name, '_');

                        int num_Particles_Tissue = stoi(line_Data[1]) - stoi(line_Data[0]) + 1;

                        cout << "Sequences migrating: " << num_Particles_Tissue << endl;

                        current_cell_load_per_Tissue[tissue] = current_cell_load_per_Tissue[tissue] + num_Particles_Tissue;

                        functions.config_Folder(host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation), "Tissue " + to_string(tissue) + " Generation " + to_string(current_Generation));

                        if (!filesystem::exists(output_Node_location + "/cancer_Host"))
                        {
                            functions.config_Folder(output_Node_location + "/cancer_Host", "Cancer host node");
                            functions.create_File(output_Node_location + "/cancer_Host" + "/sequence_Profiles.csv", "Sequence_ID\tTissue");
                            functions.create_File(output_Node_location + "/cancer_Host" + "/sequence_parent_Progeny_relationships.csv", "Source\tTarget");
                        }

                        filesystem::copy_file(entry.path().string(), host_Folder + "/" + to_string(tissue) + "/generation_" + to_string(current_Generation) + "/" + file_Name + ".nfasta");

                        fstream sequence_Profile;
                        sequence_Profile.open(output_Node_location + "/cancer_Host" + "/sequence_Profiles.csv", ios::app);

                        for (int sequence_Num = 0; sequence_Num < num_Particles_Tissue; sequence_Num++)
                        {
                            sequence_Profile << tissue_Names[tissue] << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue_Names[tissue] << endl;
                            last_seq_Num++;
                        }

                        sequence_Profile.close();
                    }
                }
            }
        }
    }
}

void cancer_Host::sequence_Write_Configurator(vector<string> &sequence_Write_Store_All, vector<string> &sequence_Write_Store,
                                              int &max_sequences_per_File, const string &folder_Location, int &last_seq_Num,
                                              vector<char> &seq_Status,
                                              string sequence_Profiles_Location, string tissue, int current_Generation)
{
    for (int sequence_Collect = 0; sequence_Collect < sequence_Write_Store.size(); sequence_Collect++)
    {
        sequence_Write_Store_All.push_back(sequence_Write_Store[sequence_Collect]);
    }

    sequence_Write_Store.clear();

    if (sequence_Write_Store_All.size() >= max_sequences_per_File)
    {
        int full_Write_Count = sequence_Write_Store_All.size() / max_sequences_per_File;

        for (int full = 0; full < full_Write_Count; full++)
        {

            string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + max_sequences_per_File - 1) + ".nfasta";
            fstream fasta_File;
            fasta_File.open(fasta_file_Location, ios::out);
            fstream sequence_Profile;
            //"Sequence_ID\tHost\tTissue"
            sequence_Profile.open(sequence_Profiles_Location, ios::app);

            if (fasta_File.is_open())
            {
                for (int write_Seq = (full * max_sequences_per_File); write_Seq < ((full * max_sequences_per_File) + max_sequences_per_File); write_Seq++)
                {
                    fasta_File << ">" << last_seq_Num;
                    if (seq_Status.size() == 0)
                    {
                        fasta_File << "_A";
                    }
                    else
                    {
                        fasta_File << "_" << seq_Status[write_Seq];
                    }
                    fasta_File << endl;
                    fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                    sequence_Profile << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue << endl;
                    last_seq_Num++;
                }

                fasta_File.close();
                sequence_Profile.close();
            }
            else
            {
                cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
                exit(-1);
            }
            // last_seq_Num = last_seq_Num + max_sequences_per_File;
        }

        // int parital_Write_Count = sequence_Write_Store_All.size() % max_sequences_per_File;
        vector<string> sequence_Write_Store_temp;
        vector<char> temp_seq_Status;
        for (int fill = full_Write_Count * max_sequences_per_File; fill < sequence_Write_Store_All.size(); fill++)
        {
            sequence_Write_Store_temp.push_back(sequence_Write_Store_All[fill]);
            if (seq_Status.size() > 0)
            {
                temp_seq_Status.push_back(seq_Status[fill]);
            }
        }

        sequence_Write_Store_All.clear();
        sequence_Write_Store_All = sequence_Write_Store_temp;

        if (seq_Status.size() > 0)
        {
            seq_Status.clear();
            seq_Status = temp_seq_Status;
        }
    }
}

void cancer_Host::partial_Write_Check(vector<string> &sequence_Write_Store_All,
                                            const string &folder_Location, int &last_seq_Num,
                                            vector<char> &seq_Status,
                                            string sequence_Profiles_Location, string tissue, int current_Generation)
{
    if (sequence_Write_Store_All.size() > 0)
    {
        string fasta_file_Location = folder_Location + "/" + to_string(last_seq_Num) + "_" + to_string(last_seq_Num + sequence_Write_Store_All.size() - 1) + ".nfasta";
        fstream fasta_File;
        fasta_File.open(fasta_file_Location, ios::out);
        fstream sequence_Profile;
        sequence_Profile.open(sequence_Profiles_Location, ios::app);
        if (fasta_File.is_open())
        {
            for (int write_Seq = 0; write_Seq < sequence_Write_Store_All.size(); write_Seq++)
            {
                fasta_File << ">" << last_seq_Num;
                if (seq_Status.size() == 0)
                {
                    fasta_File << "_A";
                }
                else
                {
                    fasta_File << "_" << seq_Status[write_Seq];
                }
                fasta_File << endl;
                fasta_File << sequence_Write_Store_All[write_Seq] << endl;
                sequence_Profile << tissue << "_" << current_Generation << "_" << last_seq_Num << "\t" << tissue << endl;
                last_seq_Num++;
            }

            fasta_File.close();
            sequence_Profile.close();
        }
        else
        {
            cout << "ERROR: COULD NOT CREATE NFASTA FILE: " << fasta_file_Location << endl;
            exit(-1);
        }
        sequence_Write_Store_All.clear();
    }
}

void cancer_Host::intialize_Tissues(string &host_Folder, vector<vector<string>> &tissue_Sequences, functions_library &functions, int &current_Generation)
{
    current_cell_load_per_Tissue = (int *)malloc(sizeof(int) * num_Tissues);
    dead_Particle_count = (int *)malloc(sizeof(int) * num_Tissues);
    parents_Prev_generation = (int *)malloc(sizeof(int) * num_Tissues);
    current_Generation = 0;

    set<int> init_removed_by_Transfer_Indexes;

    for (int tissue = 0; tissue < num_Tissues; tissue++)
    {
        vector<string> tissue_Sequence;
        functions.config_Folder(host_Folder + "/" + to_string(tissue), "Tissue " + to_string(tissue + 1));

        tissue_Sequences.push_back(tissue_Sequence);
        current_cell_load_per_Tissue[tissue] = 0;
        dead_Particle_count[tissue] = 0;

        parents_Prev_generation[tissue] = 0;

        removed_by_Transfer_Indexes.push_back(init_removed_by_Transfer_Indexes);
    }
}