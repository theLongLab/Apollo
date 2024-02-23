#include "bfs.cuh"

bfs::bfs(string parameter_Master_Location)
{
    cout << "Initiating Breath First Search to identify pedigree of a sequence\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Pedigree Node ID\"",
        "\"Pedigree Tissue\"",
        "\"Pedigree Generation\"",
        "\"Pedigree Sequence\"",
        "\"Nodes master profile\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    string node_main_Index = "";
    int tissue_main_Index = -1;
    string nsequence = "";
    int num_tissues_per_Node = 0;

    this->intermediate_Folder_location = Parameters.get_STRING(found_Parameters[0]);
    this->output_Folder_location = Parameters.get_STRING(found_Parameters[1]);

    cout << "\nReading target sequence\n";
    string pedigree_Sequence_loation = Parameters.get_STRING(found_Parameters[5]);
    fstream pedigree_File;
    pedigree_File.open(pedigree_Sequence_loation, ios::in);

    if (pedigree_File.is_open())
    {
        string line;
        // skip first line;
        getline(pedigree_File, line);
        // sequence line
        getline(pedigree_File, line);
        // cout << line << endl;
        nsequence = line;
        pedigree_File.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN TARGET SEQUENCE FILE: " << pedigree_Sequence_loation << endl;
        exit(-1);
    }

    for (int base = 0; base < nsequence.size(); base++)
    {
        if (nsequence[base] == 'A' || nsequence[base] == 'a')
        {
            nsequence[base] = '0';
        }
        else if (nsequence[base] == 'T' || nsequence[base] == 't')
        {
            nsequence[base] = '1';
        }
        else if (nsequence[base] == 'G' || nsequence[base] == 'g')
        {
            nsequence[base] = '2';
        }
        else if (nsequence[base] == 'C' || nsequence[base] == 'c')
        {
            nsequence[base] = '3';
        }
        else if (nsequence[base] != '0' && nsequence[base] != '1' && nsequence[base] != '2' && nsequence[base] != '3')
        {
            cout << "ERROR: UNRECOGNISED BASE IN SEQUENCE: " << nsequence[base];
            exit(-1);
        }
    }

    cout << "Target sequence loaded\n";
    // cout << nsequence << endl;

    string node_ID = Parameters.get_STRING(found_Parameters[2]);
    cout << "\nGetting index of target node: " << node_ID << "\n";

    string node_Index_file_location = intermediate_Folder_location + "/index_Data/node_Index.csv";
    fstream node_index_File;
    node_index_File.open(node_Index_file_location, ios::in);

    if (node_index_File.is_open())
    {
        string line;
        vector<string> line_Data;

        // skip first header line
        getline(node_index_File, line);

        while (getline(node_index_File, line))
        {
            function.split(line_Data, line, '\t');
            node_Indexes.push_back(make_pair(stoi(line_Data[0]), line_Data[1]));
            if (line_Data[1] == node_ID)
            {
                node_main_Index = line_Data[0];
                // break;
            }
        }
        node_index_File.close();
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN NODE INDEX FILE: " << node_Index_file_location;
        exit(-1);
    }

    if (node_main_Index != "")
    {
        cout << "Target node's index: " << node_main_Index << endl;
    }
    else
    {
        cout << "ERROR: UNABLE TO FIND THE NODE ID: " << node_ID << endl;
        exit(-1);
    }

    sort(node_Indexes.begin(), node_Indexes.end());

    string tissue_Name = Parameters.get_STRING(found_Parameters[3]);
    cout << "\nGetting index of tissue: " << tissue_Name << endl;

    string node_Master_location = Parameters.get_STRING(found_Parameters[6]);

    vector<pair<string, string>> Tissue_profiles_block_Data = Parameters.get_block_from_File(node_Master_location, "Tissue profiles");
    num_tissues_per_Node = Parameters.get_INT(Tissue_profiles_block_Data, "Number of tissues");

    if (num_tissues_per_Node > 0)
    {
        cout << "\nNumber of tissues in a node: " << num_tissues_per_Node << endl;

        for (int tissue = 0; tissue < num_tissues_per_Node; tissue++)
        {
            string check_Tissue = "Tissue " + to_string(tissue + 1) + " Name";
            tissue_Names.push_back(Parameters.get_STRING(Tissue_profiles_block_Data, check_Tissue));
            cout << check_Tissue << ": " << tissue_Names[tissue] << endl;

            if (tissue_Names[tissue] == tissue_Name)
            {
                tissue_main_Index = tissue;
            }
        }
    }
    else
    {
        cout << "ERROR: TISSUE NUMBER HAS TO BE GREATER THAN ZERO.\n\n";
    }

    if (tissue_main_Index != -1)
    {
        cout << "\nTissue index of " << tissue_Name << " : " << tissue_main_Index << endl;
    }
    else
    {
        cout << "UNABLE TO FIND TISSUE INDEX: " << tissue_Name << endl;
        exit(-1);
    }

    generation = Parameters.get_INT(found_Parameters[4]);

    cout << "\nIndentifying matching sequences from target\n";
    string sequence_Search_folder = intermediate_Folder_location + "/sequence_Data/" + node_main_Index;

    int re_tar_sequence_Folder = check_Tar_Folder(sequence_Search_folder);

    sequence_Search_folder = sequence_Search_folder + "/" + to_string(tissue_main_Index);

    int re_tar_Tissue_folder = check_Tar_Folder(sequence_Search_folder);

    sequence_Search_folder = sequence_Search_folder + "/generation_" + to_string(generation);

    int re_tar_Generation = check_Tar_Folder(sequence_Search_folder);

    for (const auto &entry : filesystem::directory_iterator(sequence_Search_folder))
    {
        if (filesystem::is_regular_file(entry) && entry.path().extension() == ".nfasta")
        {
            string check_File_location = entry.path().string();
            cout << "Checking nfasta file: " << check_File_location << endl;

            fstream nfasta_File;
            nfasta_File.open(check_File_location, ios::in);
            if (nfasta_File.is_open())
            {
                string line;
                string header;
                while (getline(nfasta_File, line))
                {
                    if (line.at(0) == '>')
                    {
                        header = line;
                    }
                    else
                    {
                        if (line == nsequence)
                        {
                            header = header.substr(1);
                            header = header.substr(0, header.find('_'));
                            search_sequence_IDs.push_back(node_ID + "_" + tissue_Name + "_" + to_string(generation) + "_" + header);
                        }
                    }
                }

                nfasta_File.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN NFASTA FILE: " << check_File_location << endl;
                exit(-1);
            }
        }
    }

    if (search_sequence_IDs.size() > 0)
    {
        cout << "\nFound " << search_sequence_IDs.size() << " matching sequence(s)\n";
        // this->current_node_ID = node_ID;
        //  cout << search_sequence_IDs[0] << endl;

        if (!filesystem::exists(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information"))
        {
            filesystem::create_directory(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information");
        }

        if (!filesystem::exists(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information/" + filesystem::path(pedigree_Sequence_loation).stem().string()))
        {
            filesystem::create_directory(this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information/" + filesystem::path(pedigree_Sequence_loation).stem().string());
        }

        pedigree_Folder_location = this->output_Folder_location + "/node_Data/" + node_ID + "/pedigree_Information/" + filesystem::path(pedigree_Sequence_loation).stem().string();
    }
    else
    {
        cout << "ERROR SEQUENCE NOT FOUND\n";
        exit(-1);
    }
}

void bfs::ingress()
{
    functions_library functions = functions_library();

    for (int sequence = 0; sequence < search_sequence_IDs.size(); sequence++)
    {
        string ID_Sequence_Original = search_sequence_IDs[sequence];
        cout << "\nIdentifying pedigree of sequence: " << ID_Sequence_Original << endl;

        vector<string> sequence_Information;
        functions.split(sequence_Information, ID_Sequence_Original, '_');
        string current_node_ID = sequence_Information[0] + "_" + sequence_Information[1];

        vector<pair<string, string>> queue;
        queue.push_back(make_pair(sequence_Information[0] + "_" + sequence_Information[1], ID_Sequence_Original));
        int queue_Count = 0;

        functions.create_File(pedigree_Folder_location + "/" + ID_Sequence_Original + "_pedigree_Relationships.csv", "Source\tTarget\tType");
        functions.create_File(pedigree_Folder_location + "/" + ID_Sequence_Original + "_sequence_Information.csv", "ID\tHost\tTissue\tGeneration");

        set<string> sequence_IDs;
        sequence_IDs.insert(ID_Sequence_Original);

        vector<string> line_Data;

        do
        {
            string ID_Sequence = queue[queue_Count].second;
            string line;

            fstream node_File;
            node_File.open(this->output_Folder_location + "/node_Data/" + queue[queue_Count].first + "/sequence_parent_Progeny_relationships.csv",
                           ios::in);

            // initialize this using node generational_Summary
            vector<pair<string, string>> progeny_Parent;
            vector<string> Type;

            int collect = -1;

            if (node_File.is_open())
            {
                cout << "\nReading node parent progeny file: " << queue[queue_Count].first << endl;

                // skip header
                getline(node_File, line);

                while (getline(node_File, line))
                {
                    functions.split(line_Data, line, '\t');
                    if (collect == 1 && line_Data[1] != ID_Sequence)
                    {
                        cout << "Captured target sequence: " << ID_Sequence << endl;
                        break;
                    }
                    progeny_Parent.push_back(make_pair(line_Data[1], line_Data[0]));
                    Type.push_back(line_Data[2]);
                    if (line_Data[1] == ID_Sequence)
                    {
                        collect = 1;
                    }
                }
                node_File.close();
            }
            else
            {
                cout << "ERROR UNABLE TO OPEN PARENT PROGENY FILE: " << this->output_Folder_location << "/node_Data/" << queue[queue_Count].first << "/sequence_parent_Progeny_relationships.csv" << endl;
                exit(-1);
            }

            cout << "Collecting parent progeny pedigree\n";

            int get_Parents = 0;
            vector<string> progeny_List;
            progeny_List.push_back(ID_Sequence);
            int track_progeny_Parent = progeny_Parent.size() - 1;

            fstream source_Target_pedigree;
            source_Target_pedigree.open(pedigree_Folder_location + "/" + ID_Sequence_Original + "_pedigree_Relationships.csv", ios::app);

            if (source_Target_pedigree.is_open())
            {
                int catch_Check = -1;
                do
                {
                    if (progeny_Parent[track_progeny_Parent].first == progeny_List[get_Parents])
                    {
                        catch_Check = 1;
                        source_Target_pedigree << progeny_Parent[track_progeny_Parent].second << "\t"
                                               << progeny_Parent[track_progeny_Parent].first << "\t"
                                               << Type[track_progeny_Parent] << "\n";

                        // if parent not present in progeny list add it to progeny list
                        auto it = sequence_IDs.find(progeny_Parent[track_progeny_Parent].second);
                        if (it == sequence_IDs.end())
                        {
                            sequence_IDs.insert(progeny_Parent[track_progeny_Parent].second);
                            // check if new host
                            functions.split(line_Data, progeny_Parent[track_progeny_Parent].second, '_');
                            string check_Node = line_Data[0] + "_" + line_Data[1];
                            if (check_Node != queue[queue_Count].first)
                            {
                                queue.push_back(make_pair(check_Node, progeny_Parent[track_progeny_Parent].second));
                            }
                            else
                            {
                                progeny_List.push_back(progeny_Parent[track_progeny_Parent].second);
                            }
                        }
                    }
                    track_progeny_Parent--;

                    if (catch_Check == 1 && progeny_Parent[track_progeny_Parent].first != progeny_List[get_Parents])
                    {
                        get_Parents++;
                        catch_Check = -1;
                    }

                } while (get_Parents < progeny_List.size() && track_progeny_Parent >= 0);

                source_Target_pedigree.close();
            }
            else
            {
                cout << "ERROR: UNABLE TO OPEN FILE: " << pedigree_Folder_location << "/" << ID_Sequence_Original << "_pedigree_Relationships.csv" << endl;
                exit(-1);
            }

            queue_Count++;
        } while (queue_Count < queue.size());

        cout << "Writing sequence information\n";

        vector<string> IDs_complete(sequence_IDs.begin(), sequence_IDs.end());
        fstream sequence_Information_File;
        sequence_Information_File.open(pedigree_Folder_location + "/" + ID_Sequence_Original + "_sequence_Information.csv", ios::app);
        if (sequence_Information_File.is_open())
        {
            for (string write_Line : IDs_complete)
            {
                sequence_Information_File << write_Line << "\t";
                functions.split(line_Data, write_Line, '_');
                sequence_Information_File << line_Data[0] << "_" << line_Data[1] << "\t"
                                          << line_Data[2] << "\t" << line_Data[3] << endl;
            }
            sequence_Information_File.close();
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN FILE: " << pedigree_Folder_location << "/" << ID_Sequence_Original << "_sequence_Information.csv\n";
            exit(-1);
        }
    }
}

int bfs::check_Tar_Folder(string location)
{
    int re_Tar = -1;

    if (filesystem::exists(location + ".tar"))
    {
        cout << "Extracting folder: " << location + ".tar\n";
        string command = "tar -xf" + location + ".tar -C .";

        int result = system(command.c_str());

        if (result == 0)
        {
            // The command executed successfully
            cout << "Successfully untarred the folder." << endl;
        }
        else
        {
            // An error occurred during the execution of the command
            cout << "Failed to untar the folder." << endl;
            exit(-1);
        }

        re_Tar = 1;
    }

    return re_Tar;
}
