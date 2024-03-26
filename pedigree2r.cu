#include "pedigree2r.cuh"

pedigree2r::pedigree2r(string parameter_Master_Location)
{
    cout << "Initializing conversion\n";
    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"Pedigree folder location\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    pedigree_Folder_location = Parameters.get_STRING(found_Parameters[0]);
    cout << "\nLooking at Pedigree folder: " << pedigree_Folder_location << endl
         << endl;

    vector<string> line_Data;

    for (const auto &entry : filesystem::directory_iterator(pedigree_Folder_location))
    {
        if (filesystem::is_regular_file(entry))
        {
            // cout << file_to_Check << endl;
            if (entry.path().extension().string() == ".csv")
            {
                // file_to_Check = entry.path().filename().stem().string();
                function.split(line_Data, entry.path().filename().stem().string(), '_');
                if (line_Data[5] == "pedigree" && line_Data[6] == "Relationships" && line_Data.size() == 7)
                {
                    cout << "Found file: " << entry.path().string() << endl;
                    this->files_to_Process.push_back(make_pair(entry.path().filename().stem().string(), entry.path().string()));
                }
            }
        }
    }

    cout << "\nFound " << files_to_Process.size() << " pedigree Realtionship file(s)\n\n";
}

void pedigree2r::ingress()
{
    functions_library functions = functions_library();
    vector<string> line_Data;

    for (int file = 0; file < files_to_Process.size(); file++)
    {
        cout << "Processing file " << file + 1 << " of " << files_to_Process.size() << endl;
        string file_Location = files_to_Process[file].second;

        cout << "File processing: " << file_Location << endl;

        fstream pedigree_File;
        pedigree_File.open(files_to_Process[file].second, ios::in);

        if (pedigree_File.is_open())
        {

            int generations = -1;
            functions.split(line_Data, files_to_Process[file].first, '_');
            generations = stoi(line_Data[3]);

            cout << "File has " << generations << " generation(s)\n";

            fstream converted_File;
            converted_File.open(pedigree_Folder_location + "/" + files_to_Process[file].first + "_r_Conversion.csv", ios::out);

            if (converted_File.is_open())
            {
                cout << "\nIntialized conversion file: " << pedigree_Folder_location << "/" + files_to_Process[file].first << "_r_Conversion.csv\n";

                for (int gen = 0; gen < generations + 1; gen++)
                {
                    converted_File << gen << "\t";
                }
                converted_File << "sequence_ID\n";
                converted_File.flush();

                cout << "\nProcessing pedigree file\n";

                string line;

                // skip first header line
                char delim = '\t';
                getline(pedigree_File, line);
                functions.split(line_Data, line, delim);
                if (line_Data.size() == 1)
                {
                    delim = ',';
                }

                vector<vector<string>> parents_Progeny;
                vector<string> parents;

                while (getline(pedigree_File, line))
                {
                    functions.split(line_Data, line, delim);

                    int present = -1;
                    for (int check = 0; check < parents.size(); check++)
                    {
                        if (parents[check] == line_Data[0])
                        {
                            present = check;
                            break;
                        }
                    }

                    if (present != -1)
                    {
                        parents_Progeny[present].push_back(line_Data[1]);
                    }
                    else
                    {
                        parents.push_back(line_Data[0]);
                        vector<string> progeny;
                        progeny.push_back(line_Data[1]);
                        parents_Progeny.push_back(progeny);
                    }
                }

                pedigree_File.close();

                cout << parents_Progeny.size() << " parents found\n\n";

                cout << "Mapping converted file\n";

                vector<string> lines;

                for (int round = parents.size() - 1; round >= 0; round--)
                {
                    vector<pair<string, string>> queue;

                    line = parents[round];
                    functions.split(line_Data, line, '_');
                    if (line_Data[3] != "0")
                    {
                        break;
                    }
                    else
                    {
                        cout << "Processing parent: " << line << endl;
                    }

                    for (int progeny = 0; progeny < parents_Progeny[round].size(); progeny++)
                    {
                        functions.split(line_Data, line, '_');
                        queue.push_back(make_pair(parents_Progeny[round][progeny], line_Data[4]));
                    }

                    int track_Queue = 0;
                    do
                    {
                        string find = queue[track_Queue].first;
                        cout << "Finding node: " << find << ": ";
                        functions.split(line_Data, find, '_');
                        queue[track_Queue].second = queue[track_Queue].second + "\t" + line_Data[4];

                        int found = -1;
                        for (int parent_check = 0; parent_check < parents.size(); parent_check++)
                        {
                            if (parents[parent_check] == find)
                            {
                                found = parent_check;
                                cout << "found\n";
                                break;
                            }
                        }

                        if (found != -1)
                        {
                            for (int progeny = 0; progeny < parents_Progeny[found].size(); progeny++)
                            {
                                queue.push_back(make_pair(parents_Progeny[found][progeny], queue[track_Queue].second));
                            }
                        }
                        else
                        {
                            cout << "Writing lines\n";
                            // lines.push_back(queue[track_Queue].second + "\t" + find);
                            converted_File << queue[track_Queue].second << "\n";
                            converted_File.flush();
                        }

                        track_Queue++;

                    } while (track_Queue < queue.size());
                }

                // cout << "Writing lines\n";
                // for (int line_Write = 0; line_Write < lines.size(); line_Write++)
                // {
                //     converted_File << lines[line_Write] << "\n";
                // }
                converted_File.close();
            }
        }
    }

    cout << "\nConversion complete for all pedigree files in folder\n\n";
}