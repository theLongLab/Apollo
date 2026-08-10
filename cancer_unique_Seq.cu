#include "cancer_unique_Seq.cuh"

cancer_unique_Seq::cancer_unique_Seq(string parameter_Master_Location)
{
    cout << "Initiating cancer Sequence collection\n"
         << endl;

    parameter_load Parameters = parameter_load();
    functions_library functions = functions_library();

    vector<string> parameters_List = {"\"Intermediate folders\"", "\"Output folders\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    intermediate_Folder_location = functions.path_Check(Parameters.get_STRING(found_Parameters[0]));
    output_Folder_location = functions.path_Check(Parameters.get_STRING(found_Parameters[1])) + "summary_Haplotypes/";

    if (filesystem::exists(intermediate_Folder_location))
    {
        cout << "Intermediate folder location: " << this->intermediate_Folder_location << endl;
    }
    else
    {
        cerr << "ERROR: UNABLE TO FIND THE INTERMEDIATES FOLDER: " << this->intermediate_Folder_location << endl;
        exit(-1);
    }

    if (!filesystem::exists(output_Folder_location))
    {
        functions.config_Folder(output_Folder_location, "output Folder");
    }
}

void cancer_unique_Seq::ingress()
{
    string tissue_Folder = intermediate_Folder_location + "sequence_Data/cancer_Host/";
    functions_library functions = functions_library();

    vector<string> line_Data;
    string line;

    if (!filesystem::exists(tissue_Folder))
    {
        cerr << "ERROR: UNABLE TO FIND CANCER HOST FOLDER: " << tissue_Folder << endl;
        exit(-1);
    }

    fstream output_File_Alive;
    fstream output_File_Dead;

    output_File_Alive.open(this->output_Folder_location + "alive_Haplotypes.tsv", ios::out);
    output_File_Dead.open(this->output_Folder_location + "dead_Haplotypes.tsv", ios::out);

    if (output_File_Alive.is_open() && output_File_Dead.is_open())
    {
        // string Replication_Prob;
        // string Gen_Death_prob;
        // string Replication_Factor;
        // string metastatic_Prob;
        // string survivability;

        string columns = "Tissue\tGeneration\tHaplotype\tCount\tReplication_Prob\tGen_Death_prob\tReplication_Factor\tmetastatic_Prob\tsurvivability\n";

        output_File_Alive << columns;
        output_File_Alive.close();
        output_File_Dead << columns;
        output_File_Dead.close();
    }
    else
    {
        cerr << "ERROR: UNABLE TO CREATE OUTPUT FILES" << this->output_Folder_location << endl;
        exit(-1);
    }

    for (const auto &tissue : filesystem::directory_iterator(tissue_Folder))
    {
        if (filesystem::is_directory(tissue))
        {
            string tissue_ID = tissue.path().filename().string();
            cout << "\nDetected Tissue: " << tissue_ID << endl;

            for (const auto &generation : filesystem::directory_iterator(tissue.path().string()))
            {
                if (filesystem::is_directory(generation.path().string()))
                {
                    functions.split(line_Data, generation.path().filename().string(), '_');
                    string generation_Index = line_Data[1];
                    cout << "\nProcessing generation: " << generation_Index << endl;

                    unordered_map<string, sequence_Details> sequence_Count_Alive;
                    unordered_map<string, sequence_Details> sequence_Count_Dead;

                    for (const auto &nfasta : filesystem::directory_iterator(generation.path().string()))
                    {
                        if (nfasta.path().extension().string() == ".nfasta")
                        {
                            cout << "Reading file: " << nfasta.path().string() << endl;

                            fstream nfasta_File;
                            nfasta_File.open(nfasta.path().string(), ios::in);

                            if (nfasta_File.is_open())
                            {
                                string dead_Or_Alive = "";

                                string Replication_Prob_temp = "";
                                string Gen_Death_prob_temp;
                                string Replication_Factor_temp;
                                string metastatic_Prob_temp;
                                string survivability_temp;

                                while (getline(nfasta_File, line))
                                {
                                    if (line.empty())
                                    {
                                        continue;
                                    }
                                    if (line.at(0) == '>')
                                    {
                                        functions.split(line_Data, line, '_');

                                        if (line_Data.size() < 8)
                                        {
                                            cerr << "ERROR NO 8 VALUES IN NFASTA FILE: " << nfasta.path().string() << endl;
                                            exit(-1);
                                        }

                                        dead_Or_Alive = line_Data[1];

                                        Replication_Prob_temp = line_Data[2];
                                        Gen_Death_prob_temp = line_Data[3];
                                        Replication_Factor_temp = line_Data[5];
                                        metastatic_Prob_temp = line_Data[6];
                                        survivability_temp = line_Data[7];
                                        // to_write_Sequence_Store_NEXT_Gen.push_back(make_pair(to_string(last_index_Seq_Written) +
                                        //                                      survival_Status +
                                        //                                      to_string(progeny_Configuration_Cancer[progeny_Index][2]) +
                                        //                                      "_" + to_string(progeny_Configuration_Cancer[progeny_Index][1]) +
                                        //                                      "_" + to_string(progeny_Elapsed[progeny_Index]) +
                                        //                                      "_" + to_string(progeny_Configuration_Cancer[progeny_Index][0]) +
                                        //                                      "_" + to_string(progeny_Configuration_Cancer[progeny_Index][3]) +
                                        //                                      "_" + to_string(progeny_Configuration_Cancer[progeny_Index][4]),
                                        //                                  sequence));
                                    }
                                    else
                                    {
                                        if (dead_Or_Alive == "A")
                                        {
                                            sequence_Details &entry = sequence_Count_Alive[line];
                                            if (entry.count == 0) // first time seeing this sequence
                                            {
                                                entry.Replication_Prob = Replication_Prob_temp;
                                                entry.Gen_Death_prob = Gen_Death_prob_temp;
                                                entry.Replication_Factor = Replication_Factor_temp;
                                                entry.metastatic_Prob = metastatic_Prob_temp;
                                                entry.survivability = survivability_temp;
                                            }
                                            entry.count++;
                                        }
                                        else if (dead_Or_Alive == "D")
                                        {
                                            sequence_Details &entry = sequence_Count_Dead[line];
                                            if (entry.count == 0) // first time seeing this sequence
                                            {
                                                entry.Replication_Prob = Replication_Prob_temp;
                                                entry.Gen_Death_prob = Gen_Death_prob_temp;
                                                entry.Replication_Factor = Replication_Factor_temp;
                                                entry.metastatic_Prob = metastatic_Prob_temp;
                                                entry.survivability = survivability_temp;
                                            }
                                            entry.count++;
                                        }
                                        else
                                        {
                                            cerr << "UNIDENTIFIED STATUS: " << dead_Or_Alive << endl;
                                            exit(-1);
                                        }
                                        dead_Or_Alive = "";
                                    }
                                }
                                nfasta_File.close();
                            }
                            else
                            {
                                cerr << "ERROR: UNABLE TO OPEN NFASTA FILE" << nfasta.path().string() << endl;
                            }
                        }
                    }
                    cout << "Total unique sequences ALIVE: " << sequence_Count_Alive.size() << endl;
                    cout << "Total unique sequences DEAD: " << sequence_Count_Dead.size() << endl;

                    output_File_Alive.open(this->output_Folder_location + "alive_Haplotypes.tsv", ios::app);

                    if (output_File_Alive.is_open())
                    {
                        for (const auto &[sequence, details] : sequence_Count_Alive)
                        {
                            output_File_Alive << tissue_ID << "\t"
                                              << generation_Index << "\t"
                                              << sequence << "\t"
                                              << details.count << "\t"
                                              << details.Replication_Prob << "\t"
                                              << details.Gen_Death_prob << "\t"
                                              << details.Replication_Factor << "\t"
                                              << details.metastatic_Prob << "\t"
                                              << details.survivability << "\n";
                        }

                        output_File_Alive.close();
                    }
                    else
                    {
                        cerr << "ERROR: UNABLE TO OPEN ALIVE OUTPUT FILE: " << this->output_Folder_location << "alive_Haplotypes.tsv";
                        exit(-1);
                    }

                    output_File_Dead.open(this->output_Folder_location + "dead_Haplotypes.tsv", ios::app);

                    if (output_File_Dead.is_open())
                    {
                        for (const auto &[sequence, details] : sequence_Count_Dead)
                        {
                            output_File_Dead << tissue_ID << "\t"
                                             << generation_Index << "\t"
                                             << sequence << "\t"
                                             << details.count << "\t"
                                             << details.Replication_Prob << "\t"
                                             << details.Gen_Death_prob << "\t"
                                             << details.Replication_Factor << "\t"
                                             << details.metastatic_Prob << "\t"
                                             << details.survivability << "\n";
                        }

                        output_File_Dead.close();
                    }
                    else
                    {
                        cerr << "ERROR: UNABLE TO OPEN DEAD OUTPUT FILE: " << this->output_Folder_location << "dead_Haplotypes.tsv";
                        exit(-1);
                    }
                }
            }
        }
    }
}