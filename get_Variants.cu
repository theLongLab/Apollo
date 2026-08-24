#include "get_Variants.cuh"

get_Variants::get_Variants(string reference_Sequence_location, string variant_File_location,
                           string output_File_location)
{
    cout << "Intialize variant calling" << endl;

    if (filesystem::exists(reference_Sequence_location))
    {
        fstream reference_File;
        reference_File.open(reference_Sequence_location, ios::in);
        if (reference_File.is_open())
        {
            string line;
            while (getline(reference_File, line))
            {
                if (line.at(0) != '>')
                {
                    reference_Sequence = line;
                }
            }
            reference_File.close();

            if (reference_Sequence.size() <= 0)
            {
                cerr << "ERROR: NO REFERENCE SEQUENCE FOUND, LENGTH IS 0" << endl;
                exit(-1);
            }
            cout << "Reference sequence length: " << reference_Sequence.size() << endl;
            cout << "Reference sequence: " << reference_Sequence << endl;
        }
        else
        {
            cerr << "ERROR: UNABLE TO OPEN REFERENCE FILE: " << reference_Sequence_location << endl;
            exit(-1);
        }
    }
    else
    {
        cerr << "ERROR: UNABLE TO FIND REFERENCE FILE: " << reference_Sequence_location << endl;
        exit(-1);
    }

    if (filesystem::exists(variant_File_location))
    {
        fstream variant_File;
        variant_File.open(variant_File_location, ios::in);

        if (variant_File.is_open())
        {
            string line;
            // skip header;
            getline(variant_File, line);
            while (getline(variant_File, line))
            {
                variant_Lines.push_back(line);
            }
            variant_File.close();

            lines_To_process = variant_Lines.size();
            cout << lines_To_process << " line(s) collected" << endl;
        }
        else
        {
            cerr << "ERROR: UNABLE TO OPEN VARIANT FILE: " << variant_File_location << endl;
            exit(-1);
        }
    }
    else
    {
        cerr << "ERROR: UNABLE TO FIND VARIANT FILE: " << variant_File_location << endl;
        exit(-1);
    }

    num_Threads = thread::hardware_concurrency();
    cout << "\nUsing " << num_Threads << " thread(s)" << endl;

    this->output_File_location = output_File_location;
}

void get_Variants::ingress()
{
    functions_library functions = functions_library();

    vector<thread> threads_vec;

    // shared_mutex g_mutex;

    int actual_thread_Use = num_Threads;
    vector<pair<int, int>> start_Stop_threads = functions.fixed_thread_start_Stop(this->num_Threads, lines_To_process, actual_thread_Use);

    num_Threads = actual_thread_Use;

    cout << "Processing lines" << endl;

    for (int thread_ID = 0; thread_ID < actual_thread_Use; thread_ID++)
    {
        threads_vec.push_back(thread{&get_Variants::process_Lines, this, thread_ID + 1, start_Stop_threads[thread_ID].first, start_Stop_threads[thread_ID].second, ref(functions)});
    }

    for (thread &t : threads_vec)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    threads_vec.clear();

    fstream write_Out;
    write_Out.open(this->output_File_location, ios::out);

    if (write_Out.is_open())
    {
        cout << "Writing to file: " << output_File_location << flush;
        write_Out << "Tissue\tGeneration\tHaplotype\tCount\tReplication_Prob\tGen_Death_prob\tReplication_Factor\tmetastatic_Prob\tsurvivability\tSequence\tMutations" << endl;

        for (int line = 0; line < variant_Lines.size(); line++)
        {
            write_Out << variant_Lines[line] << "\n";
        }
        write_Out.close();
        cout << ": DONE" << endl;
    }
    else
    {
        cerr << "ERROR: UNABLE TO CREATE OUTPUT FILE: " << output_File_location << endl;
        exit(-1);
    }
}

void get_Variants::process_Lines(int thread_ID, int start, int stop, functions_library &functions)
{
    string status = "Completed processing " + to_string(thread_ID) +
                    " of " + to_string(num_Threads);

    vector<string> line_Data;
    for (int line = start; line < stop; line++)
    {
        functions.split(line_Data, variant_Lines[line], '\t');

        string numeric_Sequence = line_Data[2];
        if (numeric_Sequence.size() != reference_Sequence.size())
        {
            cerr << "ERROR: SEQUENCE LENGTH MISMATCH: " << variant_Lines[line] << endl;
            exit(-1);
        }
        string sequence = "";

        string mutations_Captured = "";

        for (int base_Index = 0; base_Index < numeric_Sequence.size(); base_Index++)
        {
            char base = numeric_Sequence[base_Index] == '0' ? 'A' : numeric_Sequence[base_Index] == '1' ? 'T'
                                                                : numeric_Sequence[base_Index] == '2'   ? 'G'
                                                                                                        : 'C';

            sequence += base;

            if (base != reference_Sequence[base_Index])
            {
                mutations_Captured.append(reference_Sequence[base_Index] + to_string(base_Index + 1) + base + ";");
            }
        }

        if (mutations_Captured != "")
        {
            if (mutations_Captured[mutations_Captured.size() - 1] == ';')
            {
                mutations_Captured = mutations_Captured.substr(0, mutations_Captured.size() - 1);
            }
        }
        else
        {
            mutations_Captured = "NA";
        }

        variant_Lines[line].append("\t" + sequence + "\t" + mutations_Captured);
    }

    cout << status << endl;
}