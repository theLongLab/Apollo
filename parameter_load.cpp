#include "parameter_load.h"
#include "functions_library.cuh"

parameter_load::parameter_load()
{
}

parameter_load::parameter_load(string file_location)
{
    this->file_location = file_location;
}

vector<string> parameter_load::get_parameters(string file_Location, vector<string> &parameters_List)
{
    vector<pair<string, int>> parameter_List_to_Find;
    vector<string> found_Parameters;
    int number_to_Find = parameters_List.size();

    for (size_t i = 0; i < number_to_Find; i++)
    {
        parameter_List_to_Find.push_back(make_pair(parameters_List[i], i));
        found_Parameters.push_back("");
    }

    functions_library function = functions_library();

    fstream parameter_File;
    parameter_File.open(file_Location, ios::in);

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        vector<string> line_Data;

        while (getline(parameter_File, line))
        {
            if (line != "}" && line != "")
            {
                string remove = line;
                int i = 0;
                while (remove[i] == ' ')
                {
                    i++; // Skip leading spaces
                }
                remove.erase(0, i);

                if (remove.at(0) != '#')
                {
                    if (remove.at(remove.size() - 1) == ',')
                    {
                        remove = remove.substr(0, remove.length() - 1);
                    }

                    // cout << remove << endl;
                    function.split(line_Data, remove, ':');

                    for (int check_Pos = 0; check_Pos < parameter_List_to_Find.size(); check_Pos++)
                    {
                        if (line_Data[0] == parameter_List_to_Find[check_Pos].first)
                        {
                            found_Parameters[parameter_List_to_Find[check_Pos].second] = line_Data[1];
                            parameter_List_to_Find.erase(parameter_List_to_Find.begin() + check_Pos);

                            break;
                        }
                    }
                }
            }
            if (parameter_List_to_Find.size() == 0)
            {
                break;
            }
        }
        parameter_File.close();
    }

    // check if all were found
    for (size_t i = 0; i < found_Parameters.size(); i++)
    {
        // cout << found_Parameters[i] << endl;
        if (found_Parameters[i] == "")
        {
            cout << "ERROR:\nCHECK FILE " << file_Location << "\nThe following parameter is missing: " << parameters_List[i] << endl;
            exit(-1);
        }
    }

    return found_Parameters;
}

void parameter_load::get_parameters(int &CUDA_device_ID, string &parent_SEQ_folder,
                                    float &mean_rep_time, float &standard_deviation_rep_time,
                                    float &mean_days_host, float &standard_deviation_host_time)
{
    functions_library function = functions_library();

    fstream parameter_File;
    parameter_File.open(this->file_location, ios::in);

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        vector<string> line_Data;

        while (getline(parameter_File, line))
        {
            if (line != "}" && line != "")
            {
                string remove = line;
                int i = 0;
                while (remove[i] == ' ')
                {
                    i++; // Skip leading spaces
                }
                remove.erase(0, i);

                if (remove.at(0) != '#')
                {
                    if (remove.at(remove.size() - 1) == ',')
                    {
                        remove = remove.substr(0, remove.length() - 1);
                    }

                    // cout << remove << endl;
                    function.split(line_Data, remove, ':');

                    if (line_Data[0] == "\"CUDA Device ID\"")
                    {
                        // cout << line_Data[1] << endl;
                        CUDA_device_ID = get_INT(line_Data[1]);
                        // cout << CUDA_device_ID << endl;
                    }
                    else if (line_Data[0] == "\"Parent sequences folder\"")
                    {
                        // cout << line_Data[1] << endl;
                        parent_SEQ_folder = get_STRING(line_Data[1]);
                        // cout << parent_SEQ_folder << endl;
                    }
                    else if (line_Data[0] == "\"Mean replication time\"")
                    {
                        mean_rep_time = get_FLOAT(line_Data[1]);
                    }
                    else if (line_Data[0] == "\"Standard_deviation replication time\"")
                    {
                        standard_deviation_rep_time = get_FLOAT(line_Data[1]);
                    }
                    else if (line_Data[0] == "\"Mean days in host\"")
                    {
                        mean_days_host = get_FLOAT(line_Data[1]);
                    }
                    else if (line_Data[0] == "\"Standard_deviation in host\"")
                    {
                        standard_deviation_host_time = get_FLOAT(line_Data[1]);
                    }
                }
            }
        }

        parameter_File.close();
    }
}

float parameter_load::get_FLOAT(string value)
{
    return stof(value.substr(1, value.length() - 2));
}

int parameter_load::get_INT(string value)
{
    return stoi(value);
}

string parameter_load::get_STRING(string value)
{
    return (value.substr(1, value.length() - 2));
}