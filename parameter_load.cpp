#include "parameter_load.h"

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

vector<pair<string, string>> parameter_load::get_block_from_File(string &file_parameter_Location, string block_Header)
{
    vector<pair<string, string>> block_Data;
    int activate_Collection = 0;

    functions_library function = functions_library();

    fstream parameter_File;
    parameter_File.open(file_parameter_Location, ios::in);

    if (parameter_File.is_open())
    {
        string line;
        getline(parameter_File, line);

        vector<string> line_Data;

        while (getline(parameter_File, line))
        {
            line_Data = clean_Line(line, function);

            if (line_Data.size() != 0)
            {
                if (line_Data[0] == "\"" + block_Header + "\"")
                {
                    activate_Collection = 1;
                    // cout << line_Data[0] << endl;
                    break;
                }
            }
        }

        if (activate_Collection == 1)
        {
            int count_bracket = 0;

            while (getline(parameter_File, line))
            {
                line_Data = clean_Line(line, function);

                if (line_Data.size() != 0)
                {
                    if (line_Data.size() > 1)
                    {
                        // cout << line << endl;
                        // cout << line_Data[1] << endl;
                        // exit(-1);
                        if (line_Data[1] == "{")
                        {
                            count_bracket++;
                        }

                        if (count_bracket >= 0)
                        {
                            block_Data.push_back(make_pair(line_Data[0], line_Data[1]));
                        }
                        else
                        {
                            break;
                        }
                    }
                    else if (line_Data[0] == "}")
                    {
                        count_bracket--;

                        if (count_bracket >= 0)
                        {
                            block_Data.push_back(make_pair(line_Data[0], ""));
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            cout << "SYSTEM ERROR: " << block_Header << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
            exit(-1);
        }

        parameter_File.close();
    }

    return block_Data;
}

vector<pair<string, string>> parameter_load::get_block_from_block(vector<pair<string, string>> &block, string block_Header)
{
    vector<pair<string, string>> block_Data;

    int catch_Index = -1;

    for (int check = 0; check < block.size(); check++)
    {
        if (block[check].first == "\"" + block_Header + "\"")
        {
            catch_Index = check + 1;
            break;
        }
    }

    if (catch_Index != -1)
    {
        int count_bracket = 0;

        for (int check = catch_Index; check < block.size(); check++)
        {
            if (block[check].second == "{")
            {
                count_bracket++;
            }
            else if (block[check].first == "}")
            {
                count_bracket--;
            }

            if (count_bracket >= 0)
            {
                block_Data.push_back(make_pair(block[check].first, block[check].second));
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        cout << "SYSTEM ERROR: " << block_Header << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }

    return block_Data;
}

vector<pair<string, string>> parameter_load::check_block_from_block(vector<pair<string, string>> &block, string block_Header)
{
    vector<pair<string, string>> block_Data;

    int catch_Index = -1;

    for (int check = 0; check < block.size(); check++)
    {
        if (block[check].first == "\"" + block_Header + "\"")
        {
            catch_Index = check + 1;
            break;
        }
    }

    if (catch_Index != -1)
    {
        int count_bracket = 0;

        for (int check = catch_Index; check < block.size(); check++)
        {
            if (block[check].second == "{")
            {
                count_bracket++;
            }
            else if (block[check].first == "}")
            {
                count_bracket--;
            }

            if (count_bracket >= 0)
            {
                block_Data.push_back(make_pair(block[check].first, block[check].second));
            }
            else
            {
                break;
            }
        }
    }

    return block_Data;
}

int parameter_load::get_INT(vector<pair<string, string>> block, string value)
{
    int integer_Value;
    int found = -1;

    for (int i = 0; i < block.size(); i++)
    {
        if (block[i].first == "\"" + value + "\"")
        {
            integer_Value = stoi(block[i].second);
            found = 1;
            break;
        }
    }

    if (found == 1)
    {
        return integer_Value;
    }
    else
    {
        cout << "SYSTEM ERROR: " << value << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }
}

float parameter_load::get_FLOAT(vector<pair<string, string>> block, string value)
{
    float return_Value;
    int found = -1;

    for (int i = 0; i < block.size(); i++)
    {
        if (block[i].first == "\"" + value + "\"")
        {
            return_Value = stof(block[i].second.substr(1, block[i].second.length() - 2));
            found = 1;
            break;
        }
    }

    if (found == 1)
    {
        return return_Value;
    }
    else
    {
        cout << "SYSTEM ERROR: " << value << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }
}

string parameter_load::get_STRING(vector<pair<string, string>> block, string value)
{
    string return_Value;
    int found = -1;

    for (int i = 0; i < block.size(); i++)
    {
        if (block[i].first == "\"" + value + "\"")
        {
            return_Value = block[i].second.substr(1, block[i].second.length() - 2);
            found = 1;
            break;
        }
    }

    if (found == 1)
    {
        return return_Value;
    }
    else
    {
        cout << "SYSTEM ERROR: " << value << " DOES NOT EXIST IN THE CURRENT BLOCK.\n";
        exit(-1);
    }
}

vector<string> parameter_load::clean_Line(string line, functions_library &function)
{
    vector<string> line_Data;

    if (line != "")
    {
        // string trim_Line = line;
        int i = 0;
        while (line[i] == ' ')
        {
            i++; // Skip leading spaces
        }

        line.erase(0, i);

        if (line.at(0) != '#')
        {
            if (line.at(line.size() - 1) == ',')
            {
                line = line.substr(0, line.length() - 1);
            }
            function.split(line_Data, line, ':');
        }
    }

    return line_Data;
}