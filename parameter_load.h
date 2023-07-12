#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <filesystem>
#include <list>
#include <set>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

class parameter_load
{
private:
    string file_location;

public:
    parameter_load(string file_Location);
    parameter_load();

    void get_parameters(int &CUDA_device_ID, string &parent_SEQ_folder,
                        float &mean_rep_time, float &standard_deviation_rep_time,
                        float &mean_days_host, float &standard_deviation_host_time);

    vector<string> get_parameters(string file_Location, vector<string> &parameters_List);

    int get_INT(string value);
    string get_STRING(string value);
    float get_FLOAT(string value);
};