#include "cancer.cuh"

cancer::cancer(string parameter_Master_Location)
{
    cout << "Initiating cancer simulation\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"CUDA Device IDs\"",
        "\"CPU cores\"",
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Multi read\"",
        "\"Nodes master profile\"",
        "\"Sequence master profile\"",
        "\"Intermediate Sequences per file\"",
        "\"Process cell rate\"",
        "\"Start date\"",
        "\"Stop after generations\"",
        "\"Enable folder management\"",
        "\"First infection\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    cout << "\nConfiguring folders:\n";

    intermediate_Folder_location = Parameters.get_STRING(found_Parameters[2]);
    output_Folder_location = Parameters.get_STRING(found_Parameters[3]);

    function.config_Folder(intermediate_Folder_location, "Intermediate");
    function.config_Folder(output_Folder_location, "Output");

    max_sequences_per_File = Parameters.get_INT(found_Parameters[7]);
    max_Cells_at_a_time = Parameters.get_INT(found_Parameters[8]);
    start_Date = Parameters.get_STRING(found_Parameters[9]);
    first_Infection = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[12]));

    if (max_sequences_per_File <= 0)
    {
        cout << "ERROR: Intermediate Sequences per file PARAMETER MUST BE GREATER THAN ZERO.\n";
        exit(-1);
    }

    if (max_Cells_at_a_time <= 0)
    {
        cout << "ERROR: Process cell rate PARAMETER MUST BE GREATER THAN ZERO.\n";
        exit(-1);
    }

    cout << "\nFolder management: ";

    if (function.to_Upper_Case(Parameters.get_STRING(found_Parameters[11])) == "YES")
    {
        enable_Folder_management = "YES";

        vector<string> folder_Management = {"\"Compress folders\""};
        vector<string> folder_management_Parameters = Parameters.get_parameters(parameter_Master_Location, folder_Management);

        if (function.to_Upper_Case(Parameters.get_STRING(folder_Management[0])) == "YES")
        {
            enable_Compression = "YES";
        }
    }
    cout << enable_Folder_management << endl;
    cout << "Folder compression: " << enable_Compression << endl;

    cout << "\nConfiguring node master profiles:\n";
    node_Master_location = Parameters.get_STRING(found_Parameters[5]);

    cout << "\nConfiguring sequence master profiles:\n";
    sequence_Master_location = Parameters.get_STRING(found_Parameters[6]);

    stop_after_generations = function.to_Upper_Case(Parameters.get_STRING(found_Parameters[10]));

    if (stop_after_generations == "YES")
    {
        cout << "\nConfiguring simulation termination by overall generations run\n";

        vector<string> parameters_List_stop_Gen_Mode = {"\"Mode to stop\""};
        vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

        if (function.to_Upper_Case(Parameters.get_STRING(stop_Gen_Parameters[0])) == "GENERATIONS")
        {
            stop_gen_Mode = 0;
            cout << "Stop after given number of generations: ";
            parameters_List_stop_Gen_Mode.clear();
            stop_Gen_Parameters.clear();

            vector<string> parameters_List_stop_Gen_Mode = {"\"Number of generations\""};
            vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

            stop_generations_Count = Parameters.get_INT(stop_Gen_Parameters[0]);
            cout << stop_generations_Count << endl;
        }
        else
        {
            cout << "Stop after given date: ";
            stop_gen_Mode = 1;
            parameters_List_stop_Gen_Mode.clear();
            stop_Gen_Parameters.clear();

            vector<string> parameters_List_stop_Gen_Mode = {"\"End date\""};
            vector<string> stop_Gen_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List_stop_Gen_Mode);

            string stop_Date_String = Parameters.get_STRING(stop_Gen_Parameters[0]);

            cout << stop_Date_String << endl;
            vector<string> split_Date;
            function.split(split_Date, stop_Date_String, '-');

            stop_Date = function.date_to_Decimal(stoi(split_Date[0]), stoi(split_Date[1]), stoi(split_Date[2]));
            cout << "Decimal date: " << stop_Date << endl;
        }
    }
    else
    {
        cout << "\nSimulation termination by overall generations run is deactivated\n";
    }

    cout << "\nConfiguring hardware resources:\n\n";
    this->CPU_cores = Parameters.get_INT(found_Parameters[1]);
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    cout << "Maximum cells processed at a time: " << this->max_Cells_at_a_time << endl
         << endl;

    cout << "Maximum sequences per file (Intermediary): " << this->max_sequences_per_File << endl
         << endl;
    this->multi_Read = Parameters.get_STRING(found_Parameters[4]);
    transform(multi_Read.begin(), multi_Read.end(), multi_Read.begin(), ::toupper);
    cout << "Multiple read and write: " << this->multi_Read << endl
         << endl;

    string cuda_IDs_String = Parameters.get_STRING(found_Parameters[0]);
    vector<string> cuda_IDs;
    function.split(cuda_IDs, cuda_IDs_String, ',');
    if (cuda_IDs.size() > 0)
    {
        this->num_Cuda_devices = cuda_IDs.size();
        CUDA_device_IDs = (int *)malloc(sizeof(int) * num_Cuda_devices);
        tot_Blocks = (int *)malloc(sizeof(int) * num_Cuda_devices);
        tot_ThreadsperBlock = (int *)malloc(sizeof(int) * num_Cuda_devices);
        function.print_Cuda_devices(cuda_IDs, this->CUDA_device_IDs, num_Cuda_devices, this->tot_Blocks, this->tot_ThreadsperBlock);
    }
    else
    {
        cout << "ERROR: THERE HAS TO BE AT LEAST ONE CUDA DEVICE SELECTED\n";
        exit(-1);
    }
}

void cancer::ingress()
{
    cout << "\nConfiguring parameters\n";
}