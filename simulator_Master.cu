#include "simulator_Master.cuh"
#include "functions_library.cuh"

simulator_Master::simulator_Master(string parameter_Master_Location)
{
    cout << "Intializing Simulator (based on CATE's engine)\n";

    parameter_load Parameters = parameter_load();
    functions_library function = functions_library();

    vector<string> parameters_List = {
        "\"CUDA Device ID\"",
        "\"CPU cores\"",
        "\"GPU max units\"",
        "\"Intermediate folders\"",
        "\"Output folders\"",
        "\"Multi read\"",
        "\"Network profile\""};

    vector<string> found_Parameters = Parameters.get_parameters(parameter_Master_Location, parameters_List);

    cout << "Configuring folders:\n";

    output_Folder_location = Parameters.get_STRING(found_Parameters[4]);
    intermediate_Folder_location = Parameters.get_STRING(found_Parameters[3]);

    function.config_Folder(intermediate_Folder_location, "Intermediate");
    function.config_Folder(output_Folder_location, "Output");

    cout << "\nConfiguring hardware resources:\n\n";
    this->CPU_cores = Parameters.get_INT(found_Parameters[1]);
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    this->multi_Read = Parameters.get_STRING(found_Parameters[5]);
    cout << "Multiple read and write: " << this->multi_Read << endl
         << endl;

    this->CUDA_device_number = Parameters.get_INT(found_Parameters[0]);
    function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

    this->gpu_Limit = Parameters.get_INT(found_Parameters[2]);

    cout << "Per round GPU max unit: " << this->gpu_Limit << endl
         << endl;

    configure_Network_Profile(Parameters.get_STRING(found_Parameters[6]), Parameters);
}

void simulator_Master::configure_Network_Profile(string network_Profile_File, parameter_load &Parameters)
{
    cout << "Configuring network profile: " << network_Profile_File << endl;
}
