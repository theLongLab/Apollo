#include "within_host_test_2.cuh"
#include "functions_library.cuh"
#include "parameter_load.h"

within_host_test_2::within_host_test_2(int CUDA_device_number, int CPU_cores, int gpu_Limit, string output_Folder, string reference_genome_Location, string replication_profile_file, mt19937 gen, string sequence_Profile_file, string write_Progeny_parents, string write_Sequences, string intermediate_Folders, string multi_READ, int at_a_Time_cells)
{
    functions_library function = functions_library();
    this->gen = gen;

    this->CPU_cores = CPU_cores;
    cout << "Available CPU cores: " << this->CPU_cores << endl
         << endl;

    this->multi_READ = multi_READ;
    cout << "Multiple read and write: " << this->multi_READ << endl
         << endl;

    this->CUDA_device_number = CUDA_device_number;
    function.print_Cuda_device(this->CUDA_device_number, this->tot_Blocks, this->tot_ThreadsperBlock);

    this->gpu_Limit = gpu_Limit;

    cout << "Per round GPU max unit: " << this->gpu_Limit << endl
         << endl;

    this->at_a_Time_cells = at_a_Time_cells;

    cout << "Number of cells simulated at a time: " << this->at_a_Time_cells << endl
         << endl;

    this->output_Folder = output_Folder + "/within_host";
    function.config_Folder(this->output_Folder, "Within host results");
    this->intermediate_Folders = intermediate_Folders;
    this->intermediate_sequence_Store = this->intermediate_Folders + "/sequences";
    function.config_Folder(this->intermediate_sequence_Store, "Intermediate sequences");
    this->intermediate_profile_Store = this->intermediate_Folders + "/sequence_profiles";
    function.config_Folder(this->intermediate_profile_Store, "Intermediate sequence profiles");
    // cout << "Output results folder: " << this->output_Folder << endl;

    if (write_Progeny_parents == "YES")
    {
        this->write_Progeny_parent = "YES";
        this->progeny_Parent_folder = this->output_Folder + "/progeny_parent_Data";
        function.config_Folder(this->progeny_Parent_folder, "Progeny parent relationship");
        this->progeny_File = progeny_Parent_folder + "/progeny_parent_relationships.csv";
        this->progeny_Recombination_File = progeny_Parent_folder + "/progeny_Recombination_details.csv";
        this->sequence_Profiles = progeny_Parent_folder + "/sequence_Profiles.csv";
        function.config_File_delete_create(this->progeny_File, "Source\tTarget\tType");
    }

    if (write_Sequences == "YES")
    {
        this->write_Sequences = "YES";
        this->progeny_sequences_Folder = this->output_Folder + "/sequences";
        function.config_Folder(this->progeny_sequences_Folder, "Progeny sequences");
    }

    this->references_Folder = reference_genome_Location;
    cout << "Reference sequences folder: " << this->references_Folder << endl;

    this->replication_profile_file = replication_profile_file;
    cout << "Replication profile: " << this->replication_profile_file << endl;

    this->sequence_profile_file = sequence_Profile_file;
    cout << "Sequence profile: " << this->sequence_profile_file << endl;
    cout << endl;
}

__global__ void cuda_convert_back_Sequence(int genome_Size, int **sequences, char *cuda_Sequence, int sequence_Index)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Size)
    {

        int base_Numeral = sequences[sequence_Index][tid];

        if (base_Numeral == 0)
        {
            cuda_Sequence[tid] = 'A';
        }
        else if (base_Numeral == 1)
        {
            cuda_Sequence[tid] = 'T';
        }
        else if (base_Numeral == 2)
        {
            cuda_Sequence[tid] = 'G';
        }
        else if (base_Numeral == 3)
        {
            cuda_Sequence[tid] = 'C';
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_fill_Master_sequences(int genome_Size, int *master_Sequence, char *cuda_reference)
{

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Size)
    {
        // A = 0
        // T = 1
        // G = 2
        // C = 3

        char base = cuda_reference[tid];

        if (base == 'A' || base == 'a')
        {
            master_Sequence[tid] = 0;
        }
        else if (base == 'T' || base == 't')
        {
            master_Sequence[tid] = 1;
        }
        else if (base == 'G' || base == 'g')
        {
            master_Sequence[tid] = 2;
        }
        else if (base == 'C' || base == 'c')
        {
            master_Sequence[tid] = 3;
        }

        tid += blockDim.x * gridDim.x;
    }
}

void within_host_test_2::ingress(float rep_time, float host_days, string mode)
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    // int size_Array = 2;
    // int *test_Sum = (int *)malloc(sizeof(int) * size_Array);

    // for (int i = 0; i < size_Array; i++)
    // {
    //     test_Sum[i] = 1;
    // }

    // int *cuda_test_Sum;
    // cudaMallocManaged(&cuda_test_Sum, size_Array * sizeof(int));

    // cudaMemcpy(cuda_test_Sum, test_Sum, size_Array * sizeof(int), cudaMemcpyHostToDevice);

    // int sum = function.sum_CUDA(cuda_test_Sum, size_Array);

    // cout << sum << endl;

    // exit(-1);

    cout << "Time in host (days)        : " << host_days << endl;
    cout << "Time per replication (days): " << rep_time << endl;
    int generations = ceil(host_days / rep_time);
    cout << "Generations (upper round)  : " << generations << endl
         << endl;

    cout << "Sequences found: " << endl;

    vector<string> reference_Sequence_list = function.get_Files(references_Folder, ".fasta");

    int parents_in_current_generation = reference_Sequence_list.size();

    cout << "\nTotal reference sequences: " << parents_in_current_generation << endl
         << endl;

    vector<string> parent_Sequence_Headers;

    // int parents_in_current_generation = total_ref_Sequences;

    int *parents = (int *)malloc(parents_in_current_generation * sizeof(int));
    // vector<string> parent_IDs;

    // for (size_t i = 0; i < parents_in_current_generation; i++)
    // {
    //     parents[i] = i;
    //     if (this->write_Progeny_parent == "YES")
    //     {
    //         parent_IDs.push_back("0_" + to_string(i));
    //     }
    // }

    string parent_Sequences_Store = this->intermediate_sequence_Store + "/generation_0";
    function.config_Folder(parent_Sequences_Store, "Generation 0 sequence store");

    for (int num_Sequences_ref = 0; num_Sequences_ref < parents_in_current_generation; num_Sequences_ref++)
    {
        int *cuda_parent_Sequences;
        parents[num_Sequences_ref] = num_Sequences_ref;
        // parent_IDs.push_back("0_" + to_string(num_Sequences_ref));

        int genome_bp;
        string header;
        string reference_Seq = "";

        cout << "Processing sequence " << num_Sequences_ref + 1 << endl;

        reference_Seq = function.read_Reference(reference_Sequence_list[num_Sequences_ref], header, genome_bp);

        // cout << reference_Seq << endl;

        parent_Sequence_Headers.push_back(header.substr(1));

        if (num_Sequences_ref == 0)
        {
            this->genome_SIZE = genome_bp;
            function.genome_Size = genome_bp;

            // cudaMallocManaged(&cuda_parent_Sequences, (genome_SIZE + 1) * reference_Sequence_list.size() * sizeof(int));
            // int **tmp = (int **)malloc(reference_Sequence_list.size() * sizeof(tmp[0]));
            // for (int i = 0; i < reference_Sequence_list.size(); i++)
            // {
            //     cudaMalloc((void **)&tmp[i], (genome_SIZE + 1) * sizeof(tmp[0][0]));
            // }
            // cudaMemcpy(cuda_parent_Sequences, tmp, reference_Sequence_list.size() * sizeof(int *), cudaMemcpyHostToDevice);

            // free(tmp);
        }
        else
        {
            if (genome_bp != genome_SIZE)
            {
                cout << "ERROR Sequence sizes do not match: " << reference_Sequence_list[num_Sequences_ref] << endl;
                exit(-1);
            }
        }

        cudaMalloc(&cuda_parent_Sequences, genome_SIZE * sizeof(int));

        char *reference_full, *cuda_reference;
        reference_full = (char *)malloc((reference_Seq.size() + 1) * sizeof(char));
        cudaMallocManaged(&cuda_reference, (reference_Seq.size() + 1) * sizeof(char));
        strcpy(reference_full, reference_Seq.c_str());
        cudaMemcpy(cuda_reference, reference_full, (reference_Seq.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

        free(reference_full);

        // cuda_fill_Master_sequences(int genome_Size, int **master_Sequence, char *cuda_reference, int ref_Num)
        cuda_fill_Master_sequences<<<tot_Blocks, tot_ThreadsperBlock>>>(this->genome_SIZE, cuda_parent_Sequences, cuda_reference);
        cudaDeviceSynchronize();

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error (Master fill sequences): %s\n", cudaGetErrorString(err));
        }

        int *parent_Sequences = (int *)malloc(genome_SIZE * sizeof(int));
        cudaMemcpy(parent_Sequences, cuda_parent_Sequences, genome_SIZE * sizeof(int), cudaMemcpyDeviceToHost);

        string n_FASTA_location = parent_Sequences_Store + "/0_" + to_string(num_Sequences_ref) + ".nfasta";

        function.config_File_delete_create(n_FASTA_location);
        fstream nfasta_Write;
        nfasta_Write.open(n_FASTA_location, ios::out);

        if (nfasta_Write.is_open())
        {
            for (int base_Write = 0; base_Write < genome_SIZE; base_Write++)
            {
                nfasta_Write << to_string(parent_Sequences[base_Write]);
            }
            nfasta_Write.close();
        }

        free(parent_Sequences);

        cudaFree(cuda_reference);
        cudaFree(cuda_parent_Sequences);
    }

    // function.find_Unique_values();

    // exit(-1);

    // Load sequence profiles and mutation profiles

    // MUTATION PROFILE LOAD IN PROGRESS
    parameter_load Parameters = parameter_load();

    vector<string> activate_List = {"\"Mutation activate\"",
                                    "\"Recombination activate\"",
                                    "\"Proof reading activate\""};

    vector<string> activate_Parameters = Parameters.get_parameters(replication_profile_file, activate_List);
    string mutation_Activate = Parameters.get_STRING(activate_Parameters[0]);
    string recombination_Activate = Parameters.get_STRING(activate_Parameters[1]);
    string proof_Reading = Parameters.get_STRING(activate_Parameters[2]);

    // int rows_seqeunce_mutation_Tracker = 0;

    transform(mutation_Activate.begin(), mutation_Activate.end(), mutation_Activate.begin(), ::toupper);
    if (mutation_Activate == "YES")
    {
        this->mutation_Activate_parent = 1;
        // rows_seqeunce_mutation_Tracker++;
    }

    transform(recombination_Activate.begin(), recombination_Activate.end(), recombination_Activate.begin(), ::toupper);
    if (recombination_Activate == "YES")
    {
        this->recombination_Activate_parent = 1;
        // rows_seqeunce_mutation_Tracker = rows_seqeunce_mutation_Tracker + 2;
    }

    transform(proof_Reading.begin(), proof_Reading.end(), proof_Reading.begin(), ::toupper);
    if (proof_Reading == "YES")
    {
        this->proof_reading_Activate_parent = 1;
    }

    sequence_Mutation_tracker = function.create_Fill_2D_array(7, this->genome_SIZE, -1);

    // for (size_t row = 0; row < 3; row++)
    // {
    //     for (size_t col = 0; col < this->genome_SIZE; col++)
    //     {
    //         cout << sequence_Mutation_tracker[row][col] << " ";
    //     }
    //     // cin.ignore();
    //     // cout << row + 1 << endl;
    // }

    // exit(-1);

    configure_Mutation_Profiles(generations);
    fitness_Profiles();
    survivability_Profiles();
    configure_Recombination_Profiles();
    configure_proof_reading_Profiles();

    // cout << sequence_Mutation_tracker[4][149] << "\n";

    // for (size_t r = 0; r < 5; r++)
    // {
    //     for (size_t i = 0; i < this->genome_SIZE; i++)
    //     {
    //         cout << sequence_Mutation_tracker[r][i] << " ";
    //     }
    //     cout << endl;
    // }

    configure_sequence_Profile(parents_in_current_generation, parent_Sequence_Headers);
    // exit(-1);
    //  for (size_t r = 0; r < 2; r++)
    //  {
    //      for (size_t c = 0; c < 10; c++)
    //      {
    //          cout << current_gen_Parent_data[r][c] << " ";
    //      }
    //      cout << endl;
    //  }

    // exit(-1);

    // Proof checking mechanism

    cout << "Processing replication phase data:" << endl;
    allocate_Phases(generations);

    cout << endl;

    cout << "Getting cell distribution data for viral units:" << endl;
    string cell_distribution_Type;
    get_Cell_distribution(cell_distribution_Type);
    if (cell_Limit != -1)
    {
        cout << "Number of cells available for infecetion is limited\n";
        cout << "Available cell distribution shape: " << total_Cells_shape << endl
             << "Available cell distribution scale: " << total_Cells_scale << endl;

        gamma_distribution<double> gamma_Cells(total_Cells_shape, total_Cells_scale);

        cell_Limit = (int)gamma_Cells(gen);
        cout << "Cells available for infection: " << cell_Limit << endl;
    }

    cout << "Distribution type: " << cell_distribution_Type << endl;

    if (cell_distribution_Type == "Gamma")
    {
        cout << "Shape: " << this->cell_shape << endl;
        cout << "Scale: " << this->cell_scale << endl;
    }
    else if (cell_distribution_Type == "Negative binomial")
    {
        cout << "R: " << this->cell_r << endl;
        cout << "Probability: " << this->cell_prob << endl;
    }

    cout << endl;

    // exit(-1);

    cout << "Getting progeny distribution data for viral units:" << endl;
    // string progeny_distribution_Type;
    get_Progeny_distribution(function);
    cout << "Distribution type: " << function.progeny_distribution_Type << endl;

    if (function.progeny_distribution_Type == "Gamma")
    {
        cout << "Shape: " << function.progeny_shape << endl;
        cout << "Scale: " << function.progeny_scale << endl;
    }
    else if (function.progeny_distribution_Type == "Negative binomial")
    {
        cout << "Mean: " << function.progeny_mean << endl;
        cout << "Dispersion: " << function.progeny_dispersion << endl;

        function.progeny_prob = function.progeny_dispersion / (function.progeny_mean + function.progeny_dispersion);
        function.progeny_r = round((function.progeny_mean * function.progeny_prob) / (1 - function.progeny_prob));

        cout << "R: " << function.progeny_r << endl;
        cout << "Probability: " << function.progeny_prob << endl;
    }

    cout << endl;

    // exit(-1);

    // rows
    cout << "Intializing within host engine\n"
         << endl;

    // DELETE
    // parents_in_current_generation = 15;

    // int **parent_child;
    //= function.create_INT_2D_arrays(sequences_in_current_generation, 2);

    // for (size_t i = 0; i < sequences_in_current_generation; i++)
    // {
    //     cout << name_ID[parent_child[i][0]] << "\t" << name_ID[parent_child[i][1]] << endl;
    // }

    string parent_Profiles_Store = this->intermediate_profile_Store + "/generation_0";
    function.config_Folder(parent_Profiles_Store, "Generation 0 profile store");

    // //! Remove
    // float **CUDA_current_gen_Parent_data = function.float_2D_Array_load_to_CUDA(this->current_gen_Parent_data, parents_in_current_generation, 1 + (3 * recombination_hotspots));
    // //! Remove
    // float *CUDA_parent_Proof_reading_probability;
    // if (proof_reading_Activate_parent != 0)
    // {
    //     CUDA_parent_Proof_reading_probability = function.copy_1D_to_CUDA_FLOAT(sequences_Proof_reading_probability, parents_in_current_generation);
    // }

    // int **CUDA_recombination_hotspots_start_stop;
    int *stride_Array;

    if (recombination_Activate_parent != 0)
    {
        function.recombination_hotspots = this->recombination_hotspots;
        function.CUDA_recombination_hotspots_start_stop = function.int_2D_Array_load_to_CUDA(recombination_hotspots_start_stop, recombination_hotspots, 2);

        float **A_0_Recombination = function.create_FLOAT_2D_arrays(rows_Fitness + rows_Prob + rows_Selectivity + rows_Survivability, 6);
        float **T_1_Recombination = function.create_FLOAT_2D_arrays(rows_Fitness + rows_Prob + rows_Selectivity + rows_Survivability, 6);
        float **G_2_Recombination = function.create_FLOAT_2D_arrays(rows_Fitness + rows_Prob + rows_Selectivity + rows_Survivability, 6);
        float **C_3_Recombination = function.create_FLOAT_2D_arrays(rows_Fitness + rows_Prob + rows_Selectivity + rows_Survivability, 6);

        stride_Array = (int *)malloc(5 * sizeof(int));

        stride_Array[0] = 0;
        stride_Array[1] = stride_Array[0] + rows_Prob;
        stride_Array[2] = stride_Array[1] + rows_Selectivity;
        stride_Array[3] = stride_Array[2] + rows_Fitness;
        stride_Array[4] = stride_Array[3] + rows_Survivability;

        function.CUDA_stride_Array = function.copy_1D_to_CUDA_INT(stride_Array, 5);

        for (int row_Prob = 0; row_Prob < rows_Prob; row_Prob++)
        {
            for (int col = 0; col < 6; col++)
            {
                A_0_Recombination[stride_Array[0] + row_Prob][col] = A_0_probability_Recombination[row_Prob][col];
                T_1_Recombination[stride_Array[0] + row_Prob][col] = T_1_probability_Recombination[row_Prob][col];
                G_2_Recombination[stride_Array[0] + row_Prob][col] = G_2_probability_Recombination[row_Prob][col];
                C_3_Recombination[stride_Array[0] + row_Prob][col] = C_3_probability_Recombination[row_Prob][col];
            }
        }

        for (int row_Selectivity = 0; row_Selectivity < rows_Selectivity; row_Selectivity++)
        {
            for (int col = 0; col < 6; col++)
            {
                A_0_Recombination[stride_Array[1] + row_Selectivity][col] = A_0_selectivity_Recombination[row_Selectivity][col];
                T_1_Recombination[stride_Array[1] + row_Selectivity][col] = T_1_selectivity_Recombination[row_Selectivity][col];
                G_2_Recombination[stride_Array[1] + row_Selectivity][col] = G_2_selectivity_Recombination[row_Selectivity][col];
                C_3_Recombination[stride_Array[1] + row_Selectivity][col] = C_3_selectivity_Recombination[row_Selectivity][col];
            }
        }

        for (int row_Fitness = 0; row_Fitness < rows_Fitness; row_Fitness++)
        {
            for (int col = 0; col < 6; col++)
            {
                A_0_Recombination[stride_Array[2] + row_Fitness][col] = A_0_fitness_Recombination[row_Fitness][col];
                T_1_Recombination[stride_Array[2] + row_Fitness][col] = T_1_fitness_Recombination[row_Fitness][col];
                G_2_Recombination[stride_Array[2] + row_Fitness][col] = G_2_fitness_Recombination[row_Fitness][col];
                C_3_Recombination[stride_Array[2] + row_Fitness][col] = C_3_fitness_Recombination[row_Fitness][col];
            }
        }

        for (int row_Survivability = 0; row_Survivability < rows_Survivability; row_Survivability++)
        {
            for (int col = 0; col < 6; col++)
            {
                A_0_Recombination[stride_Array[3] + row_Survivability][col] = A_0_survivability_Recombination[row_Survivability][col];
                T_1_Recombination[stride_Array[3] + row_Survivability][col] = T_1_survivability_Recombination[row_Survivability][col];
                G_2_Recombination[stride_Array[3] + row_Survivability][col] = G_2_survivability_Recombination[row_Survivability][col];
                C_3_Recombination[stride_Array[3] + row_Survivability][col] = C_3_survivability_Recombination[row_Survivability][col];
            }
        }

        // for (int type = 0; type < 3; type++)
        // {
        //     for (int row = stride_Array[type]; row < stride_Array[type + 1]; row++)
        //     {
        //         for (int col = 0; col < 6; col++)
        //         {
        //             cout << A_0_Recombination[row][col] << " ";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }

        function.CUDA_A_0_Recombination = function.float_2D_Array_load_to_CUDA(A_0_Recombination, (rows_Prob + rows_Fitness + rows_Selectivity + rows_Survivability), 6);
        function.CUDA_T_1_Recombination = function.float_2D_Array_load_to_CUDA(T_1_Recombination, (rows_Prob + rows_Fitness + rows_Selectivity + rows_Survivability), 6);
        function.CUDA_G_2_Recombination = function.float_2D_Array_load_to_CUDA(G_2_Recombination, (rows_Prob + rows_Fitness + rows_Selectivity + rows_Survivability), 6);
        function.CUDA_C_3_Recombination = function.float_2D_Array_load_to_CUDA(C_3_Recombination, (rows_Prob + rows_Fitness + rows_Selectivity + rows_Survivability), 6);

        free(A_0_Recombination);
        free(T_1_Recombination);
        free(G_2_Recombination);
        free(C_3_Recombination);

        free(A_0_fitness_Recombination);
        free(T_1_fitness_Recombination);
        free(G_2_fitness_Recombination);
        free(C_3_fitness_Recombination);

        free(A_0_selectivity_Recombination);
        free(T_1_selectivity_Recombination);
        free(G_2_selectivity_Recombination);
        free(C_3_selectivity_Recombination);

        free(A_0_probability_Recombination);
        free(T_1_probability_Recombination);
        free(G_2_probability_Recombination);
        free(C_3_probability_Recombination);

        free(A_0_survivability_Recombination);
        free(T_1_survivability_Recombination);
        free(G_2_survivability_Recombination);
        free(C_3_survivability_Recombination);

        free(recombination_hotspots_start_stop);
        free(stride_Array);
    }

    if (mutation_Activate_parent != 0)
    {
        function.mutation_hotspots = this->mutation_hotspots;
        function.CUDA_mutation_rates_Hotspot_generation = function.float_2D_Array_load_to_CUDA(this->mutation_rates_Hotspot_generation, this->mutation_hotspots, generations);
        function.CUDA_mutation_Regions_start_stop = function.int_2D_Array_load_to_CUDA(this->mutation_Regions_start_stop, mutation_hotspots, 2);

        free(mutation_rates_Hotspot_generation);
        free(mutation_Regions_start_stop);

        function.CUDA_A_0_mutation = function.float_2D_Array_load_to_CUDA(A_0_mutation, this->mutation_hotspots, 4);
        function.CUDA_T_1_mutation = function.float_2D_Array_load_to_CUDA(T_1_mutation, this->mutation_hotspots, 4);
        function.CUDA_G_2_mutation = function.float_2D_Array_load_to_CUDA(G_2_mutation, this->mutation_hotspots, 4);
        function.CUDA_C_3_mutation = function.float_2D_Array_load_to_CUDA(C_3_mutation, this->mutation_hotspots, 4);

        free(A_0_mutation);
        free(T_1_mutation);
        free(G_2_mutation);
        free(C_3_mutation);
    }

    if (this->fitness_points != -1)
    {
        function.CUDA_A_0_fitness = function.float_2D_Array_load_to_CUDA(A_0_fitness, this->fitness_points, 4);
        function.CUDA_T_1_fitness = function.float_2D_Array_load_to_CUDA(T_1_fitness, this->fitness_points, 4);
        function.CUDA_G_2_fitness = function.float_2D_Array_load_to_CUDA(G_2_fitness, this->fitness_points, 4);
        function.CUDA_C_3_fitness = function.float_2D_Array_load_to_CUDA(C_3_fitness, this->fitness_points, 4);

        free(A_0_fitness);
        free(T_1_fitness);
        free(G_2_fitness);
        free(C_3_fitness);
    }

    if (this->survivability_points != -1)
    {
        function.CUDA_A_0_survivability = function.float_2D_Array_load_to_CUDA(A_0_survivability, this->fitness_points, 4);
        function.CUDA_T_1_survivability = function.float_2D_Array_load_to_CUDA(T_1_survivability, this->fitness_points, 4);
        function.CUDA_G_2_survivability = function.float_2D_Array_load_to_CUDA(G_2_survivability, this->fitness_points, 4);
        function.CUDA_C_3_survivability = function.float_2D_Array_load_to_CUDA(C_3_survivability, this->fitness_points, 4);

        free(A_0_survivability);
        free(T_1_survivability);
        free(G_2_survivability);
        free(C_3_survivability);
    }

    if (this->proof_Reading_mutations != -1)
    {

        function.CUDA_A_0_probability_Proof_reading = function.float_2D_Array_load_to_CUDA(A_0_probability_Proof_reading, this->proof_Reading_mutations, 4);
        function.CUDA_T_1_probability_Proof_reading = function.float_2D_Array_load_to_CUDA(T_1_probability_Proof_reading, this->proof_Reading_mutations, 4);
        function.CUDA_G_2_probability_Proof_reading = function.float_2D_Array_load_to_CUDA(G_2_probability_Proof_reading, this->proof_Reading_mutations, 4);
        function.CUDA_C_3_probability_Proof_reading = function.float_2D_Array_load_to_CUDA(C_3_probability_Proof_reading, this->proof_Reading_mutations, 4);

        free(A_0_probability_Proof_reading);
        free(T_1_probability_Proof_reading);
        free(G_2_probability_Proof_reading);
        free(C_3_probability_Proof_reading);
    }

    if (this->write_Progeny_parent == "YES")
    {
        string header_recombination = "ID\tParent_ID";
        string header_profiles = "ID\tGeneration\tSurvivability\tFitness";

        if (proof_reading_Activate_parent != 0)
        {
            header_profiles = header_profiles + "\tProof_reading_accuracy";
        }

        int recomb_Columns = 0;

        if (recombination_hotspots != -1)
        {
            for (size_t i = 0; i < recombination_hotspots; i++)
            {
                header_recombination = header_recombination + "\t" + "Recombination_hotspot_" + to_string(i + 1);
                header_profiles = header_profiles +
                                  "\t" + "Recombination_hotspot_" + to_string(i + 1) + "_Selectivity" +
                                  "\t" + "Recombination_hotspot_" + to_string(i + 1) + "_Fitness" +
                                  "\t" + "Recombination_hotspot_" + to_string(i + 1) + "_Ratio_Progeny" +
                                  "\t" + "Recombination_hotspot_" + to_string(i + 1) + "_Survivability";
            }
            recomb_Columns = recombination_hotspots;
        }

        header_profiles = header_profiles + "\tSurvive_or_Not";

        function.config_File_delete_create(this->progeny_Recombination_File, header_recombination);
        function.config_File_delete_create(this->sequence_Profiles, header_profiles);

        fstream profile_File;
        profile_File.open(sequence_Profiles, ios::app);

        if (profile_File.is_open())
        {
            for (size_t i = 0; i < parents_in_current_generation; i++)
            {
                fstream parent_Profile_Files;
                fstream parent_Probability_Files;
                fstream survivability_Files;

                parent_Profile_Files.open(parent_Profiles_Store + "/0_" + to_string(parents[i]) + ".profile", ios::out);
                survivability_Files.open(parent_Profiles_Store + "/0_" + to_string(parents[i]) + ".profile_surv", ios::out);
                if (proof_reading_Activate_parent != 0)
                {
                    function.proof_reading_Activate_parent = 1;
                    parent_Probability_Files.open(parent_Profiles_Store + "/0_" + to_string(parents[i]) + ".profile_prob", ios::out);
                }

                parent_Profile_Files << to_string(current_gen_Parent_data[i][0]);
                if (proof_reading_Activate_parent != 0)
                {
                    parent_Probability_Files << to_string(sequences_Proof_reading_probability[i]);
                    parent_Probability_Files.close();
                }

                survivability_Files << to_string(sequences_Survivability[i][0]);
                for (int hotspot_surv = 0; hotspot_surv < recomb_Columns; hotspot_surv++)
                {
                    survivability_Files << "\t" << to_string(sequences_Survivability[i][hotspot_surv + 1]);
                }
                survivability_Files << "\n";
                survivability_Files.close();

                profile_File << "0_" << to_string(parents[i]) << "\t"
                             << "0\t" << to_string(sequences_Survivability[i][0])
                             << "\t" << to_string(current_gen_Parent_data[i][0]) << "\t" << to_string(sequences_Proof_reading_probability[i]);

                int track = 0;
                for (int hotspots = 0; hotspots < (recomb_Columns * 3); hotspots++)
                {
                    profile_File << "\t" << to_string(current_gen_Parent_data[i][hotspots + 1]);
                    if (hotspots == ((track * 3) + 2))
                    {
                        profile_File << "\t" << to_string(sequences_Survivability[i][track + 1]);
                        track++;
                    }
                    parent_Profile_Files << "\t" << to_string(current_gen_Parent_data[i][hotspots + 1]);
                }
                profile_File << "\tYes";
                profile_File << "\n";

                parent_Profile_Files.close();
            }
            profile_File.close();
        }
    }

    // exit(-1);

    function.CUDA_sequence_Mutation_tracker = function.int_2D_Array_load_to_CUDA(this->sequence_Mutation_tracker, 7, this->genome_SIZE);

    // Assign remaining variables to functions

    // for (size_t i = 0; i < generations; i++)
    // {
    //     cout << "**************************" << endl;
    //     cout << "Generation: " << i << endl;
    //     cout << "generation_modes: " << generation_modes[i] << endl;
    //     cout << "phase_Modes: " << phase_Modes[generation_modes[i]] << endl;
    //     cout << "phase_paramters: " << phase_parameters[generation_modes[i]][0] << " " << phase_parameters[generation_modes[i]][1] << " " << phase_parameters[generation_modes[i]][2] << endl;
    //     cout << "**************************" << endl;
    // }

    // exit(-1);

    string generation_Summary_Folder = this->output_Folder + "/generations_Summary";
    function.config_Folder(generation_Summary_Folder, "Generational summary data");
    string generation_Summary_File = generation_Summary_Folder + "/generation_Summary_File.csv";
    function.config_File_delete_create(generation_Summary_File, "Generation\tType\tParents_number\tProgeny_number");

    function.mode = mode;
    function.cells_of_parents = this->progeny_Parent_folder + "/cells_of_parents.csv";
    function.config_File_delete_create(function.cells_of_parents, "ID\tParent_Cell_ID");
    function.cells_of_progeny = this->progeny_Parent_folder + "/cells_of_progeny.csv";
    function.config_File_delete_create(function.cells_of_progeny, "ID\tProgeny_Cell_ID");

    for (int generation = 0; generation < generations; generation++)
    {
        // parents_in_current_generation = 50;
        int gen_Real = generation + 1;
        cout << "Processing " << gen_Real << " of " << generations << " generations" << endl
             << endl;

        // Assigning viral particles per cell

        string generation_line = to_string(generation + 1) + "\t";

        if (phase_Modes[generation_modes[generation]] == 0)
        {
            generation_line = generation_line + "Growth\t";
        }
        else if (phase_Modes[generation_modes[generation]] == 1)
        {
            generation_line = generation_line + "Stationary\t";
        }
        else
        {
            generation_line = generation_line + "Depriciation\t";
        }

        if (parents_in_current_generation > 0)
        {
            // parents_in_current_generation = 50;
            // free(parents);
            // parents = (int *)malloc(parents_in_current_generation * sizeof(int));
            // for (size_t i = 0; i < parents_in_current_generation; i++)
            // {
            //     parents[i] = i;
            // }

            cout << "Potential parents in generation " << generation + 1 << ": " << parents_in_current_generation << endl;

            generation_line = generation_line + to_string(parents_in_current_generation);

            // cout << "Calculating parental fitness" << endl;
            // float *cuda_Parental_Fitness = parent_Fitness(parents_in_current_generation, CUDA_current_gen_Parent_data);

            string progeny_Sequences_Store = this->intermediate_sequence_Store + "/generation_" + to_string(generation + 1);
            function.config_Folder(progeny_Sequences_Store, "Generation " + to_string(generation + 1) + " sequence store");
            string progeny_Profile_Store = this->intermediate_profile_Store + "/generation_" + to_string(generation + 1);
            function.config_Folder(progeny_Profile_Store, "Generation " + to_string(generation + 1) + " sequence profile store");

            // int *parents_Cells = (int *)malloc(parents_in_current_generation * sizeof(int));

            vector<int> cells_Start_Stop;

            // int cell_ID = 0;
            cells_Start_Stop.push_back(0);

            //  int max_Count = 0;

            cout << "Attaching viral unit(s) to cell(s)" << endl;

            int assigned_parents = 0;

            while (assigned_parents < parents_in_current_generation)
            {
                int num_particles;
                if (cell_distribution_Type == "Gamma")
                {
                    gamma_distribution<float> dist(cell_shape, cell_scale);
                    num_particles = (int)round(dist(gen));
                    // cout << cell_Count << endl;
                }
                else if (cell_distribution_Type == "Negative binomial")
                {
                    negative_binomial_distribution<int> dist(cell_r, cell_prob);
                    num_particles = dist(gen);
                }
                // cout << "p: " << num_particles << endl;

                int current_Count = 0;

                if (num_particles > 0)
                {
                    // unique_Cells++;
                    for (int i = 0; i < num_particles; i++)
                    {
                        if (assigned_parents >= parents_in_current_generation)
                        {
                            break;
                        }
                        else
                        {
                            // cout << "a: " << assigned_parents << endl;
                            // parents_Cells[assigned_parents] = cell_ID;
                            assigned_parents++;
                            current_Count++;
                        }
                    }

                    // if (current_Count > max_Count)
                    // {
                    //     max_Count = current_Count;
                    // }

                    // cell_ID++;
                    cells_Start_Stop.push_back(cells_Start_Stop[cells_Start_Stop.size() - 1] + current_Count);
                }

                if (cell_Limit != -1)
                {
                    if ((cells_Start_Stop.size() - 1) >= cell_Limit)
                    {
                        // cout << "cell_Limit reached: " << cell_Limit << endl;
                        // cout << "cells infected: " << cells_Start_Stop.size() - 1;
                        // exit(-1);
                        break;
                    }
                }
            }

            // scramble parents array then assign them to the cells.

            // int array_Size = sizeof(parents) / sizeof(parents[0]);

            // for (size_t i = 0; i < cells_Start_Stop.size() - 1; i++)
            // {
            //     for (int print = cells_Start_Stop[i]; print < cells_Start_Stop[i + 1]; print++)
            //     {
            //         cout << parents[print] << " ";
            //     }
            //     cout << endl;
            // }

            random_shuffle(&parents[0], &parents[parents_in_current_generation]);
            // int temp = parents[0];
            // parents[0] = parents[1];
            // parents[1] = temp;

            int num_Unique_cells = cells_Start_Stop.size() - 1;

            cout << "Number of infected cell(s): " << num_Unique_cells << endl
                 << endl;

            vector<pair<int, int>> start_Stop_Rounds;

            int full_Rounds = num_Unique_cells / this->at_a_Time_cells;
            int partial_Rounds = num_Unique_cells % this->at_a_Time_cells;

            for (int full = 0; full < full_Rounds; full++)
            {
                int start = full * this->at_a_Time_cells;
                int stop = start + this->at_a_Time_cells;
                start_Stop_Rounds.push_back(make_pair(start, stop));
            }

            if (partial_Rounds != 0)
            {
                int start = num_Unique_cells - partial_Rounds;
                start_Stop_Rounds.push_back(make_pair(start, num_Unique_cells));
            }

            int sum_Progeny_in_Generation = 0;
            int cells_Processed = 0;

            vector<int> surviving_Progeny;

            for (int cells_Process = 0; cells_Process < start_Stop_Rounds.size(); cells_Process++)
            {
                int num_of_Cells = start_Stop_Rounds[cells_Process].second - start_Stop_Rounds[cells_Process].first;
                cout << "Processing round " << cells_Process + 1 << " of " << start_Stop_Rounds.size() << ": " << num_of_Cells << " cell(s)" << endl;

                function.process_Cells(this->multi_READ, generation, sum_Progeny_in_Generation,
                                       num_of_Cells, start_Stop_Rounds[cells_Process].first, start_Stop_Rounds[cells_Process].second,
                                       cells_Start_Stop,
                                       parent_Profiles_Store, parents,
                                       parent_Sequences_Store,
                                       progeny_File, progeny_Recombination_File, sequence_Profiles,
                                       progeny_Sequences_Store, progeny_Profile_Store,
                                       cells_Processed, surviving_Progeny);
                // cout << endl;
            }

            cout << "Completed generation via " << start_Stop_Rounds.size() << " rounds" << endl;
            cout << "Simulated progeny: " << sum_Progeny_in_Generation << endl;
            cout << "Progeny that survived till parenthood: " << surviving_Progeny.size() << endl;
            cout << "Progeny that perished before becoming parents: " << sum_Progeny_in_Generation - surviving_Progeny.size() << endl;

            // exit(-1);

            generation_line = generation_line + "\t" + to_string(sum_Progeny_in_Generation) + "\n";

            cout << "Purging parent profile intermediaries" << endl;
            if (filesystem::exists(parent_Profiles_Store))
            {
                filesystem::remove_all(parent_Profiles_Store);

                // string sourceFolder = parent_Sequences_Store;
                string tar_Folder = parent_Sequences_Store + ".tar";

                string command_Tar = "tar -cf " + tar_Folder + " " + parent_Sequences_Store + " && rm -R " + parent_Sequences_Store;

                int result = system(command_Tar.c_str());

                if (result == 0)
                {
                    cout << "Tar successful" << endl;
                }
                else
                {
                    cout << "Failed to tar the sequence folder: " << parent_Sequences_Store << endl;
                    exit(-1);
                }
            }

            parent_Profiles_Store = progeny_Profile_Store;
            parent_Sequences_Store = progeny_Sequences_Store;

            cout << endl;

            free(parents);

            // if (generation == 1)
            // {
            //     cout << "Stop check" << endl;
            //     exit(-1);
            // }

            if (phase_Modes[generation_modes[generation]] == 0)
            {
                cout << "Growth phase" << endl;
                // free(parents);
                parents_in_current_generation = surviving_Progeny.size();
                parents = (int *)malloc(parents_in_current_generation * sizeof(int));

                for (int fill = 0; fill < parents_in_current_generation; fill++)
                {
                    parents[fill] = surviving_Progeny[fill];
                    //  cout << parents[fill] << endl;
                }
                // exit(-1);
            }
            else if (phase_Modes[generation_modes[generation]] == 1)
            {
                cout << "Stationary phase" << endl;
                if (phase_parameters[generation_modes[generation]][2] == -1)
                {
                    phase_parameters[generation_modes[generation]][2] = parents_in_current_generation;
                }

                normal_distribution<float> dist(phase_parameters[generation_modes[generation]][2], phase_parameters[generation_modes[generation]][0]);
                int parents_in_Next_gen = round(dist(gen));

                if (parents_in_Next_gen >= surviving_Progeny.size())
                {

                    parents_in_current_generation = surviving_Progeny.size();
                    parents = (int *)malloc(parents_in_current_generation * sizeof(int));

                    for (int fill = 0; fill < parents_in_current_generation; fill++)
                    {
                        parents[fill] = surviving_Progeny[fill];
                    }
                }
                else
                {
                    cout << "Parents qualifying to the next generation: " << parents_in_Next_gen << endl;
                    priority_queue<pair<float, int>, vector<pair<float, int>>, CustomComparator> pq;

                    // we reduce to only the fittest
                    if (multi_READ == "YES")
                    {
                        cout << "Multi read of files to get parents with highest fitness" << endl;

                        int full_Rounds = surviving_Progeny.size() / this->gpu_Limit;
                        int partial_Rounds = surviving_Progeny.size() % this->gpu_Limit;

                        vector<pair<int, int>> start_stops;

                        for (int full = 0; full < full_Rounds; full++)
                        {
                            int start = full * this->gpu_Limit;
                            int stop = start + this->gpu_Limit;
                            start_stops.push_back(make_pair(start, stop));
                        }

                        if (partial_Rounds != 0)
                        {
                            int start = surviving_Progeny.size() - partial_Rounds;
                            start_stops.push_back(make_pair(start, surviving_Progeny.size()));
                        }

                        for (size_t i = 0; i < start_stops.size(); i++)
                        {
                            vector<thread> threads_vec;
                            int num_of_values_current = start_stops[i].second - start_stops[i].first;
                            int num_per_Core = num_of_values_current / this->CPU_cores;
                            int remainder = num_of_values_current % this->CPU_cores;

                            for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
                            {
                                int start_Cell = core_ID * num_per_Core;
                                int stop_Cell = start_Cell + num_per_Core;

                                threads_vec.push_back(thread{&within_host_test_2::read_Profiles_multi_Thread, this, start_Cell, stop_Cell, generation, parent_Profiles_Store, start_stops[i].first, surviving_Progeny});
                            }

                            if (remainder != 0)
                            {
                                int start_Cell = num_of_values_current - remainder;
                                int stop_Cell = num_of_values_current;

                                threads_vec.push_back(thread{&within_host_test_2::read_Profiles_multi_Thread, this, start_Cell, stop_Cell, generation, parent_Profiles_Store, start_stops[i].first, surviving_Progeny});
                            }

                            for (thread &t : threads_vec)
                            {
                                if (t.joinable())
                                {
                                    t.join();
                                }
                            }

                            threads_vec.clear();
                            // cout << index_Fitness.size() << endl;
                            for (int i = 0; i < index_Fitness.size(); i++)
                            {
                                pq.push(index_Fitness[i]);
                                // remove dive by 2
                                if (pq.size() > parents_in_Next_gen)
                                {
                                    // If the priority queue size exceeds 10, remove the smallest value
                                    pq.pop();
                                }
                            }
                            index_Fitness.clear();
                        }
                    }
                    else
                    {
                        cout << "Single read of files to get parents with highest fitness" << endl;
                        for (int parent = 0; parent < surviving_Progeny.size(); parent++)
                        {
                            fstream parent_File;
                            parent_File.open(parent_Profiles_Store + "/" + to_string(generation + 1) + "_" + to_string(surviving_Progeny[parent]) + ".profile", ios::in);
                            if (parent_File.is_open())
                            {
                                string line;
                                getline(parent_File, line);
                                vector<string> line_Data;
                                function.split(line_Data, line, '\t');

                                vector<pair<float, int>> pairedVector;

                                float overall_Fitness = stof(line_Data[0]);

                                for (int hotspot = 0; hotspot < function.recombination_hotspots; hotspot++)
                                {
                                    overall_Fitness = overall_Fitness * stof(line_Data[(hotspot * 3) + 2]);
                                }

                                // get product of fitness
                                pairedVector.push_back(make_pair(overall_Fitness, parent));

                                pq.push(pairedVector[0]);
                                // remove dive by 2
                                if (pq.size() > parents_in_Next_gen)
                                {
                                    // If the priority queue size exceeds 10, remove the smallest value
                                    pq.pop();
                                }

                                parent_File.close();
                            }
                        }
                    }

                    cout << "Configuring fittest parents" << endl;
                    parents_in_current_generation = parents_in_Next_gen;
                    parents = (int *)malloc(parents_in_current_generation * sizeof(int));

                    int index = 0;

                    while (!pq.empty())
                    {
                        pair<float, int> value = pq.top();
                        pq.pop();
                        // cout << "Value: " << value.first << ", Parent: " << value.second << endl;
                        parents[index] = value.second;
                        index++;
                    }
                }
            }
            else
            {
                cout << "Depriciation phase" << endl;
                if (phase_parameters[generation_modes[generation]][2] == -1)
                {
                    phase_parameters[generation_modes[generation]][2] = parents_in_current_generation;
                }

                float reduction_Ratio = function.beta_Distribution(phase_parameters[generation_modes[generation]][0], phase_parameters[generation_modes[generation]][1], gen);
                cout << "Reduction ratio: " << reduction_Ratio << endl;

                float num_To_remove = phase_parameters[generation_modes[generation]][2] * reduction_Ratio;

                int parents_in_Next_gen = (int)round(phase_parameters[generation_modes[generation]][2] - num_To_remove);

                cout << "Parents moving to next generation: " << parents_in_Next_gen << endl;

                phase_parameters[generation_modes[generation]][2] = parents_in_Next_gen;

                if (parents_in_Next_gen >= surviving_Progeny.size())
                {

                    parents_in_current_generation = surviving_Progeny.size();
                    parents = (int *)malloc(parents_in_current_generation * sizeof(int));

                    for (int fill = 0; fill < parents_in_current_generation; fill++)
                    {
                        parents[fill] = surviving_Progeny[fill];
                    }
                }
                else if (parents_in_Next_gen > 0)
                {
                    // cout << "Parents qualifying to the next generation: " << parents_in_Next_gen << endl;
                    priority_queue<pair<float, int>, vector<pair<float, int>>, CustomComparator> pq;

                    if (multi_READ == "YES")
                    {
                        cout << "Multi read of files to get parents with highest fitness" << endl;

                        int full_Rounds = surviving_Progeny.size() / this->gpu_Limit;
                        int partial_Rounds = surviving_Progeny.size() % this->gpu_Limit;

                        vector<pair<int, int>> start_stops;

                        for (int full = 0; full < full_Rounds; full++)
                        {
                            int start = full * this->gpu_Limit;
                            int stop = start + this->gpu_Limit;
                            start_stops.push_back(make_pair(start, stop));
                        }

                        if (partial_Rounds != 0)
                        {
                            int start = surviving_Progeny.size() - partial_Rounds;
                            start_stops.push_back(make_pair(start, surviving_Progeny.size()));
                        }

                        for (size_t i = 0; i < start_stops.size(); i++)
                        {
                            vector<thread> threads_vec;
                            int num_of_values_current = start_stops[i].second - start_stops[i].first;
                            int num_per_Core = num_of_values_current / this->CPU_cores;
                            int remainder = num_of_values_current % this->CPU_cores;

                            for (int core_ID = 0; core_ID < this->CPU_cores; core_ID++)
                            {
                                int start_Cell = core_ID * num_per_Core;
                                int stop_Cell = start_Cell + num_per_Core;

                                threads_vec.push_back(thread{&within_host_test_2::read_Profiles_multi_Thread, this, start_Cell, stop_Cell, generation, parent_Profiles_Store, start_stops[i].first, surviving_Progeny});
                            }

                            if (remainder != 0)
                            {
                                int start_Cell = num_of_values_current - remainder;
                                int stop_Cell = num_of_values_current;

                                threads_vec.push_back(thread{&within_host_test_2::read_Profiles_multi_Thread, this, start_Cell, stop_Cell, generation, parent_Profiles_Store, start_stops[i].first, surviving_Progeny});
                            }

                            for (thread &t : threads_vec)
                            {
                                if (t.joinable())
                                {
                                    t.join();
                                }
                            }

                            threads_vec.clear();
                            // cout << index_Fitness.size() << endl;
                            for (int i = 0; i < index_Fitness.size(); i++)
                            {
                                pq.push(index_Fitness[i]);
                                // remove dive by 2
                                if (pq.size() > parents_in_Next_gen)
                                {
                                    // If the priority queue size exceeds 10, remove the smallest value
                                    pq.pop();
                                }
                            }
                            index_Fitness.clear();
                        }
                    }
                    else
                    {
                        cout << "Single read of files to get parents with highest fitness" << endl;
                        for (int parent = 0; parent < surviving_Progeny.size(); parent++)
                        {
                            fstream parent_File;
                            parent_File.open(parent_Profiles_Store + "/" + to_string(generation + 1) + "_" + to_string(surviving_Progeny[parent]) + ".profile", ios::in);
                            if (parent_File.is_open())
                            {
                                string line;
                                getline(parent_File, line);
                                vector<string> line_Data;
                                function.split(line_Data, line, '\t');

                                vector<pair<float, int>> pairedVector;

                                float overall_Fitness = stof(line_Data[0]);

                                for (int hotspot = 0; hotspot < function.recombination_hotspots; hotspot++)
                                {
                                    overall_Fitness = overall_Fitness * stof(line_Data[(hotspot * 3) + 2]);
                                }

                                // get product of fitness
                                pairedVector.push_back(make_pair(overall_Fitness, parent));

                                pq.push(pairedVector[0]);
                                // remove dive by 2
                                if (pq.size() > parents_in_Next_gen)
                                {
                                    // If the priority queue size exceeds 10, remove the smallest value
                                    pq.pop();
                                }

                                parent_File.close();
                            }
                        }
                    }
                    cout << "Configuring fittest parents" << endl;
                    parents_in_current_generation = parents_in_Next_gen;
                    parents = (int *)malloc(parents_in_current_generation * sizeof(int));

                    int index = 0;

                    while (!pq.empty())
                    {
                        pair<float, int> value = pq.top();
                        pq.pop();
                        // cout << "Value: " << value.first << ", Parent: " << value.second << endl;
                        parents[index] = value.second;
                        index++;
                    }
                }
                else
                {
                    parents_in_current_generation = 0;
                }

                // exit(-1);
            }

            // cout << "max:" << max_Count << endl;

            // assign viruses to cells

            // int **cell_and_their_viruses = function.create_INT_2D_arrays(num_Unique_cells, max_Count);
            // int *per_Cell_max_viruses = (int *)malloc(num_Unique_cells * sizeof(int));

            // int *parent_and_their_cell = (int *)malloc(parents_in_current_generation * sizeof(int));

            // int viral_Count_per_Cell = 0;

            // int cell_ID_current = parents_Cells[0];

            // for (int assign_Cells = 0; assign_Cells < parents_in_current_generation; assign_Cells++)
            // {
            //     if (cell_ID_current != parents_Cells[assign_Cells])
            //     {
            //         per_Cell_max_viruses[cell_ID_current] = viral_Count_per_Cell;
            //         cell_ID_current = parents_Cells[assign_Cells];
            //         viral_Count_per_Cell = 0;
            //     }

            //     cell_and_their_viruses[parents_Cells[assign_Cells]][viral_Count_per_Cell] = parents[assign_Cells];
            //     parent_and_their_cell[parents[assign_Cells]] = parents_Cells[assign_Cells];
            //     viral_Count_per_Cell++;
            // }

            // free(parents_Cells);

            // for (size_t i = 0; i < parents_in_current_generation; i++)
            // {
            //    cout<< parent_and_their_cell[i] << " ";
            // }
            // cout << endl;

            // per_Cell_max_viruses[cell_ID_current] = viral_Count_per_Cell;

            // for (int row = 0; row < num_Unique_cells; row++)
            // {
            //     for (int col = 0; col < per_Cell_max_viruses[row]; col++)
            //     {
            //         cout << cell_and_their_viruses[row][col] << " ";
            //     }
            //     cout << endl;
            // }

            // determine the number of progeny per viral parent
            // int **cuda_Progeny_numbers;
            //= function.progeny_distribution_CUDA(this->progeny_distribution_Type, parents_in_current_generation,
            //                                                                 this->progeny_shape, this->progeny_scale,
            //                                                                 this->progeny_mean, this->progeny_dispersion,
            //                                                                 cuda_Parental_Fitness, CUDA_current_gen_Parent_data, this->recombination_hotspots);

            // cudaFree(cuda_Parental_Fitness);

            // int sum_Progeny = function.sum_CUDA(cuda_Progeny_numbers, parents_in_current_generation);
            // cout << "Sum total progeny: " << sum_Progeny << endl;

            // int *stride_Progeny_Index = (int *)malloc((parents_in_current_generation + 1) * sizeof(int));
            // int **progeny_recom_Index_Cuda = function.create_Progeny_Array(parents_in_current_generation, stride_Progeny_Index, sum_Progeny, recombination_hotspots, cuda_Progeny_numbers, -1);
            // cout << endl;

            // cudaFree(cuda_Progeny_numbers);

            // if (recombination_hotspots != -1)
            // {
            //     int **cell_and_their_viruses_CUDA = function.int_2D_Array_load_to_CUDA(cell_and_their_viruses, num_Unique_cells, max_Count);
            //     int *parent_and_their_cell_CUDA = function.copy_1D_to_CUDA_INT(parent_and_their_cell, parents_in_current_generation);
            //     int *per_Cell_max_viruses_CUDA = function.copy_1D_to_CUDA_INT(per_Cell_max_viruses, num_Unique_cells);

            //     function.progeny_Recombination_parents_array(progeny_recom_Index_Cuda, sum_Progeny, this->recombination_hotspots,
            //                                                  parent_and_their_cell_CUDA, cell_and_their_viruses_CUDA, per_Cell_max_viruses_CUDA,
            //                                                  CUDA_current_gen_Parent_data, max_Count, num_Unique_cells);

            //     int *stride_Progeny_Index_CUDA = function.copy_1D_to_CUDA_INT(stride_Progeny_Index, parents_in_current_generation + 1);
            //     function.progeny_Shuffle(progeny_recom_Index_Cuda, this->recombination_hotspots, parents_in_current_generation, stride_Progeny_Index_CUDA);

            //     cudaFree(cell_and_their_viruses_CUDA);
            //     cudaFree(parent_and_their_cell_CUDA);
            //     cudaFree(per_Cell_max_viruses_CUDA);

            //     cout << endl;
            // }

            // // **create sequences
            // // Configure progeny profiles
            // float *CUDA_progeny_Proof_reading_probability;
            // cudaMallocManaged(&CUDA_progeny_Proof_reading_probability, sum_Progeny * sizeof(float));

            // float **CUDA_current_gen_Progeny_data = function.create_current_Progeny_data(sum_Progeny, this->recombination_hotspots, progeny_recom_Index_Cuda, CUDA_current_gen_Parent_data, CUDA_parent_Proof_reading_probability, CUDA_progeny_Proof_reading_probability);

            // // ! Change
            // vector<string> progeny_IDs;
            // function.hard_Load_progeny(generation, parent_Sequences_Store, genome_SIZE,
            //                            parent_IDs,
            //                            sum_Progeny,
            //                            progeny_recom_Index_Cuda, recombination_hotspots, CUDA_recombination_hotspots_start_stop,
            //                            multi_READ,
            //                            progeny_Sequences_Store,
            //                            this->mutation_hotspots,
            //                            CUDA_current_gen_Progeny_data, CUDA_progeny_Proof_reading_probability, CUDA_mutation_rates_Hotspot_generation, CUDA_mutation_Regions_start_stop, CUDA_sequence_Mutation_tracker, proof_reading_Activate_parent,
            //                            progeny_IDs);

            // cout << "Copying progeny parent data to host for next generation" << endl;

            // float **current_gen_Progeny_data = function.create_FLOAT_2D_arrays(sum_Progeny, (1 + (3 * this->recombination_hotspots)));
            // float *progeny_Proof_reading_probability = (float *)malloc(sum_Progeny * sizeof(float));
            // cudaMemcpy(progeny_Proof_reading_probability, CUDA_progeny_Proof_reading_probability, sum_Progeny * sizeof(float), cudaMemcpyDeviceToHost);

            // for (int progeny = 0; progeny < sum_Progeny; progeny++)
            // {
            //     // cudaMemcpy(progeny_recom_Index[progeny], progeny_recom_Index_Cuda[progeny], (columns + 1) * sizeof(progeny_recom_Index_Cuda[0][0]), cudaMemcpyDeviceToHost);
            //     cudaMemcpy(current_gen_Progeny_data[progeny], CUDA_current_gen_Progeny_data[progeny], ((1 + (3 * this->recombination_hotspots)) + 1) * sizeof(CUDA_current_gen_Progeny_data[0][0]), cudaMemcpyDeviceToHost);
            // }

            // if (this->write_Progeny_parent == "YES")
            // {
            //     cout << "\nWriting progeny parent relationship data: " << endl;
            //     int columns = 1;
            //     if (this->recombination_hotspots != -1)
            //     {
            //         columns = columns + recombination_hotspots;
            //     }

            //     int **progeny_recom_Index = function.create_INT_2D_arrays(sum_Progeny, columns);
            //     cout << "Copying progeny parent data to host" << endl;

            //     for (int progeny = 0; progeny < sum_Progeny; progeny++)
            //     {
            //         cudaMemcpy(progeny_recom_Index[progeny], progeny_recom_Index_Cuda[progeny], (columns + 1) * sizeof(progeny_recom_Index_Cuda[0][0]), cudaMemcpyDeviceToHost);
            //         // cudaMemcpy(current_gen_Progeny_data[progeny], CUDA_current_gen_Progeny_data[progeny], ((1 + (3 * this->recombination_hotspots)) + 1) * sizeof(CUDA_current_gen_Progeny_data[0][0]), cudaMemcpyDeviceToHost);
            //     }

            //     cout << "Writing progeny parent data: " << this->progeny_File << endl;

            //     fstream progeny_File_write;
            //     fstream progney_Recombination_write;
            //     fstream sequence_Profile_write;

            //     progeny_File_write.open(this->progeny_File, ios::app);
            //     progney_Recombination_write.open(this->progeny_Recombination_File, ios::app);
            //     sequence_Profile_write.open(this->sequence_Profiles, ios::app);

            //     // string sequence_File;
            //     // fstream sequence_Write;
            //     // if (this->write_Sequences == "YES")
            //     // {
            //     //     sequence_File = this->progeny_sequences_Folder + "/generation_" + to_string(generation + 1) + ".fasta";
            //     //     function.config_File_delete_create(sequence_File);
            //     //     sequence_Write.open(sequence_File, ios::app);
            //     // }

            //     if (progeny_File_write.is_open())
            //     {
            //         for (int progeny_Index = 0; progeny_Index < sum_Progeny; progeny_Index++)
            //         {

            //             progeny_IDs.push_back(to_string(generation + 1) + "_" + to_string(progeny_Index));

            //             progeny_File_write << parent_IDs[progeny_recom_Index[progeny_Index][0]] << "\t" << progeny_IDs[progeny_Index] << "\tprimary_parent"
            //                                << "\n";

            //             progney_Recombination_write << progeny_IDs[progeny_Index] << "\t" << parent_IDs[progeny_recom_Index[progeny_Index][0]];

            //             sequence_Profile_write << progeny_IDs[progeny_Index] << "\t" << to_string(generation + 1) << "\t" << to_string(current_gen_Progeny_data[progeny_Index][0]) << "\t" << to_string(progeny_Proof_reading_probability[progeny_Index]);

            //             for (int column = 1; column < columns; column++)
            //             {
            //                 int recom_Parent_Index = progeny_recom_Index[progeny_Index][column];
            //                 if (recom_Parent_Index != -1)
            //                 {
            //                     progeny_File_write << parent_IDs[recom_Parent_Index] << "\t" << progeny_IDs[progeny_Index] << "\trecombinant_" << column << "_parent"
            //                                        << "\n";
            //                     progney_Recombination_write << "\t" << parent_IDs[recom_Parent_Index];
            //                 }
            //                 else
            //                 {
            //                     progney_Recombination_write << "\t" << parent_IDs[progeny_recom_Index[progeny_Index][0]];
            //                 }

            //                 int recom_start = ((column - 1) * 3) + 1;
            //                 sequence_Profile_write << "\t" << to_string(current_gen_Progeny_data[progeny_Index][recom_start]) << "\t" << to_string(current_gen_Progeny_data[progeny_Index][recom_start + 1]) << "\t" << to_string(current_gen_Progeny_data[progeny_Index][recom_start + 2]);
            //             }

            //             sequence_Profile_write << "\n";
            //             progney_Recombination_write << "\n";

            //             // if (this->write_Sequences == "YES")
            //             // {
            //             //     sequence_Write << ">" << progeny_IDs[progeny_Index] << "\n";

            //             //     // convert cuda sequences int to char
            //             //     char *sequence, *CUDA_sequence;

            //             //     sequence = (char *)malloc((genome_SIZE + 1) * sizeof(char));
            //             //     cudaMallocManaged(&CUDA_sequence, (genome_SIZE + 1) * sizeof(char));

            //             //     // cuda_convert_back_Sequence(int genome_Size, int **sequences, char *cuda_Sequence, int sequence_Index)
            //             //     cuda_convert_back_Sequence<<<tot_Blocks, tot_ThreadsperBlock>>>(genome_SIZE, cuda_progeny_Sequences, CUDA_sequence, progeny_Index);

            //             //     cudaError_t err = cudaGetLastError();
            //             //     if (err != cudaSuccess)
            //             //     {
            //             //         printf("CUDA Error: %s\n", cudaGetErrorString(err));

            //             //         // Possibly: exit(-1) if program cannot continue....
            //             //     }
            //             //     cudaDeviceSynchronize();

            //             //     cudaMemcpy(sequence, CUDA_sequence, (genome_SIZE + 1) * sizeof(char), cudaMemcpyDeviceToHost);
            //             //     string sequence_string = sequence;
            //             //     sequence_string = sequence_string.substr(0, genome_SIZE);
            //             //     sequence_Write << sequence_string << "\n";

            //             //     sequence_Write.flush();

            //             //     free(sequence);
            //             //     cudaFree(CUDA_sequence);
            //             // }
            //         }
            //         progney_Recombination_write.close();
            //         progeny_File_write.close();
            //         sequence_Profile_write.close();

            //         // if (this->write_Sequences == "YES")
            //         // {
            //         //     sequence_Write.close();
            //         // }
            //     }

            //     free(progeny_recom_Index);

            //     free(progeny_Proof_reading_probability);
            // }

            // cudaFree(progeny_recom_Index_Cuda);

            // // RELOAD to next generation
            // if (generation_modes[generation] == 0)
            // {
            //     cout << "Growth generation" << endl
            //          << endl;
            //     parents_in_current_generation = sum_Progeny;
            //     cudaFree(CUDA_current_gen_Parent_data);

            //     parent_Sequences_Store = progeny_Sequences_Store;

            //     CUDA_current_gen_Parent_data = function.float_2D_Array_load_to_CUDA(current_gen_Progeny_data, parents_in_current_generation, 1 + (3 * this->recombination_hotspots));

            //     // cout << "Growth generation" << endl;

            //     // sum_Progeny, (1 + (3 * this->recombination_hotspots))

            //     // cudaMemcpy(CUDA_current_gen_Parent_data, CUDA_current_gen_Progeny_data, sum_Progeny * (1 + (3 * this->recombination_hotspots)) * sizeof(float), cudaMemcpyDeviceToDevice);
            //     //  CUDA_current_gen_Parent_data = CUDA_current_gen_Progeny_data;

            //     free(parents);
            //     parents = (int *)malloc(parents_in_current_generation * sizeof(int));
            //     parent_IDs.clear();
            //     parent_IDs = progeny_IDs;
            //     progeny_IDs.clear();

            //     for (int i = 0; i < parents_in_current_generation; i++)
            //     {
            //         parents[i] = i;
            //     }

            //     cudaFree(CUDA_parent_Proof_reading_probability);
            //     cudaMalloc((void **)&CUDA_parent_Proof_reading_probability, parents_in_current_generation * sizeof(float));
            //     cudaMemcpy(CUDA_parent_Proof_reading_probability, progeny_Proof_reading_probability, parents_in_current_generation * sizeof(float), cudaMemcpyHostToDevice);
            // }
            // else
            // {
            //     cout << "OK" << endl;
            //     exit(-1);
            // }

            // free(current_gen_Progeny_data);

            // cudaFree(CUDA_current_gen_Progeny_data);
            // cudaFree(CUDA_progeny_Proof_reading_probability);

            // // cudaFree(cuda_progeny_Sequences);
        }
        else
        {
            cout << "WARNING: No parents in generation to simulate." << endl;
            generation_line = generation_line + "0\t0\n";
        }

        fstream generation_File_Write;
        generation_File_Write.open(generation_Summary_File, ios::app);

        if (generation_File_Write.is_open())
        {
            generation_File_Write << generation_line;
            generation_File_Write.close();
        }

        // REMOVE
        // if (generation == 2)
        // {
        //     break;
        // }
    }

    cout << "Purging parent profile intermediaries" << endl;
    if (filesystem::exists(parent_Profiles_Store))
    {
        filesystem::remove_all(parent_Profiles_Store);

        // string sourceFolder = parent_Sequences_Store;
        string tar_Folder = parent_Sequences_Store + ".tar";

        string command_Tar = "tar -cf " + tar_Folder + " " + parent_Sequences_Store + " && rm -R " + parent_Sequences_Store;

        int result = system(command_Tar.c_str());

        if (result == 0)
        {
            cout << "Tar successful" << endl;
        }
        else
        {
            cout << "Failed to tar the sequence folder: " << parent_Sequences_Store << endl;
            exit(-1);
        }
    }

    // cudaFree(cuda_parent_Sequences);
    free(generation_modes);
    free(phase_Modes);
}

__global__ void CUDA_parental_Fitness(int num_Values, int start_Index, float **CUDA_current_gen_Parent_data, int hotspot_Count, float *cuda_Parental_Fitness)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < num_Values)
    {
        int parent_ID = tid + start_Index;

        float fitness_Parent = CUDA_current_gen_Parent_data[parent_ID][0];

        for (int hotspot = 0; hotspot < hotspot_Count; hotspot++)
        {
            fitness_Parent = fitness_Parent * CUDA_current_gen_Parent_data[parent_ID][(hotspot * 3) + 2];
        }

        cuda_Parental_Fitness[parent_ID] = fitness_Parent;

        tid += blockDim.x * gridDim.x;
    }
}

void within_host_test_2::read_Profiles_multi_Thread(int start, int stop,
                                                    int generation_Current, string parent_Profiles_Store, int start_Index, vector<int> surviving_Progeny)
{
    vector<pair<float, int>> pairedVector;

    functions_library function = functions_library();

    int recom_Hotspot = this->recombination_hotspots;
    if (recom_Hotspot == -1)
    {
        recom_Hotspot = 0;
    }

    for (int i = start; i < stop; i++)
    {
        fstream parent_File;
        parent_File.open(parent_Profiles_Store + "/" + to_string(generation_Current + 1) + "_" + to_string(surviving_Progeny[i + start_Index]) + ".profile", ios::in);
        // cout << parent_Profiles_Store + "/" + to_string(generation_Current + 1) + "_" + to_string(i + start_Index) + ".profile" << endl;

        if (parent_File.is_open())
        {
            string line;
            getline(parent_File, line);
            vector<string> line_Data;
            function.split(line_Data, line, '\t');

            // vector<pair<float, int>> pairedVector;

            float overall_Fitness = stof(line_Data[0]);

            for (int hotspot = 0; hotspot < recom_Hotspot; hotspot++)
            {
                overall_Fitness = overall_Fitness * stof(line_Data[(hotspot * 3) + 2]);
            }

            // get product of fitness
            pairedVector.push_back(make_pair(overall_Fitness, i + start_Index));

            parent_File.close();
        }
    }

    unique_lock<shared_mutex> ul(g_mutex);
    for (int i = 0; i < pairedVector.size(); i++)
    {
        index_Fitness.push_back(make_pair(pairedVector[i].first, pairedVector[i].second));
    }
}

float *within_host_test_2::parent_Fitness(int &parents_in_current_generation, float **CUDA_current_gen_Parent_data)
{
    float *cuda_Parental_Fitness;
    cudaMallocManaged(&cuda_Parental_Fitness, parents_in_current_generation * sizeof(float));

    int full_Rounds = parents_in_current_generation / this->gpu_Limit;
    int partial_Rounds = parents_in_current_generation % this->gpu_Limit;

    vector<pair<int, int>> start_stops;

    for (int full = 0; full < full_Rounds; full++)
    {
        int start = full * this->gpu_Limit;
        int stop = start + this->gpu_Limit;
        start_stops.push_back(make_pair(start, stop));
    }

    if (partial_Rounds != 0)
    {
        int start = parents_in_current_generation - partial_Rounds;
        start_stops.push_back(make_pair(start, parents_in_current_generation));
    }

    for (size_t i = 0; i < start_stops.size(); i++)
    {
        int num_of_values_current = start_stops[i].second - start_stops[i].first;
        // CUDA_parental_Fitness_(int num_Values, int start_Index, float **CUDA_current_gen_Parent_data, int hotspot_Count, float *cuda_Parental_Fitness)

        CUDA_parental_Fitness<<<tot_Blocks, tot_ThreadsperBlock>>>(num_of_values_current, start_stops[i].first, CUDA_current_gen_Parent_data, recombination_hotspots, cuda_Parental_Fitness);

        cudaError_t err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            printf("CUDA Error: %s\n", cudaGetErrorString(err));

            // Possibly: exit(-1) if program cannot continue....
        }
        cudaDeviceSynchronize();
    }

    cout << "Completed fitness generation via " << start_stops.size() << " GPU rounds" << endl;

    return cuda_Parental_Fitness;
}

// int **within_host_test_2::standard_Progeny_Numbers(int &parents_in_current_generation, float *cuda_Parental_Fitness, functions_library &function, float **CUDA_current_gen_Parent_data)
// {
//     cout << "Determining progeny numbers" << endl;

//     int rows = this->recombination_hotspots;
//     if (rows == -1)
//     {
//         rows = 0;
//     }

//     // cudaMallocManaged(&cuda_Progeny_numbers, parents_in_current_generation * sizeof(int));

//     // progeny_distribution_CUDA(string &distribution_Type, int &num_of_parents, float &shape, float &scale, float &mean, float &dispersion_Parameter, int *cuda_Progeny_numbers, float *cuda__Parent_finess)
//     int **cuda_Progeny_numbers = function.progeny_distribution_CUDA(this->progeny_distribution_Type, parents_in_current_generation,
//                                                                     this->progeny_shape, this->progeny_scale,
//                                                                     this->progeny_mean, this->progeny_dispersion,
//                                                                     cuda_Parental_Fitness, CUDA_current_gen_Parent_data, rows);

//     return cuda_Progeny_numbers;
//     cudaFree(cuda_Progeny_numbers);
// }

void within_host_test_2::configure_proof_reading_Profiles()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);
    parameter_load Parameters = parameter_load();

    if (this->proof_reading_Activate_parent == 1 && this->mutation_Activate_parent == 1)
    {
        cout << "Configuring Proof reading profiles\n"
             << endl;

        fstream replication_Profile;
        replication_Profile.open(this->replication_profile_file, ios::in);
        if (replication_Profile.is_open())
        {
            string line;
            vector<string> line_Data;

            int proof_BLOCK = 0;

            while (getline(replication_Profile, line))
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

                        function.split(line_Data, remove, ':');

                        if (line_Data[0] == "\"Proof Reading\"")
                        {
                            proof_BLOCK = 1;
                            break;
                        }
                    }
                }
            }

            if (proof_BLOCK == 1)
            {
                cout << "Proof reading mutational effects" << endl;

                // int mutation_Number = 0;

                while (getline(replication_Profile, line))
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

                            function.split(line_Data, remove, ':');

                            if (line_Data[0] == "}")
                            {
                                break;
                            }
                            else
                            {
                                if (line_Data[0] == "\"Mutation effect types\"")
                                {
                                    proof_Reading_mutations = Parameters.get_INT(line_Data[1]);
                                    break;
                                }
                            }
                        }
                    }
                }

                if (this->proof_Reading_mutations != -1)
                {
                    if (this->proof_Reading_mutations == 0)
                    {
                        cout << "Since mutation effects is 0, deactivating Proof Reading mutations." << endl;
                        this->proof_Reading_mutations = -1;
                    }
                    else
                    {
                        int mutational_effects_Block = -1;
                        while (getline(replication_Profile, line))
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

                                    function.split(line_Data, remove, ':');

                                    if (line_Data[0] == "}")
                                    {
                                        break;
                                    }
                                    else
                                    {

                                        if (line_Data[0] == "\"Mutation effects\"")
                                        {
                                            mutational_effects_Block = 1;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if (mutational_effects_Block == 1)
                        {
                            cout << "Number of proof reading mutations: " << this->proof_Reading_mutations << endl;
                            string keyword = "Effect ";

                            A_0_probability_Proof_reading = function.create_Fill_2D_array_FLOAT(proof_Reading_mutations, 4, 0);
                            T_1_probability_Proof_reading = function.create_Fill_2D_array_FLOAT(proof_Reading_mutations, 4, 0);
                            G_2_probability_Proof_reading = function.create_Fill_2D_array_FLOAT(proof_Reading_mutations, 4, 0);
                            C_3_probability_Proof_reading = function.create_Fill_2D_array_FLOAT(proof_Reading_mutations, 4, 0);

                            this->position_Proof_reading_mutations = (int *)malloc(proof_Reading_mutations * sizeof(int));

                            for (int effect = 0; effect < proof_Reading_mutations; effect++)
                            {
                                string check_Key = keyword + to_string(effect + 1);

                                int effects_Block = -1;

                                while (getline(replication_Profile, line))
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

                                            function.split(line_Data, remove, ':');

                                            if (line_Data[0] == "}")
                                            {
                                                break;
                                            }
                                            else
                                            {

                                                if (Parameters.get_STRING(line_Data[0]) == check_Key)
                                                {
                                                    cout << "Configuring " << check_Key << endl;
                                                    effects_Block = 1;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (effects_Block == 1)
                                {
                                    while (getline(replication_Profile, line))
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

                                                function.split(line_Data, remove, ':');

                                                if (line_Data[0] == "}")
                                                {
                                                    break;
                                                }
                                                else
                                                {
                                                    if (line_Data[0] == "\"Position\"")
                                                    {
                                                        this->position_Proof_reading_mutations[effect] = Parameters.get_INT(line_Data[1]);
                                                        this->sequence_Mutation_tracker[4][position_Proof_reading_mutations[effect] - 1] = effect;
                                                        cout << "Position: " << position_Proof_reading_mutations[effect] << endl;
                                                        continue;
                                                    }
                                                    else
                                                    {
                                                        char base;
                                                        int mutation;
                                                        function.get_base_mutation(Parameters.get_STRING(line_Data[0]), base, mutation);

                                                        // cout << line_Data[0] << endl;

                                                        if (base == 'A')
                                                        {
                                                            A_0_probability_Proof_reading[effect][mutation] = Parameters.get_FLOAT(line_Data[1]);
                                                        }
                                                        else if (base == 'T')
                                                        {
                                                            T_1_probability_Proof_reading[effect][mutation] = Parameters.get_FLOAT(line_Data[1]);
                                                        }
                                                        else if (base == 'G')
                                                        {
                                                            G_2_probability_Proof_reading[effect][mutation] = Parameters.get_FLOAT(line_Data[1]);
                                                        }
                                                        else if (base == 'C')
                                                        {
                                                            C_3_probability_Proof_reading[effect][mutation] = Parameters.get_FLOAT(line_Data[1]);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    cout << "ERROR: " << check_Key << " MISSING." << endl;
                                    exit(-1);
                                }
                            }
                            // for (size_t i = 0; i < 4; i++)
                            // {
                            //     cout << T_1_probability_Proof_reading[0][i] << endl;
                            // }
                        }
                        else
                        {
                            cout << "ERROR: \"Mutation effects\" PARAMETER WAS NOT FOUND." << endl;
                            exit(-1);
                        }
                    }
                }
            }
            else
            {
                cout << "ERROR PROOF READING MUTATION BLOCK ABSENT. EVEN if NO MUTATIONAL EFFECTS IT SHOULD BE SET TO ZERO." << endl;
                exit(-1);
            }

            replication_Profile.close();
        }
    }
    else
    {
        cout << "Proof reading mechanism is deactivated" << endl;
    }
}

void within_host_test_2::configure_sequence_Profile(int parents_in_current_generation, vector<string> parent_headers)
{
    functions_library function = functions_library();
    parameter_load Parameters = parameter_load();

    cout << "Configuring parent sequences" << endl
         << endl;

    current_gen_Parent_data = function.create_FLOAT_2D_arrays(parents_in_current_generation, 1 + (3 * recombination_hotspots));

    if (proof_reading_Activate_parent != 0)
    {
        sequences_Proof_reading_probability = (float *)malloc(parents_in_current_generation * sizeof(float));
    }

    // survivability = (float *)malloc(parents_in_current_generation * sizeof(float));
    int cols_survivability;
    if (recombination_hotspots == -1)
    {
        cols_survivability = 0;
    }
    else
    {
        cols_survivability = recombination_hotspots;
    }
    sequences_Survivability = function.create_FLOAT_2D_arrays(parents_in_current_generation, cols_survivability + 1);

    fstream sequence_Profile;
    sequence_Profile.open(this->sequence_profile_file, ios::in);

    if (sequence_Profile.is_open())
    {
        string line;
        vector<string> line_Data;

        int num_Sequences_Present = -1;

        while (getline(sequence_Profile, line))
        {
            // cout << line << endl;
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

                    if (line_Data[0] == "}")
                    {
                        break;
                    }
                    else
                    {
                        if (line_Data[0] == "\"Number of sequences\"")
                        {
                            num_Sequences_Present = Parameters.get_INT(line_Data[1]);
                            break;
                        }
                    }
                }
            }
        }
        if (num_Sequences_Present == parents_in_current_generation)
        {
            cout << "Processing " << num_Sequences_Present << " sequences\n"
                 << endl;
            string keyword = "Sequence ";

            for (int sequence = 0; sequence < num_Sequences_Present; sequence++)
            {
                string check_Sequence = keyword + to_string(sequence + 1);
                // cout << "Configuring sequences" << endl;
                string sequence_Name;

                int sequence_Block = -1;
                while (getline(sequence_Profile, line))
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

                            // if (line_Data[0] == "}")
                            // {
                            //     break;
                            // }
                            // else
                            //{
                            if (Parameters.get_STRING(line_Data[0]) == check_Sequence)
                            {
                                sequence_Block = 1;
                                break;
                            }
                            //}
                        }
                    }
                }
                if (sequence_Block != -1)
                {
                    cout << "Configuring " << check_Sequence << endl;
                    int seq_ID_check, fitness_check, proof_check = -1;
                    while (getline(sequence_Profile, line))
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

                                if (line_Data[0] == "}")
                                {
                                    break;
                                }
                                else
                                {
                                    if (line_Data[0] == "\"Sequence ID\"")
                                    {
                                        sequence_Name = Parameters.get_STRING(line_Data[1]);
                                        // cout << sequence_Name << endl;
                                        for (int headers = 0; headers < parent_headers.size(); headers++)
                                        {
                                            // cout << parent_headers[headers] << endl;
                                            if (parent_headers[headers] == sequence_Name)
                                            {
                                                // cout << "Found" << endl;
                                                cout << "Sequence ID: " << sequence_Name << endl;
                                                seq_ID_check = headers;
                                                break;
                                            }
                                        }
                                        if (seq_ID_check == -1)
                                        {
                                            cout << "ERROR: " << sequence_Name << " SEQUENCE NOT FOUND." << endl;
                                            exit(-1);
                                        }
                                    }
                                    else if (line_Data[0] == "\"Fitness\"")
                                    {
                                        current_gen_Parent_data[seq_ID_check][0] = Parameters.get_FLOAT(line_Data[1]);
                                        cout << "Fitness: " << current_gen_Parent_data[seq_ID_check][0] << endl;
                                        fitness_check = 1;
                                    }
                                    else if (line_Data[0] == "\"Survivability\"")
                                    {
                                        sequences_Survivability[seq_ID_check][0] = Parameters.get_FLOAT(line_Data[1]);
                                        cout << "Survivability probability: " << sequences_Survivability[seq_ID_check][0] << endl;
                                    }
                                    else
                                    {
                                        if (proof_reading_Activate_parent != 0)
                                        {
                                            if (line_Data[0] == "\"Proof reading accuracy\"")
                                            {
                                                sequences_Proof_reading_probability[seq_ID_check] = Parameters.get_FLOAT(line_Data[1]);
                                                cout << "Proof reading accuracy: " << sequences_Proof_reading_probability[seq_ID_check] << endl;
                                                proof_check = 1;
                                            }
                                        }
                                    }
                                    if (proof_reading_Activate_parent != 0)
                                    {
                                        if ((seq_ID_check != -1) && (fitness_check != -1) && (proof_check != -1))
                                        {
                                            break;
                                        }
                                    }
                                    else if ((seq_ID_check != -1) && (fitness_check != -1))
                                    {
                                        proof_check = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if (seq_ID_check != -1 && fitness_check != -1 && proof_check != -1)
                    {
                        int recombination_Block = -1;
                        while (getline(sequence_Profile, line))
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

                                    if (line_Data[0] == "}")
                                    {
                                        break;
                                    }
                                    else
                                    {
                                        if (line_Data[0] == "\"Recombination hotspots\"")
                                        {
                                            recombination_Block = 1;
                                            cout << "Collecting recombination data" << endl;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if (recombination_Block != -1)
                        {
                            int hotspot_Number = -1;
                            while (getline(sequence_Profile, line))
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

                                        if (line_Data[0] == "}")
                                        {
                                            break;
                                        }
                                        else
                                        {
                                            if (line_Data[0] == "\"Number of hotspots\"")
                                            {
                                                hotspot_Number = Parameters.get_INT(line_Data[1]);
                                                cout << "Number of recombination hotspots: " << hotspot_Number << endl;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            if (hotspot_Number == this->recombination_hotspots)
                            {
                                string hotspot_Keyword = "Hotspot ";

                                for (int recom_Hotspot = 0; recom_Hotspot < hotspot_Number; recom_Hotspot++)
                                {
                                    string check_hotspot_Recomb = hotspot_Keyword + to_string(recom_Hotspot + 1);

                                    int hotspot_Block = -1;

                                    while (getline(sequence_Profile, line))
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

                                                if (line_Data[0] == "}")
                                                {
                                                    break;
                                                }
                                                else
                                                {
                                                    if (Parameters.get_STRING(line_Data[0]) == check_hotspot_Recomb)
                                                    {
                                                        cout << "Collecting data from " << check_hotspot_Recomb << endl;
                                                        hotspot_Block = 1;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (hotspot_Block != -1)
                                    {
                                        int column_Selectivity = (recom_Hotspot * 3) + 1;
                                        int column_Fitness = column_Selectivity + 1;
                                        int column_Progeny_ratio = column_Fitness + 1;

                                        while (getline(sequence_Profile, line))
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

                                                    if (line_Data[0] == "}")
                                                    {
                                                        break;
                                                    }
                                                    else
                                                    {
                                                        if (line_Data[0] == "\"Ratio of recombinant progeny\"")
                                                        {
                                                            current_gen_Parent_data[seq_ID_check][column_Progeny_ratio] = Parameters.get_FLOAT(line_Data[1]);
                                                            cout << "Ratio of recombinant progeny: " << current_gen_Parent_data[seq_ID_check][column_Progeny_ratio] << endl;
                                                        }
                                                        else if (line_Data[0] == "\"Selectivity\"")
                                                        {
                                                            current_gen_Parent_data[seq_ID_check][column_Selectivity] = Parameters.get_FLOAT(line_Data[1]);
                                                            cout << "Selectivity: " << current_gen_Parent_data[seq_ID_check][column_Selectivity] << endl;
                                                        }
                                                        else if (line_Data[0] == "\"Fitness\"")
                                                        {
                                                            current_gen_Parent_data[seq_ID_check][column_Fitness] = Parameters.get_FLOAT(line_Data[1]);
                                                            cout << "Fitness: " << current_gen_Parent_data[seq_ID_check][column_Fitness] << endl;
                                                        }
                                                        else if (line_Data[0] == "\"Survivability\"")
                                                        {
                                                            sequences_Survivability[seq_ID_check][recom_Hotspot + 1] = Parameters.get_FLOAT(line_Data[1]);
                                                            cout << "Survivability" << sequences_Survivability[seq_ID_check][recom_Hotspot + 1] << endl;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        cout << "ERROR " << check_hotspot_Recomb << "MISSING FROM CHROMOSOME " << sequence << endl;
                                        exit(-1);
                                    }
                                }
                            }
                            else
                            {
                                cout << "NUMBER OF HOTSPOTS OF " << hotspot_Number << " DOES NOT MATCH THAT OF THE REPLICATION PROFILE " << this->recombination_hotspots << endl;
                                exit(-1);
                            }
                        }
                    }
                    else
                    {
                        cout << "ERROR IN SEQUENCE ID OR FITNESS VALUES FOR SEQUENCE " << sequence + 1 << endl;
                        exit(-1);
                    }
                }
                else
                {
                    cout << "ERROR SEQUENCE " << sequence + 1 << " MISSING IN THE SEQUENCE PROFILE FILE: " << sequence_profile_file << endl;
                    exit(-1);
                }
                cout << endl;
            }
        }
        else
        {
            cout << "ERROR: NUMBER OF PARENT SEQUENCE(S) IN THE FOLDER ARE "
                 << parents_in_current_generation
                 << ",\nIT DOES NOT MATCH THE NUMBER OF SEQUENCES IN THE FILE WHICH IS: " << num_Sequences_Present << endl;
            exit(-1);
        }
        sequence_Profile.close();
    }
}

void within_host_test_2::configure_Recombination_Profiles()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);
    parameter_load Parameters = parameter_load();

    if (this->recombination_Activate_parent == 1)
    {
        cout << "Configuring recombination profiles\n"
             << endl;

        fstream replication_Profile;
        replication_Profile.open(this->replication_profile_file, ios::in);
        if (replication_Profile.is_open())
        {
            string line;
            // getline(replication_Profile, line);

            vector<string> line_Data;

            int recombination_BLOCK = 0;

            while (getline(replication_Profile, line))
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

                        if (line_Data[0] == "\"Recombination\"")
                        {
                            recombination_BLOCK = 1;
                            continue;
                        }

                        if (recombination_BLOCK == 1)
                        {
                            if (line_Data[0] == "}")
                            {
                                break;
                            }
                            else
                            {
                                if (line_Data[0] == "\"Number of recombination hotspots\"")
                                {
                                    this->recombination_hotspots = Parameters.get_INT(line_Data[1]);
                                    cout << "Number of recombination hotspots: " << this->recombination_hotspots << endl
                                         << endl;
                                }

                                if (this->recombination_hotspots != -1)
                                {
                                    string hotspot_keyword = "Hotspot ";

                                    // this->recombination_probability = (float *)malloc(recombination_hotspots * sizeof(float));
                                    this->recombination_hotspots_start_stop = function.create_INT_2D_arrays(recombination_hotspots, 2);
                                    this->recomb_effects_prob_selectivity_fitness = function.create_Fill_2D_array(recombination_hotspots, 4, 0);

                                    for (int hotspots = 0; hotspots < recombination_hotspots; hotspots++)
                                    {
                                        vector<string> effect_type;
                                        vector<int> effect_positions;
                                        vector<vector<string>> effect_base_mutations;

                                        string check_hotspot = hotspot_keyword + to_string(hotspots + 1);
                                        int hotspot_block = 0;
                                        while (getline(replication_Profile, line))
                                        {
                                            if (line != "}" && line != "")
                                            {
                                                remove = line;
                                                i = 0;
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
                                                    // cout << remove << endl;
                                                    if (Parameters.get_STRING(line_Data[0]) == check_hotspot)
                                                    {
                                                        hotspot_block = 1;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        if (hotspot_block == 1)
                                        {
                                            cout << "Processing " << check_hotspot << endl;

                                            int data_1, data_3, data_4 = 0;
                                            int mutation_effect_Count = 0;

                                            while (getline(replication_Profile, line))
                                            {
                                                if (line != "}" && line != "")
                                                {
                                                    remove = line;
                                                    i = 0;
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
                                                        function.split(line_Data, remove, ':');

                                                        if (line_Data[0] == "\"Region\"")
                                                        {
                                                            recombination_hotspots_start_stop[hotspots][0] = stoi(line_Data[1].substr(1, line_Data[1].size()));
                                                            recombination_hotspots_start_stop[hotspots][1] = stoi(line_Data[2].substr(0, line_Data[2].size() - 1));
                                                            cout << "Target region: From " << recombination_hotspots_start_stop[hotspots][0] << " to " << recombination_hotspots_start_stop[hotspots][1] << endl;
                                                            data_1++;
                                                            data_4++;
                                                        }
                                                        // else if (line_Data[0] == "\"Probability of recombination\"")
                                                        // {
                                                        //     recombination_probability[hotspots] = Parameters.get_FLOAT(line_Data[1]);
                                                        //     cout << "Probability of recombination: " << recombination_probability[hotspots] << endl;
                                                        //     data_2++;
                                                        //     data_4++;
                                                        // }
                                                        else if (line_Data[0] == "\"Mutation effect types\"")
                                                        {
                                                            // cout << remove << endl;
                                                            mutation_effect_Count = Parameters.get_INT(line_Data[1]);
                                                            data_3++;
                                                            data_4++;
                                                        }
                                                    }
                                                }
                                                if (data_4 == 2)
                                                {
                                                    break;
                                                }
                                            }
                                            if (data_4 != 2)
                                            {
                                                cout << "ERROR IN " << check_hotspot << " DATA MISSING: \n";
                                                // cout << data_1 + data_2 + data_3 << endl;
                                                if (data_1 == 0)
                                                {
                                                    cout << "Check if \"Region\" parameter is present" << endl;
                                                }
                                                // if (data_2 == 0)
                                                // {
                                                //     cout << "Check if \"Probability of recombination\" parameter is present" << endl;
                                                // }
                                                if (data_3 == 0)
                                                {
                                                    cout << "Check if \"Mutation effect types\" parameter is present" << endl;
                                                }
                                                exit(-1);
                                            }
                                            if (mutation_effect_Count != 0)
                                            {
                                                int mutation_effect_Block = 0;
                                                cout << mutation_effect_Count << " mutation effect(s) present" << endl;
                                                while (getline(replication_Profile, line))
                                                {
                                                    if (line != "}" && line != "")
                                                    {
                                                        remove = line;
                                                        i = 0;
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
                                                            function.split(line_Data, remove, ':');

                                                            if (line_Data[0] == "\"Mutation effects\"")
                                                            {
                                                                mutation_effect_Block = 1;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                                if (mutation_effect_Block == 1)
                                                {
                                                    string keyword_Effect = "Effect ";

                                                    for (int effect = 0; effect < mutation_effect_Count; effect++)
                                                    {
                                                        int effect_detail_block = 0;
                                                        string check_Effect = keyword_Effect + to_string(effect + 1);

                                                        while (getline(replication_Profile, line))
                                                        {
                                                            if (line != "}" && line != "")
                                                            {
                                                                remove = line;
                                                                i = 0;
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
                                                                    function.split(line_Data, remove, ':');

                                                                    if (Parameters.get_STRING(line_Data[0]) == check_Effect)
                                                                    {
                                                                        effect_detail_block = 1;
                                                                        break;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        if (effect_detail_block == 1)
                                                        {
                                                            cout << "\nConfiguring " << check_Effect << endl;
                                                            vector<string> base_mutations;

                                                            while (getline(replication_Profile, line))
                                                            {
                                                                if (line != "}" && line != "")
                                                                {
                                                                    remove = line;
                                                                    i = 0;
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
                                                                        function.split(line_Data, remove, ':');

                                                                        if (line_Data[0] == "}")
                                                                        {
                                                                            break;
                                                                        }
                                                                        else
                                                                        {
                                                                            // vector<string> effect_type;
                                                                            // vector<int> effect_positions;
                                                                            // vector<vector<string>> effect_base_mutations;
                                                                            if (line_Data[0] == "\"Position\"")
                                                                            {
                                                                                int position_Effect = Parameters.get_INT(line_Data[1]);
                                                                                effect_positions.push_back(position_Effect);
                                                                                cout << "Position: " << position_Effect << endl;
                                                                            }
                                                                            else if (line_Data[0] == "\"Effect type\"")
                                                                            {
                                                                                string effect_Type = Parameters.get_STRING(line_Data[1]);
                                                                                effect_type.push_back(effect_Type);
                                                                                cout << "Effect type: " << effect_Type << endl;
                                                                            }
                                                                            else
                                                                            {
                                                                                base_mutations.push_back(remove);
                                                                            }
                                                                            // cout << remove << endl;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            effect_base_mutations.push_back(base_mutations);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        effect_type_All_hotspots.push_back(effect_type);
                                        effect_positions_All_hotspots.push_back(effect_positions);
                                        effect_base_mutations_All_hotspots.push_back(effect_base_mutations);
                                        cout << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            replication_Profile.close();
        }

        cout << "Assembly of recombination effects into arrays:\n"
             << endl;

        int recomb_positions = effect_positions_All_hotspots.size();
        // cout << recomb_positions << endl;

        vector<pair<int, int>> prob_hotspot_effect;
        vector<pair<int, int>> selectivity_hotspot_effect;
        vector<pair<int, int>> fitness_hotspot_effect;
        vector<pair<int, int>> survivability_hotspot_effect;

        for (int hotspot = 0; hotspot < recomb_positions; hotspot++)
        {
            cout << "Hotspot " << hotspot + 1;

            vector<string> effect_types = effect_type_All_hotspots[hotspot];
            vector<int> effect_positions = effect_positions_All_hotspots[hotspot];

            if (effect_types.size() > 0)
            {
                cout << ": Recombination hotspot effect types" << endl;
            }
            else
            {
                cout << ": No mutation effects" << endl;
            }

            for (int effects = 0; effects < effect_types.size(); effects++)
            {
                if (effect_types[effects] == "Probability")
                {
                    cout << effect_types[effects];
                    recomb_effects_prob_selectivity_fitness[hotspot][0] = recomb_effects_prob_selectivity_fitness[hotspot][0] + 1;
                    if (sequence_Mutation_tracker[1][effect_positions[effects] - 1] == -1)
                    {
                        sequence_Mutation_tracker[1][effect_positions[effects] - 1] = 0;
                    }
                    sequence_Mutation_tracker[1][effect_positions[effects] - 1]++;
                    prob_hotspot_effect.push_back(make_pair(hotspot, effects));
                    rows_Prob++;
                }
                else if (effect_types[effects] == "Selectivity")
                {
                    cout << effect_types[effects];
                    recomb_effects_prob_selectivity_fitness[hotspot][1] = recomb_effects_prob_selectivity_fitness[hotspot][1] + 1;
                    if (sequence_Mutation_tracker[2][effect_positions[effects] - 1] == -1)
                    {
                        sequence_Mutation_tracker[2][effect_positions[effects] - 1] = 0;
                    }
                    sequence_Mutation_tracker[2][effect_positions[effects] - 1]++;
                    selectivity_hotspot_effect.push_back(make_pair(hotspot, effects));
                    rows_Selectivity++;
                }
                else if (effect_types[effects] == "Fitness")
                {
                    cout << effect_types[effects];
                    recomb_effects_prob_selectivity_fitness[hotspot][2] = recomb_effects_prob_selectivity_fitness[hotspot][2] + 1;
                    if (sequence_Mutation_tracker[3][effect_positions[effects] - 1] == -1)
                    {
                        sequence_Mutation_tracker[3][effect_positions[effects] - 1] = 0;
                    }
                    sequence_Mutation_tracker[3][effect_positions[effects] - 1]++;
                    fitness_hotspot_effect.push_back(make_pair(hotspot, effects));
                    rows_Fitness++;
                }
                else if (effect_types[effects] == "Survivability")
                {
                    cout << effect_types[effects];
                    recomb_effects_prob_selectivity_fitness[hotspot][3] = recomb_effects_prob_selectivity_fitness[hotspot][3] + 1;
                    if (sequence_Mutation_tracker[6][effect_positions[effects] - 1] == -1)
                    {
                        sequence_Mutation_tracker[6][effect_positions[effects] - 1] = 0;
                    }
                    sequence_Mutation_tracker[6][effect_positions[effects] - 1]++;
                    survivability_hotspot_effect.push_back(make_pair(hotspot, effects));
                    rows_Survivability++;
                }
                cout << endl;
            }

            // cout << endl;
        }

        // cout << rows_Prob << endl;
        // cout << rows_Selectivity << endl;
        // cout << rows_Fitness << endl;

        cout << "\nConfiguration of mutation bases" << endl;

        A_0_probability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Prob, 6, 0);
        T_1_probability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Prob, 6, 0);
        G_2_probability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Prob, 6, 0);
        C_3_probability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Prob, 6, 0);

        vector<string> line_Data;

        for (size_t prob = 0; prob < prob_hotspot_effect.size(); prob++)
        {
            vector<vector<string>> hotspot = effect_base_mutations_All_hotspots[prob_hotspot_effect[prob].first];
            vector<string> base_mutations = hotspot[prob_hotspot_effect[prob].second];

            A_0_probability_Recombination[prob][0] = (float)prob_hotspot_effect[prob].first;
            A_0_probability_Recombination[prob][1] = (float)effect_positions_All_hotspots[prob_hotspot_effect[prob].first][prob_hotspot_effect[prob].second];

            T_1_probability_Recombination[prob][0] = (float)prob_hotspot_effect[prob].first;
            T_1_probability_Recombination[prob][1] = (float)effect_positions_All_hotspots[prob_hotspot_effect[prob].first][prob_hotspot_effect[prob].second];

            G_2_probability_Recombination[prob][0] = (float)prob_hotspot_effect[prob].first;
            G_2_probability_Recombination[prob][1] = (float)effect_positions_All_hotspots[prob_hotspot_effect[prob].first][prob_hotspot_effect[prob].second];

            C_3_probability_Recombination[prob][0] = (float)prob_hotspot_effect[prob].first;
            C_3_probability_Recombination[prob][1] = (float)effect_positions_All_hotspots[prob_hotspot_effect[prob].first][prob_hotspot_effect[prob].second];

            for (int bases = 0; bases < base_mutations.size(); bases++)
            {
                // cout << base_mutations[bases] << endl;
                function.split(line_Data, base_mutations[bases], ':');
                string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                float Base_Change_Fitness = Parameters.get_FLOAT(line_Data[1]);

                char base;
                int mutation_Index;

                function.get_base_mutation(query_Base_Change, base, mutation_Index);

                if (base == 'A')
                {
                    A_0_probability_Recombination[prob][mutation_Index + 2] = Base_Change_Fitness;
                    // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                }
                else if (base == 'T')
                {
                    T_1_probability_Recombination[prob][mutation_Index + 2] = Base_Change_Fitness;
                }
                else if (base == 'G')
                {
                    G_2_probability_Recombination[prob][mutation_Index + 2] = Base_Change_Fitness;
                }
                else if (base == 'C')
                {
                    C_3_probability_Recombination[prob][mutation_Index + 2] = Base_Change_Fitness;
                }
                else
                {
                    cout << "ERROR IN RECOMBINATION POINT HOTSPOT: " << prob_hotspot_effect[prob].first << " EFFECT: " << prob_hotspot_effect[prob].second << "\nMUTATION BASE " << query_Base_Change << endl;
                    exit(-1);
                }
            }
        }

        A_0_selectivity_Recombination = function.create_Fill_2D_array_FLOAT(rows_Selectivity, 6, 1);
        T_1_selectivity_Recombination = function.create_Fill_2D_array_FLOAT(rows_Selectivity, 6, 1);
        G_2_selectivity_Recombination = function.create_Fill_2D_array_FLOAT(rows_Selectivity, 6, 1);
        C_3_selectivity_Recombination = function.create_Fill_2D_array_FLOAT(rows_Selectivity, 6, 1);

        // vector<string> line_Data;

        for (size_t select = 0; select < selectivity_hotspot_effect.size(); select++)
        {
            vector<vector<string>> hotspot = effect_base_mutations_All_hotspots[selectivity_hotspot_effect[select].first];
            vector<string> base_mutations = hotspot[selectivity_hotspot_effect[select].second];

            A_0_selectivity_Recombination[select][0] = (float)selectivity_hotspot_effect[select].first;
            A_0_selectivity_Recombination[select][1] = (float)effect_positions_All_hotspots[selectivity_hotspot_effect[select].first][selectivity_hotspot_effect[select].second];

            T_1_selectivity_Recombination[select][0] = (float)selectivity_hotspot_effect[select].first;
            T_1_selectivity_Recombination[select][1] = (float)effect_positions_All_hotspots[selectivity_hotspot_effect[select].first][selectivity_hotspot_effect[select].second];

            G_2_selectivity_Recombination[select][0] = (float)selectivity_hotspot_effect[select].first;
            G_2_selectivity_Recombination[select][1] = (float)effect_positions_All_hotspots[selectivity_hotspot_effect[select].first][selectivity_hotspot_effect[select].second];

            C_3_selectivity_Recombination[select][0] = (float)selectivity_hotspot_effect[select].first;
            C_3_selectivity_Recombination[select][1] = (float)effect_positions_All_hotspots[selectivity_hotspot_effect[select].first][selectivity_hotspot_effect[select].second];

            for (int bases = 0; bases < base_mutations.size(); bases++)
            {
                // cout << base_mutations[bases] << endl;
                function.split(line_Data, base_mutations[bases], ':');
                string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                float Base_Change_Fitness = Parameters.get_FLOAT(line_Data[1]);

                char base;
                int mutation_Index;

                function.get_base_mutation(query_Base_Change, base, mutation_Index);

                if (base == 'A')
                {
                    A_0_selectivity_Recombination[select][mutation_Index + 2] = Base_Change_Fitness;
                    // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                }
                else if (base == 'T')
                {
                    T_1_selectivity_Recombination[select][mutation_Index + 2] = Base_Change_Fitness;
                }
                else if (base == 'G')
                {
                    G_2_selectivity_Recombination[select][mutation_Index + 2] = Base_Change_Fitness;
                }
                else if (base == 'C')
                {
                    C_3_selectivity_Recombination[select][mutation_Index + 2] = Base_Change_Fitness;
                }
                else
                {
                    cout << "ERROR IN RECOMBINATION POINT HOTSPOT: " << selectivity_hotspot_effect[select].first << " EFFECT: " << selectivity_hotspot_effect[select].second << "\nMUTATION BASE " << query_Base_Change << endl;
                    exit(-1);
                }
            }
        }

        A_0_fitness_Recombination = function.create_Fill_2D_array_FLOAT(rows_Fitness, 6, 1);
        T_1_fitness_Recombination = function.create_Fill_2D_array_FLOAT(rows_Fitness, 6, 1);
        G_2_fitness_Recombination = function.create_Fill_2D_array_FLOAT(rows_Fitness, 6, 1);
        C_3_fitness_Recombination = function.create_Fill_2D_array_FLOAT(rows_Fitness, 6, 1);

        // vector<string> line_Data;

        for (size_t fit = 0; fit < fitness_hotspot_effect.size(); fit++)
        {
            vector<vector<string>> hotspot = effect_base_mutations_All_hotspots[fitness_hotspot_effect[fit].first];
            vector<string> base_mutations = hotspot[fitness_hotspot_effect[fit].second];

            A_0_fitness_Recombination[fit][0] = (float)fitness_hotspot_effect[fit].first;
            A_0_fitness_Recombination[fit][1] = (float)effect_positions_All_hotspots[fitness_hotspot_effect[fit].first][fitness_hotspot_effect[fit].second];

            T_1_fitness_Recombination[fit][0] = (float)fitness_hotspot_effect[fit].first;
            T_1_fitness_Recombination[fit][1] = (float)effect_positions_All_hotspots[fitness_hotspot_effect[fit].first][fitness_hotspot_effect[fit].second];

            G_2_fitness_Recombination[fit][0] = (float)fitness_hotspot_effect[fit].first;
            G_2_fitness_Recombination[fit][1] = (float)effect_positions_All_hotspots[fitness_hotspot_effect[fit].first][fitness_hotspot_effect[fit].second];

            C_3_fitness_Recombination[fit][0] = (float)fitness_hotspot_effect[fit].first;
            C_3_fitness_Recombination[fit][1] = (float)effect_positions_All_hotspots[fitness_hotspot_effect[fit].first][fitness_hotspot_effect[fit].second];

            for (int bases = 0; bases < base_mutations.size(); bases++)
            {
                // cout << base_mutations[bases] << endl;
                function.split(line_Data, base_mutations[bases], ':');
                string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                float Base_Change_Fitness = Parameters.get_FLOAT(line_Data[1]);

                char base;
                int mutation_Index;

                function.get_base_mutation(query_Base_Change, base, mutation_Index);

                if (base == 'A')
                {
                    A_0_fitness_Recombination[fit][mutation_Index + 2] = Base_Change_Fitness;
                    // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                }
                else if (base == 'T')
                {
                    T_1_fitness_Recombination[fit][mutation_Index + 2] = Base_Change_Fitness;
                }
                else if (base == 'G')
                {
                    G_2_fitness_Recombination[fit][mutation_Index + 2] = Base_Change_Fitness;
                }
                else if (base == 'C')
                {
                    C_3_fitness_Recombination[fit][mutation_Index + 2] = Base_Change_Fitness;
                }
                else
                {
                    cout << "ERROR IN RECOMBINATION POINT HOTSPOT: " << fitness_hotspot_effect[fit].first << " EFFECT: " << fitness_hotspot_effect[fit].second << "\nMUTATION BASE " << query_Base_Change << endl;
                    exit(-1);
                }
            }
        }

        A_0_survivability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Survivability, 6, 1);
        T_1_survivability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Survivability, 6, 1);
        G_2_survivability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Survivability, 6, 1);
        C_3_survivability_Recombination = function.create_Fill_2D_array_FLOAT(rows_Survivability, 6, 1);

        for (size_t surv = 0; surv < survivability_hotspot_effect.size(); surv++)
        {
            vector<vector<string>> hotspot = effect_base_mutations_All_hotspots[survivability_hotspot_effect[surv].first];
            vector<string> base_mutations = hotspot[survivability_hotspot_effect[surv].second];

            A_0_survivability_Recombination[surv][0] = (float)survivability_hotspot_effect[surv].first;
            A_0_survivability_Recombination[surv][1] = (float)effect_positions_All_hotspots[survivability_hotspot_effect[surv].first][survivability_hotspot_effect[surv].second];

            T_1_survivability_Recombination[surv][0] = (float)survivability_hotspot_effect[surv].first;
            T_1_survivability_Recombination[surv][1] = (float)effect_positions_All_hotspots[survivability_hotspot_effect[surv].first][survivability_hotspot_effect[surv].second];

            G_2_survivability_Recombination[surv][0] = (float)survivability_hotspot_effect[surv].first;
            G_2_survivability_Recombination[surv][1] = (float)effect_positions_All_hotspots[survivability_hotspot_effect[surv].first][survivability_hotspot_effect[surv].second];

            C_3_survivability_Recombination[surv][0] = (float)survivability_hotspot_effect[surv].first;
            C_3_survivability_Recombination[surv][1] = (float)effect_positions_All_hotspots[survivability_hotspot_effect[surv].first][survivability_hotspot_effect[surv].second];

            for (int bases = 0; bases < base_mutations.size(); bases++)
            {
                function.split(line_Data, base_mutations[bases], ':');
                string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                float Base_Change_Survivability = Parameters.get_FLOAT(line_Data[1]);

                char base;
                int mutation_Index;

                function.get_base_mutation(query_Base_Change, base, mutation_Index);
                if (base == 'A')
                {
                    A_0_survivability_Recombination[surv][mutation_Index + 2] = Base_Change_Survivability;
                    // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                }
                else if (base == 'T')
                {
                    T_1_survivability_Recombination[surv][mutation_Index + 2] = Base_Change_Survivability;
                }
                else if (base == 'G')
                {
                    G_2_survivability_Recombination[surv][mutation_Index + 2] = Base_Change_Survivability;
                }
                else if (base == 'C')
                {
                    C_3_survivability_Recombination[surv][mutation_Index + 2] = Base_Change_Survivability;
                }
                else
                {
                    cout << "ERROR IN RECOMBINATION POINT HOTSPOT: " << survivability_hotspot_effect[surv].first << " EFFECT: " << survivability_hotspot_effect[surv].second << "\nMUTATION BASE " << query_Base_Change << endl;
                    exit(-1);
                }
            }
        }

        //  //TEST

        // for (size_t i = 0; i < 6; i++)
        // {
        //     cout << A_0_fitness_Recombination[0][i] << " ";
        // }
        // cout << endl;

        // for (size_t r = 0; r < recombination_hotspots; r++)
        // {
        //     for (size_t i = 0; i < 3; i++)
        //     {
        //         cout << recomb_effects_prob_selectivity_fitness[r][i] << " ";
        //     }
        //     cout << endl;
        // }
    }
    else
    {
        cout << "Recombinations are deactivated" << endl;
    }
    cout << endl;
}

void within_host_test_2::survivability_Profiles()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);
    if (this->mutation_Activate_parent == 1)
    {
        cout << "Configuring mutation survivability profiles\n"
             << endl;

        fstream replication_Profile;
        replication_Profile.open(this->replication_profile_file, ios::in);
        if (replication_Profile.is_open())
        {
            parameter_load Parameters = parameter_load();

            string line;
            getline(replication_Profile, line);

            vector<string> line_Data;

            int survivability_Profile_BLOCK = 0;

            while (getline(replication_Profile, line))
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

                        if (line_Data[0] == "\"Mutation survivability effects\"")
                        {
                            survivability_Profile_BLOCK = 1;
                            break;
                        }
                    }
                }
            }
            if (survivability_Profile_BLOCK == 1)
            {
                // cout << "Block" << endl;
                while (getline(replication_Profile, line))
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

                            if (line_Data[0] == "}")
                            {
                                break;
                            }
                            else
                            {
                                if (line_Data[0] == "\"Number of survivability points\"")
                                {
                                    this->survivability_points = Parameters.get_INT(line_Data[1]);
                                    cout << "Number of survivability points: " << this->survivability_points << endl
                                         << endl;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (this->survivability_points != -1)
                {
                    string survivability_keyword = "Point ";
                    survivability_point_Locations = (int *)malloc(survivability_points * sizeof(int));

                    this->A_0_survivability = function.create_FLOAT_2D_arrays(this->survivability_points, 4);
                    this->T_1_survivability = function.create_FLOAT_2D_arrays(this->survivability_points, 4);
                    this->G_2_survivability = function.create_FLOAT_2D_arrays(this->survivability_points, 4);
                    this->C_3_survivability = function.create_FLOAT_2D_arrays(this->survivability_points, 4);

                    for (int surv_Point = 0; surv_Point < survivability_points; surv_Point++)
                    {
                        A_0_survivability[surv_Point][0] = -1;
                        T_1_survivability[surv_Point][1] = -1;
                        G_2_survivability[surv_Point][2] = -1;
                        C_3_survivability[surv_Point][3] = -1;

                        string check_survivability_Point = survivability_keyword + to_string(surv_Point + 1);
                        int point_Block = -1;

                        while (getline(replication_Profile, line))
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
                                    function.split(line_Data, remove, ':');
                                    if (Parameters.get_STRING(line_Data[0]) == check_survivability_Point)
                                    {
                                        cout << "Configuring " << check_survivability_Point << endl;
                                        point_Block = 1;
                                        break;
                                    }
                                }
                            }
                        }
                        if (point_Block == 1)
                        {
                            int pos_Found = -1;
                            while (getline(replication_Profile, line))
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
                                        function.split(line_Data, remove, ':');

                                        if (line_Data[0] == "}")
                                        {
                                            break;
                                        }
                                        else
                                        {
                                            if (line_Data[0] == "\"Position\"")
                                            {
                                                survivability_point_Locations[surv_Point] = stoi(line_Data[1]);
                                                sequence_Mutation_tracker[5][fitness_point_Locations[surv_Point] - 1] = surv_Point;
                                                cout << "Target point: " << survivability_point_Locations[surv_Point] << endl;
                                                pos_Found = 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            if (pos_Found == 1)
                            {
                                cout << "Configuring survivability effects\n"
                                     << endl;
                                while (getline(replication_Profile, line))
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

                                            if (line_Data[0] == "}")
                                            {
                                                break;
                                            }
                                            else
                                            {
                                                string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                                                float Base_Change_Fitness = Parameters.get_FLOAT(line_Data[1]);

                                                char base;
                                                int mutation_Index;

                                                function.get_base_mutation(query_Base_Change, base, mutation_Index);

                                                if (base == 'A')
                                                {
                                                    A_0_survivability[surv_Point][mutation_Index] = Base_Change_Fitness;
                                                    // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                                                }
                                                else if (base == 'T')
                                                {
                                                    T_1_survivability[surv_Point][mutation_Index] = Base_Change_Fitness;
                                                }
                                                else if (base == 'G')
                                                {
                                                    G_2_survivability[surv_Point][mutation_Index] = Base_Change_Fitness;
                                                }
                                                else if (base == 'C')
                                                {
                                                    C_3_survivability[surv_Point][mutation_Index] = Base_Change_Fitness;
                                                }
                                                else
                                                {
                                                    cout << "ERROR IN SURVIVABILITY POINT " << surv_Point << " MUTATION BASE " << query_Base_Change << endl;
                                                    exit(-1);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else
                            {
                                cout << "ERROR: " << check_survivability_Point << " PARAMETER \"Position\" IS MISSING." << endl;
                                exit(-1);
                            }
                        }
                        else
                        {
                            cout << "ERROR: " << check_survivability_Point << " IS MISSING." << endl;
                            exit(-1);
                        }
                    }
                }
                else
                {
                    cout << "ERROR: PARAMETER \"Number of survivability points\" IS MISSING." << endl;
                    exit(-1);
                }
            }
            else
            {
                cout << "ERROR: SURVIVABILITY BLOCK IS MISSING" << endl;
                exit(-1);
            }

            replication_Profile.close();
        }
    }
    else
    {
        cout << "Mutations are deactivated" << endl;
    }

    // CHECK
    // for (size_t i = 0; i < 4; i++)
    // {
    //     cout << A_0_survivability[1][i] << "\t";
    // }
    // cout << endl;
}

void within_host_test_2::fitness_Profiles()
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    if (this->mutation_Activate_parent == 1)
    {
        cout << "Configuring mutation fitness profiles\n"
             << endl;

        fstream replication_Profile;
        replication_Profile.open(this->replication_profile_file, ios::in);
        if (replication_Profile.is_open())
        {
            parameter_load Parameters = parameter_load();

            string line;
            getline(replication_Profile, line);

            vector<string> line_Data;

            int fitness_Profile_BLOCK = 0;

            // int one_Check = 0;

            while (getline(replication_Profile, line))
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

                        if (line_Data[0] == "\"Mutation fitess effects\"")
                        {
                            fitness_Profile_BLOCK = 1;
                            break;
                        }
                    }
                }
            }

            if (fitness_Profile_BLOCK == 1)
            {

                while (getline(replication_Profile, line))
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

                            if (line_Data[0] == "}")
                            {
                                break;
                            }
                            else
                            {

                                if (line_Data[0] == "\"Number of fitess points\"")
                                {
                                    this->fitness_points = Parameters.get_INT(line_Data[1]);
                                    cout << "Number of fitess points: " << this->fitness_points << endl
                                         << endl;
                                    break;
                                }
                            }
                        }
                    }
                }

                if (this->fitness_points != -1)
                {
                    string fitness_keyword = "Point ";

                    fitness_point_Locations = (int *)malloc(fitness_points * sizeof(int));

                    // create a sequence length and fill. Saves search time

                    this->A_0_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);
                    this->T_1_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);
                    this->G_2_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);
                    this->C_3_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);

                    for (int fit_Point = 0; fit_Point < fitness_points; fit_Point++)
                    {
                        A_0_fitness[fit_Point][0] = -1;
                        T_1_fitness[fit_Point][1] = -1;
                        G_2_fitness[fit_Point][2] = -1;
                        C_3_fitness[fit_Point][3] = -1;

                        string check_fitness_Point = fitness_keyword + to_string(fit_Point + 1);
                        int point_Block = -1;

                        while (getline(replication_Profile, line))
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

                                    // if (line_Data[0] == "}")
                                    // {
                                    //     break;
                                    // }
                                    // else
                                    //{

                                    if (Parameters.get_STRING(line_Data[0]) == check_fitness_Point)
                                    {
                                        cout << "Configuring " << check_fitness_Point << endl;
                                        point_Block = 1;
                                        break;
                                    }
                                    //  }
                                }
                            }
                        }
                        if (point_Block == 1)
                        {
                            int pos_Found = -1;
                            while (getline(replication_Profile, line))
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

                                        if (line_Data[0] == "}")
                                        {
                                            break;
                                        }
                                        else
                                        {
                                            if (line_Data[0] == "\"Position\"")
                                            {
                                                fitness_point_Locations[fit_Point] = stoi(line_Data[1]);
                                                sequence_Mutation_tracker[0][fitness_point_Locations[fit_Point] - 1] = fit_Point;
                                                cout << "Target point: " << fitness_point_Locations[fit_Point] << endl;
                                                pos_Found = 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            if (pos_Found == 1)
                            {
                                cout << "Configuring fitness effects\n"
                                     << endl;
                                while (getline(replication_Profile, line))
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

                                            if (line_Data[0] == "}")
                                            {
                                                break;
                                            }
                                            else
                                            {
                                                string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                                                float Base_Change_Fitness = Parameters.get_FLOAT(line_Data[1]);

                                                char base;
                                                int mutation_Index;

                                                function.get_base_mutation(query_Base_Change, base, mutation_Index);

                                                if (base == 'A')
                                                {
                                                    A_0_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                                                    // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                                                }
                                                else if (base == 'T')
                                                {
                                                    T_1_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                                                }
                                                else if (base == 'G')
                                                {
                                                    G_2_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                                                }
                                                else if (base == 'C')
                                                {
                                                    C_3_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                                                }
                                                else
                                                {
                                                    cout << "ERROR IN FITNESS POINT " << fit_Point << " MUTATION BASE " << query_Base_Change << endl;
                                                    exit(-1);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else
                            {
                                cout << "ERROR: " << check_fitness_Point << " PARAMETER \"Position\" IS MISSING." << endl;
                                exit(-1);
                            }
                        }
                        else
                        {
                            cout << "ERROR: " << check_fitness_Point << " IS MISSING." << endl;
                            exit(-1);
                        }
                    }
                }
                else
                {
                    cout << "ERROR: PARAMETER \"Number of fitess points\" IS MISSING." << endl;
                    exit(-1);
                }

                // if (line_Data[0] == "}")
                // {
                //     break;
                // }

                // if (line_Data[0] == "\"Number of fitess points\"")
                // {
                //     this->fitness_points = Parameters.get_INT(line_Data[1]);
                //     cout << "Number of fitess points: " << this->fitness_points << endl
                //          << endl;
                // }
                // else
                // {
                //     if (this->fitness_points != -1)
                //     {
                //         string fitness_keyword = "Point ";

                //         fitness_point_Locations = (int *)malloc(fitness_points * sizeof(int));

                //         // create a sequence length and fill. Saves search time

                //         this->A_0_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);
                //         this->T_1_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);
                //         this->G_2_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);
                //         this->C_3_fitness = function.create_FLOAT_2D_arrays(this->fitness_points, 4);

                //         for (int fit_Point = 0; fit_Point < fitness_points; fit_Point++)
                //         {
                //             A_0_fitness[fit_Point][0] = -1;
                //             T_1_fitness[fit_Point][1] = -1;
                //             G_2_fitness[fit_Point][2] = -1;
                //             C_3_fitness[fit_Point][3] = -1;

                //             string check_fitness_Point = fitness_keyword + to_string(fit_Point + 1);

                //             if (Parameters.get_STRING(line_Data[0]) == check_fitness_Point)
                //             {
                //                 cout << "Configuring " << check_fitness_Point << ":\n";

                //                 string fitness_Line;
                //                 while (getline(replication_Profile, fitness_Line))
                //                 {
                //                     if (fitness_Line != "}" && fitness_Line != "")
                //                     {
                //                         string remov_Fintess_Line = fitness_Line;
                //                         int i_2 = 0;
                //                         while (remov_Fintess_Line[i_2] == ' ')
                //                         {
                //                             i_2++; // Skip leading spaces
                //                         }
                //                         remov_Fintess_Line.erase(0, i_2);

                //                         if (remov_Fintess_Line.at(0) != '#')
                //                         {
                //                             if (remov_Fintess_Line.at(remov_Fintess_Line.size() - 1) == ',')
                //                             {
                //                                 remov_Fintess_Line = remov_Fintess_Line.substr(0, remov_Fintess_Line.length() - 1);
                //                             }

                //                             function.split(line_Data, remov_Fintess_Line, ':');

                //                             if (line_Data[0] == "}")
                //                             {
                //                                 one_Check = 0;
                //                                 break;
                //                             }
                //                             else
                //                             {
                //                                 // cout << remov_Fintess_Line << endl;
                //                                 if (line_Data[0] == "\"Position\"")
                //                                 {
                //                                     fitness_point_Locations[fit_Point] = stoi(line_Data[1]);
                //                                     sequence_Mutation_tracker[0][fitness_point_Locations[fit_Point] - 1] = fitness_point_Locations[fit_Point];
                //                                     cout << "Target point: " << fitness_point_Locations[fit_Point] << endl;
                //                                     continue;
                //                                 }
                //                                 else
                //                                 {
                //                                     if (one_Check == 0)
                //                                     {
                //                                         cout << "Configuring fitness effects\n"
                //                                              << endl;
                //                                         one_Check = 1;
                //                                     }
                //                                     string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                //                                     float Base_Change_Fitness = Parameters.get_FLOAT(line_Data[1]);

                //                                     char base;
                //                                     int mutation_Index;

                //                                     function.get_base_mutation(query_Base_Change, base, mutation_Index);

                //                                     if (base == 'A')
                //                                     {
                //                                         A_0_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                //                                         // cout << remov_Fintess_Line << "\t" << A_0_fitness[fit_Point][mutation_Index] << endl;
                //                                     }
                //                                     else if (base == 'T')
                //                                     {
                //                                         T_1_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                //                                     }
                //                                     else if (base == 'G')
                //                                     {
                //                                         G_2_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                //                                     }
                //                                     else if (base == 'C')
                //                                     {
                //                                         C_3_fitness[fit_Point][mutation_Index] = Base_Change_Fitness;
                //                                     }
                //                                     else
                //                                     {
                //                                         cout << "ERROR IN FITNESS POINT " << fit_Point << " MUTATION BASE " << query_Base_Change << endl;
                //                                         exit(-1);
                //                                     }
                //                                 }
                //                             }
                //                         }
                //                     }
                //                 }
                //             }
                //             // REMOVE
                //             // break;
                //         }
                //     }
                // }
            }
            else
            {
                cout << "ERROR: FITNESS BLOCK IS MISSING" << endl;
                exit(-1);
            }

            replication_Profile.close();
        }
    }
    else
    {
        cout << "Mutations are deactivated" << endl;
    }
    //  cout << endl;
}

void within_host_test_2::configure_Mutation_Profiles(int num_Generations)
{
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);
    // parameter_load Parameters = parameter_load();

    if (this->mutation_Activate_parent == 1)
    {
        cout << "Configuring mutation profiles\n"
             << endl;

        fstream replication_Profile;
        replication_Profile.open(this->replication_profile_file, ios::in);
        if (replication_Profile.is_open())
        {
            parameter_load Parameters = parameter_load();

            string line;
            getline(replication_Profile, line);

            vector<string> line_Data;

            int mutation_Profile_BLOCK = 0;
            // int block_Open = 0;

            while (getline(replication_Profile, line))
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

                        if (line_Data[0] == "\"Mutation rate\"")
                        {
                            mutation_Profile_BLOCK = 1;
                            break;
                        }
                    }
                }
            }
            if (mutation_Profile_BLOCK == 1)
            {

                while (getline(replication_Profile, line))
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

                            if (line_Data[0] == "}")
                            {
                                break;
                            }
                            else
                            {

                                if (line_Data[0] == "\"Number of mutation hotspots\"")
                                {
                                    this->mutation_hotspots = Parameters.get_INT(line_Data[1]);
                                    cout << "Number of mutation hotspots: " << this->mutation_hotspots << endl
                                         << endl;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (this->mutation_hotspots != -1)
                {
                    string keyword = "Hotspot ";

                    this->mutation_rates_Hotspot_generation = function.create_FLOAT_2D_arrays(this->mutation_hotspots, num_Generations);
                    this->mutation_Regions_start_stop = function.create_INT_2D_arrays(this->mutation_hotspots, 2);

                    this->A_0_mutation = function.create_FLOAT_2D_arrays(this->mutation_hotspots, 4);
                    this->T_1_mutation = function.create_FLOAT_2D_arrays(this->mutation_hotspots, 4);
                    this->G_2_mutation = function.create_FLOAT_2D_arrays(this->mutation_hotspots, 4);
                    this->C_3_mutation = function.create_FLOAT_2D_arrays(this->mutation_hotspots, 4);

                    for (int hotspot = 0; hotspot < this->mutation_hotspots; hotspot++)
                    {
                        int hotspot_Block = -1;
                        string check_keyword = keyword + to_string(hotspot + 1);
                        while (getline(replication_Profile, line))
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

                                    // if (line_Data[0] == "}")
                                    // {
                                    //     break;
                                    // }
                                    // else
                                    //{

                                    if (Parameters.get_STRING(line_Data[0]) == check_keyword)
                                    {
                                        cout << "Configuring " << check_keyword << endl;
                                        // this->mutation_hotspots = Parameters.get_INT(line_Data[1]);
                                        hotspot_Block = 1;
                                        break;
                                    }
                                    //}
                                }
                            }
                        }
                        if (hotspot_Block == 1)
                        {
                            int region_Check = -1;
                            while (getline(replication_Profile, line))
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

                                        if (line_Data[0] == "}")
                                        {
                                            break;
                                        }
                                        else
                                        {

                                            if (line_Data[0] == "\"Region\"")
                                            {
                                                mutation_Regions_start_stop[hotspot][0] = stoi(line_Data[1].substr(1, line_Data[1].size()));
                                                mutation_Regions_start_stop[hotspot][1] = stoi(line_Data[2].substr(0, line_Data[2].size() - 1));
                                                cout << "Target region: From " << mutation_Regions_start_stop[hotspot][0] << " to " << mutation_Regions_start_stop[hotspot][1] << endl;
                                                region_Check = 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            if (region_Check == 1)
                            {
                                int clock_Model = -1;
                                while (getline(replication_Profile, line))
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

                                            if (line_Data[0] == "}")
                                            {
                                                break;
                                            }
                                            else
                                            {

                                                if (line_Data[0] == "\"Clock model\"")
                                                {
                                                    cout << "Configuring clock model" << endl;
                                                    clock_Model = 1;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (clock_Model == 1)
                                {
                                    string clock_Model_Type = "";
                                    int clock_Check = -1;
                                    while (getline(replication_Profile, line))
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

                                                if (line_Data[0] == "}")
                                                {
                                                    break;
                                                }
                                                else
                                                {

                                                    if (line_Data[0] == "\"Type\"")
                                                    {
                                                        clock_Model_Type = Parameters.get_STRING(line_Data[1]);
                                                        clock_Check = 1;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (clock_Check == 1)
                                    {
                                        if (clock_Model_Type == "Relaxed gamma")
                                        {
                                            cout << "Clock model: Relaxed Gamma" << endl;

                                            float shape_Param = -1;
                                            float scale_Param = -1;

                                            int line_Count = 0;

                                            while (getline(replication_Profile, line))
                                            {
                                                if (line_Count == 2)
                                                {
                                                    break;
                                                }
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

                                                        if (line_Data[0] == "}")
                                                        {
                                                            break;
                                                        }
                                                        else
                                                        {

                                                            if (line_Data[0] == "\"Shape of mutation rate\"")
                                                            {
                                                                shape_Param = Parameters.get_FLOAT(line_Data[1]);
                                                                line_Count++;
                                                                continue;
                                                            }
                                                            else if (line_Data[0] == "\"Scale of mutation rate\"")
                                                            {
                                                                scale_Param = Parameters.get_FLOAT(line_Data[1]);
                                                                line_Count++;
                                                                continue;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (line_Count == 2)
                                            {
                                                float *mutation_rates_Generations = function.gamma_distribution_CUDA(num_Generations, shape_Param, scale_Param);
                                                function.array_Copy_Float(mutation_rates_Hotspot_generation, mutation_rates_Generations, hotspot, num_Generations);
                                            }
                                            else
                                            {
                                                cout << "ERROR IN Relaxed gamma CLOCK MODEL PARAMTERS." << endl;
                                                exit(-1);
                                            }
                                        }
                                        else if (clock_Model_Type == "Strict fixed")
                                        {
                                            cout << "Clock model: Strict Fixed" << endl;

                                            float mutation_rate_Generations = -1;
                                            int line_Count = 0;

                                            while (getline(replication_Profile, line))
                                            {
                                                if (line_Count == 1)
                                                {
                                                    break;
                                                }
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

                                                        if (line_Data[0] == "}")
                                                        {
                                                            break;
                                                        }
                                                        else
                                                        {

                                                            if (line_Data[0] == "\"Fixed mutation rate\"")
                                                            {
                                                                mutation_rate_Generations = Parameters.get_FLOAT(line_Data[1]);
                                                                line_Count++;
                                                                continue;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (line_Count == 1)
                                            {
                                                cout << "Fixed mutation rate: " << mutation_rate_Generations << endl;

                                                for (int fill_Gen = 0; fill_Gen < num_Generations; fill_Gen++)
                                                {
                                                    mutation_rates_Hotspot_generation[hotspot][fill_Gen] = mutation_rate_Generations;
                                                    // cout << mutation_rates_Hotspot_generation[hotspot][fill_Gen] << endl;
                                                }
                                            }
                                            else
                                            {
                                                cout << "ERROR IN Strict fixed CLOCK MODEL PARAMTER." << endl;
                                                exit(-1);
                                            }
                                        }
                                        else if (clock_Model_Type == "Strict Gamma")
                                        {
                                            cout << "Clock model: Strict Gamma" << endl;

                                            float shape_Param = -1;
                                            float scale_Param = -1;

                                            int line_Count = 0;

                                            while (getline(replication_Profile, line))
                                            {
                                                if (line_Count == 2)
                                                {
                                                    break;
                                                }
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

                                                        if (line_Data[0] == "}")
                                                        {
                                                            break;
                                                        }
                                                        else
                                                        {

                                                            if (line_Data[0] == "\"Shape of mutation rate\"")
                                                            {
                                                                shape_Param = Parameters.get_FLOAT(line_Data[1]);
                                                                line_Count++;
                                                                continue;
                                                            }
                                                            else if (line_Data[0] == "\"Scale of mutation rate\"")
                                                            {
                                                                scale_Param = Parameters.get_FLOAT(line_Data[1]);
                                                                line_Count++;
                                                                continue;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if (line_Count == 2)
                                            {
                                                cout << "Shape: " << shape_Param << endl
                                                     << "Scale: " << scale_Param << endl;

                                                gamma_distribution<float> gamma(shape_Param, scale_Param);

                                                float mutation_rate_Generations = gamma(gen);

                                                for (int fill_Gen = 0; fill_Gen < num_Generations; fill_Gen++)
                                                {
                                                    mutation_rates_Hotspot_generation[hotspot][fill_Gen] = mutation_rate_Generations;
                                                    // cout << mutation_rates_Hotspot_generation[hotspot][fill_Gen] << endl;
                                                }
                                            }
                                            else
                                            {
                                                cout << "ERROR IN Strict gamma CLOCK MODEL PARAMTERS." << endl;
                                                exit(-1);
                                            }
                                        }
                                        else
                                        {
                                            cout << check_keyword << " UNKNOWN CLOCK TYPE: " << clock_Model_Type << endl;
                                            exit(-1);
                                        }
                                    }
                                    else
                                    {
                                        cout << check_keyword << " \"Type\" PARAMETER IS MISSING FOR CLOCK MODEL." << endl;
                                        exit(-1);
                                    }

                                    // Probability config
                                    cout << "Configuring mutation probabilities\n"
                                         << endl;

                                    int mutation_Prob_block = -1;

                                    while (getline(replication_Profile, line))
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

                                                if (line_Data[0] == "}")
                                                {
                                                    break;
                                                }
                                                else
                                                {

                                                    if (line_Data[0] == "\"Probability\"")
                                                    {
                                                        mutation_Prob_block = 1;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (mutation_Prob_block == 1)
                                    {
                                        while (getline(replication_Profile, line))
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

                                                    if (line_Data[0] == "}")
                                                    {
                                                        break;
                                                    }
                                                    else
                                                    {
                                                        string query_Base_Change = Parameters.get_STRING(line_Data[0]);
                                                        float Base_Change_Probability = Parameters.get_FLOAT(line_Data[1]);

                                                        char base;
                                                        int mutation_Index;

                                                        function.get_base_mutation(query_Base_Change, base, mutation_Index);

                                                        // cout << query_Base_Change << " : " << to_string(Base_Change_Probability) << " : " << base << " " << mutation_Index << endl;

                                                        if (base == 'A')
                                                        {
                                                            A_0_mutation[hotspot][mutation_Index] = Base_Change_Probability;
                                                            // cout << "hotspot " << mutation_Index << endl;
                                                            // cout << remove << "\t" << A_0_mutation[hotspot][mutation_Index] << endl;
                                                        }
                                                        else if (base == 'T')
                                                        {
                                                            T_1_mutation[hotspot][mutation_Index] = Base_Change_Probability;
                                                        }
                                                        else if (base == 'G')
                                                        {
                                                            G_2_mutation[hotspot][mutation_Index] = Base_Change_Probability;
                                                        }
                                                        else if (base == 'C')
                                                        {
                                                            C_3_mutation[hotspot][mutation_Index] = Base_Change_Probability;
                                                        }
                                                        else
                                                        {
                                                            cout << "ERROR IN HOTSPOT " << hotspot << " MUTATION BASE " << query_Base_Change << endl;
                                                            exit(-1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        cout << check_keyword << " \"Probability\" PARAMETER IS MISSING." << endl;
                                        exit(-1);
                                    }
                                }
                                else
                                {
                                    cout << check_keyword << " \"Clock model\" PARAMETER IS MISSING." << endl;
                                    exit(-1);
                                }
                            }
                            else
                            {
                                cout << check_keyword << " \"Region\" PARAMETER IS MISSING." << endl;
                                exit(-1);
                            }
                        }
                        else
                        {
                            cout << check_keyword << " BLOCK IS MISSING." << endl;
                            exit(-1);
                        }
                    }
                }
                else
                {
                    cout << "PARAMETER \"Number of mutation hotspots\" MISSING." << endl;
                    exit(-1);
                }
            }

            replication_Profile.close();
        }
    }
    else
    {
        cout << "Mutations are deactivated\n"
             << endl;
    }
    // cout << endl;
}

void within_host_test_2::get_Progeny_distribution(functions_library &function)
{
    // functions_library function = functions_library();
    parameter_load Parameters = parameter_load();

    fstream replication_Profile;
    replication_Profile.open(this->replication_profile_file, ios::in);

    if (replication_Profile.is_open())
    {
        string line;
        getline(replication_Profile, line);

        vector<string> line_Data;

        int progeny_block_FOUND = 0;

        while (getline(replication_Profile, line))
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

                    if (line_Data[0] == "\"Progeny distribution\"")
                    {
                        progeny_block_FOUND = 1;
                    }

                    if (progeny_block_FOUND == 1 && line_Data[0] != "\"Progeny distribution\"")
                    {
                        if (line_Data[0] == "}")
                        {
                            break;
                        }
                        else
                        {
                            if (line_Data[0] == "\"Distribution type progeny\"")
                            {
                                function.progeny_distribution_Type = Parameters.get_STRING(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Shape progeny per_unit\"")
                            {
                                function.progeny_shape = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Scale progeny per_unit\"")
                            {
                                function.progeny_scale = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Mean progeny per_unit\"")
                            {
                                function.progeny_mean = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Dispersion progeny per_unit\"")
                            {
                                function.progeny_dispersion = Parameters.get_FLOAT(line_Data[1]);
                            }
                        }
                    }
                }
            }
        }
        replication_Profile.close();
    }

    if (function.progeny_distribution_Type == "Gamma")
    {
        if (function.progeny_shape == 0)
        {
            cout << "ERROR IN PARAMETER \"Shape progeny per_unit\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
        if (function.progeny_scale == 0)
        {
            cout << "ERROR IN PARAMETER \"Scale progeny per_unit\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
    }
    else if (function.progeny_distribution_Type == "Negative binomial")
    {
        if (function.progeny_mean == 0)
        {
            cout << "ERROR IN PARAMETER \"R progeny per_unit\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
        if (function.progeny_dispersion == 0)
        {
            cout << "ERROR IN PARAMETER \"Probability progeny per_unit\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
    }
    else
    {
        cout << "ERROR IN PARAMETER \"Distribution type progeny\": CHECK FILE: " << this->replication_profile_file << endl;
    }
}

void within_host_test_2::get_Cell_distribution(string &distribution_Type)
{
    functions_library function = functions_library();
    parameter_load Parameters = parameter_load();

    fstream replication_Profile;
    replication_Profile.open(this->replication_profile_file, ios::in);

    if (replication_Profile.is_open())
    {
        string line;
        getline(replication_Profile, line);

        vector<string> line_Data;

        int cell_block_FOUND = 0;

        while (getline(replication_Profile, line))
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

                    if (line_Data[0] == "\"Viral distribution per cell\"")
                    {
                        cell_block_FOUND = 1;
                    }

                    if (cell_block_FOUND == 1 && line_Data[0] != "\"Viral distribution per cell\"")
                    {
                        if (line_Data[0] == "}")
                        {
                            break;
                        }
                        else
                        {
                            if (line_Data[0] == "\"Distribution type cell\"")
                            {
                                distribution_Type = Parameters.get_STRING(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Shape units per_cell\"")
                            {
                                this->cell_shape = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Scale units per_cell\"")
                            {
                                this->cell_scale = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"R units per_cell\"")
                            {
                                this->cell_r = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Probability units per_cell\"")
                            {
                                this->cell_prob = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Cell limits\"")
                            {
                                if (Parameters.get_STRING(line_Data[1]) == "Yes")
                                {
                                    this->cell_Limit = 1;
                                }
                            }
                            else if (line_Data[0] == "\"Shape cells\"")
                            {
                                this->total_Cells_shape = Parameters.get_FLOAT(line_Data[1]);
                            }
                            else if (line_Data[0] == "\"Scale cells\"")
                            {
                                this->total_Cells_scale = Parameters.get_FLOAT(line_Data[1]);
                            }
                        }
                    }
                }
            }
        }
        replication_Profile.close();
    }

    if (distribution_Type == "Gamma")
    {
        if (cell_shape == 0)
        {
            cout << "ERROR IN PARAMETER \"Shape units per_cell\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
        if (cell_scale == 0)
        {
            cout << "ERROR IN PARAMETER \"Scale units per_cell\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
    }
    else if (distribution_Type == "Negative binomial")
    {
        if (cell_r == 0)
        {
            cout << "ERROR IN PARAMETER \"R units per_cell\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
        if (cell_prob == 0)
        {
            cout << "ERROR IN PARAMETER \"Probability units per_cell\" (cannot be zero):\nCHECK FILE: " << this->replication_profile_file << endl;
        }
    }
    else
    {
        cout << "ERROR IN PARAMETER \"Distribution type cell\": CHECK FILE: " << this->replication_profile_file << endl;
    }
}

void within_host_test_2::allocate_Phases(int num_Generations)
{
    // num_Generations = 55;
    functions_library function = functions_library(tot_Blocks, tot_ThreadsperBlock, gpu_Limit, CPU_cores);

    generation_modes = (int *)malloc(num_Generations * sizeof(int));
    float *times_of_Phases;

    fstream replication_Profile;
    replication_Profile.open(this->replication_profile_file, ios::in);

    int number_of_Phases = -1;

    if (replication_Profile.is_open())
    {
        string line;
        getline(replication_Profile, line);

        parameter_load Parameters = parameter_load();

        vector<string> line_Data;

        int phase_block_FOUND = 0;

        while (getline(replication_Profile, line))
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

                    if (line_Data[0] == "\"Phases\"")
                    {
                        phase_block_FOUND = 1;
                        break;
                    }
                }
            }
        }
        if (phase_block_FOUND == 1)
        {
            while (getline(replication_Profile, line))
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

                        if (line_Data[0] == "}")
                        {
                            break;
                        }
                        else
                        {
                            if (line_Data[0] == "\"Number of phases\"")
                            {
                                number_of_Phases = Parameters.get_INT(line_Data[1]);
                                cout << "Number of Phases: " << number_of_Phases << endl
                                     << endl;
                                break;
                            }
                        }
                    }
                }
            }
            if (number_of_Phases != -1)
            {
                string keyword = "Phase ";

                this->phase_Modes = (int *)malloc(number_of_Phases * sizeof(int));
                this->phase_parameters = function.create_Fill_2D_array_FLOAT(number_of_Phases, 3, -1);
                times_of_Phases = (float *)malloc(number_of_Phases * sizeof(float));

                for (int phase = 0; phase < number_of_Phases; phase++)
                {
                    int phase_Block = -1;
                    times_of_Phases[phase] = -1;
                    string check_Keyword = keyword + to_string(phase + 1);

                    while (getline(replication_Profile, line))
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
                                function.split(line_Data, remove, ':');

                                if (Parameters.get_STRING(line_Data[0]) == check_Keyword)
                                {
                                    cout << "Configuring " << check_Keyword << endl;
                                    // this->mutation_hotspots = Parameters.get_INT(line_Data[1]);
                                    phase_Block = 1;
                                    break;
                                }
                            }
                        }
                    }
                    if (phase_Block == 1)
                    {
                        while (getline(replication_Profile, line))
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

                                    if (line_Data[0] == "}")
                                    {
                                        break;
                                    }
                                    else
                                    {
                                        if (line_Data[0] == "\"Time\"")
                                        {
                                            times_of_Phases[phase] = Parameters.get_FLOAT(line_Data[1]);
                                            cout << "Time ratio: " << times_of_Phases[phase] << endl;
                                        }
                                        else if (line_Data[0] == "\"Mode\"")
                                        {
                                            string mode_Phase = Parameters.get_STRING(line_Data[1]);
                                            transform(mode_Phase.begin(), mode_Phase.end(), mode_Phase.begin(), ::toupper);
                                            cout << "Mode: " << mode_Phase << endl;
                                            if (mode_Phase == "GROWTH")
                                            {
                                                phase_Modes[phase] = 0;
                                            }
                                            else if (mode_Phase == "STATIONARY")
                                            {
                                                phase_Modes[phase] = 1;
                                            }
                                            else if (mode_Phase == "DEPRICIATION")
                                            {
                                                phase_Modes[phase] = 2;
                                            }
                                            else
                                            {
                                                cout << "ERROR IN PHASE " << phase + 1 << " MODE: " << mode_Phase << endl;
                                                exit(-1);
                                            }
                                        }
                                        else if (line_Data[0] == "\"Variance\"")
                                        {
                                            phase_parameters[phase][0] = Parameters.get_FLOAT(line_Data[1]);
                                            cout << "Variance for Stationary: " << phase_parameters[phase][0] << endl;
                                        }
                                        else if (line_Data[0] == "\"Alpha\"")
                                        {
                                            phase_parameters[phase][0] = Parameters.get_FLOAT(line_Data[1]);
                                            cout << "Alpha parameter for Deprication: " << phase_parameters[phase][0] << endl;
                                        }
                                        else if (line_Data[0] == "\"Beta\"")
                                        {
                                            phase_parameters[phase][1] = Parameters.get_FLOAT(line_Data[1]);
                                            cout << "Beta parameter for Deprication: " << phase_parameters[phase][1] << endl;
                                        }
                                    }
                                }
                            }
                        }
                        // check configuration
                        if (phase_Modes[phase] == 1)
                        {
                            if (phase_parameters[phase][0] == -1)
                            {
                                cout << "PHASE " << phase + 1 << " STATIONARY MODE VARIANCE NOT CONFIGURED." << endl;
                                exit(-1);
                            }
                        }
                        else if (phase_Modes[phase] == 2)
                        {
                            string error_1 = "";
                            string error_2 = "";
                            if (phase_parameters[phase][0] == -1)
                            {
                                error_1 = "PHASE " + to_string(phase + 1) + " DEPRICIATION MODE ALPHA NOT CONFIGURED.\n";
                            }
                            if (phase_parameters[phase][1] == -1)
                            {
                                error_2 = "PHASE " + to_string(phase + 1) + " DEPRICIATION MODE BETA NOT CONFIGURED.\n";
                            }

                            if (error_1 != "" || error_2 != "")
                            {
                                cout << error_1;
                                cout << error_2;
                                exit(-1);
                            }
                        }
                        if (times_of_Phases[phase] == -1)
                        {
                            cout << "PHASE " << phase + 1 << " TIME NOT CONFIGURED." << endl;
                            exit(-1);
                        }
                    }
                    else
                    {
                        cout << "MISSING \"Phase " << phase + 1 << "\" BLOCK." << endl;
                        exit(-1);
                    }
                    cout << endl;
                }
            }
        }
        else
        {
            cout << "PARAMETER \"Phases\" MISSING." << endl;
            exit(-1);
        }
        replication_Profile.close();
    }

    int total = 0;

    for (size_t i = 0; i < number_of_Phases; i++)
    {
        // cout << phase_Times[i] << "\t" << phase_Modes[i] << endl;
        int phase_Fill = round(times_of_Phases[i] * (float)num_Generations) + total;

        for (int fill = total; fill < phase_Fill; fill++)
        {
            if (fill >= num_Generations)
            {
                break;
            }
            else
            {
                generation_modes[fill] = i;
                total++;
            }
        }
    }

    if (total < num_Generations)
    {
        for (size_t i = total; i < num_Generations; i++)
        {
            generation_modes[i] = number_of_Phases - 1;
        }
    }

    int common = 0;
    int pre_gen = 0;

    for (size_t i = 0; i < num_Generations; i++)
    {
        if (i == 0)
        {
            common = generation_modes[i];
        }
        else
        {
            if (common != generation_modes[i])
            {
                cout << "From generation " << pre_gen + 1 << " to " << i << ": Phase: " << generation_modes[i - 1] + 1 << endl;
                pre_gen = i;
                common = generation_modes[i];
            }
        }
        // cout << "Generation: " << i + 1 << "\t"
        //      << "Phase: " << generation_modes[i] + 1 << endl;
    }

    cout << "From generation " << pre_gen + 1 << " to " << num_Generations << ": Phase: " << generation_modes[num_Generations - 1] + 1 << endl;

    // exit(-1);

    // vector<string> phase_Data;
    // vector<string> phase_HEADERS;

    // if (replication_Profile.is_open())
    // {
    //     string line;
    //     getline(replication_Profile, line);

    //     vector<string> line_Data;

    //     int phase_block_FOUND = 0;

    //     while (getline(replication_Profile, line))
    //     {
    //         if (line != "}" && line != "")
    //         {
    //             string remove = line;
    //             int i = 0;
    //             while (remove[i] == ' ')
    //             {
    //                 i++; // Skip leading spaces
    //             }
    //             remove.erase(0, i);

    //             if (remove.at(0) != '#')
    //             {
    //                 if (remove.at(remove.size() - 1) == ',')
    //                 {
    //                     remove = remove.substr(0, remove.length() - 1);
    //                 }

    //                 // cout << remove << endl;
    //                 function.split(line_Data, remove, ':');

    //                 if (line_Data[0] == "\"Phases\"")
    //                 {
    //                     phase_block_FOUND = 1;
    //                 }

    //                 if (phase_block_FOUND == 1 && line_Data[0] != "\"Phases\"")
    //                 {
    //                     if (line_Data[0] == "}")
    //                     {
    //                         break;
    //                     }
    //                     else
    //                     {
    //                         if (line_Data[0] == "\"Number of phases\"")
    //                         {
    //                             number_of_Phases = stoi(line_Data[1]);
    //                         }
    //                         else
    //                         {
    //                             phase_HEADERS.push_back(line_Data[0]);
    //                             phase_Data.push_back(line_Data[1]);
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     replication_Profile.close();
    // }

    // cout << "Identified " << number_of_Phases << " phases:" << endl;

    // // phase modes
    // // 0 = growth
    // // 1 = stationary
    // // 2 = depreciation

    // float *phase_Times = (float *)malloc(number_of_Phases * sizeof(float));
    // phase_Modes = (int *)malloc(number_of_Phases * sizeof(int));

    // parameter_load Parameters = parameter_load();

    // for (size_t i = 0; i < phase_HEADERS.size(); i++)
    // {
    //     string end_line = phase_HEADERS[i].substr(phase_HEADERS[i].size() - 5);

    //     if (end_line == "time\"")
    //     {
    //         for (int phase_time_value = 0; phase_time_value < number_of_Phases; phase_time_value++)
    //         {
    //             int phase_Value = phase_time_value + 1;
    //             string check_for = "\"Phase " + to_string(phase_Value) + " time\"";

    //             if (phase_HEADERS[i] == check_for)
    //             {
    //                 phase_Times[phase_time_value] = Parameters.get_FLOAT(phase_Data[i]);
    //                 break;
    //             }
    //         }
    //     }
    //     else if (end_line == "mode\"")
    //     {
    //         for (int phase_time_value = 0; phase_time_value < number_of_Phases; phase_time_value++)
    //         {
    //             int phase_Value = phase_time_value + 1;
    //             string check_for = "\"Phase " + to_string(phase_Value) + " mode\"";

    //             if (phase_HEADERS[i] == check_for)
    //             {
    //                 if (phase_Data[i] == "\"Growth\"")
    //                 {
    //                     phase_Modes[phase_time_value] = 0;
    //                 }
    //                 else if (phase_Data[i] == "\"Stationary\"")
    //                 {
    //                     phase_Modes[phase_time_value] = 1;
    //                 }
    //                 else if (phase_Data[i] == "\"Depreciation\"")
    //                 {
    //                     phase_Modes[phase_time_value] = 2;
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }

    // int total = 0;

    // for (size_t i = 0; i < number_of_Phases; i++)
    // {
    //     // cout << phase_Times[i] << "\t" << phase_Modes[i] << endl;
    //     int phase_Fill = round(phase_Times[i] * (float)num_Generations) + total;

    //     for (int fill = total; fill < phase_Fill; fill++)
    //     {
    //         if (fill >= num_Generations)
    //         {
    //             break;
    //         }
    //         else
    //         {
    //             generation_modes[fill] = i;
    //             total++;
    //         }
    //     }
    // }

    // if (total < num_Generations)
    // {
    //     for (size_t i = total; i < num_Generations; i++)
    //     {
    //         generation_modes[i] = number_of_Phases - 1;
    //     }
    // }

    // int common = 0;
    // int pre_gen = 0;

    // for (size_t i = 0; i < num_Generations; i++)
    // {
    //     if (i == 0)
    //     {
    //         common = generation_modes[i];
    //     }
    //     else
    //     {
    //         if (common != generation_modes[i])
    //         {
    //             cout << "From generation " << pre_gen + 1 << " to " << i << ": Phase: " << generation_modes[i - 1] + 1 << endl;
    //             pre_gen = i;
    //             common = generation_modes[i];
    //         }
    //     }
    //     // cout << "Generation: " << i + 1 << "\t"
    //     //      << "Phase: " << generation_modes[i] + 1 << endl;
    // }

    // cout << "From generation " << pre_gen + 1 << " to " << num_Generations << ": Phase: " << generation_modes[num_Generations - 1] + 1 << endl;
}