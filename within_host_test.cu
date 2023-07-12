#include "within_host_test.cuh"
#include "functions_library.cuh"

within_host_test::within_host_test(int cuda_Device_ID, string reference_Genome_location)
{
    cout << "TEST PHASE: WITHIN HOST EVOLUTION\n"
         << endl;

    cudaSetDevice(cuda_Device_ID);

    functions_library functions = functions_library();
    functions.print_Cuda_device(cuda_Device_ID, this->tot_Blocks, this->tot_ThreadsperBlock);
    this->full_Reference_genome = functions.read_Reference(reference_Genome_location, this->reference_Genome_header, this->genome_Size);
}

__global__ void cuda_fill_Master_sequence(int genome_Size, int **master_Sequence, char *cuda_reference)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Size)
    {
        // A = 0
        // T = 1
        // G = 2
        // C = 3

        char base = cuda_reference[tid];

        if (base == 'A')
        {
            master_Sequence[0][tid] = 0;
        }
        else if (base == 'T')
        {
            master_Sequence[0][tid] = 1;
        }
        else if (base == 'G')
        {
            master_Sequence[0][tid] = 2;
        }
        else if (base == 'C')
        {
            master_Sequence[0][tid] = 3;
        }

        // master_Sequence[0][tid] = 0;
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void setup_distribution(curandState_t *d_states, int progeny_Sequences)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < progeny_Sequences)
    {
        curand_init(tid + 7, tid, 0, &d_states[tid]);
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void normal_dist(curandState_t *d_states, double *values, int progeny_Sequences)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < progeny_Sequences)
    {
        curandState_t localState = d_states[tid];
        double z = (double)(curand_poisson(&localState, 2));
        values[tid] = z;
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void uniform_distribution(int progeny_Sequences, int *values, int value)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < progeny_Sequences)
    {
        values[tid] = value;
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_progeny_Sequence_fill(int genome_Size, int **parent_child, int **parent_sequences, int **progeny_Sequences, int progeny_in_generation, int ID_first_in_generation)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < progeny_in_generation)
    {
        int current_progeny_ID = tid + ID_first_in_generation - progeny_in_generation;
        int parent = parent_child[current_progeny_ID][0];

        for (int i = 0; i < genome_Size; i++)
        {
            // replication error rate
                // CODE for error checking if present
            progeny_Sequences[tid][i] = parent_sequences[parent][i];
        }

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_progeny_life(int progeny_in_generation, int ID_first_in_generation, int *distribution_Progeny_number, int *number_of_progeny)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < progeny_in_generation)
    {
        // int current_progeny_ID = tid + ID_first_in_generation;

        // CODE to mutate progeny sequence


        // code to multiply
        number_of_progeny[tid] = distribution_Progeny_number[tid];

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_old_to_new_parent_child(int **old_parent_child, int **new_parent_child, int old_progeny_count)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < old_progeny_count)
    {
        new_parent_child[tid][0] = old_parent_child[tid][0];
        new_parent_child[tid][1] = old_parent_child[tid][1];

        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_old_to_new_parent_sequences(int **old_sequences, int **new_parent_sequences, int number_of_old_parents, int genome_Size)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Size)
    {
        for (int i = 0; i < number_of_old_parents; i++)
        {
            new_parent_sequences[i][tid] = old_sequences[i][tid];
        }
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_fill_new_parents(int number_of_old_parents, int number_of_new_progeny, int **progeny_Sequences, int **new_parent_sequences, int genome_Size)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < number_of_new_progeny)
    {
        int index = number_of_old_parents + tid;
        for (int i = 0; i < genome_Size; i++)
        {
            new_parent_sequences[index][i] = progeny_Sequences[tid][i];
        }
        tid += blockDim.x * gridDim.x;
    }
}

__global__ void cuda_Sequence_swap(int genome_Size, int total_virues, int **temp_sequences, int **parent_sequences)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < genome_Size)
    {
        for (int i = 0; i < total_virues; i++)
        {
            parent_sequences[i][tid] = temp_sequences[i][tid];
        }

        tid += blockDim.x * gridDim.x;
    }
}

void within_host_test::ingress()
{
    cout << "Simulating within host\n\n";

    int total_viruses = 1;
    int number_of_generations = 3;

    // int progeny_current_Generation = 0;

    int **cuda_parent_Sequences;
    cudaMallocManaged(&cuda_parent_Sequences, genome_Size * 1 * sizeof(int));
    int **tmp = (int **)malloc(1 * sizeof(tmp[0]));
    for (int i = 0; i < 1; i++)
    {
        cudaMalloc((void **)&tmp[i], genome_Size * sizeof(tmp[0][0]));
    }
    cudaMemcpy(cuda_parent_Sequences, tmp, 1 * sizeof(int *), cudaMemcpyHostToDevice);

    free(tmp);

    // load reference sequence to be converted
    char *reference_full, *cuda_reference;
    reference_full = (char *)malloc((full_Reference_genome.size() + 1) * sizeof(char));
    cudaMallocManaged(&cuda_reference, (full_Reference_genome.size() + 1) * sizeof(char));
    strcpy(reference_full, full_Reference_genome.c_str());
    cudaMemcpy(cuda_reference, reference_full, (full_Reference_genome.size() + 1) * sizeof(char), cudaMemcpyHostToDevice);

    free(reference_full);

    cuda_fill_Master_sequence<<<tot_Blocks, tot_ThreadsperBlock>>>(this->genome_Size, cuda_parent_Sequences, cuda_reference);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        printf("CUDA Error 1: %s\n", cudaGetErrorString(err));

        // Possibly: exit(-1) if program cannot continue....
    }

    // int **parent_sequence = (int **)malloc(1 * sizeof(int *));
    // for (int i = 0; i < 1; i++)
    // {
    //     parent_sequence[i] = (int *)malloc(genome_Size * sizeof(int));
    // }

    // for (int i = 0; i < 1; i++)
    // {
    //     cudaMemcpy(parent_sequence[i], cuda_parent_Sequences[i], genome_Size * sizeof(cuda_parent_Sequences[0][0]), cudaMemcpyDeviceToHost);
    // }

    // for (size_t i = 0; i < genome_Size; i++)
    // {
    //     cout << parent_sequence[0][i];
    // }
    // cout << endl;

    // exit(0);

    // change later to distribution draw
    int new_progeny_master = 2;
    total_viruses = total_viruses + new_progeny_master;

    // actually one less, but makes loop easier
    // int ID_Last_child_in_generation = ID_First_child_in_generation + new_progeny_master;

    vector<pair<int, int>> parent_child_Full;

    // rows
    int **parent_child, **cuda_parent_child;
    parent_child = (int **)malloc(total_viruses * sizeof(int *));

    // columns
    for (int i = 0; i < total_viruses; i++)
    {
        parent_child[i] = (int *)malloc(2 * sizeof(int));
    }

    for (int row = 0; row < total_viruses; row++)
    {
        parent_child_Full.push_back(make_pair(0, row));
        parent_child[row][0] = 0;
        parent_child[row][1] = row;
    }

    cudaMallocManaged(&cuda_parent_child, total_viruses * 2 * sizeof(int));
    int **tmp_2 = (int **)malloc(total_viruses * sizeof(tmp_2[0]));
    for (size_t i = 0; i < total_viruses; i++)
    {
        cudaMalloc((void **)&tmp_2[i], 2 * sizeof(tmp_2[0][0]));
    }
    cudaMemcpy(cuda_parent_child, tmp_2, total_viruses * sizeof(int *), cudaMemcpyHostToDevice);

    for (size_t i = 0; i < total_viruses; i++)
    {
        cudaMemcpy(tmp_2[i], parent_child[i], 2 * sizeof(cuda_parent_child[0][0]), cudaMemcpyHostToDevice);
    }

    free(tmp_2);

    int parents = 1;

    char distribution_Type = 'N';
    distribution_Type = 'N';

    random_device rd;
    mt19937 gen(rd());

    for (int generation = 0; generation < number_of_generations; generation++)
    {

        cout << "Processing generation: " << generation + 1 << " of " << number_of_generations << endl;

        int ID_First_child_in_generation = total_viruses;

        // int *new_progeny_child, *cuda_new_progeny_child;
        // cudaMallocManaged(&cuda_new_progeny_child, current_progeny * sizeof(int));

        // store all new sequences progeny
        int **cuda_progeny_Sequences;
        cudaMallocManaged(&cuda_progeny_Sequences, new_progeny_master * this->genome_Size * sizeof(int));

        int **tmp_3 = (int **)malloc(new_progeny_master * sizeof(tmp_3[0]));
        for (int i = 0; i < new_progeny_master; i++)
        {
            cudaMalloc((void **)&tmp_3[i], this->genome_Size * sizeof(tmp_3[0][0]));
        }

        cudaMemcpy(cuda_progeny_Sequences, tmp_3, new_progeny_master * sizeof(int *), cudaMemcpyHostToDevice);
        free(tmp_3);

        cout << "Filling progeny sequences" << endl;
        // cuda_progeny_Sequence_fill(int genome_Size, int **parent_child, int **parent_sequences, int **progeny_Sequences, int progeny_in_generation, int ID_first_in_generation)
        cuda_progeny_Sequence_fill<<<tot_Blocks, tot_ThreadsperBlock>>>(genome_Size, cuda_parent_child, cuda_parent_Sequences, cuda_progeny_Sequences, new_progeny_master, ID_First_child_in_generation);
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            printf("CUDA Error 2: %s\n", cudaGetErrorString(err));
            // Possibly: exit(-1) if program cannot continue....
        }

        int *cuda_progeny_values, *cuda_progeny_distribution;
        cudaMallocManaged(&cuda_progeny_distribution, new_progeny_master * sizeof(int));
        cudaMallocManaged(&cuda_progeny_values, new_progeny_master * sizeof(int));

        // now we see how many kids the progeny in turn produce

        cout << "Generating distributions" << endl;

        if (distribution_Type == 'U')
        {
            cout << "Uniform Distribution" << endl;
            uniform_distribution<<<tot_Blocks, tot_ThreadsperBlock>>>(new_progeny_master, cuda_progeny_distribution, 2);
            cudaDeviceSynchronize();
            // cuda_progeny_life(int progeny_in_generation, int ID_first_in_generation, int *distribution_Progeny_number, int *number_of_progeny)
            err = cudaGetLastError();
            if (err != cudaSuccess)
            {
                printf("CUDA Error 4: %s\n", cudaGetErrorString(err));
                // Possibly: exit(-1) if program cannot continue....
            }
        }
        else if (distribution_Type == 'N')
        {
            cout << "Negative Binomial Distribution" << endl;
            negative_binomial_distribution<int> d(5, 0.75);

            int *temp__progeny_distribution;
            temp__progeny_distribution = (int *)malloc(new_progeny_master * sizeof(int));
            for (int i = 0; i < new_progeny_master; i++)
            {
                temp__progeny_distribution[i] = d(gen);
                cout << temp__progeny_distribution[i] << endl;
            }
            cudaMemcpy(cuda_progeny_distribution, temp__progeny_distribution, new_progeny_master * sizeof(int), cudaMemcpyHostToDevice);
            free(temp__progeny_distribution);
        }

        cout << "Simulating viral life cycle" << endl;
        cuda_progeny_life<<<tot_Blocks, tot_ThreadsperBlock>>>(new_progeny_master, ID_First_child_in_generation, cuda_progeny_distribution, cuda_progeny_values);
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            printf("CUDA Error 5: %s\n", cudaGetErrorString(err));
            // Possibly: exit(-1) if program cannot continue....
        }

        // parent sequences
        int **temp_cuda_parent_Sequences;
        cudaMallocManaged(&temp_cuda_parent_Sequences, genome_Size * total_viruses * sizeof(int));
        int **tmp_4 = (int **)malloc(total_viruses * sizeof(tmp_4[0]));
        for (int i = 0; i < total_viruses; i++)
        {
            cudaMalloc((void **)&tmp_4[i], genome_Size * sizeof(tmp_4[0][0]));
        }
        cudaMemcpy(temp_cuda_parent_Sequences, tmp_4, total_viruses * sizeof(int *), cudaMemcpyHostToDevice);

        free(tmp_4);

        // fill_old parents sequences
        // cuda_old_to_new_parent_sequences(int **old_sequences, int **new_parent_sequences, int number_of_old_parents, int genome_Size)
        cout << "Converting progeny sequences to parent sequences 1" << endl;
        cuda_old_to_new_parent_sequences<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_parent_Sequences, temp_cuda_parent_Sequences, parents, genome_Size);
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            printf("CUDA Error 6: %s\n", cudaGetErrorString(err));
            // Possibly: exit(-1) if program cannot continue....
        }

        // cuda_fill_new_parents(int number_of_old_parents, int number_of_new_progeny, int **progeny_Sequences, int **new_parent_sequences, int genome_Size)
        cout << "Converting progeny sequences to parent sequences 2" << endl;
        cuda_fill_new_parents<<<tot_Blocks, tot_ThreadsperBlock>>>(parents, new_progeny_master, cuda_progeny_Sequences, temp_cuda_parent_Sequences, genome_Size);
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            printf("CUDA Error 7: %s\n", cudaGetErrorString(err));
            // Possibly: exit(-1) if program cannot continue....
        }

        parents = parents + new_progeny_master;

        cudaFree(cuda_parent_Sequences);

        cudaMallocManaged(&cuda_parent_Sequences, genome_Size * total_viruses * sizeof(int));
        int **tmp_5 = (int **)malloc(total_viruses * sizeof(tmp_5[0]));
        for (int i = 0; i < total_viruses; i++)
        {
            cudaMalloc((void **)&tmp_5[i], genome_Size * sizeof(tmp_5[0][0]));
        }
        cudaMemcpy(cuda_parent_Sequences, tmp_5, total_viruses * sizeof(int *), cudaMemcpyHostToDevice);

        free(tmp_5);

        // cuda_Sequence_swap(int genome_Size, int total_virues, int **temp_sequences, int **parent_sequences)
        cout << "Converting progeny sequences to parent sequences 3" << endl;
        cuda_Sequence_swap<<<tot_Blocks, tot_ThreadsperBlock>>>(genome_Size, total_viruses, temp_cuda_parent_Sequences, cuda_parent_Sequences);
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess)
        {
            printf("CUDA Error 8: %s\n", cudaGetErrorString(err));
            // Possibly: exit(-1) if program cannot continue....
        }
        // parent sequence end

        cudaFree(temp_cuda_parent_Sequences);

        int *progeny_values = (int *)malloc(new_progeny_master * sizeof(int));
        cudaMemcpy(progeny_values, cuda_progeny_values, new_progeny_master * sizeof(int), cudaMemcpyDeviceToHost);

        // new parent child first
        int temp_new_progeny_master = 0;
        for (int progeny = 0; progeny < new_progeny_master; progeny++)
        {
            // cout << progeny_values[progeny] << endl;
            temp_new_progeny_master = temp_new_progeny_master + progeny_values[progeny];
        }

        int temp_total_viruses = total_viruses + temp_new_progeny_master;

        // int **temp_cuda_parent_child;

        // cudaMallocManaged(&temp_cuda_parent_child, temp_total_viruses * 2 * sizeof(int));
        // int **tmp_6 = (int **)malloc(temp_total_viruses * sizeof(tmp_6[0]));
        // for (size_t i = 0; i < temp_total_viruses; i++)
        // {
        //     cudaMalloc((void **)&tmp_6[i], 2 * sizeof(tmp_6[0][0]));
        // }
        // cudaMemcpy(temp_cuda_parent_child, tmp_6, temp_total_viruses * sizeof(int *), cudaMemcpyHostToDevice);

        // free(tmp_6);

        // cuda_old_to_new_parent_child(int **old_parent_child, int **new_parent_child, int old_progeny_count)
        cout << "Configuring parent child relationships" << endl;
        // cuda_old_to_new_parent_child<<<tot_Blocks, tot_ThreadsperBlock>>>(cuda_parent_child, temp_cuda_parent_child, total_viruses);
        // cudaDeviceSynchronize();

        // err = cudaGetLastError();
        // if (err != cudaSuccess)
        // {
        //     printf("CUDA Error 9: %s\n", cudaGetErrorString(err));
        //     // Possibly: exit(-1) if program cannot continue....
        // }

        // int **temp_parent_child = (int **)malloc(total_viruses * sizeof(int *));

        // // columns
        // for (int i = 0; i < total_viruses; i++)
        // {
        //     temp_parent_child[i] = (int *)malloc(2 * sizeof(int));
        // }

        // for (int i = 0; i < total_viruses; i++)
        // {
        //     temp_parent_child[i][0] = parent_child[i][0];
        //     temp_parent_child[i][1] = parent_child[i][1];
        // }

        free(parent_child);
        parent_child = (int **)malloc(temp_total_viruses * sizeof(int *));
        for (int i = 0; i < temp_total_viruses; i++)
        {
            parent_child[i] = (int *)malloc(2 * sizeof(int));
        }

        for (int i = 0; i < total_viruses; i++)
        {
            cudaMemcpy(parent_child[i], cuda_parent_child[i], 2 * sizeof(cuda_parent_child[0][0]), cudaMemcpyDeviceToHost);
        }

        // free(temp_parent_child);

        int child_ID = 0;
        int start_parent_child = total_viruses;
        for (int progeny = 0; progeny < new_progeny_master; progeny++)
        {
            int children = progeny_values[progeny];
            int progeny_ID = ID_First_child_in_generation - new_progeny_master + progeny;

            for (int iterate = 0; iterate < children; iterate++)
            {
                parent_child[start_parent_child][0] = progeny_ID;
                parent_child[start_parent_child][1] = ID_First_child_in_generation + child_ID;
                start_parent_child++;
                child_ID++;
            }
        }

        cudaFree(cuda_parent_child);
        // cudaFree(temp_cuda_parent_child);

        cudaMallocManaged(&cuda_parent_child, temp_total_viruses * 2 * sizeof(int));
        int **tmp_7 = (int **)malloc(temp_total_viruses * sizeof(tmp_7[0]));
        for (size_t i = 0; i < temp_total_viruses; i++)
        {
            cudaMalloc((void **)&tmp_7[i], 2 * sizeof(tmp_7[0][0]));
        }
        cudaMemcpy(cuda_parent_child, tmp_7, temp_total_viruses * sizeof(int *), cudaMemcpyHostToDevice);

        for (size_t i = 0; i < temp_total_viruses; i++)
        {
            cudaMemcpy(tmp_7[i], parent_child[i], 2 * sizeof(cuda_parent_child[0][0]), cudaMemcpyHostToDevice);
        }

        free(tmp_7);

        cudaFree(cuda_progeny_Sequences);
        cudaFree(cuda_progeny_values);
        cudaFree(cuda_progeny_distribution);

        new_progeny_master = temp_new_progeny_master;
        total_viruses = temp_total_viruses;

        cout << endl;
    }

    for (size_t row = 0; row < total_viruses; row++)
    {
        for (size_t column = 0; column < 2; column++)
        {
            cout << parent_child[row][column] << "\t";
        }
        cout << "\n";
    }

    cudaFree(cuda_parent_Sequences);
}