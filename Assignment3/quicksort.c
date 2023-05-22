#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

double find_median (double *data, int len);
double find_mean (double *data, int len);
int compare (const void *a, const void *b);
double* merge(double* a, double* b, int length_a, int length_b);
double pivot_selection(int strategy, int group_rank, int group_size, double *local_array, int local_length, MPI_Comm comm);

int main(int argc, char **argv)
{
    // check if the number of arguments is correct
    if (argc != 4) {
        printf("Usage: %s <input file> <output file> <pivot strategies>\n", argv[0]);
        return -1;
    }

    char *inputfile = argv[1];
    char *outputfile = argv[2];
    int pivot_strategy = atoi(argv[3]);

    FILE *fp;

    int rank, size, N;
    double *data, *local_array;
    double *final_array;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // world rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request request_final;
    MPI_Status status_final;
    MPI_Request request1;
    MPI_Request request2;

    // rank 0 read input file
    if (rank == 0) {
        fp = fopen(inputfile, "r");
        if (fp == NULL) {
            perror("Error opening input file");
            return -1;
        }

        // read the number of elements to be sorted
        if (EOF == fscanf(fp, "%d", &N)) {
            perror("Error when reading number of elements from input file");
            return -1;
        }

        data = (double *)malloc(N * sizeof(double));
        // read the data from input file
        for (int i = 0; i < N; i++) {
            if (EOF == fscanf(fp, "%lf", &data[i])) {
                perror("Error when reading data from input file");
                return -1;
            }
        }
    }


    double max_time, sort_time;
    // start timer
    double start_time = MPI_Wtime();

    // Algorithm step 1
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast N to all processes
    int local_length = N / size;  // not necessarily divisible
    local_array = (double *)malloc(local_length * sizeof(double));
    // scatter the array to all processes
    MPI_Scatter(data, local_length, MPI_DOUBLE, local_array, local_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    const int last_length = local_length + N % size;
    if (last_length > local_length){
        if (rank == 0){
            MPI_Send(data+(N-last_length), last_length, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD);
        }
        if (rank == size-1){
            free(local_array);
            local_length = last_length;
            local_array = (double *)malloc(local_length * sizeof(double));
            MPI_Recv(local_array, local_length, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);      
        }
    }

    // Algorithm step 2
    double time_for_oneprocessor = MPI_Wtime();
    qsort(local_array, local_length, sizeof(double), compare); // sort the local array
    double runtime_for_oneprocessor = MPI_Wtime() - time_for_oneprocessor;
    
    // Algorithm step 3 and 4
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm newcomm;

    // At first this is world rank and size
    int group_rank = rank;
    int group_size = size;
    double pivot;
    
    if (size > 1) {
        double sort_start = MPI_Wtime();

        while (group_size > 1) {
            // Algorithm step 3.1
            // Select pivot element within each process set
            pivot = pivot_selection(pivot_strategy, group_rank, group_size, local_array, local_length, comm);
            MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);
            
            // Algorithm step 3.2
            // locally divide the data into two sets
            // find the smallest integer 'division_index' such that local_array[division_index] >= pivot
            // in case all elements in local_array are smaller than pivot, the division_index will be local_length
            int division_index = local_length / 2;
            for (;(division_index > 0) && local_array[division_index]>=pivot; division_index--){}
            for (;(division_index < local_length) && local_array[division_index]<pivot; division_index++){}
            // there is an option to specify two groups in advance and divide the local_array into two sets called sent_data and kept_data
            // dont need to malloc memory for the low and high part of the array
            int length_small = division_index;
            int length_large = local_length - division_index;

            // Algorithm step 3.3
            // exchange data pairwise
            // we explicitly split the processes in the final step
            // we split the processes into two groups from the middle point
            int pair_big, pair_small;
            int received_low, received_high;
            if (group_rank < group_size / 2){
                pair_big = group_rank + group_size / 2;
            }
            else {
                pair_small = group_rank - group_size / 2;
            }
            double* recv_array;
            
            if (group_rank < group_size/2){
                // for group_rank smaller then the half of the group_size, they will send their higher part data to the pair_big.
                // use MPI_Get_count to check the amount of elements we sent
                // actually we dont need to specify how many elements we are gonna send
                MPI_Isend(local_array+length_small, length_large, MPI_DOUBLE, pair_big, 0, comm, &request1);
                // probe for a incoming message
                MPI_Probe(pair_big, 1, comm, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &received_low);
                recv_array = (double *)malloc(received_low * sizeof(double));
                // MPI_Recv(recv_array, received_low, MPI_DOUBLE, pair_big, 1, comm, &status);
                MPI_Irecv(recv_array, received_low, MPI_DOUBLE, pair_big, 1, comm, &request2);
            } else {
                MPI_Isend(local_array, length_small, MPI_DOUBLE, pair_small, 1, comm, &request2);
                // probe for a incoming message
                MPI_Probe(pair_small, 0, comm, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &received_high);
                recv_array = (double *)malloc(received_high * sizeof(double));
                // MPI_Recv(recv_array, received_high, MPI_DOUBLE, pair_small, 0, comm, &status);
                MPI_Irecv(recv_array, received_high, MPI_DOUBLE, pair_small, 0, comm, &request1);
            }

            // Algorithm step 3.4
            // locally merge two sets of numbers
            if (group_rank < group_size/2){ 
                local_length = length_small + received_low;  // update the local_length
                local_array = merge(local_array, recv_array, length_small, received_low);
            } else {
                local_length = length_large + received_high; // update the local_length
                local_array = merge(local_array + length_small, recv_array, length_large, received_high);
            }
            
            MPI_Comm_split(comm, group_rank<group_size/2, group_rank, &newcomm);
            MPI_Comm_rank(newcomm, &group_rank);
            MPI_Comm_size(newcomm, &group_size);
            comm = newcomm;
        }


        // Global merge
        // if rank != 0, the process sends its local array to process 0.
        if (rank != 0) {
            MPI_Isend(local_array, local_length, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_final);
        }
        else {
            // if rank == 0, receive message from other processes and merge it together. 
            final_array = (double *)malloc(N * sizeof(double));
            int buffer_size = 0;
            int current_loc = 0;
            MPI_Isend(local_array, local_length, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_final);
            for (int i=0; i<size; i++) {
                MPI_Probe(i, i, MPI_COMM_WORLD, &status_final);
                MPI_Get_count(&status_final, MPI_DOUBLE, &buffer_size);
                MPI_Recv(final_array + current_loc, buffer_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status_final);
                current_loc += buffer_size;
            }

        }
        

        double run_time = MPI_Wtime() - start_time;
        MPI_Reduce(&run_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    } 
    
    else {
        final_array = local_array;
    }

    if (size > 1 && rank == 0) printf("%lf\n", max_time);
    if (size == 1 && rank == 0) printf("%lf\n", runtime_for_oneprocessor);

    // write the result and print the time
    if (rank == 0) {
        // check if the result is correct
        int check = 1;
        for (int i=0; i<N-1; i++) {
            if (final_array[i] > final_array[i+1]) {
                printf("Final_array[%d] is %lf, final_array[%d] is %lf\n", i, final_array[i], i+1, final_array[i+1]);
                check = 0;
            }
        }
        if (check == 1) {
            printf("The result is correct!\n");
        }
        
        free(data);

        // need to write data to the output file.
        if (size>1) free(final_array);
    }
    
    MPI_Finalize();
    return 0;
}

double pivot_selection(int strategy, int group_rank, int group_size, double *local_array, int local_length, MPI_Comm comm) {
    double pivot;
    switch(strategy){
        // Select the median in one processor in each group of processors.
        case 1:
            if (group_rank == 0){
                pivot = find_median(local_array, local_length);
            }
        break;
        // Select the median of all medians in each processor group.
        case 2:
        {
            double local_pivot = find_median(local_array, local_length);
            double all_pivots[group_size];
            MPI_Gather(&local_pivot, 1, MPI_DOUBLE, all_pivots, 1, MPI_DOUBLE, 0, comm);
            if (group_rank == 0){
                qsort(all_pivots, group_size, sizeof(double), compare);
                pivot = find_median(all_pivots, group_size);
            }
        }
        break;
        // Select the mean value of all medians in each processor group.
        case 3: 
        {
            double local_pivot = find_median(local_array, local_length);
            double all_pivots[group_size];
            MPI_Gather(&local_pivot, 1, MPI_DOUBLE, all_pivots, 1, MPI_DOUBLE, 0, comm);
            if (group_rank == 0) pivot = find_mean(all_pivots, group_size);
        }
        break;
    }
    return pivot;
}

double* merge(double* a, double* b, int length_a, int length_b){
    int len = length_a + length_b;
    double* array = (double *)malloc(len * sizeof(double));
    for (int i = 0, j = 0, k = 0; k < len;) {
        if (i < length_a && (j >= length_b || a[i] < b[j])) {
            array[k++] = a[i++];
        } else {
            array[k++] = b[j++];
        }
    }
    return array;
}

/* define the compare function */
int compare (const void *a, const void *b){
    // the sort is in ascending order
    return ( *(double*)a - *(double*)b );
}

/* define pivot strategies */ 
double find_median (double *data, int len)
{
    // is there a possibility that length equals to 0?

    if (len == 1) {
        return data[0];
    }
    // if length is odd then return the middle element
    if (len % 2 == 1) {
        return data[len / 2];
    }
    // if length is even then return the average of the two middle elements
    return (data[len / 2 - 1] + data[len / 2]) / 2;
}

double find_mean (double *data, int len)
{
    if (len == 0) {
        return 0;
    }
    double sum = 0, mean = 0;
    for (int i = 0; i < len; i++) {
        sum += data[i];
    }
    mean = sum / (double)len;
    return mean;
}