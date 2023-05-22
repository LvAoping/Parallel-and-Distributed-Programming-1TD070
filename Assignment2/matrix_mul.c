#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {

    if (argc != 3) {
        printf("Usage: %s <inputfile> <outputfile>\n", argv[0]);
        return 1;
    }

    char *inputfile = argv[1];
    char *outputfile = argv[2];

    FILE *fp;

    int rank, size, n;
    double *A, *B, *C, *local_A, *local_B, *local_C, *local_B_recv;

    MPI_Status status;
    MPI_Datatype coltype, col;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // rank 0 read input file
    if (rank == 0) {
        fp = fopen(inputfile, "r");
        if (fp == NULL) {
            perror("Error opening input file");
            return 1;
        }
        
        // read matrix size
        fscanf(fp, "%d", &n);

        A = (double *)malloc(n * n * sizeof(double));
        B = (double *)malloc(n * n * sizeof(double));
        C = (double *)malloc(n * n * sizeof(double));

        // read matrix A
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fscanf(fp, "%lf", &A[i * n + j]);
            }
        }

        // read matrix B
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fscanf(fp, "%lf", &B[i * n + j]);
            }
        }

        MPI_Type_vector(n, 1, n, MPI_DOUBLE, &coltype);
        MPI_Type_commit(&coltype);
        MPI_Type_create_resized(coltype, 0, 1*sizeof(double), &col);
        MPI_Type_commit(&col);
    }

    /* everyone calls bcast, data is taken from root and ends up in everyone's n */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // assume n is divisible by number of processes
    int chunk_size = n / size;

    // allocate memory for local matrix
    local_A = (double *)malloc(chunk_size * n * sizeof(double));
    local_B = (double *)malloc(chunk_size * n * sizeof(double));
    local_B_recv = (double *)malloc(chunk_size * n * sizeof(double));
    local_C = (double *)calloc(chunk_size * n, sizeof(double));

    // wait for all processes to allocate memory
    MPI_Barrier(MPI_COMM_WORLD);

    // start timer
    double starttime = MPI_Wtime();

    // scatter matrix A
    MPI_Scatter(A, chunk_size * n, MPI_DOUBLE, local_A, chunk_size * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // scatter matrix B
    MPI_Scatter(B, chunk_size, col, local_B_recv, chunk_size * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // define how to communicate
    // current process send to right process
    int right = (rank + 1) % size;
    // current process receive from left process
    int left = (rank - 1 + size) % size;
    for (int step = 0; step < size; step++) {
        int index = (rank - step + size) % size;

        // calculate local C
        for (int i = 0; i < chunk_size; i++) {
            for (int j = 0; j < chunk_size; j++) {
                for (int k = 0; k < n; k++) {
                    local_C[i*n+(index*chunk_size+j)] += local_A[i*n+k] * local_B_recv[j*n+k];
                }
            }
        }

        // send local B to right process
        // MPI_Sendrecv_replace(local_B_recv, chunk_size * n, MPI_DOUBLE, right, 0, left, 0, MPI_COMM_WORLD, &status);
        for (int i = 0; i < chunk_size * n; i++) {
            local_B[i] = local_B_recv[i];
        }

        MPI_Sendrecv(local_B, chunk_size * n, MPI_DOUBLE, right, step, local_B_recv, chunk_size * n, MPI_DOUBLE, left, step, MPI_COMM_WORLD, &status);
    }

    // gather local C
    MPI_Gather(local_C, chunk_size * n, MPI_DOUBLE, C, chunk_size * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double maxtime;
    double runtime = MPI_Wtime() - starttime;
    MPI_Reduce(&runtime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // print runtime
    if (rank == 0) {
        printf("Runtime: %f\n", maxtime);
    }

    // print output file
    if (rank == 0) {
        FILE *output = NULL;
        output = fopen(outputfile, "w+");
        if (output == NULL) {
            perror("Error opening output file");
            return 1;
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fprintf(output, "%lf ", C[i*n+j]);
            }
        }

        fclose(fp);
        fclose(output);
        free(A);
        free(B);
        free(C);
        MPI_Type_free(&coltype);
        MPI_Type_free(&col);
    }

    free(local_A);
    free(local_B);
    free(local_B_recv);
    free(local_C);

    MPI_Finalize();
    return 0;
}