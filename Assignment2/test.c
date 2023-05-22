#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

void multiply_blocks(int *A, int *B, int *C, int N, int block_size) {
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            for (int k = 0; k < block_size; k++) {
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int N = 4;
    int block_size = sqrt(world_size);
    int submatrix_size = N / block_size;

    srand(time(0));

    int *A = (int *) malloc(N * N * sizeof(int));
    int *B = (int *) malloc(N * N * sizeof(int));
    int *C = (int *) malloc(N * N * sizeof(int));

    if (world_rank == 0) {
        for (int i = 0; i < N * N; i++) {
            A[i] = rand() % 100;
            B[i] = rand() % 100;
            C[i] = 0;
        }
    }

    int *local_A = (int *) malloc(submatrix_size * submatrix_size * sizeof(int));
    int *local_B = (int *) malloc(submatrix_size * submatrix_size * sizeof(int));
    int *local_C = (int *) malloc(submatrix_size * submatrix_size * sizeof(int));

    MPI_Datatype blocktype;
    MPI_Type_vector(submatrix_size, submatrix_size, N, MPI_INT, &blocktype);
    MPI_Type_commit(&blocktype);

    int sizes[2] = {N, N};
    int subsizes[2] = {submatrix_size, submatrix_size};
    int starts[2] = {0, 0};
    MPI_Datatype subarrtype;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarrtype);
    MPI_Type_commit(&subarrtype);

    int *sendcounts = NULL;
    int *displs = NULL;

    if (world_rank == 0) {
        sendcounts = (int *) malloc(world_size * sizeof(int));
        displs = (int *) malloc(world_size * sizeof(int));

        for (int i = 0; i < world_size; i++) {
            sendcounts[i] = 1;
            displs[i] = i;
        }
    }

    MPI_Scatterv(A, sendcounts, displs, subarrtype, local_A, submatrix_size * submatrix_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, sendcounts, displs, subarrtype, local_B, submatrix_size * submatrix_size, MPI_INT, 0, MPI_COMM_WORLD);

    multiply_blocks(local_A, local_B, local_C, N, submatrix_size);

        MPI_Gatherv(local_C, submatrix_size * submatrix_size, MPI_INT, C, sendcounts, displs, subarrtype, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        printf("Matrix C:\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%d ", C[i * N + j]);
            }
            printf("\n");
        }

        free(sendcounts);
        free(displs);
    }

    free(A);
    free(B);
    free(C);
    free(local_A);
    free(local_B);
    free(local_C);

    MPI_Type_free(&blocktype);
    MPI_Type_free(&subarrtype);

    MPI_Finalize();
    return 0;
}

