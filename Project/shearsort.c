#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

// Comparator functions for qsort
int asc_comparator(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

int desc_comparator(const void* a, const void* b) {
    return (*(int*)b - *(int*)a);
}

// Function for sorting rows
void sort_row(int* row, int n, int direction) {
    if (direction == 1)
        qsort(row, n, sizeof(int), asc_comparator);
    else
        qsort(row, n, sizeof(int), desc_comparator);
}

// Function for sorting columns
void sort_column(int* column, int n) {
    qsort(column, n, sizeof(int), asc_comparator);
}

int main(int argc, char** argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = sqrt(size);  // assuming size is a perfect square
    int* matrix = (int*)malloc(n * n * sizeof(int));

    // Assume the matrix is populated with data here




    // Distribute rows among processes
    int* row = (int*)malloc(n * sizeof(int));
    MPI_Scatter(matrix, n, MPI_INT, row, n, MPI_INT, 0, MPI_COMM_WORLD);

    // Sort rows in alternating order
    sort_row(row, n, rank % 2 == 0 ? 1 : -1);

    // Gather rows back to root process
    MPI_Gather(row, n, MPI_INT, matrix, n, MPI_INT, 0, MPI_COMM_WORLD);

    // Distribute columns among processes
    int* column = (int*)malloc(n * sizeof(int));
    if (rank == 0) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                column[j] = matrix[j * n + i];
            }
            MPI_Send(column, n, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Recv(column, n, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Sort columns
    sort_column(column, n);

    // Gather columns back to root process
    if (rank == 0) {
        for (int i = 0; i < n; ++i) {
            MPI_Recv(column, n, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < n; ++j) {
                matrix[j * n + i] = column[j];
            }
        }
    }
    else {
        MPI_Send(column, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    free(row);
    free(column);
    free(matrix);

    MPI_Finalize();
    return 0;
}


// https://github.com/pottu/shearsort/tree/main
