#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

// Function to exchange columns between root and other processes
void exchange_column(int* matrix, int* column, int n, int rank) {
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
}

// Function to populate the matrix
void populate_matrix(int* matrix, int n) {
    srand(time(NULL));
    for (int i = 0; i < n * n; ++i) {
        matrix[i] = rand() % 100; // random values between 0 and 99
    }
}

// Function to print the matrix
void print_matrix(int* matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%d ", matrix[i * n + j]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if the number of processes is a perfect square
    int n = sqrt(size);
    if (n * n != size) {
        if (rank == 0) printf("Number of MPI processes must be a perfect square.\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int* matrix = (int*)malloc(n * n * sizeof(int));

    // Populate the matrix with data
    if (rank == 0) {
        populate_matrix(matrix, n);

        // Print the initial matrix
        printf("Initial matrix:\n");
        print_matrix(matrix, n);
    }

    // Distribute rows among processes
    int* row = (int*)malloc(n * sizeof(int));
    int* column = (int*)malloc(n * sizeof(int));
    
    for (int iter = 0; iter < n; iter++) {
        MPI_Scatter(matrix, n, MPI_INT, row, n, MPI_INT, 0, MPI_COMM_WORLD);

        // Sort rows in alternating order
        sort_row(row, n, (iter + rank) % 2 == 0 ? 1 : -1);

        // Gather rows back to root process
        MPI_Gather(row, n, MPI_INT, matrix, n, MPI_INT, 0, MPI_COMM_WORLD);

        // Distribute columns among processes
        distribute_columns(matrix, column, n, rank);

        // Sort columns
        sort_column(column, n);

        // Gather columns back to root process
        exchange_column(matrix, column, n, rank);
    }

    if (rank == 0) {
        // Print the sorted matrix
        printf("Sorted matrix:\n");
        print_matrix(matrix, n);
    }

    free(row);
    free(column);
    free(matrix);

    MPI_Finalize();
    return 0;
}
