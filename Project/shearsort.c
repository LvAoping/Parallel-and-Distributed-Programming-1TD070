#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>


/* Helper Functions */

int ascending(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

int descending(const void* a, const void* b) {
    return (*(int*)b - *(int*)a);
}

void sort_row(int width, int height, int *matrix)
{
  for (int row = 0; row < height; row++) {
    if ((row & 1) == 0) {
      // Even row, sort ascending
      qsort(&matrix[row*width], width, sizeof(int), ascending);
    } else {
      // Odd row, sort descending
      qsort(&matrix[row*width], width, sizeof(int), descending);
    }
  }
}

void sort_column(int width, int height, int* matrix) {
  for (int col = 0; col < width; col++) {
    int* column_data = (int*) malloc(height * sizeof(int));
    if (!column_data) {
        perror("Couldn't allocate memory for the column data");
    }
    // Extract the column data
    for (int i = 0; i < height; i++)
        column_data[i] = matrix[i * width + col];
    // Sort the column data
    qsort(column_data, height, sizeof(int), ascending);
    // Store back the sorted column data
    for (int i = 0; i < height; i++)
        matrix[i * width + col] = column_data[i];
    free(column_data);
  }
}

void merge(const int *array1, int size1, const int *array2, int size2, int *result, bool ascending) 
{
  int i = 0; // index into array1
  int j = 0; // index into array2
  int k = 0; // index into result
  while (i < size1 && j < size2) {
    if (ascending ? array1[i] < array2[j] : array1[i] >= array2[j]) {
      result[k++] = array1[i++];
    } else {
      result[k++] = array2[j++];
    }
  }
  if (i == size1) {
    while (j < size2) {
      result[k++] = array2[j++];
    }
  } else {
    while (i < size1) {
      result[k++] = array1[i++];
    }
  }
}

void block_exchange(int target, int rank, int width, int height, int *block, 
     int *target_block, MPI_Datatype row_send_type, MPI_Datatype row_recv_type)
{
  
  int offset = ((rank & 1) == 0) ? width : 0;

  // Exchange rows between the current and partner process
  MPI_Sendrecv( block + offset, 1, row_send_type, target, rank, target_block, 1, row_recv_type, target, target, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  int partner_row = 0;
  int merged[width*2];

  for (int row = ((rank & 1) == 0) ? 0 : 1, partner_row = 0; row < height; row += 2, partner_row++) {
      size_t copysize = width * sizeof(int);

      // Merge the two sequences
      bool ascending = (rank & 1) == 0;
      merge(&block[row * width], width, &target_block[partner_row * width], width, merged, ascending);

      // Copy the merged sequence back into the block
      if (rank < target) {
          memcpy(&block[row * width], merged, copysize);
          memcpy(&target_block[partner_row * width], &merged[width], copysize);
      } else {
          memcpy(&block[row * width], &merged[width], copysize);
          memcpy(&target_block[partner_row * width], merged, copysize);
      }
  }

  // Replace the send data with the received data
  MPI_Sendrecv_replace( target_block, width * ceil((double)height/2), MPI_INT, target, rank, target, target, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  partner_row = 0;
  for (int row = ((rank & 1) == 0) ? 1 : 0; row < height; row += 2) {
    memcpy(&block[row * width], &target_block[partner_row * width], width * sizeof(int));
    partner_row += 1;
  }
}

void odd_even_transpositon_sort(int width, int height, int *block, int *target_block, int rank, 
     int size, MPI_Datatype row_send_type, MPI_Datatype row_recv_type)
{
  for (int i = 0; i < size; i++) {
    int target = -1;
    if ((i & 1) == 0) {
      target = ((rank & 1) == 0) ? rank + 1 : rank - 1;
    } else {
      target = ((rank & 1) == 0) ? rank - 1 : rank + 1;
    }
    if (target >= 0 && target < size) {
      block_exchange(target, rank, width, height, block, target_block, row_send_type, row_recv_type);
    }
  }
}


/* Input and Output */

int* read_input(int* out_n, char* input_file)
{
  FILE* file = fopen(input_file, "r");
  if (!file) {
    perror("Couldn't open input file");
    return NULL;
  }
  if (fscanf(file, "%d", out_n) != 1) {
    perror("Couldn't read element count from input file");
    fclose(file);
    return NULL;
  }
  int n = *out_n;
  int* matrix = (int*) malloc(n * n * sizeof(int));
  if (!matrix) {
    perror("Couldn't allocate memory for the matrix");
    fclose(file);
    return NULL;
  }
  for (int i = 0; i < n*n; ++i) {
    if (feof(file) || fscanf(file, "%d", &matrix[i]) != 1) {
      perror("Couldn't read elements from input file");
      free(matrix);
      fclose(file);
      return NULL;
    }
  }
  fclose(file);
  return matrix;
}


void print_matrix(int n, int* matrix)
{
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      printf("%d ", matrix[row*n + col]);
    }
    printf("\n");
  }
}

void write_output(int n, int* matrix, char* output_file)
{
  if (!matrix || !output_file) {
      perror("Null matrix or output file name provided");
  }
  FILE *file = fopen(output_file, "w");
  if (!file) {
      perror("Couldn't open output file");
  }
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      if (fprintf(file, "%d ", matrix[row*n + col]) < 0) {
        perror("Failed to write to output file");
        fclose(file);
      }
      }
      if (fprintf(file, "\n") < 0) {
        perror("Failed to write to output file");
        fclose(file);
      }
  }
  if (fclose(file) == EOF) {
      perror("Failed to close output file");
  }
}

/*  Check Functions  */

bool check_integrity(int n, int *matrix, char *input_file) 
{
  int *initial = read_input(&n, input_file);
  int *result  = (int*) malloc(n * n * sizeof(int));
  memcpy(result, matrix, n * n * sizeof(int));
  // Sort both matrices and compare them to check if any element missing
  qsort(initial, n * n, sizeof(int), ascending); 
  qsort(result, n * n, sizeof(int), ascending);
  for (int i = 0; i < n * n; i++) {
    if (result[i] != initial[i]) {
      free(initial);
      free(result);
      return false;
    }
  }
  free(initial);
  free(result);
  return true;
}


bool check_sorted(int n, int *matrix)
{
  int prev = matrix[0];
  
  for (int row = 0; row < n; row++) 
  {
    int direction = (row & 1) == 0 ? 1 : -1;  // Determine direction based on row parity
    int col = direction > 0 ? 0 : n - 1;  // Calculate starting column based on direction
    
    while(col >= 0 && col < n)  // Traverse through the columns based on the direction
    {
      if (prev > matrix[row * n + col]) 
      {
        return false;
      }
      prev = matrix[row * n + col];
      col += direction; // Move column index according to direction
    }
  }
  
  return true;
}

// Main checker function
bool checker(int n, int* matrix, char *input_file) 
{
  return check_sorted(n, matrix) && check_integrity(n, matrix, input_file);
}

/* Main Function */
int main(int argc, char **argv) 
{
  MPI_Init(&argc, &argv);

  // Find rank and size
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Argument Check
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <input_file> <output_file> <output_type>\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(EXIT_FAILURE);
  }

  char* input_file = argv[1];
  char* output_file = argv[2];
  bool check_solution = true; // Set to true to check if solution is correct
  int output_type = atoi(argv[3]); // 0 for no output, 1 for write to file, 2 for print to stdout


  // Initialize matrix
  int n; // Size of matrix
  int* matrix = NULL; // Input matrix
  if (rank == 0) {
    matrix = read_input(&n, input_file);
    if (n % size != 0) {
      perror("Matrix size should be divisible by the number of processes.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(EXIT_FAILURE);
    }
  }

  // Broadcast the size of input to all processes
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Width and height of matrix blocks
  const int width = n / size;
  const int height = n;

  // Types for sending and receiving columns of matrix and blocks
  MPI_Datatype col_matrix_type, col_block_type;
  MPI_Type_vector(height, 1, height, MPI_INT, &col_matrix_type);
  MPI_Type_create_resized(col_matrix_type, 0, sizeof(int), &col_matrix_type);
  MPI_Type_commit(&col_matrix_type);
  MPI_Type_vector(height, 1, width, MPI_INT, &col_block_type);
  MPI_Type_create_resized(col_block_type, 0, sizeof(int), &col_block_type);
  MPI_Type_commit(&col_block_type);

  // Types for sending and receiving row blocks
  MPI_Datatype row_send_type, row_recv_type;
  MPI_Type_vector(ceil((double)height/2), width, width*2, MPI_INT, &row_send_type);
  MPI_Type_vector(ceil((double)height/2), width, width, MPI_INT, &row_recv_type);
  MPI_Type_commit(&row_send_type);
  MPI_Type_commit(&row_recv_type);

  // Each process' individual block of columns
  int *block = malloc(width * height * sizeof(*block));
  // Used for odd-even transposition sort.
  // Allocating here allows reusage.
  int* target_block = malloc(width * ceil((double)height/2) * sizeof(*target_block));

  // Start timer
  const double start = MPI_Wtime();

  // Scatter matrix to all processes
  MPI_Scatter(matrix, width, col_matrix_type, block, width, col_block_type, 0, MPI_COMM_WORLD);

  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Sort rows locally.
    sort_row(width, height, block);
    // Sort rows globally.
    odd_even_transpositon_sort( width, height, block, target_block, rank, size, row_send_type, row_recv_type);
    
    if (step < num_steps-1) {
      // Sort columns locally
      sort_column(width, height, block);
    }
  }
  
  // Gather blocks in temporary array
  int* matrix_tmp = NULL;
  if (rank == 0) {
    matrix_tmp = malloc(n*n * sizeof(*matrix_tmp));
  }

  // Note the gathered blocks are not in the correct order
  MPI_Gather(block, width*height, MPI_INT, matrix_tmp, width*height, MPI_INT, 0, MPI_COMM_WORLD);

  // Rearange the blocks to correct ordering
  if (rank == 0) {
    for (int rank = 0; rank < size; rank++) {
      for (int row = 0; row < height; row++) {
        memcpy(&matrix[(rank * width) + (row * n)], &matrix_tmp[(rank * width * height) + (row * width)], width * sizeof(int));
      }
    }
  }


  // Stop timer
  double my_execution_time = MPI_Wtime() - start;
	double maxtime;
	MPI_Reduce(&my_execution_time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Print execution time
  if (rank == 0) {
    printf("%lf\n", maxtime);
  }

  if (rank == 0) {
    if (check_solution) {
      // Check if correctly sorted
      if (checker(n, matrix, input_file)) {
        printf("Correct!\n");
      } else {
        printf("Incorrect!\n");
      }
    }
    if (output_type != 0) {
      // Print sorted matrix
      if (output_type == 1) {
        write_output(n, matrix, output_file);
      } 
      else if (output_type == 2) {
        print_matrix(n, matrix);
      } 
      else {
        fprintf(stderr, "Usage: <input-file> <output-file> <flag, 0 for no output, 1 for write to file, 2 for print to stdout>\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(EXIT_FAILURE);
      }
    } 
  }

  // Free memory and MPI type then finalize MPI
  free(matrix);
  free(matrix_tmp);
  free(block);
  free(target_block);
  MPI_Type_free(&col_matrix_type);
  MPI_Type_free(&col_block_type);
  MPI_Type_free(&row_send_type);
  MPI_Type_free(&row_recv_type);
  MPI_Finalize();
  return 0;
}