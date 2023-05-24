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
      qsort(&matrix[row*width], width, sizeof(int), ascending);
    } else {
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

void merge(const int *v1, int n1, const int *v2, int n2, int *result, bool ascending) 
{
  int i = 0;
  int j = 0;
  int k = 0;
  while (i < n1 && j < n2) {
    if (ascending ? v1[i] < v2[j] : v1[i] >= v2[j]) {
      result[k++] = v1[i++];
    } else {
      result[k++] = v2[j++];
    }
  }
  if (i == n1) {
    while (j < n2) {
      result[k++] = v2[j++];
    }
  } else {
    while (i < n1) {
      result[k++] = v1[i++];
    }
  }
}

void exchange_and_merge(int partner, int rank, int width, int height, int *block, 
     int *partner_block, MPI_Datatype row_send_type, MPI_Datatype row_recv_type)
{
  int offset = ((rank & 1) == 0) ? width : 0;

  MPI_Sendrecv( block + offset, 1, row_send_type, partner, rank, partner_block, 1, row_recv_type, partner, partner, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  int partner_row = 0;
  int merged[width*2];

  for (int row = ((rank & 1) == 0) ? 0 : 1, partner_row = 0; row < height; row += 2, partner_row++) {
      size_t copysize = width * sizeof(int);

      // Merge two sorted sequences
      bool ascending = (rank & 1) == 0;
      merge(&block[row * width], width, &partner_block[partner_row * width], width, merged, ascending);

      // Copy the merged sequences back into the correct position
      if (rank < partner) {
          memcpy(&block[row * width], merged, copysize);
          memcpy(&partner_block[partner_row * width], &merged[width], copysize);
      } else {
          memcpy(&block[row * width], &merged[width], copysize);
          memcpy(&partner_block[partner_row * width], merged, copysize);
      }
  }

  MPI_Sendrecv_replace( partner_block, width * ceil((double)height/2), MPI_INT, partner, rank, partner, partner, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  partner_row = 0;
  for (int row = ((rank & 1) == 0) ? 1 : 0; row < height; row += 2) {
    memcpy(&block[row * width], &partner_block[partner_row * width], width * sizeof(int));
    partner_row += 1;
  }
}

void parity_sort(int width, int height, int *block, int *partner_block, int rank, 
     int size, MPI_Datatype row_send_type, MPI_Datatype row_recv_type)
{
  for (int i = 0; i < size; i++) {
    int partner = -1;
    if ((i & 1) == 0) {
      partner = ((rank & 1) == 0) ? rank + 1 : rank - 1;
    } else {
      partner = ((rank & 1) == 0) ? rank - 1 : rank + 1;
    }
    if (partner >= 0 && partner < size) {
      exchange_and_merge(partner, rank, width, height, block, partner_block, row_send_type, row_recv_type);
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

bool check_same_elements(int n, int *matrix, char *input_file) 
{
  int *initial = read_input(&n, input_file);
  int *result  = (int*) malloc(n * n * sizeof(int));
  memcpy(result, matrix, n * n * sizeof(int));
  qsort(initial, n * n, sizeof(int), ascending); // Sort both matrices
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
  return check_sorted(n, matrix) && check_same_elements(n, matrix, input_file);
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

  // Broadcast matrix size to all PEs.
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Width and height of matrix blocks.
  const int width = n / size;
  const int height = n;

  // Types for sending and receiving columns.
  MPI_Datatype col_matrix_type, col_block_type;

  MPI_Type_vector(height, 1, height, MPI_INT, &col_matrix_type);
  MPI_Type_create_resized(col_matrix_type, 0, sizeof(int), &col_matrix_type);
  MPI_Type_commit(&col_matrix_type);
  MPI_Type_vector(height, 1, width, MPI_INT, &col_block_type);
  MPI_Type_create_resized(col_block_type, 0, sizeof(int), &col_block_type);
  MPI_Type_commit(&col_block_type);

  // Types for sending and receiving row blocks.
  MPI_Datatype row_send_type, row_recv_type;
  MPI_Type_vector(ceil((double)height/2), width, width*2, MPI_INT, &row_send_type);
  MPI_Type_vector(ceil((double)height/2), width, width, MPI_INT, &row_recv_type);
  MPI_Type_commit(&row_send_type);
  MPI_Type_commit(&row_recv_type);

  // Each process' individual block of columns.
  int *block = malloc(width * height * sizeof(*block));
  // Used for odd-even transposition sort.
  // Allocating here allows reusage.
  int* partner_block = malloc(width * ceil((double)height/2) * sizeof(*partner_block));

  // Start timer.
  const double start = MPI_Wtime();

  // Scatter initial matrix column blocks.
  MPI_Scatter(matrix, width, col_matrix_type, block, width, col_block_type, 0, MPI_COMM_WORLD);

  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Sort rows locally.
    sort_row(width, height, block);
    // Sort rows globally.
    parity_sort(width, height, block, partner_block, rank, size, 
                  row_send_type, row_recv_type);
    
    if (step < num_steps-1) {
      // Sort columns locally.
      sort_column(width, height, block);
    }
  }
  
  // Gather blocks in temporary array.
  int* matrix_tmp = NULL;
  if (rank == 0) {
    matrix_tmp = malloc(n*n * sizeof(*matrix_tmp));
  }

  MPI_Gather(block, width*height, MPI_INT, matrix_tmp, width*height, MPI_INT, 0, MPI_COMM_WORLD);

  // Copy blocks to correct ordering.
  if (rank == 0) {
    for (int p = 0; p < size; p++) {
      for (int row = 0; row < height; row++) {
        memcpy(&matrix[(p * width) + (row * n)], &matrix_tmp[(p * width * height) + (row * width)], width * sizeof(int));
      }
    }
  }





  // Stop timer.
  double my_execution_time = MPI_Wtime() - start;
	double maxtime;
	MPI_Reduce(&my_execution_time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Print execution time.
  if (rank == 0) {
    printf("%lf\n", maxtime);
  }

  if (rank == 0) {
    if (check_solution) {
      // Check if correctly sorted.
      if (checker(n, matrix, input_file)) {
        printf("Correct!\n");
      } else {
        printf("Incorrect!\n");
      }
    }
    if (output_type != 0) {
      // Print sorted matrix.
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
  free(partner_block);
  MPI_Type_free(&col_matrix_type);
  MPI_Type_free(&col_block_type);
  MPI_Type_free(&row_send_type);
  MPI_Type_free(&row_recv_type);
  MPI_Finalize();
  return 0;
}