#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>



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





int ascending(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

int descending(const void* a, const void* b) {
    return (*(int*)b - *(int*)a);
}



// bool check_same_elements(int n, int *matrix, char *input_file)
// {
//   int *initial = read_input(&n, input_file);
//   int *result  = malloc(n * n * sizeof(*result));
//   memcpy(result, matrix, n*n * sizeof(int));

//   qsort(initial, n*n, sizeof(int), ascending);
//   qsort(result, n*n, sizeof(int), ascending);

//   bool ret = true;
//   for (int i = 0; i < n*n; i++) {
//     if (result[i] != initial[i]) {
//       ret = false;
//       break;
//     }
//   }
//   free(initial);
//   free(result);
//   return ret;
// }

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

// bool check_sorted(int n, int *matrix)
// {
//   int prev = matrix[0];
//   for (int row = 0; row < n; row++) {
//     if ((row & 1) == 0) { // Even row: go left to right
//       for (int col = 0; col < n; col++) {
//         if (prev > matrix[row*n + col]) {
//           return false;
//         }
//         prev = matrix[row*n + col]; 
//       }
//     }
//     else { // Odd row: go right to left.
//       for (int col = n-1; col >= 0; col--) {
//         if (prev > matrix[row*n + col]) {
//           return false;
//         }
//         prev = matrix[row*n + col]; 
//       }
//     }
//   }
//   return true;
// }

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

void sort_row(int w, int h, int *matrix)
{
  for (int row = 0; row < h; row++) {
    if ((row & 1) == 0) {
      qsort(&matrix[row*w], w, sizeof(int), ascending);
    } else {
      qsort(&matrix[row*w], w, sizeof(int), descending);
    }
  }
}

void sort_column(int w, int h, int* matrix) {
  for (int col = 0; col < w; col++) {
    int* column_data = (int*) malloc(h * sizeof(int));
    if (!column_data) {
        perror("Couldn't allocate memory for the column data");
    }
    // Extract the column data
    for (int i = 0; i < h; i++)
        column_data[i] = matrix[i * w + col];
    // Sort the column data
    qsort(column_data, h, sizeof(int), ascending);
    // Store back the sorted column data
    for (int i = 0; i < h; i++)
        matrix[i * w + col] = column_data[i];
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



void exchange_and_merge(int partner, int rank, int w, int h, int *block, 
     int *partner_block, MPI_Datatype TYPE_ROW_SEND, MPI_Datatype TYPE_ROW_RECV)
{
  int offset = ((rank & 1) == 0) ? w : 0;

  MPI_Sendrecv(
      block + offset,   // const void *sendbuf
      1,                // int sendcount
      TYPE_ROW_SEND,    // MPI_Datatype sendtype
      partner,          // int dest 
      rank,             // int sendtag
      partner_block,    // void *recvbuf
      1,                // int recvcount
      TYPE_ROW_RECV,    // MPI_Datatype recvtype
      partner,          // int source
      partner,          // int recvtag
      MPI_COMM_WORLD,   // MPI_Comm comm
      MPI_STATUS_IGNORE // MPI_Status *status
  );

  int partner_row = 0;
  int merged[w*2];
  // for (int row = ((rank & 1) == 0) ? 0 : 1; row < h; row += 2) {
  //   size_t s = w * sizeof(int);

  //   if ((rank & 1) == 0) {
  //     merge(&block[row * w], w, &partner_block[partner_row * w], w, merged, true);
  //   } else {
  //     merge(&block[row * w], w, &partner_block[partner_row * w], w, merged, false);
  //   }

  //   if (rank < partner) {
  //     memcpy(&block[row * w], merged, s);
  //     memcpy(&partner_block[partner_row * w], &merged[w], s);
  //   } else {
  //     memcpy(&block[row * w], &merged[w], s);
  //     memcpy(&partner_block[partner_row * w], merged, s);
  //   }
  //   partner_row += 1;
  // }

  for (int row = ((rank & 1) == 0) ? 0 : 1, partner_row = 0; row < h; row += 2, partner_row++) {
      size_t copysize = w * sizeof(int);

      // Merge two sorted sequences
      bool ascending = (rank & 1) == 0;
      merge(&block[row * w], w, &partner_block[partner_row * w], w, merged, ascending);

      // Copy the merged sequences back into the correct position
      if (rank < partner) {
          memcpy(&block[row * w], merged, copysize);
          memcpy(&partner_block[partner_row * w], &merged[w], copysize);
      } else {
          memcpy(&block[row * w], &merged[w], copysize);
          memcpy(&partner_block[partner_row * w], merged, copysize);
      }
  }

  MPI_Sendrecv_replace(
      partner_block,         // void *buf, 
      w * ceil((double)h/2), // int count, 
      MPI_INT,           // MPI_Datatype datatype,
      partner,               // int dest, 
      rank,                  // int sendtag, 
      partner,               // int source, 
      partner,               // int recvtag,
      MPI_COMM_WORLD,        // MPI_Comm comm, 
      MPI_STATUS_IGNORE      // MPI_Status *status
  );

  partner_row = 0;
  for (int row = ((rank & 1) == 0) ? 1 : 0; row < h; row += 2) {
    memcpy(&block[row * w], &partner_block[partner_row * w], w * sizeof(int));
    partner_row += 1;
  }
}



void odd_even_sort(int w, int h, int *block, int *partner_block, int rank, 
     int size, MPI_Datatype TYPE_ROW_SEND, MPI_Datatype TYPE_ROW_RECV)
{
  for (int i = 0; i < size; i++) {
    int partner = -1;
    if ((i & 1) == 0) {
      partner = ((rank & 1) == 0) ? rank + 1 : rank - 1;
    } else {
      partner = ((rank & 1) == 0) ? rank - 1 : rank + 1;
    }
    if (partner >= 0 && partner < size) {
      exchange_and_merge(partner, rank, w, h, block, partner_block, TYPE_ROW_SEND, TYPE_ROW_RECV);
    }
  }
}

/* Main Function */
int main(int argc, char **argv) 
{
  MPI_Init(&argc, &argv);
  
  // Find own rank and number of processes
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Parse arguments.
  char* input_file = NULL;
  char* output_file = NULL;
  bool check_solution = false;
  bool suppress_output = false;

  if (rank == 0) {
    int opt;
    while ((opt = getopt(argc, argv, "cso:")) != -1) {
      switch (opt) {
      case 'c': 
        check_solution = true; 
        break;
      case 's': 
        suppress_output = true; 
        break;
      case 'o': 
        output_file = optarg;
        break;
      default:
        fprintf(stderr, "Usage: %s [-cs] [-o <output-file>] <input-file> \n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(EXIT_FAILURE);
      }
    }

    if (optind >= argc) {
      fprintf(stderr, "Expected input file.\n");
      fprintf(stderr, "Usage: %s [-cs] [-o <output-file>] <input-file> \n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(EXIT_FAILURE);
    }

    input_file = argv[optind];
  }

  
  int n; // Size of matrix
  int* matrix = NULL; // Input matrix

  // Root reads input matrix.
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
  const int w = n / size;
  const int h = n;

  // Each process' individual block of columns.
  int *block = malloc(w * h * sizeof(*block));
  // Used for odd-even transposition sort.
  // Allocating here allows reusage.
  int* partner_block = malloc(w * ceil((double)h/2) * sizeof(*partner_block));
  
  // Types for sending and receiving columns.
  MPI_Datatype TYPE_TMP, TYPE_TMP2, TYPE_COL_MATRIX, TYPE_COL_BLOCK;
  MPI_Type_vector(h, 1, h, MPI_INT, &TYPE_TMP);
  MPI_Type_create_resized(TYPE_TMP, 0, sizeof(int), &TYPE_COL_MATRIX);
  MPI_Type_commit(&TYPE_COL_MATRIX);
  
  MPI_Type_vector(h, 1, w, MPI_INT, &TYPE_TMP2);
  MPI_Type_create_resized(TYPE_TMP2, 0, sizeof(int), &TYPE_COL_BLOCK);
  MPI_Type_commit(&TYPE_COL_BLOCK);
  
  // Types for sending and receiving row blocks.
  MPI_Datatype TYPE_ROW_SEND, TYPE_ROW_RECV;
  MPI_Type_vector(ceil((double)h/2), w, w*2, MPI_INT, &TYPE_ROW_SEND);
  MPI_Type_vector(ceil((double)h/2), w, w, MPI_INT, &TYPE_ROW_RECV);
  MPI_Type_commit(&TYPE_ROW_SEND);
  MPI_Type_commit(&TYPE_ROW_RECV);



  // Start timer.
  const double start = MPI_Wtime();


  // Scatter initial matrix column blocks.
  MPI_Scatter(matrix, w, TYPE_COL_MATRIX, block, w, TYPE_COL_BLOCK, 0, MPI_COMM_WORLD);

  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Sort rows locally.
    sort_row(w, h, block);
    // Sort rows globally.
    odd_even_sort(w, h, block, partner_block, rank, size, 
                  TYPE_ROW_SEND, TYPE_ROW_RECV);
    
    if (step < num_steps-1) {
      // Sort columns locally.
      sort_column(w, h, block);
    }
  }

  // Gather sorted blocks.
  // NOTE: Not working on UPPMAX. 
  // Use simpler gather and copy elements afterwards as below.
  
  //MPI_Gather(
  //    block,           // const void *sendbuf,
  //    w,               // int sendcount,
  //    TYPE_COL_BLOCK,  // MPI_Datatype sendtype,
  //    M,               // void *recvbuf,
  //    w,               // int recvcount,
  //    TYPE_COL_MATRIX, // MPI_Datatype recvtype,
  //    0,            // int root,
  //    MPI_COMM_WORLD   // MPI_Comm comm
  //);

  
  // Gather blocks in temporary array.
  int* tmp = NULL;
  if (rank == 0) {
    tmp = malloc(n*n * sizeof(*tmp));
  }

  MPI_Gather(block, w*h, MPI_INT, tmp, w*h, MPI_INT, 0, MPI_COMM_WORLD);

  // Copy blocks to correct ordering.
  if (rank == 0) {
    for (int p = 0; p < size; p++) {
      for (int row = 0; row < h; row++) {
        memcpy(&matrix[(p * w) + (row * n)], &tmp[(p * w * h) + (row * w)], w * sizeof(int));
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
    if (!suppress_output) {
      // Print sorted matrix.
      if (output_file) {
        write_output(n, matrix, output_file);
      } else {
        print_matrix(n, matrix);
      }
    }
    if (check_solution) {
      // Check if correctly sorted.
      if (checker(n, matrix, input_file)) {
        printf("Correct!\n");
      } else {
        printf("Incorrect!\n");
      }
    }
  }

  // Clean up.
  free(matrix);
  free(tmp);
  free(block);
  free(partner_block);
  MPI_Type_free(&TYPE_COL_MATRIX);
  MPI_Type_free(&TYPE_COL_BLOCK);
  MPI_Type_free(&TYPE_ROW_SEND);
  MPI_Type_free(&TYPE_ROW_RECV);
  MPI_Finalize();
  return 0;
}