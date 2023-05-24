#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>



// ---- INPUT / OUTPUT -------------------------------------
int* read_input(int* out_n, char* input_file)
{
  FILE* file = fopen(input_file, "r");
  if (!file) {
    perror("Couldn't open input file");
    exit(-1);
  }
  if (fscanf(file, "%d", out_n) != 1) {
    perror("Couldn't read element count from input file");
    exit(-1);
  }
  int n = *out_n;
  int* matrix = calloc(n*n, sizeof(int));
  for (int i = 0; i < n*n; ++i) {
    if (fscanf(file, "%d", &matrix[i]) != 1) {
      perror("Couldn't read elements from input file");
			exit(-1);
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

void print_matrix_file(int n, int* matrix, char* output_file)
{
  FILE *outfile = fopen(output_file, "w");
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      fprintf(outfile, "%d ", matrix[row*n + col]);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
}

// ---- HELPERS --------------------------------------------
bool even(int n)
{
  return n % 2 == 0;
}

int ascending(const void* a, const void* b)
{
  int x = *((int*)a);
  int y = *((int*)b);
  return (x > y) - (x < y);
}

int descending(const void* a, const void* b)
{
  int x = *((int*)a);
  int y = *((int*)b);
  return (x < y) - (x > y);
}

bool less(int a, int b)
{
  return a < b;
}

bool greater_equal(int a, int b)
{
  return a >= b;
}

void swap(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

// ---- FUNCTIONALITY --------------------------------------
bool check_same_elements(int n, int *matrix, char *input_file)
{
  int *initial = read_input(&n, input_file);
  int *result  = calloc(n*n, sizeof(*result));
  memcpy(result, matrix, n*n * sizeof(int));

  qsort(initial, n*n, sizeof(int), ascending);
  qsort(result, n*n, sizeof(int), ascending);

  bool ret = true;
  for (int i = 0; i < n*n; i++) {
    if (result[i] != initial[i]) {
      ret = false;
      break;
    }
  }
  free(initial);
  free(result);
  return ret;
}

bool check_sorted(int n, int *matrix)
{
  int prev = matrix[0];
  for (int row = 0; row < n; row++) {
    if (even(row)) { // Even row: go left to right
      for (int col = 0; col < n; col++) {
        if (prev > matrix[row*n + col]) {
          return false;
        }
        prev = matrix[row*n + col]; 
      }
    }
    else { // Odd row: go right to left.
      for (int col = n-1; col >= 0; col--) {
        if (prev > matrix[row*n + col]) {
          return false;
        }
        prev = matrix[row*n + col]; 
      }
    }
  }
  return true;
}

bool checker(int n, int* matrix, char *input_file) 
{
  return check_sorted(n, matrix) && check_same_elements(n, matrix, input_file);
}

long partition(int *data, int n, int col, long left, long right)
{
  long pivot_index = left + (right - left) / 2;
  const int pivot = data[pivot_index*n+col];
  swap(&data[pivot_index*n+col], &data[right*n+col]);

  for (long i = left; i < right; i++) {
    if (data[i*n+col] <= pivot) {
      swap(&data[i*n+col], &data[left*n+col]);
      left++;
    }
  }
  swap(&data[left*n+col], &data[right*n+col]);
  return left;
}

void quicksort(int *data, int n, int col, long left, long right)
{
  if (right > left) {
    long pivot_index = partition(data, n, col, left, right);
    quicksort(data, n, col, left, pivot_index - 1);
    quicksort(data, n, col, pivot_index + 1, right);
  }
}

void sort_rows(int w, int h, int *M)
{
  for (int row = 0; row < h; row++) {
    if (even(row)) {
      qsort(&M[row*w], w, sizeof(int), ascending);
    } else {
      qsort(&M[row*w], w, sizeof(int), descending);
    }
  }
}

void sort_columns(int w, int h, int* M)
{
  for (int col = 0; col < w; col++) {
    quicksort(M, w, col, 0, h-1);
  }
}

void merge(const int *v1, int n1, const int *v2, int n2, 
           int *result, bool (*compare)(int, int)) 
{
  int i = 0;
  int j = 0;
  int k = 0;

  while (i < n1 && j < n2) {
    if (compare(v1[i], v2[j])) {
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
  int offset = even(rank) ? w : 0;

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
  for (int row = even(rank) ? 0 : 1; row < h; row += 2) {
    size_t s = w * sizeof(int);

    if (even(rank)) {
      merge(&block[row * w], w, &partner_block[partner_row * w], w, merged, less);
    } else {
      merge(&block[row * w], w, &partner_block[partner_row * w], w, merged, greater_equal);
    }

    if (rank < partner) {
      memcpy(&block[row * w], merged, s);
      memcpy(&partner_block[partner_row * w], &merged[w], s);
    } else {
      memcpy(&block[row * w], &merged[w], s);
      memcpy(&partner_block[partner_row * w], merged, s);
    }
    partner_row += 1;
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
  for (int row = even(rank) ? 1 : 0; row < h; row += 2) {
    memcpy(&block[row * w], &partner_block[partner_row * w], w * sizeof(int));
    partner_row += 1;
  }
}


void odd_even_sort(int w, int h, int *block, int *partner_block, int rank, 
     int num_PEs, MPI_Datatype TYPE_ROW_SEND, MPI_Datatype TYPE_ROW_RECV)
{
  for (int i = 0; i < num_PEs; i++) {
    int partner = -1;
    if (i % 2 == 0) {
      partner = even(rank) ? rank + 1 : rank - 1;
    } else {
      partner = even(rank) ? rank - 1 : rank + 1;
    }
    if (partner >= 0 && partner < num_PEs) {
      exchange_and_merge(partner, rank, w, h, block, partner_block, TYPE_ROW_SEND, TYPE_ROW_RECV);
    }
  }
}

// ---- MAIN -----------------------------------------------
int main(int argc, char **argv) 
{
  MPI_Init(&argc, &argv);
  
  // Find own rank and number of PEs.
  int num_PEs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_PEs);
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

    if (n % num_PEs != 0) {
      fprintf(stderr, "Matrix size should be divisible by the number of processes.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(EXIT_FAILURE);
    }
  }

  // Broadcast matrix size to all PEs.
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Width and height of matrix blocks.
  const int w = n / num_PEs;
  const int h = n;

  // Each process' individual block of columns.
  int *block = calloc(w * h, sizeof(*block));
  // Used for odd-even transposition sort.
  // Allocating here allows reusage.
  int* partner_block = calloc(w * ceil((double)h/2), sizeof(*partner_block));
  
  // // Types for sending and receiving columns.
  MPI_Datatype TYPE_TMP, TYPE_TMP2, TYPE_COL_MATRIX, TYPE_COL_SLICE;
  MPI_Type_vector(h, 1, h, MPI_INT, &TYPE_TMP);
  MPI_Type_create_resized(TYPE_TMP, 0, sizeof(int), &TYPE_COL_MATRIX);
  MPI_Type_commit(&TYPE_COL_MATRIX);
  
  MPI_Type_vector(h, 1, w, MPI_INT, &TYPE_TMP2);
  MPI_Type_create_resized(TYPE_TMP2, 0, sizeof(int), &TYPE_COL_SLICE);
  MPI_Type_commit(&TYPE_COL_SLICE);
  
  // Types for sending and receiving row blocks.
  MPI_Datatype TYPE_ROW_SEND, TYPE_ROW_RECV;
  MPI_Type_vector(ceil((double)h/2), w, w*2, MPI_INT, &TYPE_ROW_SEND);
  MPI_Type_vector(ceil((double)h/2), w, w, MPI_INT, &TYPE_ROW_RECV);
  MPI_Type_commit(&TYPE_ROW_SEND);
  MPI_Type_commit(&TYPE_ROW_RECV);



  // Start timer.
  const double start = MPI_Wtime();

  // ---- Shearsort ----------------------------------------
  // Scatter initial matrix column blocks.
  MPI_Scatter(
      matrix,          // const void *sendbuf
      w,               // int sendcount
      TYPE_COL_MATRIX, // MPI_Datatype sendtype
      block,           // void *recvbuf
      w,               // int recvcount
      TYPE_COL_SLICE,  // MPI_Datatype recvtype
      0,            // int root
      MPI_COMM_WORLD   // MPI_Comm comm
  );

  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Sort rows locally.
    sort_rows(w, h, block);
    // Sort rows globally.
    odd_even_sort(w, h, block, partner_block, rank, num_PEs, 
                  TYPE_ROW_SEND, TYPE_ROW_RECV);
    
    if (step < num_steps-1) {
      // Sort columns locally.
      sort_columns(w, h, block);
    }
  }

  // Gather sorted blocks.
  // NOTE: Not working on UPPMAX. 
  // Use simpler gather and copy elements afterwards as below.
  
  //MPI_Gather(
  //    block,           // const void *sendbuf,
  //    w,               // int sendcount,
  //    TYPE_COL_SLICE,  // MPI_Datatype sendtype,
  //    M,               // void *recvbuf,
  //    w,               // int recvcount,
  //    TYPE_COL_MATRIX, // MPI_Datatype recvtype,
  //    0,            // int root,
  //    MPI_COMM_WORLD   // MPI_Comm comm
  //);

  
  // Gather blocks in temporary array.
  int* tmp = NULL;
  if (rank == 0) {
    tmp = malloc(n * n * sizeof(int));
  }

  MPI_Gather(
      block,         // const void *sendbuf,
      w*h,           // int sendcount,
      MPI_INT,   // MPI_Datatype sendtype,
      tmp,           // void *recvbuf,
      w*h,           // int recvcount,
      MPI_INT,   // MPI_Datatype recvtype,
      0,          // int root,
      MPI_COMM_WORLD // MPI_Comm comm
  );

  // Copy blocks to correct ordering.
  if (rank == 0) {
    for (int p = 0; p < num_PEs; p++) {
      for (int row = 0; row < h; row++) {
        memcpy(&matrix[(p * w) + (row * n)], &tmp[(p * w * h) + (row * w)], w * sizeof(int));
      }
    }
  }
  // ---- End shearsort ------------------------------------
  
  // Stop timer.
  const double my_execution_time = MPI_Wtime() - start;

  // Find largest execution time.
  double slowest = 0;
  MPI_Reduce(
      &my_execution_time, 
      &slowest, 
      1,
      MPI_DOUBLE, 
      MPI_MAX, 
      0, 
      MPI_COMM_WORLD
  );

  // Output execution time.
  if (rank == 0) {
    printf("%lf\n", slowest);
  }

  if (rank == 0) {
    if (!suppress_output) {
      // Output sorted matrix.
      if (output_file) {
        print_matrix_file(n, matrix, output_file);
      } else {
        // print_matrix(n, M);
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
  MPI_Type_free(&TYPE_COL_SLICE);
  MPI_Type_free(&TYPE_ROW_SEND);
  MPI_Type_free(&TYPE_ROW_RECV);
  MPI_Finalize();
  return 0;
}