#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#define ROOT 0

// ---- INPUT / OUTPUT -------------------------------------
void bad_input()
{
  printf("Bad input file.\n");
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(EXIT_FAILURE);
}

int32_t* read_input(int* out_n, char* input_file)
{
  FILE* infile = fopen(input_file, "r");
  if (!infile) {
    bad_input();
  }
  if (fscanf(infile, "%d", out_n) != 1) {
    bad_input();
  }
  int n = *out_n;
  int32_t* M = calloc(n*n, sizeof(int32_t));
  for (int i = 0; i < n*n; ++i) {
    if (fscanf(infile, "%d", &M[i]) != 1) {
      bad_input();
    }
  }
  fclose(infile);
  return M;
}

void print_matrix(int n, int32_t* M)
{
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      printf("%d ", M[row*n + col]);
    }
    printf("\n");
  }
}

void print_matrix_file(int n, int32_t* M, char* output_file)
{
  FILE *outfile = fopen(output_file, "w");
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      fprintf(outfile, "%d ", M[row*n + col]);
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
  int32_t x = *((int32_t*)a);
  int32_t y = *((int32_t*)b);
  return (x > y) - (x < y);
}

int descending(const void* a, const void* b)
{
  int32_t x = *((int32_t*)a);
  int32_t y = *((int32_t*)b);
  return (x < y) - (x > y);
}

bool less(int32_t a, int32_t b)
{
  return a < b;
}

bool greater_equal(int32_t a, int32_t b)
{
  return a >= b;
}

void swap(int32_t *a, int32_t *b)
{
  int32_t tmp = *a;
  *a = *b;
  *b = tmp;
}

// ---- FUNCTIONALITY --------------------------------------
bool check_same_elements(int n, int32_t *M, char *input_file)
{
  int32_t *initial = read_input(&n, input_file);
  int32_t *result  = calloc(n*n, sizeof(*result));
  memcpy(result, M, n*n * sizeof(int32_t));

  qsort(initial, n*n, sizeof(int32_t), ascending);
  qsort(result, n*n, sizeof(int32_t), ascending);

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

bool check_sorted(int n, int32_t *M)
{
  int32_t prev = M[0];
  for (int row = 0; row < n; row++) {
    if (even(row)) { // Even row: go left to right
      for (int col = 0; col < n; col++) {
        if (prev > M[row*n + col]) {
          return false;
        }
        prev = M[row*n + col]; 
      }
    }
    else { // Odd row: go right to left.
      for (int col = n-1; col >= 0; col--) {
        if (prev > M[row*n + col]) {
          return false;
        }
        prev = M[row*n + col]; 
      }
    }
  }
  return true;
}

bool checker(int n, int32_t* M, char *input_file) 
{
  return check_sorted(n, M) && check_same_elements(n, M, input_file);
}

long partition(int32_t *data, int n, int col, long left, long right)
{
  long pivot_index = left + (right - left) / 2;
  const int32_t pivot = data[pivot_index*n+col];
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

void quicksort(int32_t *data, int n, int col, long left, long right)
{
  if (right > left) {
    long pivot_index = partition(data, n, col, left, right);
    quicksort(data, n, col, left, pivot_index - 1);
    quicksort(data, n, col, pivot_index + 1, right);
  }
}

void sort_rows(int w, int h, int32_t *M)
{
  for (int row = 0; row < h; row++) {
    if (even(row)) {
      qsort(&M[row*w], w, sizeof(int32_t), ascending);
    } else {
      qsort(&M[row*w], w, sizeof(int32_t), descending);
    }
  }
}

void sort_columns(int w, int h, int32_t* M)
{
  for (int col = 0; col < w; col++) {
    quicksort(M, w, col, 0, h-1);
  }
}

void merge(const int32_t *v1, int n1, const int32_t *v2, int n2, 
           int32_t *result, bool (*compare)(int32_t, int32_t)) 
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


void exchange_and_merge(int partner, int rank, int w, int h, int32_t *slice, 
     int32_t *partner_slice, MPI_Datatype TYPE_ROW_SEND, MPI_Datatype TYPE_ROW_RECV)
{
  int offset = even(rank) ? w : 0;

  MPI_Sendrecv(
      slice + offset,   // const void *sendbuf
      1,                // int sendcount
      TYPE_ROW_SEND,    // MPI_Datatype sendtype
      partner,          // int dest 
      rank,             // int sendtag
      partner_slice,    // void *recvbuf
      1,                // int recvcount
      TYPE_ROW_RECV,    // MPI_Datatype recvtype
      partner,          // int source
      partner,          // int recvtag
      MPI_COMM_WORLD,   // MPI_Comm comm
      MPI_STATUS_IGNORE // MPI_Status *status
  );

  int partner_row = 0;
  int32_t merged[w*2];
  for (int row = even(rank) ? 0 : 1; row < h; row += 2) {
    size_t s = w * sizeof(int32_t);

    if (even(rank)) {
      merge(&slice[row * w], w, &partner_slice[partner_row * w], w, merged, less);
    } else {
      merge(&slice[row * w], w, &partner_slice[partner_row * w], w, merged, greater_equal);
    }

    if (rank < partner) {
      memcpy(&slice[row * w], merged, s);
      memcpy(&partner_slice[partner_row * w], &merged[w], s);
    } else {
      memcpy(&slice[row * w], &merged[w], s);
      memcpy(&partner_slice[partner_row * w], merged, s);
    }
    partner_row += 1;
  }

  MPI_Sendrecv_replace(
      partner_slice,         // void *buf, 
      w * ceil((double)h/2), // int count, 
      MPI_INT32_T,           // MPI_Datatype datatype,
      partner,               // int dest, 
      rank,                  // int sendtag, 
      partner,               // int source, 
      partner,               // int recvtag,
      MPI_COMM_WORLD,        // MPI_Comm comm, 
      MPI_STATUS_IGNORE      // MPI_Status *status
  );

  partner_row = 0;
  for (int row = even(rank) ? 1 : 0; row < h; row += 2) {
    memcpy(&slice[row * w], &partner_slice[partner_row * w], w * sizeof(int32_t));
    partner_row += 1;
  }
}


void odd_even_sort(int w, int h, int32_t *slice, int32_t *partner_slice, int rank, 
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
      exchange_and_merge(partner, rank, w, h, slice, partner_slice, TYPE_ROW_SEND, TYPE_ROW_RECV);
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

  if (rank == ROOT) {
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
  int32_t* M = NULL; // Input matrix

  // Root reads input matrix.
  if (rank == ROOT) {
    M = read_input(&n, input_file);

    if (n % num_PEs != 0) {
      fprintf(stderr, "Matrix size should be divisible by the number of processes.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(EXIT_FAILURE);
    }
  }

  // Broadcast matrix size to all PEs.
  MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  // Width and height of matrix slices.
  const int w = n / num_PEs;
  const int h = n;

  // Each process' individual slice of columns.
  int32_t *slice = calloc(w * h, sizeof(*slice));
  // Used for odd-even transposition sort.
  // Allocating here allows reusage.
  int32_t* partner_slice = calloc(w * ceil((double)h/2), sizeof(*partner_slice));
  
  // Types for sending and receiving columns.
  MPI_Datatype TYPE_TMP, TYPE_TMP2, TYPE_COL_MATRIX, TYPE_COL_SLICE;
  MPI_Type_vector(h, 1, h, MPI_INT32_T, &TYPE_TMP);
  MPI_Type_create_resized(TYPE_TMP, 0, sizeof(int32_t), &TYPE_COL_MATRIX);
  MPI_Type_commit(&TYPE_COL_MATRIX);
  
  MPI_Type_vector(h, 1, w, MPI_INT32_T, &TYPE_TMP2);
  MPI_Type_create_resized(TYPE_TMP2, 0, sizeof(int32_t), &TYPE_COL_SLICE);
  MPI_Type_commit(&TYPE_COL_SLICE);
  
  // Types for sending and receiving row slices.
  MPI_Datatype TYPE_ROW_SEND, TYPE_ROW_RECV;
  MPI_Type_vector(ceil((double)h/2), w, w*2, MPI_INT32_T, &TYPE_ROW_SEND);
  MPI_Type_vector(ceil((double)h/2), w, w, MPI_INT32_T, &TYPE_ROW_RECV);
  MPI_Type_commit(&TYPE_ROW_SEND);
  MPI_Type_commit(&TYPE_ROW_RECV);



  // Start timer.
  const double start = MPI_Wtime();

  // ---- Shearsort ----------------------------------------
  // Scatter initial matrix column slices.
  MPI_Scatter(
      M,               // const void *sendbuf
      w,               // int sendcount
      TYPE_COL_MATRIX, // MPI_Datatype sendtype
      slice,           // void *recvbuf
      w,               // int recvcount
      TYPE_COL_SLICE,  // MPI_Datatype recvtype
      ROOT,            // int root
      MPI_COMM_WORLD   // MPI_Comm comm
  );

  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Sort rows locally.
    sort_rows(w, h, slice);
    // Sort rows globally.
    odd_even_sort(w, h, slice, partner_slice, rank, num_PEs, 
                  TYPE_ROW_SEND, TYPE_ROW_RECV);
    
    if (step < num_steps-1) {
      // Sort columns locally.
      sort_columns(w, h, slice);
    }
  }

  // Gather sorted slices.
  // NOTE: Not working on UPPMAX. 
  // Use simpler gather and copy elements afterwards as below.
  
  //MPI_Gather(
  //    slice,           // const void *sendbuf,
  //    w,               // int sendcount,
  //    TYPE_COL_SLICE,  // MPI_Datatype sendtype,
  //    M,               // void *recvbuf,
  //    w,               // int recvcount,
  //    TYPE_COL_MATRIX, // MPI_Datatype recvtype,
  //    ROOT,            // int root,
  //    MPI_COMM_WORLD   // MPI_Comm comm
  //);

  
  // Gather slices in temporary array.
  int32_t* tmp = NULL;
  if (rank == ROOT) {
    tmp = calloc(n*n, sizeof(int32_t));
  }

  MPI_Gather(
      slice,         // const void *sendbuf,
      w*h,           // int sendcount,
      MPI_INT32_T,   // MPI_Datatype sendtype,
      tmp,           // void *recvbuf,
      w*h,           // int recvcount,
      MPI_INT32_T,   // MPI_Datatype recvtype,
      ROOT,          // int root,
      MPI_COMM_WORLD // MPI_Comm comm
  );

  // Copy slices to correct ordering.
  if (rank == ROOT) {
    for (int p = 0; p < num_PEs; p++) {
      for (int row = 0; row < h; row++) {
        memcpy(&M[(p * w) + (row * n)], &tmp[(p * w * h) + (row * w)], w * sizeof(int32_t));
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
  if (rank == ROOT) {
    printf("%lf\n", slowest);
  }

  if (rank == ROOT) {
    if (!suppress_output) {
      // Output sorted matrix.
      if (output_file) {
        print_matrix_file(n, M, output_file);
      } else {
        print_matrix(n, M);
      }
    }
    if (check_solution) {
      // Check if correctly sorted.
      if (checker(n, M, input_file)) {
        printf("Correct!\n");
      } else {
        printf("Incorrect!\n");
      }
    }
  }

  // Clean up.
  free(M);
  free(tmp);
  free(slice);
  free(partner_slice);
  MPI_Type_free(&TYPE_COL_MATRIX);
  MPI_Type_free(&TYPE_COL_SLICE);
  MPI_Type_free(&TYPE_ROW_SEND);
  MPI_Type_free(&TYPE_ROW_RECV);
  MPI_Finalize();
  return 0;
}