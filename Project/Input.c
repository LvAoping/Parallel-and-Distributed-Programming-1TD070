#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int read_input(const char *file_name, double **values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}


int* read_input(int* out_n, char* input_file)
{
  FILE* file = fopen(input_file, "r");
  if (!file) {
    perror("Couldn't open input file");
    return -1;
  }
  if (fscanf(file, "%d", out_n) != 1) {
    perror("Couldn't read element count from input file");
	return -1;
  }
  int n = *out_n;
  int* M = calloc(n*n, sizeof(int));
  for (int i = 0; i < n*n; ++i) {
    if (fscanf(file, "%d", &M[i]) != 1) {
      perror("Couldn't read elements from input file");
			return -1;
    }
  }
  fclose(file);
  return M;
}