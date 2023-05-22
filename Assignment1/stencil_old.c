#include "stencil.h"

int main(int argc, char **argv) {
	if (4 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}

	char *input_name = argv[1];
	char *output_name = argv[2];
	int num_steps = atoi(argv[3]);
	int num_values, chunk;

	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv); /* Initialize MPI */
  	MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

	double *input = 0;
	double *output = 0;

	if (rank==0){
		// Read input file
		if (0 > (num_values = read_input(input_name, &input))) {
			return 2;
		}
		// Allocate data for result
		if (NULL == (output = malloc(num_values * sizeof(double)))) {
			perror("Couldn't allocate memory for output");
			return 2;
		}
	}

	// Stencil values
	MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

	chunk = num_values / size;
	double *local_input = (double *)malloc((chunk + 2*EXTENT) * sizeof(double));
	double *local_output = (double *)malloc((chunk + 2*EXTENT) * sizeof(double));
	MPI_Scatter(input, chunk, MPI_DOUBLE, local_input+EXTENT, chunk, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

	// Start timer
	double start = MPI_Wtime();
	int left = rank - 1;
	if (left < 0){left = size-1;}
	int right = rank + 1;
	if (right > size -1){right=0;}

	// Repeatedly apply stencil
	for (int s=0; s<num_steps; s++) {
		// first pass data to the left then to the right
		/*
		if (size != 1){
			if (rank==0){
				MPI_Ssend(local_input+EXTENT, EXTENT, MPI_DOUBLE, size-1, 000, MPI_COMM_WORLD);
				MPI_Recv(local_input+chunk+EXTENT, EXTENT, MPI_DOUBLE, rank+1, 000, MPI_COMM_WORLD, &status);
				MPI_Ssend(local_input+chunk, EXTENT, MPI_DOUBLE, rank+1, 111, MPI_COMM_WORLD);
				MPI_Recv(local_input, EXTENT, MPI_DOUBLE, size-1, 111, MPI_COMM_WORLD, &status);			
			} else if (rank==size-1){
				MPI_Recv(local_input+chunk+EXTENT, EXTENT, MPI_DOUBLE, 0, 000, MPI_COMM_WORLD, &status);
				MPI_Ssend(local_input+EXTENT, EXTENT, MPI_DOUBLE, rank-1, 000, MPI_COMM_WORLD);
				MPI_Recv(local_input, EXTENT, MPI_DOUBLE, rank-1, 111, MPI_COMM_WORLD, &status);	
				MPI_Ssend(local_input+chunk, EXTENT, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
			} else{
				MPI_Recv(local_input+chunk+EXTENT, EXTENT, MPI_DOUBLE, rank+1, 000, MPI_COMM_WORLD, &status);
				MPI_Ssend(local_input+EXTENT, EXTENT, MPI_DOUBLE, rank-1, 000, MPI_COMM_WORLD);
				MPI_Recv(local_input, EXTENT, MPI_DOUBLE, rank-1, 111, MPI_COMM_WORLD, &status);
				MPI_Ssend(local_input+chunk, EXTENT, MPI_DOUBLE, rank+1, 111, MPI_COMM_WORLD);
			}
		} else {
			MPI_Send(local_input+EXTENT, EXTENT, MPI_DOUBLE, 0, 000, MPI_COMM_WORLD);
			MPI_Recv(local_input+chunk+EXTENT, EXTENT, MPI_DOUBLE, 0, 000, MPI_COMM_WORLD, &status);
			MPI_Send(local_input+chunk, EXTENT, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
			MPI_Recv(local_input, EXTENT, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &status);	
		}
		*/
		MPI_Sendrecv(local_input+EXTENT, EXTENT, MPI_DOUBLE, left, 000, local_input+chunk+EXTENT, EXTENT, MPI_DOUBLE, right, 000, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(local_input+chunk, EXTENT, MPI_DOUBLE, right, 111, local_input, EXTENT, MPI_DOUBLE, left, 111, MPI_COMM_WORLD, &status);

		// Apply stencil
		for (int i=EXTENT; i<chunk+EXTENT; i++) {
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * local_input[index];
			}
			local_output[i] = result;
		}
		
		// Swap local_input and local_output
		if (s < num_steps-1) {
			double *tmp = local_input;
			local_input = local_output;
			local_output = tmp;
		}
	}
	
	// Stop timer
	double my_execution_time = MPI_Wtime() - start;
	double maxtime;
	MPI_Reduce( &my_execution_time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	MPI_Gather(local_output+EXTENT, chunk, MPI_DOUBLE, output, chunk, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// MPI_Gather(&my_execution_time, 1, MPI_DOUBLE, timelst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// double maxtime=timelst[0];

	if (rank==0){
		// Write result
		// for (int i=0; i<size; i++) {
		// 	if (timelst[i]>maxtime) {
		// 		maxtime = timelst[i];
		// 	}
		// }
		printf("%f\n", maxtime);
		#ifdef PRODUCE_OUTPUT_FILE
		if (0 != write_output(output_name, output, num_values)) {
			return 2;
		}
		#endif
		
		// Clean up
		free(input);
		free(output);
	}
	
	free(local_input);
	free(local_output);
	MPI_Finalize(); /* Shut down and clean up MPI */
	return 0;
}


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


int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}