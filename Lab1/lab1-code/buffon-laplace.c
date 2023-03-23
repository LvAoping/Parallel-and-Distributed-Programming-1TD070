#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "buffon-laplace.h"

/**
 * This program estimates the value of PI using a Buffon-Laplace simulation.
 * Table [0,1]x[0,1] divided into [Nx,Ny] boxes
 * N needles of length L are dropped
 * Usage: mpirun -np nproc ./buffon-laplace Nx Ny L N
 * @author Malin Kallen 2019-2020
 */

int main(int argc, char **argv) {
	// Set up
	MPI_Init(&argc, &argv);
	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	// Program arguments
	if (argc != 5) {
		if (rank == 0) {
			printf("\nUsage: mpirun -np nproc ./buffon-laplace Nx Ny L N, where\n");
			printf(" Nx = Number of horizontal boxes\n");
			printf(" Ny = Number of vertical boxes\n");
			printf(" L  = Length of the needles\n");
			printf(" N  = Number of needles\n\n");
		}
        MPI_Finalize();
		return 0;
	}

	const int Nx = atoi(argv[1]);
	const int Ny = atoi(argv[2]);
	const double L = atof(argv[3]);
	const int globalN = atoi(argv[4]);
	const int localN = globalN/num_proc + (rank < globalN%num_proc);
    double A = 1.0/Nx;  // Table [0,1]x[0,1]
    double B = 1.0/Ny;
    
	// Simulate needle throws and count number of crossings
	int localC = run_simulation(A, B, L, localN, Nx, Ny);

	// Approximate pi
	int globalC=0;
	MPI_Reduce(&localC, &globalC, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (0 == rank) {
		double pi_approx = globalN * (2*L*(A+B)-L*L) / (globalC*A*B);
		printf("Needles  : %d\nCrossings: %d (%.1f%%)\nPI: %.16f\n", globalN, globalC, 100.0*globalC/globalN, pi_approx);
	}

	MPI_Finalize();
	return 0;
}

int run_simulation(double A, double B, double L, int N, int Nx, int Ny) {
	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	int *seeds;
	if (0 == rank) {
		seeds = malloc(num_proc * sizeof(int));
		for (int i=0; i<num_proc; i++) {
			seeds[i] = rand();
		}
	}
	int seed;
	MPI_Scatter(seeds, 1, MPI_INTEGER, &seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	if (0 == rank) {
		free(seeds);
	}
	srand(seed);
	int C = 0;
	point_t eye, point;
	for (int i=0; i<N; i++) {
		throw_needle(L, &eye, &point, Nx, Ny);
		if (crossing(eye, point, A, B, Nx, Ny)) {
			C++;
		}
	}
	return C;
}

void throw_needle(double L, point_t *eye, point_t *point, int Nx, int Ny) {
	eye->x = 1.0*rand()/RAND_MAX;
	eye->y = 1.0*rand()/RAND_MAX;
	double alfa = 2 * PI * rand()/RAND_MAX;	// Angle of needle
	point->x = eye->x + L * cos(alfa);
	point->y = eye->y + L * sin(alfa);
}

int crossing(point_t eye, point_t point, double A, double B, int Nx, int Ny) {
 
    // Determine box of eye and point
    int g1x=(int)(eye.x/A);
    int g1y=(int)(eye.y/B);
    int g2x=(int)(point.x/A);
    int g2y=(int)(point.y/B);
    
    // Check if eye and point in different boxes or outside domain
    return g1x!=g2x || g1y!=g2y || point.x<0 || point.y<0;
}

