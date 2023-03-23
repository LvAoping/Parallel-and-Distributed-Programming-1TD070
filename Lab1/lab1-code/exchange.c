/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a, b;
  MPI_Status status;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  a = 100.0 + (double) rank;  /* Different a on different processors */

  /* Exchange variable a, notice the send-recv order */
  // if (rank == 0) {
  //   MPI_Send(&a, 1, MPI_DOUBLE, 1, 111, MPI_COMM_WORLD);
  //   MPI_Recv(&b, 1, MPI_DOUBLE, 2, 333, MPI_COMM_WORLD, &status);
  //   printf("Processor 0 got %f from processor 2\n", b);
  // } else if (rank==1) {
  //   MPI_Recv(&b, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &status);
  //   MPI_Send(&a, 1, MPI_DOUBLE, 2, 222, MPI_COMM_WORLD);
  //   printf("Processor 1 got %f from processor 0\n", b);
  // } else if (rank==2) {
  //   MPI_Send(&a, 1, MPI_DOUBLE, 0, 333, MPI_COMM_WORLD);
  //   MPI_Recv(&b, 1, MPI_DOUBLE, 1, 222, MPI_COMM_WORLD, &status);
  //   printf("Processor 2 got %f from processor 1\n", b);
  // }

if (rank == 0) {
    MPI_Send(&a, 1, MPI_DOUBLE, 1, rank+1, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, size-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    printf("Processor 0 got %f from processor %i\n", b, size-1);
  } else if (rank== size-1) {
    MPI_Send(&a, 1, MPI_DOUBLE, 0, size, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, size-2, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    printf("Processor %i got %f from processor %i\n",size-1, b, size-2);
  } else {
    MPI_Recv(&b, 1, MPI_DOUBLE, rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Send(&a, 1, MPI_DOUBLE, rank+1, rank+1, MPI_COMM_WORLD);  
    printf("Processor %i got %f from processor %i\n",rank, b, rank-1);
  }

  MPI_Finalize(); 

  return 0;
}
