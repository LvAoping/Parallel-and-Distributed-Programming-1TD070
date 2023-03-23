/**********************************************************************
 * This program deadlocks beautifully using MPI/C
 *
 * Andreas K�h�ri
 * 1999-12-29
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
  
  a = (double) rank;

  if (size != 2) { /* This if-block makes sure only two processors
                    * takes part in the execution of the code, pay no
                    * attention to it */
    if (rank == 0)
      fprintf(stderr, "\aRun on two processors only!\n");

    MPI_Finalize();
    exit(0);
  }
  
  /* Here comes the important bit */

  if (rank == 0) {
    MPI_Send(&a, 1, MPI_DOUBLE, 1, 111, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, 1, 222, MPI_COMM_WORLD, &status);
    printf("Processor 0 got %f from processor 1\n", b);
  } else {
    MPI_Send(&a, 1, MPI_DOUBLE, 0, 222, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &status);
    printf("Processor 1 got %f from processor 0\n", b);
  }

  MPI_Finalize(); /* Shut down and clean up MPI */

  return 0;
}
