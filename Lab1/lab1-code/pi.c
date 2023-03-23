/**********************************************************************
 * This program calculates pi using MPI/C
 *
 * Use midpoint rule for the numerical integration
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  int rank, size;
  const long int intervals = 100000000L ; /* The sum is [globally]
					     divided into this many
					     intervals     */
  int chunk;             /* This many iterations will I do */
  int i, istart, istop;  /* Variables for the local loop   */
  double sum, dx;

  double t_begin, t_end;  /* Variables for timing */

  MPI_Init(&argc, &argv); /* Initialize MPI */
  t_begin = MPI_Wtime();  /* Start timing */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

  chunk  = intervals/size;       /* (We assume this is an integer)   */
  istart = rank*chunk;           /* Calculate start and stop indices */
  istop  = (rank + 1)*chunk - 1; /* for the local loop               */

  dx  = 1.0/(size*chunk);
  sum = 0.0;
  for (i = istart; i <= istop; i++) { /* The local loop */
    double tmp = dx*(0.5 + i);
    sum += dx*4.0/(1.0 + tmp*tmp);
  }

  if (rank == 0) {
    double globsum = sum;

    for (i = 1; i < size; i++) { /* Collect the partial sums */
      MPI_Status status;
      MPI_Recv(&sum, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
      globsum += sum;
    }
    
    printf("PI is approx. %.16f\n",  globsum);

  } else { /* Send my partial sum to the processor with its
	      rank equal to zero */
    MPI_Send(&sum, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
  }

  
  
  MPI_Finalize(); /* Shut down and clean up MPI */
  


  return 0;
}
