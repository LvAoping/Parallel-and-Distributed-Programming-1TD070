/* PINGPONG */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

double timer();

void processor_A();
void processor_B();

int main(int argc, char *argv[]) {
  int rank, size;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

  if (size != 2) { /* This if-block makes sure only two processors
                    * takes part in the execution of the code, pay no
                    * attention to it */
    if (rank == 0)
      printf("\aRun on two processors only!\n");
  }
  else {
    if ( rank == 0 )
      processor_A();
    else
      processor_B();
  }
  MPI_Finalize();
  return 0;
}

void processor_A(){

  double *message;
  double timer1;
  int maxlen = 2000000;
  const int ping=101, pong=102;
  int len,i;
  int nreps;
  MPI_Status status;

  message = malloc(maxlen*sizeof(double));

  printf("%10s %19s %18s\n","size (B)", "time (µs)", "rate (MB/s)");

  /* warm up */
  len=100;
  for (i=1; i<=10; i++) {
    MPI_Ssend(message, len, MPI_DOUBLE, 1, ping, MPI_COMM_WORLD);
    MPI_Recv(message, len, MPI_DOUBLE, 1, pong, MPI_COMM_WORLD, &status);
  }

  for (len=1; len<=maxlen; len=(int)(len*1.25)+1) {
    timer1=timer();

    /* run more transfers for small messages to get more accurate time measurement */
    if(len < 1000)
      nreps = 10000;
    else
      nreps = 100;

    for (i=1; i<=nreps; i++) {
      MPI_Ssend(message, len, MPI_DOUBLE, 1, ping, MPI_COMM_WORLD);
      MPI_Recv(message, len, MPI_DOUBLE, 1, pong, MPI_COMM_WORLD, &status);
    }

    timer1=timer()-timer1;
    printf("%10d %18.2f %18.2f\n",8*len,1000000.0*timer1/(2.0*nreps),2.0*8.0*nreps*len/timer1/1000000.0);
  }

  free(message);
}

void processor_B(){

  double *message;
  int maxlen = 2000000;
  const int ping=101, pong=102;
  int len,i;
  int nreps;
  MPI_Status status;

  message = malloc(maxlen*sizeof(double));

  /* warm up */
  len=100;
  for (i=1; i<=10; i++) {
    MPI_Recv(message, len, MPI_DOUBLE, 0, ping, MPI_COMM_WORLD, &status);
    MPI_Ssend(message, len, MPI_DOUBLE, 0, pong, MPI_COMM_WORLD);
  }

  for (len=1; len<=maxlen; len=(int)(len*1.25)+1) {

    /* run more transfers for small messages to get more accurate time measurement */
    if(len < 1000)
      nreps = 10000;
    else
      nreps = 100;

    for (i=1; i<=nreps; i++) {
      MPI_Recv(message, len, MPI_DOUBLE, 0, ping, MPI_COMM_WORLD, &status);
      MPI_Ssend(message, len, MPI_DOUBLE, 0, pong, MPI_COMM_WORLD);
    }
  }

  free(message);
}

double timer()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
