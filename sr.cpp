/* Simple send-receive example */
#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
  int rank, p;
  int N = 1000000;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  double tt = MPI_Wtime();

  for (int i = 0; i< N;i++){

  if (rank == 0) {
    int message_out = 0;
    int message_in;

    MPI_Status status;
    //MPI_Request request_out, request_in;

    MPI_Send(&message_out, 1, MPI_INT, 1, 999, MPI_COMM_WORLD);

    MPI_Recv(&message_in,  1, MPI_INT, p-1, 999, MPI_COMM_WORLD, &status);

    printf("Rank %d received %d\n", rank, message_in);

  } else if (rank == p-1) {
    int message_out = 0;
    int message_in;

    MPI_Status status;

    MPI_Recv(&message_in, 1, MPI_INT, p-2, 999, MPI_COMM_WORLD, &status);
    message_out = message_in + p-1;

    MPI_Send(&message_out, 1, MPI_INT, 0, 999, MPI_COMM_WORLD);

    //printf("The message is %d\n", message_in);
  }
  int message_out = 0;
  int message_in;

  MPI_Status status;

  MPI_Recv(&message_in, 1, MPI_INT, rank-1, 999, MPI_COMM_WORLD, &status);
  message_out = message_in + rank;

  MPI_Send(&message_out, 1, MPI_INT, rank+1, 999, MPI_COMM_WORLD);

  }

  double elapsed = MPI_Wtime() - tt;
  if (0 == mpirank) {
    printf("Time elapsed is %f seconds.\n", elapsed);
  }

  MPI_Finalize();

  return 0;
}
