#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
  int rank, p;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  int N = 1000000;
  int N_I = 10000;
  int* array = (int*) malloc(N * sizeof(int)); 
  //int* array_out = (int*) malloc(N * sizeof(int)); 
  for (int j = 0; j < N; j++) {array[j] = 12;}
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  int message_out;
  int message_in;

  double tt = MPI_Wtime();
  for (int i = 0; i< N_I;i++){
  
  if (rank != 0) {
    MPI_Status status;

    // MPI_Recv(&message_in,  1, MPI_INT, rank-1, i, MPI_COMM_WORLD, &status);
    // message_out = message_in + rank;
    MPI_Recv(array,  N, MPI_INT, rank-1, i, MPI_COMM_WORLD, &status);
    //array_out = array_in;
  } //else message_out = 0;

  //MPI_Send(&message_out, 1, MPI_INT, (rank+1)% p, i, MPI_COMM_WORLD);
  MPI_Send(array, N, MPI_INT, (rank+1)% p, i, MPI_COMM_WORLD);
  
  if (rank == 0){
  MPI_Status status;

  //MPI_Recv(&message_in, 1, MPI_INT, p-1, i, MPI_COMM_WORLD, &status);
  //if (i==N-1) printf("Rank %d in %d received %d\n", rank, p, message_in);
  MPI_Recv(array,  N, MPI_INT, p-1, i, MPI_COMM_WORLD, &status);
  
  }

  }
  
  double elapsed = MPI_Wtime() - tt;
  if (rank == 0) {
    printf("Time elapsed is %f seconds.\n", elapsed);
    printf("bandwidth: %e GB/s\n", (N_I * N* sizeof(int)*p)/elapsed/1e9);
  }

  free(array);
  //free(array_out);

  MPI_Finalize();
  //printf("asefacwecaserfeasfssef\n");

  return 0;
}
