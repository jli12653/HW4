#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(int* prefix_sum, const int* A, long n) {
  if (n == 0) return;
  prefix_sum[0] = 0;
  for (long i = 1; i < n; i++) {
    prefix_sum[i] = prefix_sum[i-1] + A[i-1];
  }
}


int main(int argc, char *argv[]) {
  int rank, p;
  

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  long N = 100000000;

  int* A = NULL;
  int* B0 = NULL;
  int* B1 = NULL;
  double tt = 0;
  
  if (rank == 0 ){
  A = (int*) malloc(N * sizeof(int));
  B0 = (int*) malloc(N * sizeof(int));
  B1 = (int*) malloc(N * sizeof(int));
  for (long i = 0; i < N; i++) A[i] = rand();
  for (long i = 0; i < N; i++) B1[i] = 0;
  
  
  tt = MPI_Wtime();
  scan_seq(B0, A, N);
  printf("sequential-scan = %fs\n", MPI_Wtime() - tt);
  
  
  }
  
  tt = MPI_Wtime();

  long n = N/p;
  int* local = (int*) malloc(n * sizeof(int));
  int* local_sum = (int*) malloc(n * sizeof(int));
  //for (long i = 0; i < n; i++) local_sum[i] = 0;
  int* correction = (int*) malloc(p * sizeof(int));

  MPI_Scatter(A, n, MPI_INT, local, n, MPI_INT, 0, MPI_COMM_WORLD);

  int s = 0;
  local_sum[0] = 0;
  for (long i = 0; i < n-1; i++) {
    s += local[i];
    local_sum[i+1] = s;
  }
  if (rank!= p-1 ) s += local[n-1];

  //MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allgather(&s, 1, MPI_INT, correction, 1, MPI_INT, MPI_COMM_WORLD) ;

  int offset = 0;
    
  for (int i = 0; i < rank; i++){
      offset += correction[i];
  }

  for (long i = 0; i < n; i++) {
    local_sum[i] += offset;
  }

  MPI_Gather(local_sum, n, MPI_INT, B1, n, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0){
  printf("parallel-scan   = %fs\n", MPI_Wtime() - tt);
  

  int err = 0;
  for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
  printf("error = %ld\n", err);

  
  }

  free(A);
  free(B0);
  free(B1);
  free(local);
  free(local_sum);
  free(correction);

  MPI_Finalize();
  return 0;
}
