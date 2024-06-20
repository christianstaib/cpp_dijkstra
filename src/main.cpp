#include <mpi.h>
#include <omp.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int thread_nr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  if (rank == 0) {
    printf("There are %d processes in total\n", num_procs);
  }

#pragma omp parallel default(shared) private(thread_nr)
  {
    // #pragma omp single
    //     {
    //       int num_threads = omp_get_num_threads();
    //       printf("There are %d threads of process %d on %s\n", num_threads,
    //       rank,
    //              processor_name);
    //     }

    printf("x\n");

    // thread_nr = omp_get_thread_num();
    // printf("Hello World from thread %d from process %d on node %s\n",
    // thread_nr,
    //        rank, processor_name);
  }

  MPI_Finalize();
}
