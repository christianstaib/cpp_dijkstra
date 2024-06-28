#include <cstdio>
#include <mpi.h>
#include <omp.h>
#include <sched.h>

int main(int argc, char *argv[]) {
  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  printf("there are %d processes\n", rank);

#pragma omp parallel
  {
    printf("hi from thread %d of %d on process %d on node %s \n",
           omp_get_thread_num(), omp_get_num_threads(), sched_getcpu(),
           processor_name);

    // dummy workload
    int i = 0;
    while (true) {
    }
    ++i;
  }

  printf("finish %d \n", omp_get_thread_num());

  MPI_Finalize();
}
