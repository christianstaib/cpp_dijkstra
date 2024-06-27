#include <cstdint>
#include <cstdio>
#include <mpi.h>
#include <omp.h>
#include <random>
#include <sched.h>

int main(int argc, char *argv[]) {
  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  printf("starting path calculations\n");
#pragma omp parallel
  {
    //
    // setup thread local random number generator
    //

    // Create a random device and use it to seed the generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution range. Max is inclusive therfore -1
    std::uniform_int_distribution<> dis(0, 99999);

    printf("hi from thread %d of %d on %d \n", omp_get_thread_num(),
           omp_get_num_threads(), sched_getcpu());

#pragma omp for
    for (int vertex = rank; vertex < 100000000; vertex += num_procs) {
      for (int reps = 0; reps < 1; ++reps) {
        uint32_t source = dis(gen);
        uint32_t target = dis(gen);

        // dummy work to keep threads running
        for (; target > 0; --target) {
        }
      }
    }
  }

  printf("finish %d \n", omp_get_thread_num());

  MPI_Finalize();
}
