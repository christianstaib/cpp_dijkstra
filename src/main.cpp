// Compile with: mpicc -fopenmp hello_hybrid.c -o hello_hybrid

#include "graph.hpp"
#include <fstream>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  if (rank == 0) {
    printf("There are %d processes in total\n", num_procs);
  }

  printf("loading file\n");
  std::ifstream file("example.json");
  nlohmann::json data = nlohmann::json::parse(file);

  auto p2 = data.template get<graph::reversibleVecGraph>();
  printf("there are %d vertices\n", p2.number_of_vertices);

#pragma omp parallel
  {
#pragma omp for
    for (int vertex = rank; vertex < p2.number_of_vertices;
         vertex += num_procs) {

      int thread_nr = omp_get_thread_num();

      printf("%d\n t:%d, p:%d, n:%s", vertex, thread_nr, rank, processor_name);
      int32_t weight = p2.dijkstra(vertex)[123];
    }
  }

  MPI_Finalize();
}
