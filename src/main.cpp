#include "graph.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

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

  int myrank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  printf("loading file\n");
  std::ifstream file("example.json");
  json data = json::parse(file);

  auto p2 = data.template get<graph::reversibleVecGraph>();

  uint32_t num_fin_vertex = 0;

  printf("start\n");
#pragma omp parallel
  {
#pragma omp for
    for (int i = myrank; i < p2.number_of_vertices; i += nproc) {
      int x = p2.dijkstra(i)[1270];

      // #pragma omp critical
      //       {
      //         num_fin_vertex += 1;
      //
      if (i % 1000 == 0) {

        float progress_percent =
            (float)i / (float)p2.number_of_vertices * 100.0;

        printf("%f%%\n", progress_percent);
      }
    }
    //    }
  }
  MPI_Finalize();

  return 0;
}
