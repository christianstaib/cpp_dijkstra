#include "graph.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int myrank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

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
