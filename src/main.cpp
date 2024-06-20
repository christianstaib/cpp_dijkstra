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

  printf("rank:%d num_procs:%d\n", rank, num_procs);

  // printf("loading file\n");
  // std::ifstream file("example.json");
  // nlohmann::json data = nlohmann::json::parse(file);

  // auto p2 = data.template get<graph::reversibleVecGraph>();
  // printf("there are %d vertices\n", p2.number_of_vertices);

  // #pragma omp parallel
  //   {
  // #pragma omp for
  //     for (int vertex = rank; vertex < p2.number_of_vertices;
  //          vertex += num_procs) {
  //
  //       // if (vertex % 1000 == 0) {
  //       printf("%f\n", (float)vertex / (float)p2.number_of_vertices * 100.0);
  //       // }
  //
  //       int thread_nr = omp_get_thread_num();
  //       int num_threads = omp_get_num_threads();
  //
  //       // printf("%d t:%d, p:%d, n:%s, num_threads:%d\n", vertex, thread_nr,
  //       // rank,
  //       //        processor_name, num_threads);
  //       int32_t weight = p2.dijkstra(vertex)[123];
  //     }
  //   }

  MPI_Finalize();
}
