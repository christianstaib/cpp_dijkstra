// Compile with: mpicc -fopenmp hello_hybrid.c -o hello_hybrid

#include "graph.hpp"
#include <cstdint>
#include <format>
#include <fstream>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <random>
#include <stdio.h>
#include <vector>
#include <zpp_bits.h>

int main(int argc, char *argv[]) {
  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  // make sure we are set up corretly
  printf("hello from rank %d on %s\n", rank, processor_name);

  // read the graph
  std::ifstream graph_file("example.json");
  nlohmann::json graph_data = nlohmann::json::parse(graph_file);
  auto graph = graph_data.template get<graph::reversibleVecGraph>();

  // vector to save all calculated paths
  std::vector<graph::path> paths = std::vector<graph::path>();

#pragma omp parallel
  {
    printf("hello from thread %d rank %d on %s\n", omp_get_thread_num(), rank,
           processor_name);
    // Create a random device and use it to seed the generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution range. Max is inclusive therfore -1
    std::uniform_int_distribution<> dis(0, graph.number_of_vertices - 1);

#pragma omp for
    for (int vertex = rank; vertex < graph.number_of_vertices;
         vertex += num_procs) {

      printf("hello from thread %d rank %d on %s\n", omp_get_thread_num(), rank,
             processor_name);

      for (int reps = 0; reps < 10; ++reps) {
        uint32_t source = dis(gen);
        uint32_t target = dis(gen);
        auto path = graph.dijkstra(source, target);

        if (path) {
          // only one write at a time, as vector is not threadsafe
#pragma omp critical
          { paths.push_back(*path); }
        }
      }
    }
  }

  nlohmann::json j = paths;
  std::ofstream out_file(std::format("paths_%s.json", processor_name));
  out_file << j;

  MPI_Finalize();
}
