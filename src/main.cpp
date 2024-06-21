// Compile with: mpicc -fopenmp hello_hybrid.c -o hello_hybrid

#include "graph.hpp"
#include <cstdint>
#include <fstream>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include <omp.h>
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

  if (rank == 1) {
    graph::path p;
    p.vertices = {1, 2, 3, 4};
    p.weight = 10;

    nlohmann::json j = p;
    auto x = j.dump();

    MPI_Send(x.c_str(), x.size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  } else if (rank == 0) {
    MPI_Status status;
    MPI_Probe(1, 0, MPI_COMM_WORLD, &status);

    int message_size;
    MPI_Get_count(&status, MPI_CHAR, &message_size);

    std::vector<char> buffer(message_size);
    MPI_Recv(buffer.data(), message_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    printf("%s", buffer.data());
  }

  // read the graph
  std::ifstream graph_file("example.json");
  nlohmann::json graph_data = nlohmann::json::parse(graph_file);
  auto graph = graph_data.template get<graph::reversibleVecGraph>();

  // vector to save all calculated paths
  std::vector<graph::path> paths = std::vector<graph::path>();

#pragma omp parallel
  {
#pragma omp for
    for (int vertex = rank; vertex < graph.number_of_vertices;
         vertex += num_procs) {

      auto path = graph.dijkstra(vertex, 1234);

      if (path) {
#pragma omp critical
        {
          // only one write at a time
          paths.push_back(*path);
        }
      }
    }
  }

  if (rank == 0) {
    std::vector<graph::path> all_paths(0);
  }

  nlohmann::json j = paths;
  std::ofstream out_file("key.json");
  out_file << j;

  MPI_Finalize();
}
