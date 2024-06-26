#include "graph.hpp"
#include <cstdint>
#include <cstdio>
#include <format>
#include <fstream>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include <random>
#include <vector>
#include <zpp_bits.h>

int main(int argc, char *argv[]) {
  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  //
  // SENDER
  //

  if (rank == 0) {
    std::vector<std::byte> data;

    // initialize a path
    graph::path paths_sender;
    paths_sender.vertices = {1, 2, 3, 4};
    paths_sender.weight = 10;

    std::vector<graph::path> v1 = {paths_sender};
    auto out_sender = zpp::bits::out(data);
    (void)out_sender(v1);

    MPI_Send(data.data(), data.size(), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
  }

  //
  // RECEIVER
  //

  if (rank == 1) {
    MPI_Status status;
    int data_size;

    // Probe for the incoming message to get the size
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_BYTE, &data_size);
    printf("data size is %d\n", data_size);

    std::vector<std::byte> data;

    std::vector<graph::path> paths_receiver;
    auto out = zpp::bits::in(data);
    (void)out(paths_receiver);

    printf("weight is %d\n", paths_receiver[0].weight);
    printf("length of vector<bytes> is %lu\n", data.size());
  }

  //
  //
  //

  // read the graph
  std::ifstream graph_file("example.json");
  nlohmann::json graph_data = nlohmann::json::parse(graph_file);
  auto graph = graph_data.template get<graph::reversibleVecGraph>();

  // vector to save all calculated paths
  std::vector<graph::path> paths = std::vector<graph::path>();

  // Create a random device and use it to seed the generator
  std::random_device rd;
  std::mt19937 gen(rd());

  // Define the distribution range. Max is inclusive therfore -1
  std::uniform_int_distribution<> dis(0, graph.number_of_vertices - 1);

  for (int vertex = rank; vertex < graph.number_of_vertices;
       vertex += num_procs) {

    for (int reps = 0; reps < 10; ++reps) {
      uint32_t source = dis(gen);
      uint32_t target = dis(gen);
      auto path = graph.dijkstra(source, target);

      if (path) {
        paths.push_back(*path);
      }
    }
  }

  nlohmann::json j = paths;
  std::ofstream out_file(std::format("paths_%s.json", processor_name));
  out_file << j;

  MPI_Finalize();
}
