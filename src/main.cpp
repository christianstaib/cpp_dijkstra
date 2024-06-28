#include "graph.hpp"
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <random>
#include <string>
#include <vector>
#include <zpp_bits.h>

int main(int argc, char *argv[]) {

  int num_procs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  std::string graph_file_path = argv[1];
  std::string graph_file_name =
      std::filesystem::path(graph_file_path).filename();
  int num_paths = std::stoi(argv[2]) / num_procs;

  if (rank == 0) {
    printf("expecting arround %d paths\n", num_paths * num_procs);
  }

  //
  // setup graph over all nodes
  //

  std::vector<std::byte> data;

  if (rank == 0) {
    // read the graph
    graph::reversibleVecGraph graph;
    std::ifstream graph_file(graph_file_path);
    nlohmann::json graph_data = nlohmann::json::parse(graph_file);
    graph = graph_data.template get<graph::reversibleVecGraph>();

    // serialize the graph
    auto out_sender = zpp::bits::out(data);
    (void)out_sender(graph);
  }

  graph::reversibleVecGraph graph;
  {
    // broadcast size of graph
    uint64_t data_size = data.size();
    MPI_Bcast(&data_size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      printf("broadcasting graph\n");
    }
    // broadcast graph
    data.resize(data_size);
    MPI_Bcast(data.data(), data_size, MPI_BYTE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      printf("deserialing graph\n");
    }
    // deserialze graph
    auto in_receiver = zpp::bits::in(data);
    (void)in_receiver(graph);
  }

  //
  // calculates random paths
  //

  // vector to save all calculated paths
  std::vector<graph::path> paths = std::vector<graph::path>();

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("starting path calculations\n");
  }

  printf("hi from process %d on cpu %d on node %s\n", rank, sched_getcpu(),
         processor_name);

#pragma omp parallel
  {
    //
    // setup thread local random number generator
    //

    // Create a random device and use it to seed the generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution range. Max is inclusive therfore -1
    std::uniform_int_distribution<> dis(0, graph.number_of_vertices - 1);

    //  printf("hi from thread %d of %d from process %d of %d. This thread is "
    //         "running on cpu %d on node %s\n",
    //         omp_get_thread_num(), omp_get_num_threads(), rank, num_procs,
    //         sched_getcpu(), processor_name);

#pragma omp for schedule(dynamic)
    for (int vertex = rank; vertex < num_paths; vertex += num_procs) {

      uint32_t source = dis(gen);
      uint32_t target = dis(gen);
      auto path = graph.dijkstra(source, target);

      if (path) {
#pragma omp critical
        paths.push_back(*path);
      }
    }
  }

  if (rank != 0) {
    //
    // send local paths to rank 0
    //
    auto out_sender = zpp::bits::out(data);
    (void)out_sender(paths);
    MPI_Send(data.data(), data.size(), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
  } else {
    //
    // receive local paths
    //
    for (int i = 0; i < num_procs - 1; ++i) {
      // probe incomming message
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

      // retrieve size of incommin message
      int data_size = 0;
      MPI_Get_count(&status, MPI_BYTE, &data_size);
      data.resize(data_size);

      // reseive incomming message
      MPI_Recv(data.data(), data_size, MPI_BYTE, MPI_ANY_SOURCE, 0,
               MPI_COMM_WORLD, &status);

      // deserialise pathes
      std::vector<graph::path> received_paths;
      auto in_receiver = zpp::bits::in(data);
      (void)in_receiver(received_paths);

      printf("there are %lu paths\n", paths.size());
      paths.insert(paths.end(), received_paths.begin(), received_paths.end());
    }

    printf("there are %lu paths\n", paths.size());
    nlohmann::json j = paths;
    std::string out_file_name = std::format("{}.paths.json", graph_file_name);
    std::ofstream out_file(out_file_name.c_str());
    out_file << j;
  }

  MPI_Finalize();
  return 0;
}

std::vector<unsigned int> divide_into_parts(unsigned int number,
                                            unsigned int parts) {
  std::vector<unsigned int> result(parts, number / parts);
  unsigned int remainder = number % parts;

  for (unsigned int i = 0; i < remainder; ++i) {
    ++result[i];
  }

  return result;
}
