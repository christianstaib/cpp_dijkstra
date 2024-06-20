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
  int rc;
  int mpi_rank_global, mpi_numPEs_global, procnamelen;
  int mpiversion, mpisubversion;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm mpi_comm_global = MPI_COMM_WORLD;

  /* Initialize MPI */
  rc = MPI_Init(&argc, &argv);

  rc = MPI_Comm_rank(mpi_comm_global, &mpi_rank_global);
  rc = MPI_Comm_size(mpi_comm_global, &mpi_numPEs_global);

  /* Get info about MPI version, Communicator size, rank, processor name */
  rc = MPI_Get_version(&mpiversion, &mpisubversion);
  if (!mpi_rank_global)
    printf("MPI version: %d.%d\n", mpiversion, mpisubversion);

  rc = MPI_Get_processor_name(processor_name, &procnamelen);
  if (rc != MPI_SUCCESS)
    rc = MPI_Abort(MPI_COMM_WORLD, 1);
  printf("Hello from %s rank %d, member of a communicator with %d PEs\n",
         processor_name, mpi_rank_global, mpi_numPEs_global);

  /* Finalize MPI */
  MPI_Finalize();

  //  MPI_Init(&argc, &argv);
  //  int myrank;
  //  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //  printf("my rank is %d", myrank);
  //
  //  std::ifstream f("example.json");
  //  json data = json::parse(f);
  //
  //  auto p2 = data.template get<graph::reversibleVecGraph>();
  //
  //  std::cout << p2.dijkstra(123)[1270] << std::endl;
  //
  //  uint32_t num_fin_vertex = 0;
  //
  //  auto start = std::chrono::high_resolution_clock::now();
  //
  // #pragma omp parallel
  //  {
  // #pragma omp for
  //    for (int i = 0; i < p2.number_of_vertices; i++) {
  //      int x = p2.dijkstra(i)[1270];
  //
  // #pragma omp critical
  //      {
  //        num_fin_vertex += 1;
  //
  //        if (num_fin_vertex % 1000 == 0) {
  //
  //          float progress_percent =
  //              (float)num_fin_vertex / (float)p2.number_of_vertices * 100.0;
  //
  //          printf("%f%%\n", progress_percent);
  //        }
  //      }
  //    }
  //  }
  return 0;
}
