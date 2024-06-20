#include "graph.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main() {

  std::ifstream f("example.json");
  json data = json::parse(f);

  auto p2 = data.template get<graph::reversibleVecGraph>();

  std::cout << p2.dijkstra(123)[1270] << std::endl;

  uint32_t num_fin_vertex = 0;

  auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < p2.number_of_vertices; i++) {
      int x = p2.dijkstra(i)[1270];

#pragma omp critical
      {
        num_fin_vertex += 1;

        if (num_fin_vertex % 1000 == 0) {

          float progress_percent =
              (float)num_fin_vertex / (float)p2.number_of_vertices * 100.0;

          printf("%f%%\n", progress_percent);
        }
      }
    }
  }
  return 0;
}
