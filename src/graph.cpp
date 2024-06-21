#include "graph.hpp"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
graph::reversibleVecGraph::dijkstra(uint32_t source) const {

  // Initialize distances to infinity
  std::vector<uint32_t> distances(number_of_vertices,
                                  std::numeric_limits<uint32_t>::max());
  std::vector<uint32_t> predecessors(number_of_vertices,
                                     std::numeric_limits<uint32_t>::max());
  distances[source] = 0;

  // Min-heap priority queue to select the edge with the smallest weight
  using Node = std::pair<uint32_t, uint32_t>; // (distance, vertex)
  std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
  pq.push({0, source});

  while (!pq.empty()) {
    uint32_t current_distance = pq.top().first;
    uint32_t current_vertex = pq.top().second;
    pq.pop();

    // If the distance is not up to date, skip it
    if (current_distance > distances[current_vertex]) {
      continue;
    }

    // Explore the neighbors
    for (const auto &edge : out_edges[current_vertex]) {
      uint32_t next_vertex = edge.head;
      uint32_t weight = edge.weight;

      // Calculate the distance to the next vertex
      uint32_t new_distance = current_distance + weight;

      // If a shorter path is found
      if (new_distance < distances[next_vertex]) {
        distances[next_vertex] = new_distance;
        predecessors[next_vertex] = current_vertex;
        pq.push({new_distance, next_vertex});
      }
    }
  }

  return {distances, predecessors};
}

std::optional<graph::path>
graph::reversibleVecGraph::dijkstra(uint32_t source, uint32_t target) const {
  auto [distances, predecessors] = dijkstra(source);

  if (distances[target] == std::numeric_limits<uint32_t>::max()) {
    return std::nullopt;
  }

  auto weight = distances[target];

  std::vector<uint32_t> vertices = std::vector<uint32_t>();
  auto current_vertex = target;
  while (current_vertex != std::numeric_limits<uint32_t>::max()) {
    vertices.push_back(current_vertex);
    current_vertex = predecessors[current_vertex];
  }

  std::reverse(vertices.begin(), vertices.end());

  auto path = graph::path{vertices, weight};

  return std::optional<graph::path>{path};
}
