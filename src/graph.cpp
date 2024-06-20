#include "graph.hpp"

#include <queue>

std::vector<uint32_t>
graph::reversibleVecGraph::dijkstra(uint32_t start) const {

  // Initialize distances to infinity
  std::vector<uint32_t> distances(number_of_vertices,
                                  std::numeric_limits<uint32_t>::max());
  distances[start] = 0;

  // Min-heap priority queue to select the edge with the smallest weight
  using Node = std::pair<uint32_t, uint32_t>; // (distance, vertex)
  std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
  pq.push({0, start});

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
        pq.push({new_distance, next_vertex});
      }
    }
  }

  return distances;
}
