#pragma once

#include <cstdint>
#include <nlohmann/json.hpp>
#include <vector>

namespace graph {
class weightedDirectedEdge {
public:
  uint32_t tail;
  uint32_t head;
  uint32_t weight;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(weightedDirectedEdge, tail, head, weight);
};

class directedTaillessWeightedEdge {
public:
  uint32_t head;
  uint32_t weight;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(directedTaillessWeightedEdge, head, weight);
};

class directedHeadlessWeightedEdge {
public:
  uint32_t tail;
  uint32_t weight;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(directedHeadlessWeightedEdge, tail, weight);
};

class reversibleVecGraph {
public:
  std::vector<std::vector<directedTaillessWeightedEdge>> out_edges;
  std::vector<std::vector<directedHeadlessWeightedEdge>> in_edges;
  uint32_t number_of_vertices;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(reversibleVecGraph, out_edges, in_edges,
                                 number_of_vertices);

  std::vector<uint32_t> dijkstra(uint32_t start) const;
};
} // namespace graph
