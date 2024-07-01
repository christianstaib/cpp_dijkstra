#pragma once

#include <cstdint>
#include <nlohmann/json.hpp>
#include <optional>
#include <vector>
#include <zpp_bits.h>

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

class shortestPathRequest {
public:
  uint32_t source;
  uint32_t target;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(shortestPathRequest, source, target);
};

class shortestPathTestCaseC {
public:
  uint32_t source;
  uint32_t target;
  uint32_t weight;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(shortestPathTestCaseC, source, target, weight);
};

class path {
public:
  std::vector<uint32_t> vertices;
  uint32_t weight;

  // Serialization function using zpp_bits
  template <typename Archive> void serialize(Archive &archive) {
    archive(vertices, weight);
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(path, vertices, weight);
};

class reversibleVecGraph {
public:
  std::vector<std::vector<directedTaillessWeightedEdge>> out_edges;
  std::vector<std::vector<directedHeadlessWeightedEdge>> in_edges;
  uint32_t number_of_vertices;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(reversibleVecGraph, out_edges, in_edges,
                                 number_of_vertices);

  std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
  dijkstra(uint32_t source) const;

  std::optional<path> dijkstra(uint32_t source, uint32_t target) const;
};
} // namespace graph
