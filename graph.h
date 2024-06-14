#pragma once

#include <cstdint>
#include <vector>

struct DirectedTaillessWeightedEdge {
  int32_t tail;
  int32_t head;
  int32_t weight;
};

struct VecGraph {
  std::vector<std::vector<DirectedTaillessWeightedEdge>> edges;
};
