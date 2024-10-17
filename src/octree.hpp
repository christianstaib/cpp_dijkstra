#pragma once

#include <array>
#include <glm/ext/vector_double3.hpp>
#include <glm/glm.hpp> // Include GLM for vec3
#include <vector>

// origin is bottom left behind
// (0 or x: left -> right), (1 or y: bottom -> top), (2 or z: behind -> front)

namespace octree {
class Cube {
public:
  glm::dvec3 center;
  double size;

  // only gives a hint in which direction the subcube is. is not guranteed that
  // pos fits
  int find_subcube(glm::dvec3 pos);
  bool contains(glm::dvec3 pos);
  Cube create_subcube(int quadrant);
  std::array<Cube, 8> subdivide();
};

class Node {
public:
  int first_child;
  Cube cube;
  int next_pre_order;

  glm::dvec3 mass_center; // Center of mass as a 3D vector
  double mass;

  Node(Cube cube, int next_pre_order);
  bool is_leaf();
  bool is_branch();
  bool is_empty();
};

class Octree {
public:
  std::vector<Node> nodes;
  std::vector<int> parents; // parent of node[i] is nodes[parents[(i - 1)/8]]

  Octree(glm::dvec3 center, double size);
  void propagate();
  glm::dvec3 acc(glm::dvec3 pos, double theta);
  void insert(glm::dvec3 pos, double mass);
  int subdivide(int node);
};

} // namespace octree
