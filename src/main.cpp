// Your First C++ Program

#include "octree.hpp"
#include "space.hpp"
#include <cstdio>
#include <fstream>
#include <glm/ext/vector_double3.hpp>
#include <glm/gtx/component_wise.hpp>
#include <vector>

int main() {
  // std::ifstream file("planets_and_moons_state_vectors.csv");
  std::ifstream file("test_data.csv");

  int i = 0;

  std::vector<space::CelestialBody> bodies;

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      // using printf() in all tests for consistency
      space::CelestialBody body(line);
      // printf("%s\n", body.name.c_str());
      bodies.push_back(body);
    }
    file.close();
  }

  glm::dvec3 min_pos(std::numeric_limits<double>::max());
  glm::dvec3 max_pos(std::numeric_limits<double>::min());

  for (auto body : bodies) {
    min_pos = min(min_pos, body.pos);
    max_pos = max(max_pos, body.pos);
  }
  printf("min %f %f %f\n", min_pos.x, min_pos.y, min_pos.z);
  printf("max %f %f %f\n", max_pos.x, max_pos.y, max_pos.z);

  glm::dvec3 center = (max_pos + min_pos) / 2.0;
  printf("center %f %f %f\n", center.x, center.y, center.z);
  double size = compMax(max_pos - min_pos);
  printf("size %f\n", size);

  octree::Octree tree(center, size);
  for (auto body : bodies) {
    tree.insert(body.pos, body.mass);
  }

  int node_idx = 0;

  std::vector<octree::Node> nodes;
  nodes.push_back(tree.nodes[0]);
  while (!nodes.empty()) {
    octree::Node node = nodes.back();
    nodes.pop_back();
    if (node.is_branch()) {
      for (int i = 0; i < 8; ++i) {
        nodes.push_back(tree.nodes[node.first_child + i]);
      }
    }
    if (node.is_leaf() && !node.is_empty()) {
      if (!node.cube.contains(node.mass_center)) {
        printf("illegal3\n");
      }
      printf("(%f, %f, %f) %f\n", node.cube.center.x, node.cube.center.y,
             node.cube.center.z, node.cube.size);
    }
  }

  return 0;
}
