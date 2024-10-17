
#include "octree.hpp"
#include <array>
#include <cstdio>
#include <glm/ext/vector_double3.hpp>
#include <ranges>
#include <vector>

bool octree::Node::is_leaf() { return first_child == 0; }

bool octree::Node::is_branch() { return first_child != 0; }

bool octree::Node::is_empty() { return mass == 0.0; }

bool octree::Cube::contains(glm::dvec3 pos) {
  const double epsilon = 1e-10; // Small tolerance for floating-point precision
  return (center.x - size - epsilon <= pos.x) &&
         (pos.x < center.x + size + epsilon) &&
         (center.y - size - epsilon <= pos.y) &&
         (pos.y < center.y + size + epsilon) &&
         (center.z - size - epsilon <= pos.z) &&
         (pos.z < center.z + size + epsilon);
}

int octree::Cube::find_subcube(glm::dvec3 pos) {
  int subcube = 0;

  if (pos.x > center.x) {
    subcube += 1;
  }

  if (pos.y > center.y) {
    subcube += 2;
  }

  if (pos.z > center.z) {
    subcube += 4;
  }

  return subcube;
}

octree::Node::Node(octree::Cube cube, int next_pre_order)
    : first_child(0), cube(cube), next_pre_order(next_pre_order),
      mass_center(glm::dvec3(0.0)), mass(0.0) {}

octree::Cube octree::Cube::create_subcube(int quadrant) {
  octree::Cube subcube = octree::Cube{center, size / 2};
  subcube.center.x += ((float)(quadrant & 0b1) - 0.5) * size;
  subcube.center.y += ((float)((quadrant >> 1) & 0b1) - 0.5) * size;
  subcube.center.z += ((float)((quadrant >> 2) & 0b1) - 0.5) * size;

  return subcube;
}

octree::Octree::Octree(glm::dvec3 center, double size) {
  octree::Cube root{center, size};
  nodes.push_back(Node(root, 0));
};

std::array<octree::Cube, 8> octree::Cube::subdivide() {
  std::array<octree::Cube, 8> subcubes;

  // Iterate over all 8 possible quadrants (0 to 7)
  for (int quadrant = 0; quadrant < 8; ++quadrant) {
    subcubes[quadrant] = create_subcube(quadrant);
  }

  return subcubes;
}

int octree::Octree::subdivide(int node) {
  parents.push_back(node);
  int first_child = nodes.size();
  nodes[node].first_child = first_child;

  int nexts_pre_order[8] = {first_child + 1, first_child + 2,
                            first_child + 3, first_child + 4,
                            first_child + 5, first_child + 6,
                            first_child + 7, nodes[node].next_pre_order};

  std::array<octree::Cube, 8> subcubes = nodes[node].cube.subdivide();
  for (int i = 0; i < 8; ++i) {
    nodes.push_back(octree::Node(subcubes[i], nexts_pre_order[i]));
  }
  return first_child;
}

void octree::Octree::insert(glm::dvec3 new_pos, double new_mass) {
  int node_idx = 0;
  while (nodes[node_idx].is_branch()) {
    node_idx = nodes[node_idx].first_child +
               nodes[node_idx].cube.find_subcube(new_pos);
  }

  // if there is no body in the node, set it
  if (nodes[node_idx].is_empty()) {
    nodes[node_idx].mass_center = new_pos;
    nodes[node_idx].mass = new_mass;
    return;
  }

  glm::dvec3 old_pos = nodes[node_idx].mass_center;
  double old_mass = nodes[node_idx].mass;

  while (true) {
    int first_child = subdivide(node_idx);

    int offset_new = nodes[node_idx].cube.find_subcube(new_pos);
    int offset_old = nodes[node_idx].cube.find_subcube(old_pos);

    if (offset_new == offset_old) {
      node_idx = first_child + offset_new;
    } else {
      octree::Node *node = &nodes[first_child + offset_new];
      node->mass_center = new_pos;
      node->mass = new_mass;
      if (!node->cube.contains(new_pos)) {
        printf("illegal1\n");
      }

      node = &nodes[first_child + offset_old];
      node->mass_center = old_pos;
      node->mass = old_mass;
      if (!node->cube.contains(nodes[node_idx].mass_center)) {
        printf("illegal2\n");
      }
      break;
    }
  }
}

// void octree::Octree::propagate() {
//   for (auto &node : std::ranges::views::reverse(this->parents)) {
//     int first_child = this->nodes[node].first_child;
//
//     for (int child_offset = 0; child_offset < 8; ++child_offset) {
//       nodes[node].mass += nodes[first_child + child_offset].mass;
//     }
//
//     nodes[node].mass_center = glm::dvec3(0.0);
//     for (int child_offset = 0; child_offset < 8; ++child_offset) {
//       nodes[node].mass_center += nodes[first_child +
//       child_offset].mass_center *
//                                  nodes[first_child + child_offset].mass;
//     }
//     nodes[node].mass_center /= nodes[node].mass;
//   }
// }
