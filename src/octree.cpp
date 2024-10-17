
#include "octree.hpp"
#include <array>
#include <glm/ext/vector_double3.hpp>
#include <ranges>

bool octree::Node::is_leaf() { return first_child == 0; }

bool octree::Node::is_branch() { return first_child != 0; }

bool octree::Node::is_empty() { return mass == 0.0; }

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
  subcube.center.x += ((float)(quadrant & 0b001) - 0.5) * subcube.size;
  subcube.center.y += ((float)(quadrant & 0b010) - 0.5) * subcube.size;
  subcube.center.z += ((float)(quadrant & 0b100) - 0.5) * subcube.size;

  return subcube;
}

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

void octree::Octree::propagate() {
  for (auto &node : std::ranges::views::reverse(this->parents)) {
    int first_child = this->nodes[node].first_child;

    for (int child_offset = 0; child_offset < 8; ++child_offset) {
      nodes[node].mass += nodes[first_child + child_offset].mass;
    }

    nodes[node].mass_center = glm::dvec3(0.0);
    for (int child_offset = 0; child_offset < 8; ++child_offset) {
      nodes[node].mass_center += nodes[first_child + child_offset].mass_center *
                                 nodes[first_child + child_offset].mass;
    }
    nodes[node].mass_center /= nodes[node].mass;
  }
}

void octree::Octree::insert(glm::dvec3 pos, double mass) {
  int node_idx = 0;
  while (nodes[node_idx].is_branch()) {
    node_idx =
        nodes[node_idx].first_child + nodes[node_idx].cube.find_subcube(pos);
  }

  // if there is no body in the node, set it
  if (nodes[node_idx].is_empty()) {
    nodes[node_idx].mass_center = pos;
    nodes[node_idx].mass = mass;
    return;
  }

  // if pos to insert is equal to mass_center.
  // TODO remove?
  if (nodes[node_idx].mass_center == pos) {
    nodes[node_idx].mass += mass;
    return;
  }

  while (true) {
    int first_child = subdivide(node_idx);

    int subcube_idx_to_insert = nodes[node_idx].cube.find_subcube(pos);
    int subcube_idx_mass_center =
        nodes[node_idx].cube.find_subcube(nodes[node_idx].mass_center);

    if (subcube_idx_to_insert == subcube_idx_mass_center) {
      node_idx = first_child + subcube_idx_to_insert;
    } else {
      octree::Node *node_to_insert =
          &nodes[first_child + subcube_idx_to_insert];
      node_to_insert->mass_center = pos;
      node_to_insert->mass = mass;

      octree::Node *node_mass_center =
          &nodes[first_child + subcube_idx_mass_center];
      node_mass_center->mass_center = nodes[node_idx].mass_center;
      node_mass_center->mass = nodes[node_idx].mass;
      break;
    }
  }
}
