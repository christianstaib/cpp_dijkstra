#include <random>
#include <utility>
#define GLM_ENABLE_EXPERIMENTAL

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <glm/ext/vector_double3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/norm.hpp>
#include <ranges>
#include <vector>

#include "constants.hpp"
#include "octree.hpp"

namespace octree {

bool Node::is_leaf() { return first_child == 0; }

bool Node::is_branch() { return first_child != 0; }

bool Node::is_empty() { return mass == 0.0; }

bool Cube::contains(glm::dvec3 pos) {
  const double epsilon = 1e-5;  // Small tolerance for floating-point precision
  return (center.x - half_edge_length - epsilon <= pos.x) && (pos.x < center.x + half_edge_length + epsilon) &&
         (center.y - half_edge_length - epsilon <= pos.y) && (pos.y < center.y + half_edge_length + epsilon) &&
         (center.z - half_edge_length - epsilon <= pos.z) && (pos.z < center.z + half_edge_length + epsilon);
}

Cube::Cube() : center(glm::dvec3(0)), half_edge_length(0.0), squared_edge_length(0.0) {}

Cube::Cube(glm::dvec3 center, double half_edge_length)
    : center(center),
      half_edge_length(half_edge_length),
      squared_edge_length((2 * half_edge_length) * (2 * half_edge_length)) {}

int Cube::find_subcube(glm::dvec3 pos) {
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

Node::Node(octree::Cube cube, int next_pre_order)
    : first_child(0), cube(cube), next_pre_order(next_pre_order), mass_center(glm::dvec3(0.0)), mass(0.0) {}

Cube Cube::create_subcube(int quadrant) {
  Cube subcube = octree::Cube(center, half_edge_length / 2);
  subcube.center.x += ((float)(quadrant & 0b1) - 0.5) * half_edge_length;
  subcube.center.y += ((float)((quadrant >> 1) & 0b1) - 0.5) * half_edge_length;
  subcube.center.z += ((float)((quadrant >> 2) & 0b1) - 0.5) * half_edge_length;

  return subcube;
}

Octree::Octree(glm::dvec3 center, double size) {
  Cube root{center, size};
  nodes.push_back(Node(root, 0));
};

void Octree::clear(glm::dvec3 center, double half_edge_length) {
  nodes.clear();
  parents.clear();

  octree::Cube root(center, half_edge_length);
  nodes.push_back(Node(root, 0));
}

std::array<Cube, 8> Cube::subdivide() {
  std::array<Cube, 8> subcubes;

  // Iterate over all 8 possible quadrants (0 to 7)
  for (int quadrant = 0; quadrant < 8; ++quadrant) {
    subcubes[quadrant] = create_subcube(quadrant);
  }

  return subcubes;
}

int Octree::subdivide(int node) {
  parents.push_back(node);
  int first_child = nodes.size();
  nodes[node].first_child = first_child;

  int nexts_pre_order[8] = {first_child + 1, first_child + 2, first_child + 3, first_child + 4,
                            first_child + 5, first_child + 6, first_child + 7, nodes[node].next_pre_order};

  std::array<Cube, 8> subcubes = nodes[node].cube.subdivide();
  for (size_t i = 0; i < 8; ++i) {
    nodes.push_back(Node(subcubes[i], nexts_pre_order[i]));
  }
  return first_child;
}

void Octree::insert(glm::dvec3 new_pos, double new_mass) {
  int node_idx = 0;
  while (nodes[node_idx].is_branch()) {
    node_idx = nodes[node_idx].first_child + nodes[node_idx].cube.find_subcube(new_pos);
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
      Node *node = &nodes[first_child + offset_new];
      node->mass_center = new_pos;
      node->mass = new_mass;

      node = &nodes[first_child + offset_old];
      node->mass_center = old_pos;
      node->mass = old_mass;

      break;
    }
  }
}

void VecTreeNode::split() {
#pragma omp declare reduction( \
        merge : std::vector<VecTreeNode> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

  std::vector<VecTreeNode> to_split;

  for (int i = 0; i < 8; ++i) {
    to_split.push_back((*children)[i]);
  }

  while (!to_split.empty()) {
    VecTreeNode node = to_split.back();
    to_split.pop_back();

    node.children = new std::array<VecTreeNode, 8>();

    std::array<Cube, 8> cubes = node.cube.subdivide();
    for (int i = 0; i < 8; ++i) {
      (*node.children)[i].cube = cubes[i];
    }

#pragma omp parallel
    {
      std::array<std::vector<std::pair<glm::dvec3, double>>, 8> local_data;
#pragma omp for schedule(static)
      for (int i = 0; i < node.data.size(); ++i) {
        size_t idx = node.cube.find_subcube(node.data[i].first);
        local_data[idx].push_back({node.data[i]});
      }

#pragma omp critical
      {
        for (int i = 0; i < 8; ++i) {
          (*node.children)[i].data.insert((*node.children)[i].data.end(), local_data[i].begin(), local_data[i].end());
        }
      }
    }

    for (int i = 0; i < 8; ++i) {
      if ((*node.children)[i].data.size() > 100) {
        to_split.push_back((*node.children)[i]);
      }
    }

    node.data.clear();
  }
}

VecTreeNode VecTreeNode::create_root(size_t num_bodies, glm::dvec3 *positions, double *masses, Cube cube) {
  VecTreeNode root{};
  root.cube = cube;
  root.children = new std::array<VecTreeNode, 8>();

  std::array<Cube, 8> cubes = root.cube.subdivide();
  for (int i = 0; i < 8; ++i) {
    (*root.children)[i].cube = cubes[i];
  }

#pragma omp parallel
  {
    std::array<std::vector<std::pair<glm::dvec3, double>>, 8> local_data;
#pragma omp for schedule(guided)
    for (int i = 0; i < num_bodies; ++i) {
      size_t idx = root.cube.find_subcube(positions[i]);
      local_data[idx].push_back({positions[i], masses[i]});
    }

#pragma omp critical
    {
      for (int i = 0; i < 8; ++i) {
        (*root.children)[i].data.insert((*root.children)[i].data.end(), local_data[i].begin(), local_data[i].end());
      }
    }
  }

  root.data.clear();

  return root;
}

void Octree::propagate() {
  for (auto &parent : std::ranges::views::reverse(parents)) {
    int first_child = nodes[parent].first_child;

    nodes[parent].mass = 0;
    nodes[parent].mass_center = glm::dvec3(0.0);
    for (size_t child_offset = 0; child_offset < 8; ++child_offset) {
      nodes[parent].mass += nodes[first_child + child_offset].mass;
      nodes[parent].mass_center +=
          nodes[first_child + child_offset].mass_center * nodes[first_child + child_offset].mass;
    }
    nodes[parent].mass_center /= nodes[parent].mass;
  }
}

glm::dvec3 Octree::get_force(glm::dvec3 pos, double squared_theta) {
  glm::dvec3 force(0.0), subforce(0.0);
  int node_idx = 0;
  double squared_length_subforce, magnitude;

  while (true) {
    Node &node = nodes[node_idx];
    squared_length_subforce = glm::distance2(node.mass_center, pos);

    // pow(s, 2) / pow(e, 2) < pow(t, 2) <=> s / e < t
    if (node.is_leaf() || (node.cube.squared_edge_length / squared_length_subforce) < squared_theta) {
      subforce = node.mass_center - pos;
      if (squared_length_subforce) {
        squared_length_subforce += +constants::squared_softening_factor;
        magnitude = node.mass / (squared_length_subforce * sqrt(squared_length_subforce));
        force += subforce * magnitude;
      }

      if (node.next_pre_order == 0) break;
      node_idx = node.next_pre_order;
    } else {
      node_idx = node.first_child;
    }
  }
  return force * constants::gravitational_constant_in_au3_per_kg_d2;
}
}  // namespace octree
