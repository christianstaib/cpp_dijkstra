// Your First C++ Program

#include "octree.hpp"
#include "space.hpp"
#include <cstdio>
#include <fstream>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <vector>

int main() {

  //
  // read bodies
  //
  std::vector<space::CelestialBody> bodies;

  std::ifstream file("planets_and_moons_state_vectors.csv");
  // std::ifstream file("test_data.csv");

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      // using printf() in all tests for consistency
      space::CelestialBody body(line);
      body.vel = glm::dvec3(0.0);
      // printf("%s\n", body.name.c_str());
      bodies.push_back(body);
    }
    file.close();
  }

  double gravitational_constant = 6.67430e-11;
  double conversion_factor = pow(86400.0, 2.0) / pow(149597870700.0, 3.0);
  gravitational_constant *= conversion_factor;
  double softening_factor = 1e-11;

  double time_step = 1.0;

  for (int x = 0; x < 400; ++x) {
    std::vector<glm::dvec3> all_acc;
    for (int i = 0; i < bodies.size(); ++i) {
      glm::dvec3 acc(0);
      for (int j = 0; j < bodies.size(); ++j) {
        if (i != j) {
          acc += bodies[j].mass *
                 ((bodies[j].pos - bodies[i].pos) /
                  pow(glm::length2(bodies[j].pos - bodies[i].pos) +
                          softening_factor,
                      (3.0 / 2.0)));
        }
      }
      acc *= gravitational_constant;
      all_acc.push_back(acc);
    }
    for (int i = 0; i < bodies.size(); ++i) {
      auto &body = bodies[i];
      auto acc = all_acc[i];

      body.pos += body.vel * (time_step / 2.0);
      body.vel += acc * (time_step / 2.0);

      if (body.name == "Earth") {
        printf("Earth %f %f %f\n", body.pos.x, body.pos.y, body.pos.z);
      }
    }
  }

  //// //
  //// // get center and size. only needed for first step.
  //// // following steps can do this while appling movement
  //// //
  //// glm::dvec3 min_pos(std::numeric_limits<double>::max());
  //// glm::dvec3 max_pos(std::numeric_limits<double>::min());

  //// for (auto const &body : bodies) {
  ////   min_pos = min(min_pos, body.pos);
  ////   max_pos = max(max_pos, body.pos);
  //// }
  //// printf("min %f %f %f\n", min_pos.x, min_pos.y, min_pos.z);
  //// printf("max %f %f %f\n", max_pos.x, max_pos.y, max_pos.z);

  //// glm::dvec3 center = (max_pos + min_pos) / 2.0;
  //// printf("center %f %f %f\n", center.x, center.y, center.z);
  //// double size = compMax(max_pos - min_pos);
  //// printf("size %f\n", size);

  //// //
  //// // i5nit velocity
  //// //
  //// // for (int i = 0; i < bodies.size(); ++i) {
  //// //   glm::dvec3 upper(0.0);
  //// //   double lower = 0;
  //// //   for (int j = 0; j < bodies.size(); ++j) {
  //// //     upper += bodies[i].mass * bodies[i].vel;
  //// //     lower += bodies[i].mass;
  //// //   }

  //// //   bodies[i].vel -= (upper / lower);
  //// // }

  //// //
  //// // build octree
  //// //
  //// octree::Octree tree(center, size);
  //// for (auto const &body : bodies) {
  ////   tree.insert(body.pos, body.mass);
  //// }

  //// printf("\n\n\n");

  //// //
  //// // simulate
  //// //
  //// double theta = 0.0001;
  //// double time_step = 1.0;

  //// for (int it = 0; it < 1; ++it) {
  ////   tree.propagate();
  ////   for (auto &body : bodies) {
  ////     glm::dvec3 acc = tree.acc(body.pos, theta);
  ////     printf("acc %s %f\n", body.name.c_str(), glm::length(acc));

  ////     glm::dvec3 vel_halfstep = body.vel + acc * time_step / 2.0;
  ////     glm::dvec3 pos_halfstep = body.pos + vel_halfstep * time_step / 2.0;

  ////     body.pos = pos_halfstep + vel_halfstep * time_step / 2.0;
  ////     acc = tree.acc(body.pos, theta);
  ////     body.vel = vel_halfstep + acc * time_step / 2.0;
  ////   }

  ////   glm::dvec3 min_pos(std::numeric_limits<double>::max());
  ////   glm::dvec3 max_pos(std::numeric_limits<double>::min());
  ////   for (auto const &body : bodies) {
  ////     min_pos = min(min_pos, body.pos);
  ////     max_pos = max(max_pos, body.pos);
  ////   }
  ////   glm::dvec3 center = (max_pos + min_pos) / 2.0;
  ////   double size = compMax(max_pos - min_pos);

  ////   tree = octree::Octree(center, size);
  ////   for (auto &body : bodies) {
  ////     if (body.name == "Earth") {
  ////       printf("%f %f %f\n", body.pos.x, body.pos.y, body.pos.z);
  ////     }
  ////     tree.insert(body.pos, body.mass);
  ////   }
  //// }

  //
  // debug stuff
  //

  // glm::dvec3 pos(.50, -.07, 0.0);
  // glm::dvec3 acc = tree.acc(pos, 1.05);
  // printf("acc %f %f %f\n", acc.x, acc.y, acc.z);

  // int node_idx = 0;

  // std::vector<octree::Node> nodes;
  // nodes.push_back(tree.nodes[0]);
  // while (!nodes.empty()) {
  //   octree::Node node = nodes.back();
  //   nodes.pop_back();
  //   if (node.is_branch()) {
  //     for (int i = 0; i < 8; ++i) {
  //       nodes.push_back(tree.nodes[node.first_child + i]);
  //     }
  //   }
  //   if (node.is_leaf()) { // && !node.is_empty()) {
  //     if (node.mass != 0.0 && !node.cube.contains(node.mass_center)) {
  //       printf("illegal3\n");
  //     }
  //     printf("(%f, %f, %f) %f\n", node.cube.center.x, node.cube.center.y,
  //            node.cube.center.z, node.cube.size);
  //   }
  // }

  return 0;
}
