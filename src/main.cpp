// Your First C++ Program

#include "constants.hpp"
#include "octree.hpp"
#include "space.hpp"
#include <cstdio>
#include <fstream>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <vector>

/// Naive approach to get gravitional force on all bodies in (Kg*AU)/d^2.
std::vector<glm::dvec3>
get_gravitational_force(const std::vector<space::CelestialBody> &bodies) {
  std::vector<glm::dvec3> force(bodies.size(), glm::dvec3(0));

  glm::dvec3 distance_vector;
  double squared_distance;

  for (size_t i = 0; i < bodies.size(); ++i) {
    for (size_t j = 0; j < bodies.size(); ++j) {
      if (i == j) {
        continue;
      }

      // Precompute distance vector
      distance_vector = bodies[j].pos - bodies[i].pos;
      squared_distance =
          glm::length2(distance_vector) + constants::squared_softening_factor;
      // x*sqrt(x) should be faster than pow(x, 3/2)
      force[i] += (bodies[j].mass * distance_vector) /
                  (squared_distance * sqrt(squared_distance));
    }

    // multipling once at the end is faster and also better for precision
    force[i] *= constants::gravitational_constant_in_au3_per_kg_d2;
  }

  return force;
}

std::vector<space::CelestialBody> read_bodies(std::string path) {
  std::vector<space::CelestialBody> bodies;

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      space::CelestialBody body =
          space::CelestialBody::from_state_vactors(line);
      bodies.push_back(body);
    }
    file.close();
  }
  return bodies;
}

int main() {
  std::vector<space::CelestialBody> bodies =
      read_bodies("planets_and_moons_state_vectors.csv");
  //
  // read bodies
  //

  double time_step = 0.25;
  std::vector<glm::dvec3> all_acc_old = get_gravitational_force(bodies);

  for (int x = 0; x < int(365 / time_step); ++x) {
    for (int i = 0; i < bodies.size(); ++i) {
      auto &body = bodies[i];
      auto acc_old = all_acc_old[i];

      // Half-step position update
      body.pos += body.vel * time_step + 0.5 * acc_old * time_step * time_step;

      if (body.name == "Earth") {
        printf("Earth %f %f %f\n", body.pos.x, body.pos.y, body.pos.z);
      }
    }

    // Compute new accelerations based on the current positions
    std::vector<glm::dvec3> all_acc_new = get_gravitational_force(bodies);

    for (int i = 0; i < bodies.size(); ++i) {
      auto &body = bodies[i];
      auto acc_old = all_acc_old[i];
      auto acc_new = all_acc_new[i]; // Corrected to use the new acceleration

      // Half-step velocity update
      body.vel += 0.5 * (acc_old + acc_new) *
                  time_step; // Correctly using both old and new accelerations
    }

    // Update the old accelerations for the next iteration
    all_acc_old = all_acc_new;
  }

  return 0;
}
