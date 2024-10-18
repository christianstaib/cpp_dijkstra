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
      space::CelestialBody body(line);
      bodies.push_back(body);
    }
    file.close();
  }

  double gravitational_constant = 6.67430e-11;
  double conversion_factor = pow(86400.0, 2.0) / pow(149597870700.0, 3.0);
  gravitational_constant *= conversion_factor;
  double softening_factor = 1e-11;
  double time_step = 1.0;

  std::vector<glm::dvec3> all_acc_old;
  // Compute accelerations based on the current positions
  for (int i = 0; i < bodies.size(); ++i) {
    glm::dvec3 acc(0);
    for (int j = 0; j < bodies.size(); ++j) {
      if (i != j) {
        acc +=
            bodies[j].mass *
            ((bodies[j].pos - bodies[i].pos) /
             pow(glm::length2(bodies[j].pos - bodies[i].pos) + softening_factor,
                 3.0 / 2.0));
      }
    }
    acc *= gravitational_constant;
    all_acc_old.push_back(acc);
  }

  for (int x = 0; x < 365; ++x) {
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
    std::vector<glm::dvec3> all_acc_new;
    for (int i = 0; i < bodies.size(); ++i) {
      glm::dvec3 acc(0);
      for (int j = 0; j < bodies.size(); ++j) {
        if (i != j) {
          acc += bodies[j].mass *
                 ((bodies[j].pos - bodies[i].pos) /
                  pow(glm::length2(bodies[j].pos - bodies[i].pos) +
                          softening_factor,
                      3.0 / 2.0));
        }
      }
      acc *= gravitational_constant;
      all_acc_new.push_back(acc);
    }

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
