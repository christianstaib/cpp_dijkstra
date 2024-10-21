// Your First C++ Program

#include "constants.hpp"
#include "space.hpp"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <glm/detail/qualifier.hpp>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <omp.h>
#include <unordered_map>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/hash.hpp"

std::vector<space::CelestialBody> read_bodies(std::string path) {
  // Set up data structures needed for reading the asteroids.
  space::CelestialBody sun = space::CelestialBody::sun();

  std::unordered_map<std::string, space::CelestialBody> name_to_body;
  name_to_body.insert({"Sun", sun});

  std::vector<std::pair<glm::dvec3, space::CelestialBody>> pos_to_body;
  pos_to_body.push_back({sun.pos, sun});

  std::vector<space::CelestialBody> bodies;
  bodies.push_back(sun);

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
      space::DataRow row = space::DataRow::parse_asteroid(line);
      space::CelestialBody body =
          row.to_body(name_to_body.size(), name_to_body);

      // printf("%s\n", line.c_str());
      // printf("%f %f %f \n", body.pos.x, body.pos.y, body.pos.z);

      if (!body.name.empty()) {
        if (name_to_body.find(body.name) == name_to_body.end()) {
          name_to_body.insert({body.name, body});
        } else {
          printf(
              "Error: A bdoy with the name %s is already known (distance %f)\n",
              body.name.c_str(),
              glm::distance(body.pos, name_to_body.at(body.name).pos));
        }
      }

      for (const auto &entry : pos_to_body) {
        if (glm::distance(entry.second.pos, body.pos) <= 6.68459e-7) {
          printf("Error: A body at the position %f %f %f is already known (it "
                 "is called %s)\n",
                 body.pos.x, body.pos.y, body.pos.z, entry.second.name.c_str());
        }
      }
      pos_to_body.push_back({body.pos, body});

      // if (pos_to_body.find(body.pos) == pos_to_body.end()) {
      //   pos_to_body.insert({body.pos, body});
      // } else {
      //   printf("Error: A body at the position %f %f %f is already known\n",
      //          body.pos.x, body.pos.y, body.pos.z);
      // }

      bodies.push_back(body);
    }
    file.close();
  }

  return bodies;
}

std::vector<glm::dvec3>
get_gravitational_force(const std::vector<space::CelestialBody> &bodies) {
  std::vector<glm::dvec3> force(bodies.size(), glm::dvec3(0));

#pragma omp parallel for
  for (size_t i = 0; i < bodies.size(); ++i) {

    glm::dvec3 distance_vector;
    double squared_distance;
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

double kinetic_energy(std::vector<space::CelestialBody> const &bodies) {
  double kinetic_energy = 0.0;

  for (auto const &body : bodies) {
    kinetic_energy += 0.5 * body.mass * glm::length2(body.vel);
  }

  return kinetic_energy;
}

double potential_energy(std::vector<space::CelestialBody> const &bodies) {
  double potential_energy = 0.0;

  for (size_t i = 0; i < bodies.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      glm::dvec3 diff = bodies[j].pos - bodies[i].pos;
      double distance = glm::length(diff) + constants::softening_factor;

      double new_val = (constants::gravitational_constant_in_au3_per_kg_d2 *
                        bodies[i].mass * bodies[j].mass) /
                       distance;

      potential_energy -= new_val;
    }
  }

  return potential_energy;
}

int main() {
  std::vector<space::CelestialBody> bodies =
      read_bodies("planets_and_moons.csv");

  for (auto const &body : bodies) {
    if (body.name == "Pluto") {
      printf("found Pluto\n");
    }
  }

  // bodies = read_asteroids("scenario1_without_planets_and_moons.csv");
  // printf("bodies length %zu\n", bodies.size());

  int day_div = 100;
  double time_step = 1.0 / day_div;
  std::vector<glm::dvec3> all_acc_old = get_gravitational_force(bodies);

  std::ofstream myfile;
  myfile.open("data.txt");

  for (int x = 0; x < int((10 * 365) / time_step); ++x) {
    // printf("%f\n", x * time_step);

    if (x % day_div == 0) {
      for (size_t i = 0; i < bodies.size(); ++i) {
        myfile << "(" << bodies[i].name.c_str() << "," << bodies[i].pos.x << ","
               << bodies[i].pos.y << ")";
        if (i != bodies.size() - 1) {
          myfile << ",";
        } else {
          myfile << "\n";
        }

        std::string to_watch = "Jupiter";
        if (bodies[i].name == to_watch) {
          //   printf("%s\n", bodies[i].to_string().c_str());
          double kinetic_energyx = kinetic_energy(bodies);
          double potential_energyx = potential_energy(bodies);
          double energy = kinetic_energyx + potential_energyx;
          double ratio = (2 * kinetic_energyx) / glm::abs(potential_energyx);
          printf("energy at day %f is %f, kinetic_energy is %f, "
                 "potential_energy is %f, ratio is %f\n",
                 x * time_step, energy, kinetic_energyx, potential_energyx,
                 ratio);
          // printf("%s at day %f: %f %f\n", to_watch.c_str(), x * time_step,
          //        bodies[i].pos.x, bodies[i].pos.y);
        }
      }
    }

    for (int i = 0; i < bodies.size(); ++i) {
      auto &body = bodies[i];
      auto acc_old = all_acc_old[i];

      // Half-step position update
      body.pos += body.vel * time_step + 0.5 * acc_old * time_step * time_step;
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

  myfile.close();

  for (auto const &body : bodies) {
  }

  return 0;
}
