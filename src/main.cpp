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
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <omp.h>
#include <unordered_map>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/hash.hpp"

std::vector<space::CelestialBody> read_planets_and_moons(std::string path) {
  std::vector<space::CelestialBody> bodies;

  space::CelestialBody sun{
      0, "Sun", "STA", constants::sun_mass, glm::dvec3(0.0), glm::dvec3(0.0)};
  std::unordered_map<std::string, space::CelestialBody> bodies_map;
  bodies_map.insert({"Sun", sun});
  bodies.push_back(sun);

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      space::DataRow row = space::DataRow::parse_planet_moon(line);
      space::CelestialBody body = row.to_body(bodies_map.size(), bodies_map);
      bodies_map.insert({body.name, body});
      bodies.push_back(body);
    }
    file.close();
  }
  return bodies;
}

std::vector<space::CelestialBody> read_asteroids(std::string path) {
  std::vector<space::CelestialBody> bodies;

  space::CelestialBody sun{
      0, "Sun", "STA", constants::sun_mass, glm::dvec3(0.0), glm::dvec3(0.0)};
  std::unordered_map<std::string, space::CelestialBody> bodies_map;
  std::unordered_map<glm::dvec3, space::CelestialBody> bodies_map2;
  bodies_map.insert({"Sun", sun});
  bodies_map2.insert({glm::dvec3(0.0), sun});
  bodies.push_back(sun);

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      space::DataRow row = space::DataRow::parse_asteroid(line);
      space::CelestialBody body = row.to_body(bodies_map.size(), bodies_map);
      bodies_map.insert({body.name, body});
      bodies_map2.insert({body.pos, body});
      bodies.push_back(body);
    }
    file.close();
  }

  printf("%zu %zu %zu\n", bodies_map.size(), bodies_map2.size(), bodies.size());

  return bodies;
}

std::vector<glm::dvec3>
get_gravitational_force(const std::vector<space::CelestialBody> &bodies) {
  std::vector<glm::dvec3> force(bodies.size(), glm::dvec3(0));

#pragma omp parallel
  {

#pragma omp for
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
      double new_val = (constants::gravitational_constant_in_au3_per_kg_d2 *
                        bodies[i].mass * bodies[j].mass) /
                       (glm::length(bodies[j].pos - bodies[i].pos));
      potential_energy -= new_val;
      // printf("i:%zu/%zu j:%zu, potential_energy is %f, new val is %f \n", i,
      //        bodies.size(), j, potential_energy, new_val);
    }
  }

  return potential_energy;
}

int main() {
  std::vector<space::CelestialBody> bodies =
      read_planets_and_moons("planets_and_moons.csv");

  std::vector<space::CelestialBody> asteroids =
      read_asteroids("scenario1_without_planets_and_moons.csv");
  bodies.insert(bodies.end(), asteroids.begin() + 1, asteroids.end());

  for (auto const &body : bodies) {
    if (body.name == "Pluto") {
      printf("found Pluto\n");
    }
  }

  // bodies = read_asteroids("scenario1_without_planets_and_moons.csv");
  // printf("bodies length %zu\n", bodies.size());

  double time_step = 1.0 / 24.0;
  std::vector<glm::dvec3> all_acc_old = get_gravitational_force(bodies);

  std::ofstream myfile;
  myfile.open("data.txt");

  for (int x = 0; x < int((10 * 365) / time_step); ++x) {
    // printf("%f\n", x * time_step);

    if (x % 24 == 0) {
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
