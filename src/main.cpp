// Your First C++ Program

#include "constants.hpp"
#include "space.hpp"
#include <cstdio>
#include <fstream>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <unordered_map>
#include <vector>

std::vector<space::CelestialBody> read_planets_and_moons(std::string path) {
  std::vector<space::CelestialBody> bodies;

  space::CelestialBody sun{
      0, "Sun", "STA", constants::sun_mass, glm::dvec3(0.0), glm::dvec3(0.0)};
  std::unordered_map<std::string, space::CelestialBody> bodies_map;
  bodies_map.insert({"Sun", sun});

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
  bodies_map.insert({"Sun", sun});

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      space::DataRow row = space::DataRow::parse_asteroid(line);
      space::CelestialBody body = row.to_body(bodies_map.size(), bodies_map);
      bodies_map.insert({body.name, body});
      bodies.push_back(body);
    }
    file.close();
  }
  return bodies;
}

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

int main() {
  std::vector<space::CelestialBody> bodies =
      read_planets_and_moons("planets_and_moons.csv");
  printf("bodies length %zu\n", bodies.size());

  // bodies = read_asteroids("scenario1_without_planets_and_moons.csv");
  // printf("bodies length %zu\n", bodies.size());

  double time_step = 0.25;
  std::vector<glm::dvec3> all_acc_old = get_gravitational_force(bodies);

  std::ofstream myfile;
  myfile.open("data.txt");

  for (int x = 0; x < int(365 / time_step); ++x) {
    printf("%f\n", x * time_step);

    if (x % int(1 / time_step) == 0) {
      for (size_t i = 0; i < bodies.size(); ++i) {
        myfile << "(" << bodies[i].pos.x << "," << bodies[i].pos.y << ")";
        if (i != bodies.size() - 1) {
          myfile << ",";
        } else {
          myfile << "\n";
        }

        if (bodies[i].name == "Earth") {
          printf("Earth %f %f\n", bodies[i].pos.x, bodies[i].pos.y);
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

  return 0;
}
