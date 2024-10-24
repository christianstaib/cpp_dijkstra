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
#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <omp.h>
#include <unordered_map>
#include <vector>

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

      if (body.name == "Earth") {
        printf("Earth %f %f %f\n", body.pos.x, body.pos.y, body.pos.z);
      }

      if (!body.name.empty()) {
        if (name_to_body.find(body.name) == name_to_body.end()) {
          name_to_body.insert({body.name, body});
        } else {
          printf("Error: A bdoy with the name %s is already known (distance "
                 "%f km)\n",
                 body.name.c_str(),
                 glm::distance(body.pos, name_to_body.at(body.name).pos) *
                     constants::meters_per_astronomical_unit / 1000.0);
        }
      }

      for (const auto &entry : pos_to_body) {
        if (glm::distance(entry.second.pos, body.pos) <=
            10000 * constants::astronomical_units_per_meter) {
          printf("Error: A body at the position %f %f %f is already known (it "
                 "is called %s)\n",
                 body.pos.x, body.pos.y, body.pos.z, entry.second.name.c_str());
        }
      }
      pos_to_body.push_back({body.pos, body});

      bodies.push_back(body);
    }
    file.close();
  }

  return bodies;
}

void update_gravitational_force(double *masses, glm::dvec3 *positions,
                                glm::dvec3 *velocities, glm::dvec3 *forces,
                                size_t num_bodies) {

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < num_bodies; ++i) {

    glm::dvec3 distance_vector;
    double squared_distance;
    for (size_t j = 0; j < num_bodies; ++j) {
      if (i == j) {
        continue;
      }

      // Precompute distance vector
      distance_vector = positions[j] - positions[i];
      squared_distance =
          glm::length2(distance_vector) + constants::squared_softening_factor;
      // x*sqrt(x) should be faster than pow(x, 3/2)
      forces[i] += (masses[j] * distance_vector) /
                   (squared_distance * sqrt(squared_distance));
    }

    forces[i] *= constants::gravitational_constant_in_au3_per_kg_d2;
  }
}

double get_kinetic_energy(double *masses, glm::dvec3 *velocities,
                          size_t num_bodies) {
  double kinetic_energy = 0.0;

  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    kinetic_energy +=
        0.5 * masses[body_idx] * glm::length2(velocities[body_idx]);
  }

  return kinetic_energy;
}

double get_potential_energy(double *masses, glm::dvec3 *positions,
                            size_t num_bodies) {
  double potential_energy = 0.0;

  for (size_t i = 0; i < num_bodies; ++i) {
    for (size_t j = 0; j < i; ++j) {
      glm::dvec3 diff = positions[j] - positions[i];
      double distance = glm::length(diff) + constants::softening_factor;

      double new_val = (constants::gravitational_constant_in_au3_per_kg_d2 *
                        masses[i] * masses[j]) /
                       distance;

      potential_energy -= new_val;
    }
  }

  return potential_energy;
}

void write_data(std::ofstream &myfile, glm::dvec3 *positions,
                std::vector<space::CelestialBody> const &bodies) {
  for (size_t i = 0; i < bodies.size(); ++i) {
    myfile << "(" << bodies[i].name.c_str() << "," << bodies[i].type.c_str()
           << "," << positions[i].x << "," << positions[i].y << ")";
    if (i != bodies.size() - 1) {
      myfile << ",";
    } else {
      myfile << "\n";
    }

    // std::string to_watch = "Sun";
    // if (bodies[i].name == to_watch) {
    //   //   printf("%s\n", bodies[i].to_string().c_str());
    //   double kinetic_energyx = kinetic_energy(bodies);
    //   double potential_energyx = potential_energy(bodies);
    //   double energy = kinetic_energyx + potential_energyx;
    //   double ratio = (2 * kinetic_energyx) / glm::abs(potential_energyx);
    //   printf("energy at day %f is %f, kinetic_energy is %f, "
    //          "potential_energy is %f, ratio is %f\n",
    //          x * time_step, energy, kinetic_energyx, potential_energyx,
    //          ratio);
    // }
  }
}

int main() {
  std::vector<space::CelestialBody> bodies =
      read_bodies("data/planets_and_moons.csv");
  size_t num_bodies = bodies.size();

  // setup
  double *masses = new double[num_bodies];
  glm::dvec3 *positions = new glm::dvec3[num_bodies];
  glm::dvec3 *velocities = new glm::dvec3[num_bodies];
  glm::dvec3 *forces = new glm::dvec3[num_bodies];
  glm::dvec3 *old_forces = new glm::dvec3[num_bodies];
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    auto body = bodies[body_idx];
    masses[body_idx] = body.mass;
    positions[body_idx] = body.pos;
    velocities[body_idx] = body.vel;
  }

  int day_div = 24;
  double time_step = 1.0 / day_div;
  int num_iterations = int((5 * 365) / time_step);
  std::ofstream myfile;
  myfile.open("data.txt");

  update_gravitational_force(masses, positions, velocities, old_forces,
                             num_bodies);

  for (int iteration = 0; iteration < num_iterations; ++iteration) {
    if (iteration % day_div == 0) {
      printf("day %f\n", iteration * time_step);
      double kinetic_energy =
          get_kinetic_energy(masses, velocities, num_bodies);
      double potential_energy =
          get_potential_energy(masses, positions, num_bodies);
      printf("kinetic_energy %f\n", kinetic_energy);
      printf("potential_energy %f\n", potential_energy);
      write_data(myfile, positions, bodies);
    }

    // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
    for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
      positions[body_idx] += velocities[body_idx] * time_step +
                             0.5 * old_forces[body_idx] * time_step * time_step;
    }

    // a_(i + 1)
    update_gravitational_force(masses, positions, velocities, forces,
                               num_bodies);

    // v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
    for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
      velocities[body_idx] +=
          0.5 * (old_forces[body_idx] + forces[body_idx]) * time_step;
    }

    old_forces = forces;
  }

  myfile.close();

  return 0;
}
