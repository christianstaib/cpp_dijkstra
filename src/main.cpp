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
        // 3,Earth,PLA,5.972185999999999813e+24,-0.1685246483858766631,0.9687833049070307956,-4.120490278477268624e-06,-0.01723415470267592939,-0.003007696701791743449,3.562616065704365793e-08
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

void get_gravitational_force(std::vector<double> const &masses,
                             std::vector<glm::dvec3> const &positions,
                             std::vector<glm::dvec3> const &velocities,
                             std::vector<glm::dvec3> &forces) {

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < positions.size(); ++i) {

    glm::dvec3 distance_vector;
    double squared_distance;
    for (size_t j = 0; j < positions.size(); ++j) {
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

void write_data(std::ofstream &myfile, std::vector<glm::dvec3> const &positions,
                std::vector<space::CelestialBody> const &bodies) {
  for (size_t i = 0; i < positions.size(); ++i) {
    myfile << "(" << bodies[i].name.c_str() << "," << bodies[i].type.c_str()
           << "," << positions[i].x << "," << positions[i].y << ")";
    if (i != positions.size() - 1) {
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
  std::vector<space::CelestialBody> bodies = read_bodies("combined.csv");

  // setup
  std::vector<double> masses;
  std::vector<glm::dvec3> positions;
  std::vector<glm::dvec3> velocities;
  std::vector<glm::dvec3> forces;
  std::vector<glm::dvec3> old_forces;
  for (auto const &body : bodies) {
    masses.push_back(body.mass);
    positions.push_back(body.pos);
    velocities.push_back(body.vel);
    forces.push_back(glm::dvec3(0));
    old_forces.push_back(glm::dvec3(0));
  }

  int day_div = 24;
  double time_step = 1.0 / day_div;
  std::ofstream myfile;
  myfile.open("data.txt");

  get_gravitational_force(masses, positions, velocities, old_forces);
  for (int x = 0; x < int((5 * 365) / time_step); ++x) {
    if (x % day_div == 0) {
      printf("day %f\n", x * time_step);
      write_data(myfile, positions, bodies);
    }

    glm::dvec3 *p = positions.data();

    // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
    for (int i = 0; i < bodies.size(); ++i) {
      positions[i] += velocities[i] * time_step +
                      0.5 * old_forces[i] * time_step * time_step;
    }

    // a_(i + 1)
    get_gravitational_force(masses, positions, velocities, forces);

    // v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
    for (int i = 0; i < bodies.size(); ++i) {
      velocities[i] += 0.5 * (old_forces[i] + forces[i]) * time_step;
    }

    old_forces = forces;
  }

  myfile.close();

  return 0;
}
