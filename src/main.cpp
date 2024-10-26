// Your First C++ Program

#include <omp.h>

#include <chrono>
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
#include <numeric>
#include <unordered_map>
#include <vector>

#include "constants.hpp"
#include "octree.hpp"
#include "space.hpp"

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
    std::getline(file, line);  // skip header

    while (std::getline(file, line)) {
      space::DataRow row = space::DataRow::parse_asteroid(line);
      space::CelestialBody body = row.to_body(name_to_body.size(), name_to_body);

      if (body.name == "Earth") {
        printf("Earth %f %f %f\n", body.pos.x, body.pos.y, body.pos.z);
      }

      if (!body.name.empty()) {
        if (name_to_body.find(body.name) == name_to_body.end()) {
          name_to_body.insert({body.name, body});
        } else {
          printf(
              "Error: A bdoy with the name %s is already known (distance "
              "%f km)\n",
              body.name.c_str(),
              glm::distance(body.pos, name_to_body.at(body.name).pos) * constants::meters_per_astronomical_unit /
                  1000.0);
        }
      }

      for (const auto &entry : pos_to_body) {
        if (glm::distance(entry.second.pos, body.pos) <= 10000 * constants::astronomical_units_per_meter) {
          printf(
              "Error: A body at the position %f %f %f is already known (it "
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

void update_forces_naive(double *masses, glm::dvec3 *positions, glm::dvec3 *velocities, glm::dvec3 *forces,
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
      squared_distance = glm::length2(distance_vector) + constants::squared_softening_factor;
      // x*sqrt(x) should be faster than pow(x, 3/2)
      forces[i] += (masses[j] * distance_vector) / (squared_distance * sqrt(squared_distance));
    }

    forces[i] *= constants::gravitational_constant_in_au3_per_kg_d2;
  }
}

double get_kinetic_energy(size_t num_bodies, double *masses, glm::dvec3 *velocities) {
  double kinetic_energy = 0.0;

  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    kinetic_energy += 0.5 * masses[body_idx] * glm::length2(velocities[body_idx]);
  }

  return kinetic_energy;
}

double get_potential_energy(size_t num_bodies, double *masses, glm::dvec3 *positions) {
  double potential_energy = 0.0;

  for (size_t i = 0; i < num_bodies; ++i) {
    for (size_t j = 0; j < i; ++j) {
      glm::dvec3 diff = positions[j] - positions[i];
      double distance = glm::length(diff) + constants::softening_factor;

      double new_val = (constants::gravitational_constant_in_au3_per_kg_d2 * masses[i] * masses[j]) / distance;

      potential_energy -= new_val;
    }
  }

  return potential_energy;
}

void write_data(std::ofstream &myfile, glm::dvec3 *positions, std::vector<space::CelestialBody> const &bodies) {
  for (size_t i = 0; i < bodies.size(); ++i) {
    myfile << "(" << bodies[i].name.c_str() << "," << bodies[i].type.c_str() << "," << positions[i].x << ","
           << positions[i].y << ")";
    if (i != bodies.size() - 1) {
      myfile << ",";
    } else {
      myfile << "\n";
    }
  }
}

void loging(std::vector<space::CelestialBody> &bodies, size_t &num_bodies, double *masses, glm::dvec3 *positions,
            glm::dvec3 *velocities, int &day_div, double &time_step, std::ofstream &myfile, int &iteration,
            std::chrono::steady_clock::time_point *begin) {
  if (iteration % day_div == 0) {
    double ke = get_kinetic_energy(num_bodies, masses, velocities);
    double pe = get_potential_energy(num_bodies, masses, positions);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double ms_per_it =
        ((double)std::chrono::duration_cast<std::chrono::microseconds>(end - *begin).count() / (double)day_div) /
        1000.0;
    *begin = std::chrono::steady_clock::now();
    printf("day %f %f %f %f ms/it\n", iteration * time_step, ke, pe, ms_per_it);
    write_data(myfile, positions, bodies);
  }
}

// x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
// Needs to be run in a parallel section
void update_positions(size_t num_bodies, glm::dvec3 *velocities, glm::dvec3 *forces, glm::dvec3 *positions,
                      double time_step, glm::dvec3 *min_edge, glm::dvec3 *max_edge) {
#pragma omp parallel
  {
#pragma omp critical
    {
      *min_edge = glm::dvec3(std::numeric_limits<double>::max());
      *max_edge = glm::dvec3(std::numeric_limits<double>::min());
    }

    glm::dvec3 local_min_edge = glm::dvec3(std::numeric_limits<double>::max());
    glm::dvec3 local_max_edge = glm::dvec3(std::numeric_limits<double>::min());

#pragma omp for simd schedule(static)
    for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
      positions[body_idx] += velocities[body_idx] * time_step + 0.5 * forces[body_idx] * time_step * time_step;

      local_min_edge = min(local_min_edge, positions[body_idx]);
      local_max_edge = max(local_max_edge, positions[body_idx]);
    }

    // better use custom reduction?
#pragma omp critical
    {
      *min_edge = min(*min_edge, local_min_edge);
      *max_edge = max(*max_edge, local_max_edge);
    }
  }
}

void rebuild_tree(size_t num_bodies, glm::dvec3 *positions, double *masses, glm::dvec3 *min_edge, glm::dvec3 *max_edge,
                  octree::Octree &test) {
  glm::dvec3 center = (*min_edge + *max_edge) * 1.05;
  glm::dvec3 diff = *max_edge - center;
  double size = std::max(std::max(diff.x, diff.y), diff.z);

  test.clear(center, size);
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    test.insert(positions[body_idx], masses[body_idx]);
  }
  test.propagate();
}

int main() {
  std::vector<space::CelestialBody> bodies = read_bodies("data/combined.csv");
  size_t num_bodies = bodies.size();
  printf("there are %zu bodies\n", num_bodies);

  // setup
  double *masses = new double[num_bodies];
  glm::dvec3 *positions = new glm::dvec3[num_bodies];
  glm::dvec3 *velocities = new glm::dvec3[num_bodies];
  glm::dvec3 *forces = new glm::dvec3[num_bodies];

  glm::dvec3 min_edge(std::numeric_limits<double>::max());
  glm::dvec3 max_edge(std::numeric_limits<double>::min());
  glm::dvec3 center(0.0);

  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    masses[body_idx] = bodies[body_idx].mass;
    positions[body_idx] = bodies[body_idx].pos;
    velocities[body_idx] = bodies[body_idx].vel;
    min_edge = min(min_edge, positions[body_idx]);
    max_edge = max(max_edge, positions[body_idx]);
  }

  double theta = 1.05;
  int day_div = 24;
  double simulation_step_size = 1.0 / day_div;
  double visualization_step_size = 1.0;
  int num_iterations = int((12 * 365) / simulation_step_size);

  std::ofstream myfile;
  myfile.open("data.txt");

  octree::Octree test(glm::dvec3(0.0), 0.0);
  rebuild_tree(num_bodies, positions, masses, &min_edge, &max_edge, test);
#pragma omp parallel for schedule(static)
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    forces[body_idx] = test.acc(positions[body_idx], theta);
  }

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int iteration = 0; iteration < num_iterations; ++iteration) {
    loging(bodies, num_bodies, masses, positions, velocities, day_div, simulation_step_size, myfile, iteration, &begin);

    // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
    update_positions(num_bodies, velocities, forces, positions, simulation_step_size, &min_edge, &max_edge);

    // v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
    rebuild_tree(num_bodies, positions, masses, &min_edge, &max_edge, test);
#pragma omp parallel
    {
      glm::dvec3 new_force(0.0);
#pragma omp for simd schedule(static)
      for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
        new_force = test.acc(positions[body_idx], theta);
        velocities[body_idx] += 0.5 * (forces[body_idx] + new_force) * simulation_step_size;
        forces[body_idx] = new_force;
      }
    }
  }

  myfile.close();

  free(masses);
  free(positions);
  free(velocities);
  free(forces);

  return 0;
}
