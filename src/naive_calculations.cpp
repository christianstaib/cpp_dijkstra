#include "naive_calculations.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <glm/fwd.hpp>
#include <glm/gtx/norm.hpp>

#include "constants.hpp"

double naive_calculations::get_kinetic_energy(size_t num_bodies, double *masses, glm::dvec3 *velocities) {
  double kinetic_energy = 0.0;

#pragma omp parallel for schedule(static) reduction(+ : kinetic_energy)
  for (size_t body_idx = 0; body_idx < num_bodies; ++body_idx) {
    kinetic_energy += 0.5 * masses[body_idx] * glm::length2(velocities[body_idx]);
  }

  return kinetic_energy;
}

double naive_calculations::get_potential_energy(size_t num_bodies, double *masses, glm::dvec3 *positions) {
  double potential_energy = 0.0;

#pragma omp parallel for schedule(guided) reduction(- : potential_energy)
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

void naive_calculations::update_acceleration(space::MpiBodySystem &system) {
#pragma omp parallel for schedule(guided)
  for (size_t dyn_to_update_idx = 0; dyn_to_update_idx < system.num_dynamic_bodies; ++dyn_to_update_idx) {
    system.acceleration_next_timestep[dyn_to_update_idx] = glm::dvec3(0.0);
    glm::dvec3 distance_vector;
    double squared_distance;

    for (size_t static_idx = 0; static_idx < system.num_static_bodies; ++static_idx) {
      size_t static_to_update_idx = dyn_to_update_idx + system.offset_dynamic_bodies;
      if (static_to_update_idx == static_idx) {
        continue;
      }

      // Precompute distance vector
      distance_vector = system.position[static_idx] - system.position[static_to_update_idx];
      squared_distance = glm::length2(distance_vector);
      squared_distance += constants::squared_softening_factor;
      // x*sqrt(x) should be faster than pow(x, 3/2)
      glm::dvec3 test = (system.mass[static_idx] * distance_vector) / (squared_distance * sqrt(squared_distance));
      system.acceleration_next_timestep[dyn_to_update_idx] += test;
    }

    system.acceleration_next_timestep[dyn_to_update_idx] *= constants::gravitational_constant_in_au3_per_kg_d2;
  }
}

// x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
void naive_calculations::update_positions(space::MpiBodySystem &system, double step_size_days) {
  system.min_edge = glm::dvec3(std::numeric_limits<double>::max());
  system.max_edge = glm::dvec3(std::numeric_limits<double>::min());

#pragma omp parallel
  {
    glm::dvec3 local_min_edge = glm::dvec3(std::numeric_limits<double>::max());
    glm::dvec3 local_max_edge = glm::dvec3(std::numeric_limits<double>::min());

#pragma omp for simd schedule(guided)
    for (size_t dyn_idx = 0; dyn_idx < system.num_dynamic_bodies; ++dyn_idx) {
      size_t static_idx = dyn_idx + system.offset_dynamic_bodies;
      system.position[static_idx] += system.velocity[dyn_idx] * step_size_days +
                                     0.5 * system.acceleration[dyn_idx] * step_size_days * step_size_days;

      local_min_edge = min(local_min_edge, system.position[static_idx]);
      local_max_edge = max(local_max_edge, system.position[static_idx]);
    }

    // TODO min reduction
#pragma omp critical
    {
      system.min_edge = min(system.min_edge, local_min_edge);
      system.max_edge = max(system.max_edge, local_max_edge);
    }
  }
}

// v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
void naive_calculations::update_velocity(space::MpiBodySystem &system, double step_size_days) {
#pragma omp parallel for schedule(guided)
  for (size_t dyn_idx = 0; dyn_idx < system.num_dynamic_bodies; ++dyn_idx) {
    system.velocity[dyn_idx] +=
        0.5 * (system.acceleration[dyn_idx] + system.acceleration_next_timestep[dyn_idx]) * step_size_days;
  }
}
