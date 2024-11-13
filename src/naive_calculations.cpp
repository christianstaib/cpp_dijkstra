#include "naive_calculations.hpp"

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

void naive_calculations::update_acceleration(space::BodySystem &local_body_system,
                                             space::BodySystem const &global_body_system) {
#pragma omp parallel for schedule(static)
  for (size_t local_idx = 0; local_idx < local_body_system.num_bodies; ++local_idx) {
    glm::dvec3 distance_vector;
    double squared_distance;
    for (size_t global_idx = 0; global_idx < global_body_system.num_bodies; ++global_idx) {
      // Precompute distance vector
      distance_vector = local_body_system.position[global_idx] - local_body_system.position[local_idx];
      squared_distance = glm::length2(distance_vector);
      if (squared_distance != 0.0) {
        squared_distance += constants::squared_softening_factor;
        // x*sqrt(x) should be faster than pow(x, 3/2)
        local_body_system.acceleration_next_timestep[local_idx] +=
            (local_body_system.mass[global_idx] * distance_vector) / (squared_distance * sqrt(squared_distance));
      }
    }

    local_body_system.acceleration_next_timestep[local_idx] *= constants::gravitational_constant_in_au3_per_kg_d2;
  }
}

/// Naive approach to get gravitional force on all bodies in (Kg*AU)/d^2.
glm::dvec3 get_gravitational_force_ld(size_t body_idx_want_force, size_t num_bodies, glm::dvec3 *positions,
                                      double *masses) {
  long double x = 0;
  long double y = 0;
  long double z = 0;

  glm::dvec3 distance_vector;
  double squared_distance;

  for (size_t j = 0; j < num_bodies; ++j) {
    // Check body_idx_want_force == j can be skiped as distance_vector will be
    // zero in this case

    // Precompute distance vector
    distance_vector = positions[j] - positions[body_idx_want_force];
    squared_distance = glm::length2(distance_vector) + constants::squared_softening_factor;
    // x*sqrt(x) should be faster than pow(x, 3/2)
    x += (masses[j] * distance_vector.x) / (squared_distance * sqrt(squared_distance));
    y += (masses[j] * distance_vector.y) / (squared_distance * sqrt(squared_distance));
    z += (masses[j] * distance_vector.z) / (squared_distance * sqrt(squared_distance));
  }

  x *= constants::gravitational_constant_in_au3_per_kg_d2;
  y *= constants::gravitational_constant_in_au3_per_kg_d2;
  z *= constants::gravitational_constant_in_au3_per_kg_d2;

  glm::dvec3 force(x, y, z);

  // multipling once at the end is faster and also better for precision
  return force;
}

// x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
// Needs to be run in a parallel section
void naive_calculations::update_positions(space::BodySystem &local_body_system, double time_step) {
  // x_{i + 1} = x_i + v_i * dt + 0.5 * a_i dt^2
#pragma omp parallel
  {
#pragma omp critical
    {
      local_body_system.min_edge = glm::dvec3(std::numeric_limits<double>::max());
      local_body_system.max_edge = glm::dvec3(std::numeric_limits<double>::min());
    }

    glm::dvec3 local_min_edge = glm::dvec3(std::numeric_limits<double>::max());
    glm::dvec3 local_max_edge = glm::dvec3(std::numeric_limits<double>::min());

#pragma omp for simd schedule(guided)
    for (size_t body_idx = 0; body_idx < local_body_system.num_bodies; ++body_idx) {
      local_body_system.position[body_idx] += local_body_system.velocity[body_idx] * time_step +
                                              0.5 * local_body_system.acceleration[body_idx] * time_step * time_step;

      local_min_edge = min(local_min_edge, local_body_system.position[body_idx]);
      local_max_edge = max(local_max_edge, local_body_system.position[body_idx]);
    }

    // better use custom reduction?
#pragma omp critical
    {
      local_body_system.min_edge = min(local_body_system.min_edge, local_min_edge);
      local_body_system.max_edge = max(local_body_system.max_edge, local_max_edge);
    }
  }
}

// v_{i + 1} = v_i + 0.5 (a_i + a_{i + 1}) * dt
void naive_calculations::update_velocity(space::BodySystem &body_system, double step_size_days) {
#pragma omp parallel for simd
  for (size_t body_idx = 0; body_idx < body_system.num_bodies; ++body_idx) {
    body_system.velocity[body_idx] +=
        0.5 * (body_system.acceleration[body_idx] + body_system.acceleration_next_timestep[body_idx]) * step_size_days;
  }
}
