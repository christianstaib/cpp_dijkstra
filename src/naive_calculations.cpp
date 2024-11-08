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

void naive_calculations::update_acceleration(space::BodySystem &body_system) {
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < body_system.num_bodies; ++i) {
    glm::dvec3 distance_vector;
    double squared_distance;
    for (size_t j = 0; j < body_system.num_bodies; ++j) {
      if (i == j) {
        continue;
      }

      // Precompute distance vector
      distance_vector = body_system.position[j] - body_system.position[i];
      squared_distance = glm::length2(distance_vector) + constants::squared_softening_factor;
      // x*sqrt(x) should be faster than pow(x, 3/2)
      body_system.acceleration[i] +=
          (body_system.mass[j] * distance_vector) / (squared_distance * sqrt(squared_distance));
    }

    body_system.acceleration[i] *= constants::gravitational_constant_in_au3_per_kg_d2;
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
