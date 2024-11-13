#pragma once

#include "space.hpp"
#define GLM_ENABLE_EXPERIMENTAL

#include <glm/fwd.hpp>
#include <glm/gtx/norm.hpp>

namespace naive_calculations {
double get_kinetic_energy(size_t num_bodies, double *masses, glm::dvec3 *velocities);

double get_potential_energy(size_t num_bodies, double *masses, glm::dvec3 *positions);

void update_acceleration(space::BodySystem &body_system, space::BodySystem const &global_body_system);

glm::dvec3 get_gravitational_force_ld(size_t body_idx_want_force, size_t num_bodies, glm::dvec3 *positions,
                                      double *masses);

void update_positions(space::BodySystem &body_system, double time_step);

void update_velocity(space::BodySystem &body_system, double step_size_days);

}  // namespace naive_calculations
