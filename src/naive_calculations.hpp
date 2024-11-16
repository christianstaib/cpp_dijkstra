#pragma once

#include <cstddef>

#include "space.hpp"
#define GLM_ENABLE_EXPERIMENTAL

#include <glm/fwd.hpp>
#include <glm/gtx/norm.hpp>

namespace naive_calculations {
double get_kinetic_energy(size_t num_bodies, double *masses, glm::dvec3 *velocities);

double get_potential_energy(size_t num_bodies, double *masses, glm::dvec3 *positions);

void update_acceleration(space::MpiBodySystem &system);

void update_positions(space::MpiBodySystem &system, double step_size_days);

void update_velocity(space::MpiBodySystem &system, double step_size_days);

}  // namespace naive_calculations
