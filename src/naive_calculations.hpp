#pragma once

#include <glm/fwd.hpp>
#include <glm/gtx/norm.hpp>

namespace naive_calculations {
double get_kinetic_energy(size_t num_bodies, double *masses, glm::dvec3 *velocities);

double get_potential_energy(size_t num_bodies, double *masses, glm::dvec3 *positions);

void update_forces(size_t num_bodies, double *masses, glm::dvec3 *positions, glm::dvec3 *velocities,
                   glm::dvec3 *forces);

}  // namespace naive_calculations
