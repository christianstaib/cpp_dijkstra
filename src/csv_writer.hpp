#pragma once

#include <fstream>
#include <glm/fwd.hpp>
#include <vector>

#include "space.hpp"

namespace csv_writer {
void write_data(std::ofstream &myfile, glm::dvec3 *positions, std::vector<space::CelestialBody> const &bodies);
}
