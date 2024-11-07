#include "csv_writer.hpp"

void csv_writer::write_data(std::ofstream &myfile, glm::dvec3 *positions,
                            std::vector<space::CelestialBody> const &bodies) {
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
