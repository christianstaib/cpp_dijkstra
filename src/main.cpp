// Your First C++ Program

#include "space.hpp"
#include <cstdio>
#include <fstream>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <vector>

std::vector<space::DataRow> read_rows(std::string path) {
  std::vector<space::DataRow> bodies;

  std::ifstream file(path);

  if (file.is_open()) {
    std::string line;
    std::getline(file, line); // skip header
    //
    while (std::getline(file, line)) {
      space::DataRow body = space::DataRow::parse_planet_moon(line);
      bodies.push_back(body);
    }
    file.close();
  }
  return bodies;
}

int main() {
  std::vector<space::DataRow> rows = read_rows("planets_and_moons.csv");

  for (int i = 0; i < 5; ++i) {
    space::CelestialBody body = rows[i].to_body();
    printf("%s pos %f %f %f\n", body.name.c_str(), body.pos.x, body.pos.y,
           body.pos.z);
  }

  // for (space::DataRow const &row : rows) {
  //   printf("%s\n", row.name.c_str());
  // }

  return 0;
}
