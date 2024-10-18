// Your First C++ Program

#include "constants.hpp"
#include "space.hpp"
#include <cstdio>
#include <fstream>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_double3.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/quaternion.hpp>
#include <unordered_map>
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
  space::CelestialBody sun{
      0, "Sun", "STA", constants::sun_mass, glm::dvec3(0.0), glm::dvec3(0.0)};

  std::unordered_map<std::string, space::CelestialBody> bodies_map;
  std::vector<space::CelestialBody> bodies;

  bodies_map.insert({"Sun", sun});
  bodies.push_back(sun);

  std::vector<space::DataRow> rows = read_rows("planets_and_moons.csv");

  for (int i = 0; i < rows.size(); ++i) {
    space::CelestialBody body = rows[i].to_body(bodies_map.size(), bodies_map);
    bodies_map.insert({body.name, body});
    bodies.push_back(body);
  }

  printf("row length %zu bodies length %zu\n", rows.size(), bodies_map.size());

  for (auto it = bodies.begin(); it != bodies.end(); ++it) {
    printf("%s\n", it->to_string().c_str());
  }

  // for (space::DataRow const &row : rows) {
  //   printf("%s\n", row.name.c_str());
  // }

  return 0;
}
