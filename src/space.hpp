#pragma once

#include <glm/ext/vector_double3.hpp>
#include <string>

namespace space {
struct CelestialBody {
public:
  int id;
  std::string name;
  std::string class_type;
  double mass;
  glm::dvec3 pos;
  glm::dvec3 vel;

  CelestialBody(const std::string &string);
};

} // namespace space
