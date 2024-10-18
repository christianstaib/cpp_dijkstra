#pragma once

#include <glm/ext/vector_double3.hpp>
#include <string>

namespace space {
struct CelestialBody {
public:
  /// Unique of body. Might be helpful to avoid calculate force upon itself.
  int id;
  std::string name;
  /// TODO
  std::string class_type;
  /// Mass of body. Unit is Kg
  double mass;
  /// Position relative to sun. Unit is AU
  glm::dvec3 pos;
  /// Velocity relative to sun. Unit is AU/d
  glm::dvec3 vel;

  static CelestialBody from_state_vactors(const std::string &string);
};

} // namespace space
