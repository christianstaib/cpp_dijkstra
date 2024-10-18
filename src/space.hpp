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

  static CelestialBody from_state_vactor_string(const std::string &string);

  // static CelestialBody from_orbital_elements(const std::string &string);
};

struct DataRow {
public:
  /// Unitless
  double eccentricity;
  /// AU
  double semi_major_axis;
  /// rad
  double inclination;
  /// rad
  double longitude_of_the_ascending_node;
  /// rad
  double argument_of_periapsis;
  // JD (julian day)
  double mean_anomaly;
  // TODO
  double epoch;

  // TODO
  double h;
  // TODO
  double albedo;
  // TODO
  double diameter;
  // TODO
  double mass;

  // TODO
  std::string class_name;
  // TODO
  std::string name;
  // TODO
  std::string central_body;

  static DataRow parse(const std::string &line);
};
} // namespace space
