#pragma once

#include <cstddef>
#include <glm/ext/vector_double3.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace space {
struct CelestialBody {
 public:
  /// Unique of body. Might be helpful to avoid calculate force upon itself.
  int id;
  std::string name;
  /// Type of body. Renamed from class as class is keyword.
  std::string type;
  /// Mass of body. Unit is Kg
  double mass;
  /// Position relative to sun. Unit is AU
  glm::dvec3 pos;
  /// Velocity relative to sun. Unit is AU/d
  glm::dvec3 vel;

  double kinetic_energy();

  std::string to_string();

  static CelestialBody from_state_vactor_string(const std::string &string);
  static CelestialBody sun();

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
  // rad
  double mean_anomaly;
  // JD
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
  std::string type;
  // TODO
  std::string name;
  // TODO
  std::string central_body;

  CelestialBody to_body(int id, std::unordered_map<std::string, space::CelestialBody> const &bodies);

  static DataRow parse_asteroid(const std::string &line);
  static DataRow parse_planet_moon(const std::string &line);
};

std::vector<space::CelestialBody> read_bodies(std::string path);

// Define a struct for body system data
struct MpiBodySystem {
  // For all bodies, the position and masses are known to all nodes as this is necessary for updating the positions
  // managed by the node.
  size_t num_static_bodies;
  double *mass;
  glm::dvec3 *position;
  glm::dvec3 min_edge;  // element-wise min over all positions
  glm::dvec3 max_edge;  // element-wise max over all positions

  //  Each node is responsible for a continous subset of bodies.
  size_t offset_dynamic_bodies;  // node_id * num_dynamic_bodies
  size_t num_dynamic_bodies;     // num_static_bodies / num_dynamic_bodies = nodes
  glm::dvec3 *velocity;
  glm::dvec3 *acceleration;
  glm::dvec3 *acceleration_next_timestep;

  // Constructor to initialize the system from a vector of bodies
  MpiBodySystem(const std::vector<CelestialBody> &bodies, size_t num_procs, size_t rank);

  // Destructor to clean up dynamically allocated memory
  ~MpiBodySystem();
};

}  // namespace space
