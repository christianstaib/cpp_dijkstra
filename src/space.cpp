#include "space.hpp"
#include <sstream>

namespace space {

// Constructor to initialize CelestialBody from a CSV string
CelestialBody CelestialBody::from_state_vactors(const std::string &line) {
  std::stringstream lineStream(line);
  std::string token;

  // Create a new CelestialBody instance
  CelestialBody body;

  // Parse the ID
  std::getline(lineStream, token, ',');
  body.id = std::stoi(token);

  // Parse the name
  std::getline(lineStream, body.name, ',');

  // Parse the class type
  std::getline(lineStream, body.class_type, ',');

  // Parse the mass
  std::getline(lineStream, token, ',');
  body.mass = std::stod(token);

  // Parse position (pos_x, pos_y, pos_z)
  std::getline(lineStream, token, ',');
  body.pos.x = std::stod(token);

  std::getline(lineStream, token, ',');
  body.pos.y = std::stod(token);

  std::getline(lineStream, token, ',');
  body.pos.z = std::stod(token);

  // Parse velocity (vel_x, vel_y, vel_z)
  std::getline(lineStream, token, ',');
  body.vel.x = std::stod(token);

  std::getline(lineStream, token, ',');
  body.vel.y = std::stod(token);

  std::getline(lineStream, token, ',');
  body.vel.z = std::stod(token);

  // Return the populated CelestialBody instance
  return body;
}

} // namespace space
