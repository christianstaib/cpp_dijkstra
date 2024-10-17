#include "space.hpp"
#include <sstream>

namespace space {

// Constructor to initialize CelestialBody from a CSV string
CelestialBody::CelestialBody(const std::string &line) {
  std::stringstream lineStream(line);
  std::string token;

  // Parse the ID
  std::getline(lineStream, token, ',');
  id = std::stoi(token);

  // Parse the name
  std::getline(lineStream, name, ',');

  // Parse the class type
  std::getline(lineStream, class_type, ',');

  // Parse the mass
  std::getline(lineStream, token, ',');
  mass = std::stod(token);

  // Parse position (pos_x, pos_y, pos_z)
  std::getline(lineStream, token, ',');
  pos.x = std::stod(token);

  std::getline(lineStream, token, ',');
  pos.y = std::stod(token);

  std::getline(lineStream, token, ',');
  pos.z = std::stod(token);

  // Parse velocity (vel_x, vel_y, vel_z)
  std::getline(lineStream, token, ',');
  vel.x = std::stod(token);

  std::getline(lineStream, token, ',');
  vel.y = std::stod(token);

  std::getline(lineStream, token, ',');
  vel.z = std::stod(token);
}

} // namespace space
