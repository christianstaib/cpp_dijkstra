#include "space.hpp"
#include <sstream>

namespace space {

// Constructor to initialize CelestialBody from a CSV string
CelestialBody CelestialBody::from_state_vactor_string(const std::string &line) {
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

// Function to parse a single row of CSV data and return a DataRow object
DataRow parseRow(const std::string &line) {
  std::stringstream ss(line);
  DataRow row;
  std::string value;

  // Parse each value from the line and assign to the corresponding field
  std::getline(ss, value, ',');
  row.eccentricity = std::stod(value);

  std::getline(ss, value, ',');
  row.semi_major_axis = std::stod(value);

  std::getline(ss, value, ',');
  row.inclination = std::stod(value);

  std::getline(ss, value, ',');
  row.longitude_of_the_ascending_node = std::stod(value);

  std::getline(ss, value, ',');
  row.argument_of_periapsis = std::stod(value);

  std::getline(ss, value, ',');
  row.mean_anomaly = std::stod(value);

  std::getline(ss, value, ',');
  row.epoch = std::stod(value);

  std::getline(ss, value, ',');
  row.h = std::stod(value);

  std::getline(ss, value, ',');
  row.albedo = std::stod(value);

  std::getline(ss, value, ',');
  row.diameter = std::stod(value);

  std::getline(ss, value, ',');
  row.mass = std::stod(value);

  std::getline(ss, value, ',');
  row.class_name = value;

  std::getline(ss, value, ',');
  row.name = value;

  std::getline(ss, value, ',');
  row.central_body = value;

  return row;
}

} // namespace space
