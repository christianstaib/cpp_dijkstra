#include "space.hpp"
#include "constants.hpp"
#include <cmath>
#include <cstdio>
#include <random>
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
DataRow DataRow::parse_asteroid(const std::string &line) {
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
  printf("value is %s\n", value.c_str());
  row.epoch = std::stod(value);

  std::getline(ss, value, ',');
  row.h = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.albedo = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.diameter = value.empty() ? 0.0 : std::stod(value);

  row.mass = 0.0;

  std::getline(ss, value, ',');
  row.class_name = value;

  std::getline(ss, value, ',');
  row.name = value;

  row.central_body = "Sun";

  if (row.mass == 0.0) {
    if (row.albedo == 0.0 &&
        constants::geometric_albedo.contains(row.class_name)) {
      auto [min, max] = constants::geometric_albedo.at(row.class_name);
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dist(min, max);
      row.albedo = dist(gen);
    }

    if (row.diameter == 0.0) {
      printf("%s hsa no diameter\n", row.name.c_str());
    }
  }

  return row;
}

// Function to parse a single row of CSV data and return a DataRow object
DataRow DataRow::parse_planet_moon(const std::string &line) {
  std::stringstream ss(line);
  DataRow row;
  std::string value;

  // Parse each value from the line and assign to the corresponding field
  std::getline(ss, value, ',');
  row.eccentricity = std::stod(value);

  std::getline(ss, value, ',');
  row.semi_major_axis = std::stod(value);

  std::getline(ss, value, ',');
  row.inclination = std::stod(value) * (M_PI / 180.0);

  std::getline(ss, value, ',');
  row.longitude_of_the_ascending_node = std::stod(value) * (M_PI / 180.0);

  std::getline(ss, value, ',');
  row.argument_of_periapsis = std::stod(value) * (M_PI / 180.0);

  std::getline(ss, value, ',');
  row.mean_anomaly = std::stod(value);

  std::getline(ss, value, ',');
  row.epoch = std::stod(value);

  std::getline(ss, value, ',');
  row.h = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.albedo = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.diameter = value.empty() ? 0.0 : std::stod(value);

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

CelestialBody DataRow::to_body() {}

} // namespace space
