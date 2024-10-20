#include "space.hpp"
#include "constants.hpp"
#include <cmath>
#include <cstdio>
#include <glm/fwd.hpp>
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>

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
  std::getline(lineStream, body.type, ',');

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

std::string CelestialBody::to_string() {
  std::stringstream ss;
  ss << std::setprecision(std::numeric_limits<double>::max_digits10) << id
     << "," << name << "," << type << "," << mass << "," << pos.x << ","
     << pos.y << "," << pos.z << "," << vel.x << "," << vel.y << "," << vel.z;
  return ss.str();
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
  row.epoch = std::stod(value);

  std::getline(ss, value, ',');
  row.h = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.albedo = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.diameter = value.empty() ? 0.0 : std::stod(value);

  row.mass = 0.0;

  std::getline(ss, value, ',');
  row.type = value;

  std::getline(ss, value, ',');
  row.name = value;

  row.central_body = "Sun";

  if (row.albedo == 0.0) {
    printf("%s has no albedo\n", row.name.c_str());
    auto [min, max] = constants::geometric_albedo.at(row.type);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(min, max);
    row.albedo = dist(gen);
  }

  if (row.diameter == 0.0) {
    printf("%s has no diameter\n", row.name.c_str());
    row.diameter = 1329 * (1 / sqrt(row.albedo)) * pow(10.0, -0.2 * row.h);
  }

  if (row.mass == 0.0) {
    double p = 0.0;

    if (row.albedo < 0.1) {
      p = 1.38;
    } else if (row.albedo <= 0.2) {
      p = 2.71;
    } else {
      p = 5.32;
    }

    row.mass =
        (4.0 / 3.0) * M_PI * pow((row.diameter * 10e5) / 2, 3) * p / 10e5;
    printf("%s has no mass but now its %f\n", row.name.c_str(), row.mass);
    // printf("albedo %f\n", row.albedo);
    // printf("diameter %f\n", row.diameter);
    // printf("\n");
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
  row.mean_anomaly = std::stod(value) * (M_PI / 180.0);

  std::getline(ss, value, ',');
  row.epoch = std::stod(value);

  std::getline(ss, value, ',');
  row.h = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.albedo = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.diameter = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.mass = value.empty() ? 0.0 : std::stod(value);

  std::getline(ss, value, ',');
  row.type = value;

  std::getline(ss, value, ',');
  row.name = value;

  std::getline(ss, value, ',');
  row.central_body = value;

  return row;
}

CelestialBody DataRow::to_body(
    int id,
    std::unordered_map<std::string, space::CelestialBody> const &bodies) {
  CelestialBody body;

  space::CelestialBody central_body_body = bodies.at(central_body);

  // 1) Calculate mean anomaly m_t
  double m_t =
      mean_anomaly +
      (constants::sun_reference_epoch - epoch) *
          std::sqrt((constants::gravitational_constant_in_au3_per_kg_d2 *
                     central_body_body.mass) /
                    pow(semi_major_axis, 3));
  // printf("m_t %f\n", m_t);

  // 2) Solve Kepler’s equation for eccentric anomaly e_t
  double e_t = m_t;
  for (int i = 0; i < 30; ++i) {
    double f = e_t - eccentricity * sin(e_t) - m_t;
    double f_prime = 1 - eccentricity * cos(e_t);
    e_t -= f / f_prime;
  }
  // 3) Calculate the true anomaly
  double true_anomaly = 2.0 * atan2(sqrt(1 + eccentricity) * sin(e_t / 2),
                                    sqrt(1 - eccentricity) * cos(e_t / 2));
  // printf("true_anomaly %f\n", true_anomaly);

  // 4) Calculate the distance to the central to_body
  double distance_to_the_central_body =
      semi_major_axis * (1 - eccentricity * cos(e_t));
  // printf("distance_to_the_central_body %f\n", distance_to_the_central_body);

  // 5) Calculate the position ~o(t) and velocity ˙~o(t) vectors in the orbital
  // frame
  glm::dvec3 pos(cos(true_anomaly), sin(true_anomaly), 0);
  pos *= distance_to_the_central_body;
  glm::dvec3 vel(-sin(e_t), sqrt(1 - eccentricity * eccentricity) * cos(e_t),
                 0);
  vel *= sqrt(constants::gravitational_constant_in_au3_per_kg_d2 *
              central_body_body.mass * semi_major_axis) /
         distance_to_the_central_body;

  // this is correct

  // 6)
  double lo = argument_of_periapsis;           // lower omega
  double uo = longitude_of_the_ascending_node; // uper omega
  double i = inclination;
  double r11 = cos(lo) * cos(uo) - sin(lo) * cos(i) * sin(uo);
  double r12 = -(sin(lo) * cos(uo) + cos(lo) * cos(i) * sin(uo));
  double r21 = cos(lo) * sin(uo) + sin(lo) * cos(i) * cos(uo);
  double r22 = cos(lo) * cos(i) * cos(uo) - sin(lo) * sin(uo);
  double r31 = sin(lo) * sin(i);
  double r32 = cos(lo) * sin(i);
  glm::dmat3 r(glm::vec3(r11, r12, 0.0f), glm::vec3(r21, r22, 0.0f),
               glm::vec3(r31, r32, 0.0f));
  pos = pos * r; // WHY matrix on right side????? TODO
  vel = vel * r;

  // 7.)
  pos += central_body_body.pos;
  vel += central_body_body.vel;

  body.id = id;
  body.mass = mass;
  body.name = name;
  body.type = type;
  body.pos = pos;
  body.vel = vel;

  return body;
}

} // namespace space
