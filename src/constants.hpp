#pragma once

#include <string>
#include <unordered_map>

namespace constants {
/// Number of seconds in a day.
inline constexpr double seconds_in_day = 86400.0;

/// Number of meters in an astronomical unit.
inline constexpr double meters_in_astronomical_unit = 149597870700.0;

/// Conversion factor from (m^3)/(Kg*s^2) to (AU^3)/(Kg*d^2)
inline constexpr double conversion_factor_m3_per_kg_s2_to_au3_per_kg_d2 =
    (seconds_in_day * seconds_in_day) /
    (meters_in_astronomical_unit * meters_in_astronomical_unit *
     meters_in_astronomical_unit);

/// Gravitational constant in (m^3)/(Kg*s^2).
inline constexpr double gravitational_constant_in_m3_per_kg_s2 = 6.67430e-11;

/// Gravitational constant in (AU^3)/(Kg*d^2).
inline constexpr double gravitational_constant_in_au3_per_kg_d2 =
    conversion_factor_m3_per_kg_s2_to_au3_per_kg_d2 *
    gravitational_constant_in_m3_per_kg_s2;

/// Softening factor used in the force calculation to avoid numerical
/// instabilities when two points are very close to each other.
inline constexpr double softening_factor = 1e-11;

/// Squared softening factor. See softening_factor for more information.
inline constexpr double squared_softening_factor =
    softening_factor * softening_factor;

/// TODO
inline constexpr double sun_refrence_epoch = 2451544.5;

const std::unordered_map<std::string, std::tuple<double, double>>
    geometric_albedo = {{"AMO", {0.450, 0.550}}, {"APO", {0.450, 0.550}},
                        {"ATE", {0.450, 0.550}}, {"IEO", {0.450, 0.550}},
                        {"MCA", {0.450, 0.550}}, {"IMB", {0.030, 0.103}},
                        {"MBA", {0.097, 0.203}}, {"OMB", {0.197, 0.500}},
                        {"CEN", {0.450, 0.750}}, {"TJN", {0.188, 0.124}},
                        {"TNO", {0.022, 0.130}}, {"AST", {0.450, 0.550}},
                        {"PAA", {0.450, 0.550}}, {"HYA", {0.450, 0.550}}};
} // namespace constants
