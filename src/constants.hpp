#pragma once

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
} // namespace constants
