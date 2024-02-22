#pragma once

#include <string>
#include <vector>

struct EnumeratedOption {
  std::string name;
  std::string desc;
  std::vector<std::vector<std::string>> values_and_desc;

  [[nodiscard]] std::string_view sz() const;
  [[nodiscard]] std::string head() const;
  [[nodiscard]] std::string tail() const;
  [[nodiscard]] std::string impl() const;
  [[nodiscard]] std::string docs() const;
};

const std::vector<EnumeratedOption>& internal_options();

// #include "enums.h"

// namespace Options {
// /** Keep time options available to switch over them */
// ENUMCLASS(
//     TimeStep, char, hour, hours, h, minute, minutes, min, second, seconds, s)

// /** Possible line shape coefficients */
// ENUMCLASS(LineShapeCoeff, char, X0, X1, X2, X3)

// /** Possible Jacobian for Wind and Magnetic field */
// ENUMCLASS(WindMagJacobian, char, u, v, w, strength)

// /** Possible Jacobian for basic line parameters */
// ENUMCLASS(BasicCatParamJacobian, char, LineStrength, LineCenter)

// /** Possible HITRAN types of catalogs */
// ENUMCLASS(
//     HitranType,
//     char,
//     Pre2004,   // 2004 version changed the .par-length
//     Post2004,  // New par length
//     Online  // Onine expects a modern .par line followed by Upper then Lower quantum numbers
// )

// /** Possible AddLines Speedups */
// ENUMCLASS(LblSpeedup, char, None, QuadraticIndependent, LinearIndependent)

// ENUMCLASS(SortingOption, char, ByFrequency, ByEinstein)

// /** Options for setting planets */
// ENUMCLASS(planetOption,
//           char,
//           Earth,
//           Io,
//           Jupiter,
//           Mars,
//           Venus)
// }  // namespace Options
