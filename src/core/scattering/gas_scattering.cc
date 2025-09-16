/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   gas_scattering.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-03-15

  \brief  Implementation file for functions related to gas scattering.
*/

#include "gas_scattering.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>

#include <cmath>

Vector calc_rayleighPhaMat(const Numeric& theta_rad) {
  using Constant::pi;
  using Math::pow2;

  ARTS_USER_ERROR_IF(
      theta_rad != std::clamp<Numeric>(theta_rad, 0.0, pi),
      "Error in calc_rayleighPhaMat: Scattering angle must be in the range [0, pi], is {}",
      theta_rad);

  Vector pha_mat_int(6, 0.0);

  pha_mat_int[0] = 0.5 * (1 + pow2(cos(theta_rad)));  // F11
  pha_mat_int[1] = -0.5 * pow2(sin(theta_rad));       // F12
  pha_mat_int[2] = 0.5 * (1 + pow2(cos(theta_rad)));  // F22
  pha_mat_int[3] = cos(theta_rad);                    // F33
  pha_mat_int[4] = 0.0;                               // F34
  pha_mat_int[5] = cos(theta_rad);                    // F44

  pha_mat_int *= 1.5;
  return pha_mat_int;
}
