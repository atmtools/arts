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

#include "arts_conversions.h"
#include "gas_scattering.h"
#include "matpack_data.h"
#include "check_input.h"
#include <cmath>

Vector calc_rayleighPhaMat(const Numeric& theta_rad) {

  using Math::pow2;
  using Constant::pi;

  chk_if_in_range("Scattering angle", theta_rad, 0, pi);
  
  Vector pha_mat_int(6, 0.0);

  pha_mat_int[0] = 0.5 * (1 + pow2(cos(theta_rad))); // F11
  pha_mat_int[1] = -0.5 * pow2(sin(theta_rad));      // F12
  pha_mat_int[2] = 0.5 * (1 + pow2(cos(theta_rad))); // F22
  pha_mat_int[3] = cos(theta_rad);                      // F33
  pha_mat_int[4] = 0.0;                                 // F34
  pha_mat_int[5] = cos(theta_rad);                      // F44

  pha_mat_int *= 1.5;
  return pha_mat_int;
}
