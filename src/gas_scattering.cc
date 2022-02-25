/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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

#include "constants.h"
#include "gas_scattering.h"
#include "matpack.h"
#include <cmath>

Vector calc_rayleighPhaMat(const Numeric& theta_rad,
                            const Index& stokes_dim) {

  using Constant::pow2;
  using Constant::pi;

  chk_if_in_range("Scattering angle", theta_rad, 0, pi);
  
  Vector pha_mat_int(6, 0.0);

  switch (stokes_dim) {
    case 4:
      pha_mat_int[4] = 0.0;                               //F34
      pha_mat_int[5] = cos(theta_rad);                    //F44 
      [[fallthrough]];
    case 3:
      pha_mat_int[3] = cos(theta_rad);                    //F33 
      [[fallthrough]];
    case 2:
      pha_mat_int[1] = -0.5 * pow2(sin(theta_rad));       //F12 
      pha_mat_int[2] = 0.5 * (1 + pow2(cos(theta_rad)));  //F22 
      [[fallthrough]];
    case 1:
      pha_mat_int[0] = 0.5 * (1 + pow2(cos(theta_rad)));  //F11
  }
  
  pha_mat_int *= 1.5;
  return pha_mat_int;
}
