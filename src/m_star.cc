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
#include "messages.h"
#include "physics_funcs.h"


/*!
  \file   m_star.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>
  \date   2018-06-22

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

extern const Numeric PI;
extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void starBlackbodySimple(ArrayOfMatrix &star_spectrum,
                         ArrayOfVector &star_pos,
                         Index &star_do,
                         // Inputs:
                         const Vector &f_grid,
                         const Index &stokes_dim,
                         const Numeric &star_radius,
                         const Numeric &star_distance,
                         const Numeric &star_temperature,
                         const Numeric &star_latitude,
                         const Numeric &star_longitude,
                         const Verbosity &) {

  // spectrum
  const Numeric atan1 = std::atan(star_radius / star_distance);
  Matrix star_spec(f_grid.nelem(), stokes_dim,0. );

  planck(star_spec(joker,0), f_grid, star_temperature);
  star_spec *= PI * Constant::pow2(std::sin(atan1)); // calc flux of incoming rad at TOA
  star_spectrum.push_back(star_spec);

  // set flag
  star_do = 1;

  // set the position
  Vector pos{star_distance, star_latitude, star_longitude};
  star_pos.push_back(pos);
}

void starOff(Index &star_do,
             ArrayOfMatrix &star_spectrum,
             ArrayOfVector &star_pos,
             const Verbosity &){

  // set flag to False (default)
  star_do = 0;

  // create empty Array of Matrix for the star_spectrum
  star_spectrum.resize(0);

  // create empty Array of Vector for the star_pos
  star_pos.resize(0);
}

//  CREATE_OUT1;
//  out1 << "A star is rising!\n";
//  out1 << star_temperature << "\n";