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
#include "arts.h"
#include "auto_md.h"
#include "star.h"


/*!
  \file   m_star.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-02-08

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
                         ArrayOfVector &star_radius,
                         Index &star_do,
                         // Inputs:
                         const Vector &f_grid,
                         const Index &stokes_dim,
                         const Numeric &radius,
                         const Numeric &distance,
                         const Numeric &temperature,
                         const Numeric &latitude,
                         const Numeric &longitude,
                         const Verbosity &) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the star (",distance," m) \n"
                     " is smaller than the radius of the star (", radius," m )")


  // spectrum
  const Numeric atan1 = std::atan(radius / distance);
  Matrix star_spec(f_grid.nelem(), stokes_dim,0. );

  planck(star_spec(joker,0), f_grid, temperature);
  star_spec *= PI * Constant::pow2(std::sin(atan1)); // calc flux of incoming rad at TOA
  star_spectrum.push_back(star_spec);

  // set flag
  star_do = 1;

  // set the position
  Vector pos{distance, latitude, longitude};
  star_pos.push_back(pos);

  // set the radius
  star_radius.push_back(Vector(1,radius));

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



/* Workspace method: Doxygen documentation will be auto-generated */
void CosmicMicrowaveAndStarBackground(Matrix &iy,
                                      const Vector &f_grid,
                                      const Vector &rtp_pos,
                                      const Vector &rtp_los,
                                      const ArrayOfMatrix &star_spectrum,
                                      const ArrayOfVector &star_pos,
                                      const ArrayOfVector &star_radius,
                                      const Vector &refellipsoid,
                                      const Index &star_do,
                                      const Index &stokes_dim,
                                      const Index &atmosphere_dim,
                                      const Verbosity &verbosity) {

  // Cosmic microwave background
  MatrixCBR(iy, stokes_dim, f_grid, verbosity);

  // Star background
  if (star_do) {
    for (Index i_star = 0; i_star < star_pos.nelem(); i_star++) {
      get_star_background(iy,
                          star_pos[i_star],
                          star_radius[i_star][0],
                          star_spectrum[i_star],
                          rtp_pos,
                          rtp_los,
                          refellipsoid,
                          atmosphere_dim);
    }
  }
}