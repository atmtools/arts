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
void starBlackbodySimple(ArrayOfStar &star,
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

  Star star_temp;

  // spectrum
  star_temp.spectrum=Matrix(f_grid.nelem(), stokes_dim,0. );

  planck(star_temp.spectrum(joker,0), f_grid, temperature);
  star_temp.spectrum *= PI ; // outgoing flux at the surface of the star.


  star_temp.description = "Blackbody star" ;
  star_temp.radius = radius;
  star_temp.distance = distance;
  star_temp.latitude = latitude;
  star_temp.longitude = longitude;

  // set flag
  star_do = 1;

  //append
  star.push_back(star_temp);

}

/* Workspace method: Doxygen documentation will be auto-generated */
void starFromGrid(ArrayOfStar &star,
                         Index &star_do,
                         // Inputs:
                         const Vector &f_grid,
                         const Index &stokes_dim,
                         const GriddedField2& star_spectrum_raw,
                         const Numeric &radius,
                         const Numeric &distance,
                         const Numeric &temperature,
                         const Numeric &latitude,
                         const Numeric &longitude,
                         const String &description,
                         const Verbosity &verbosity) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the star (",distance," m) \n"
                     " is smaller than the radius of the star (", radius," m )")

  // interpolate field
  Matrix int_data = regrid_star_spectrum(star_spectrum_raw, f_grid, stokes_dim, temperature, verbosity);

  // create star
  Star star_temp;
  
  star_temp.spectrum = int_data; // set spectrum
  star_temp.spectrum *= PI; // outgoing flux at the surface of the star.

  star_temp.description = description;
  star_temp.radius = radius;
  star_temp.distance = distance;
  star_temp.latitude = latitude;
  star_temp.longitude = longitude;

  // set flag
  star_do = 1;

  //append
  star.push_back(star_temp);
}

void starOff(Index &star_do,
             ArrayOfStar &star,
             const Verbosity &){

  // set flag to False (default)
  star_do = 0;

  // create empty Array of Matrix for the star_spectrum
  star.resize(0);

}
