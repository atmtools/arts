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
#include "auto_md.h"

/*!
  \file   m_fluxes.cc
  \author Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2018-06-22

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

//extern const Numeric PI;
//extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void AngularGridsSetFluxCalc2(Vector& za_grid,
                             Vector& aa_grid,
                             Vector& za_grid_weights,
                             // Keywords:
                             const Index& N_za_grid,
                             const Index& N_aa_grid,
                             const String& za_grid_type,
                             const Verbosity&) {
    // Azimuth angle grid
}