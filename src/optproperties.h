/* Copyright (C) 2003 Claudia Emde <claudia@sat.physik.uni-bremen.de>

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
  \file   optproperties.h
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   2003-03-06
  
  \brief  Scattering database structure and functions.
  
   This file contains the definition of the SingleScatteringData structure
   and the functions in optproperties.cc that are of interest elsewhere.
*/


#ifndef optproperties_h
#define optproperties_h



//! An attribute to classify the particle type in a SingleScatteringData
//  structure.
/*! 
  GENERAL      General case
  MACROS_ISO   Macroscopically isotropic and mirror-symmetric scattering media
  HORIZ_AL     Horizonatally aligned plates and columns
  SPHERICAL    Spherical particles

  A detailed description of the differnent cases can be found in AUG.

*/
typedef enum{
  PTYPE_GENERAL = 10,
  PTYPE_MACROS_ISO = 20,
  PTYPE_HORIZ_AL = 30,
  PTYPE_SPHERICAL = 40
} PType;



/*===========================================================================
  === The SingleScatteringData structure
  ===========================================================================*/

//! Structure which describes the single scattering properties of a 
//  particle or a particle distribution.
/*! 
   The fields of the structure are described in the ARTS user guide (AUG).
   It is listed as a sub-entry to "data structures".  
*/
struct SingleScatteringData {
  PType     ptype;
  String    description;
  Vector    za_grid_sca;
  Vector    aa_grid_sca;
  Vector    za_grid_inc;
  Vector    aa_grid_inc;
  Tensor6   pha_mat_data;
  Tensor4   ext_mat_data;
  Tensor4   abs_vec_data;
};

#endif //optproperties_h
