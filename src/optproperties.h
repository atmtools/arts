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

  A detailed description of the different cases can be found in AUG.

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
  Vector    f_grid;
  Vector    T_grid;
  Vector    za_grid;
  Vector    aa_grid;
  Tensor6   pha_mat_data;
  Tensor4   ext_mat_data;
  Tensor4   abs_vec_data;
};

typedef Array<SingleScatteringData> ArrayOfSingleScatteringData;

ostream& operator<< (ostream &os, const SingleScatteringData &ssd);
ostream& operator<< (ostream &os, const ArrayOfSingleScatteringData &assd);


// General functions:
// =============================================================

void abs_vecTransform(//Output and Input
                      VectorView abs_vec_lab,
                      //Input
                      const Tensor3View abs_vec_data,
                      const VectorView za_datagrid,
                      const VectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca);


void ext_matTransform(//Output and Input
                      MatrixView ext_mat_lab,
                      //Input
                      const Tensor3View ext_mat_data,
                      const VectorView za_datagrid,
                      const VectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca);
 

void pha_matTransform(//Output
                      MatrixView pha_mat_lab,
                      //Input
                      const Tensor5View pha_mat_data,
                      const VectorView za_datagrid,
                      const VectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Numeric& za_inc,
                      const Numeric& aa_inc);


// Functions for the case: Randomly oriented particles: 
// ========================================================

void interpolate_scat_angle(//Output:
                            VectorView pha_mat_int,
                            Numeric& theta_rad,
                            //Input:
                            const Tensor5View pha_mat_data,
                            const VectorView za_datagrid,
                            const Numeric& za_sca_rad,
                            const Numeric& aa_sca_rad,
                            const Numeric& za_inc_rad,
                            const Numeric& aa_inc_rad);

void pha_mat_labCalc(//Output:
                      MatrixView pha_mat_lab,
                      //Input:
                      const VectorView pha_mat_int,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Numeric& za_inc,
                      const Numeric& aa_inc,
                      const Numeric& theta_rad);


#endif //optproperties_h
