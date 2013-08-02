/* Copyright (C) 2003-2012 Claudia Emde <claudia.emde@dlr.de>

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
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   2003-03-06
  
  \brief  Scattering database structure and functions.
  
   This file contains the definition of the SingleScatteringData structure
   and the functions in optproperties.cc that are of interest elsewhere.
*/


#ifndef optproperties_h
#define optproperties_h

#include "matpackVII.h"
#include "mystring.h"
#include "messages.h"
#include "gridded_fields.h"


//! An attribute to classify the particle type in a SingleScatteringData
//  structure.
/*! 
  GENERAL      General case
  MACROS_ISO   Macroscopically isotropic and mirror-symmetric scattering media
  HORIZ_AL     Horizonatally aligned plates and columns
  SPHERICAL    Spherical particles

  A detailed description of the different cases can be found in AUG.

*/
enum ParticleType {
  PARTICLE_TYPE_GENERAL = 10,
  PARTICLE_TYPE_MACROS_ISO = 20,
  PARTICLE_TYPE_HORIZ_AL = 30,
  PARTICLE_TYPE_SPHERICAL = 40
};



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
  ParticleType ptype;
  String       description;
  Vector       f_grid;
  Vector       T_grid;
  Vector       za_grid;
  Vector       aa_grid;
  Tensor7      pha_mat_data;
  Tensor5      ext_mat_data;
  Tensor5      abs_vec_data;
};

typedef Array<SingleScatteringData> ArrayOfSingleScatteringData;

ostream& operator<< (ostream& os, const SingleScatteringData& ssd);
ostream& operator<< (ostream& os, const ArrayOfSingleScatteringData& assd);


/*===========================================================================
  === The ScatteringMetaData structure
  ===========================================================================*/
/*!
   Structure which holds the meta scattering properties of a 
   particle or a particle distribution, like (size parameters, shape, etc.). 
*/
struct ScatteringMetaData {
  String    description;
  String    type;
  String    shape;
  String    p_type;
  Numeric   density;
  Numeric   d_max;
  Numeric   V;
  Numeric   A_projec;
  Numeric   asratio;
  Vector    f_grid;
  Vector    T_grid;
  Tensor3   ref_index;
};

typedef Array<ScatteringMetaData> ArrayOfScatteringMetaData;

ostream& operator<< (ostream& os, const ScatteringMetaData& ssd);
ostream& operator<< (ostream& os, const ArrayOfScatteringMetaData& assd);



// General functions:
// =============================================================

void abs_vecTransform(//Output and Input
                      VectorView abs_vec_lab,
                      //Input
                      ConstTensor3View abs_vec_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid,
                      const ParticleType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Verbosity& verbosity);


void ext_matTransform(//Output and Input
                      MatrixView ext_mat_lab,
                      //Input
                      ConstTensor3View ext_mat_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid,
                      const ParticleType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Verbosity& verbosity);
 

void pha_matTransform(//Output
                      MatrixView pha_mat_lab,
                      //Input
                      ConstTensor5View pha_mat_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid,
                      const ParticleType& ptype,
                      const Index& za_sca_idx,
                      const Index& aa_sca_idx,
                      const Index& za_inc_idx,
                      const Index& aa_inc_idx,
                      ConstVectorView scat_za_grid,
                      ConstVectorView scat_aa_grid,
                      const Verbosity& verbosity);


void ext_matFromabs_vec(//Output:
                        MatrixView ext_mat,
                        //Input:
                        ConstVectorView abs_vec,
                        const Index& stokes_dim);

// Functions for the case: Randomly oriented particles: 
// ========================================================

void interpolate_scat_angle(//Output:
                            VectorView pha_mat_int,
                            Numeric& theta_rad,
                            //Input:
                            ConstTensor5View pha_mat_data,
                            ConstVectorView za_datagrid,
                            const Numeric& za_sca,
                            const Numeric& aa_sca,
                            const Numeric& za_inc,
                            const Numeric& aa_inc);


void pha_mat_labCalc(//Output:
                      MatrixView pha_mat_lab,
                      //Input:
                      ConstVectorView pha_mat_int,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Numeric& za_inc,
                      const Numeric& aa_inc,
                      const Numeric& theta_rad);


// Get ext_mat and abs_vec from propmat_clearsky:
// ========================================================

void opt_prop_sum_propmat_clearsky(//Output:
                                      Tensor3&         ext_mat,
                                      Matrix&          abs_vec,
                                      //Input:
                                      const Tensor4    propmat_clearsky);

#endif //optproperties_h
