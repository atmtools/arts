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


//! An attribute to classify the particle type (ptype) of a SingleScatteringData
//structure (a scattering element).
/*! 
  PTYPE_GENERAL      General case
  PTYPE_TOTAL_RND    Totally randomly oriented particles
  PTYPE_AZIMUTH_RND  Azimuthally randomly oriented particles

  A detailed description of the different cases can be found in AUG.

*/
enum PType {
  PTYPE_GENERAL = 10,
  PTYPE_TOTAL_RND = 20,
  PTYPE_AZIMUTH_RND = 30
};


//! An attribute to classify the method to be used for SingleScatteringData
/*!
  NONE         Dummy value
  TMATRIX      T-Matrix method
*/
enum ParticleSSDMethod {
  PARTICLE_SSDMETHOD_NONE = 0,
  PARTICLE_SSDMETHOD_TMATRIX = 1,
};



/*===========================================================================
  === The SingleScatteringData structure
  ===========================================================================*/

/*!
   Structure which describes the single scattering properties of a scattering
   element (a single particle or a particle bulk).

   The fields of the structure are described in the ARTS user guide (AUG).
   It is listed as a sub-entry to "data structures".  
*/
struct SingleScatteringData {
  PType ptype;
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
typedef Array<Array<SingleScatteringData> > ArrayOfArrayOfSingleScatteringData;

ostream& operator<< (ostream& os, const SingleScatteringData& ssd);
ostream& operator<< (ostream& os, const ArrayOfSingleScatteringData& assd);


/*===========================================================================
  === The ScatteringMetaData structure
  ===========================================================================*/
/*!
   Structure which holds the meta data of a scattering element (a single
   particle or an ensemble of particles), mainly microphysical parameter
   necessary for calculating size and shape distributions.

  For a description of the structure members see built-in documentation of
  scat_meta_single in workspace.cc 
*/
struct ScatteringMetaData {
  String    description;
  String    source;
  String    refr_index;
  Numeric   mass;
  Numeric   diameter_max;
  Numeric   diameter_volume_equ;
  Numeric   diameter_area_equ_aerodynamical;
};

typedef Array<ScatteringMetaData> ArrayOfScatteringMetaData;
typedef Array<Array<ScatteringMetaData> > ArrayOfArrayOfScatteringMetaData;

ostream& operator<< (ostream& os, const ScatteringMetaData& ssd);
ostream& operator<< (ostream& os, const ArrayOfScatteringMetaData& assd);



// General functions:
// =============================================================

void scat_data_monoExtract(//Output
                           ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                           //Input
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Index& f_index );

void abs_vecTransform(//Output and Input
                      VectorView abs_vec_lab,
                      //Input
                      ConstTensor3View abs_vec_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Verbosity& verbosity);


void ext_matTransform(//Output and Input
                      MatrixView ext_mat_lab,
                      //Input
                      ConstTensor3View ext_mat_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Verbosity& verbosity);
 

void pha_matTransform(//Output
                      MatrixView pha_mat_lab,
                      //Input
                      ConstTensor5View pha_mat_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid,
                      const PType& ptype,
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

PType PTypeFromString(const String& ptype_string);
PType PType2FromString(const String& ptype_string);

String PTypeToString(const PType& ptype);

void ConvertAzimuthallyRandomSingleScatteringData(SingleScatteringData& ssd);

ParticleSSDMethod ParticleSSDMethodFromString(const String& particle_ssdmethod_string);

String ParticleSSDMethodToString(const ParticleSSDMethod& particle_ssdmethod_type);

#endif //optproperties_h
