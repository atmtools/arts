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

#include "gridded_fields.h"
#include "matpack_data.h"
#include "mystring.h"
#include <matpack.h>
#include <rtepack.h>

//! An attribute to classify the particle type (ptype) of a SingleScatteringData
//structure (a scattering element).
/*! 
  PTYPE_GENERAL      General case
  PTYPE_TOTAL_RND    Totally randomly oriented particles
  PTYPE_AZIMUTH_RND  Azimuthally randomly oriented particles

  A detailed description of the different cases can be found in AUG.

*/
enum PType : Index {
  PTYPE_GENERAL = 300,
  PTYPE_AZIMUTH_RND = 200,
  PTYPE_TOTAL_RND = 100,
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
  String description;
  Vector f_grid;
  Vector T_grid;
  Vector za_grid;
  Vector aa_grid;
  Tensor7 pha_mat_data;
  Tensor5 ext_mat_data;
  Tensor5 abs_vec_data;

  friend std::ostream& operator<<(std::ostream& os, const SingleScatteringData& ssd);
};

typedef Array<SingleScatteringData> ArrayOfSingleScatteringData;
typedef Array<Array<SingleScatteringData> > ArrayOfArrayOfSingleScatteringData;

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
  String description;
  String source;
  String refr_index;
  Numeric mass;
  Numeric diameter_max;
  Numeric diameter_volume_equ;
  Numeric diameter_area_equ_aerodynamical;

  friend std::ostream& operator<<(std::ostream& os, const ScatteringMetaData& ssd);
};

typedef Array<ScatteringMetaData> ArrayOfScatteringMetaData;
typedef Array<Array<ScatteringMetaData> > ArrayOfArrayOfScatteringMetaData;

// General functions:
// =============================================================

void opt_prop_Bulk(  //Output
    Tensor5& ext_mat,
    Tensor4& abs_vec,
    Index& ptype,
    //Input
    const ArrayOfTensor5& ext_mat_ss,
    const ArrayOfTensor4& abs_vec_ss,
    const ArrayOfIndex& ptypes_ss);

void opt_prop_ScatSpecBulk(  //Output
    ArrayOfTensor5& ext_mat,
    ArrayOfTensor4& abs_vec,
    ArrayOfIndex& ptype,
    //Input
    const ArrayOfArrayOfTensor5& ext_mat_se,
    const ArrayOfArrayOfTensor4& abs_vec_se,
    const ArrayOfArrayOfIndex& ptypes_se,
    ConstMatrixView pnds,
    ConstMatrixView t_ok);

void opt_prop_NScatElems(  //Output
    ArrayOfArrayOfTensor5& ext_mat,
    ArrayOfArrayOfTensor4& abs_vec,
    ArrayOfArrayOfIndex& ptypes,
    Matrix& t_ok,
    //Input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& T_array,
    const Matrix& dir_array,
    const Index& f_index,
    const Index& t_interp_order = 1);

void opt_prop_1ScatElem(  //Output
    Tensor5View ext_mat,
    Tensor4View abs_vec,
    Index& ptype,
    VectorView t_ok,
    //Input
    const SingleScatteringData& ssd,
    const Vector& T_array,
    const Matrix& dir_array,
    const Index& f_index,
    const Index& t_interp_order = 1);

void ext_mat_SSD2Stokes(  //Output
    MatrixView ext_mat_stokes,
    //Input
    ConstVectorView ext_mat_ssd,
    const Index& ptype);

void abs_vec_SSD2Stokes(  //Output
    VectorView abs_vec_stokes,
    //Input
    ConstVectorView abs_vec_ssd,
    const Index& ptype);

void pha_mat_Bulk(  //Output
    Tensor6& pha_mat,
    Index& ptype,
    //Input
    const ArrayOfTensor6& pha_mat_ss,
    const ArrayOfIndex& ptypes_ss);

void pha_mat_ScatSpecBulk(  //Output
    ArrayOfTensor6& pha_mat,
    ArrayOfIndex& ptype,
    //Input
    const ArrayOfArrayOfTensor6& pha_mat_se,
    const ArrayOfArrayOfIndex& ptypes_se,
    ConstMatrixView pnds,
    ConstMatrixView t_ok);

void pha_mat_NScatElems(  //Output
    ArrayOfArrayOfTensor6& pha_mat,
    ArrayOfArrayOfIndex& ptypes,
    Matrix& t_ok,
    //Input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& T_array,
    const Matrix& pdir_array,
    const Matrix& idir_array,
    const Index& f_index,
    const Index& t_interp_order = 1);

void pha_mat_1ScatElem(  //Output
    Tensor6View pha_mat,
    Index& ptype,
    VectorView t_ok,
    //Input
    const SingleScatteringData& ssd,
    const Vector& T_array,
    const Matrix& pdir_array,
    const Matrix& idir_array,
    const Index& f_start,
    const Index& t_interp_order = 1);

void abs_vecTransform(  //Output and Input
    StokvecVector& abs_vec_lab,
    //Input
    ConstTensor3View abs_vec_data,
    ConstVectorView za_datagrid,
    ConstVectorView aa_datagrid,
    const PType& ptype,
    const Numeric& za_sca,
    const Numeric& aa_sca);

void ext_matTransform(  //Output and Input
    PropmatVector& ext_mat_lab,
    //Input
    ConstTensor3View ext_mat_data,
    ConstVectorView za_datagrid,
    ConstVectorView aa_datagrid,
    const PType& ptype,
    const Numeric& za_sca,
    const Numeric& aa_sca);

void pha_matTransform(  //Output
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
    ConstVectorView za_grid,
    ConstVectorView aa_grid);

void ext_matFromabs_vec(  //Output
    MatrixView ext_mat,
    //Input
    ConstVectorView abs_vec);

// Functions for the case: Randomly oriented particles:
// ========================================================

Numeric scat_angle(const Numeric& za_sca,
                   const Numeric& aa_sca,
                   const Numeric& za_inc,
                   const Numeric& aa_inc);

void interpolate_scat_angle(  //Output:
    VectorView pha_mat_int,
    //Input:
    ConstTensor5View pha_mat_data,
    ConstVectorView za_datagrid,
    const Numeric theta);

void pha_mat_labCalc(  //Output:
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

void opt_prop_sum_propmat_clearsky(  //Output:
    PropmatVector& ext_mat,
    StokvecVector& abs_vec,
    //Input:
    const PropmatVector& propmat_clearsky);

PType PTypeFromString(const String& ptype_string);
PType PType2FromString(const String& ptype_string);

String PTypeToString(const PType& ptype);

void ConvertAzimuthallyRandomSingleScatteringData(SingleScatteringData& ssd);

ParticleSSDMethod ParticleSSDMethodFromString(
    const String& particle_ssdmethod_string);

String ParticleSSDMethodToString(
    const ParticleSSDMethod& particle_ssdmethod_type);

void ext_abs_pfun_from_tro(MatrixView ext_data,
                           MatrixView abs_data,
                           Tensor3View pfun_data,
                           const ArrayOfSingleScatteringData& scat_data,
                           const Index& iss,
                           ConstMatrixView pnd_data,
                           ArrayOfIndex& cloudbox_limits,
                           ConstVectorView T_grid,
                           ConstVectorView sa_grid,
                           const Index f_index = -1);

#endif  //optproperties_h
