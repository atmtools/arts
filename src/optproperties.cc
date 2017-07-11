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
  === File description
  ===========================================================================*/

/*!
  \file   optproperties.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Thu Mar  6 11:29:59 2003
  
  \brief  This file contains definitions and functions related to the
          optical properties of particles.
  
  
 */


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "check_input.h"
#include "interpolation.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "optproperties.h"
#include "xml_io.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;


#define F11 pha_mat_int[0]
#define F12 pha_mat_int[1]
#define F22 pha_mat_int[2]
#define F33 pha_mat_int[3]
#define F34 pha_mat_int[4]
#define F44 pha_mat_int[5]

#define PART_TYPE scat_data[i_ss][i_se].ptype
#define F_DATAGRID scat_data[i_ss][i_se].f_grid
#define T_DATAGRID scat_data[i_ss][i_se].T_grid
#define ZA_DATAGRID scat_data[i_ss][i_se].za_grid
#define AA_DATAGRID scat_data[i_ss][i_se].aa_grid
#define PHA_MAT_DATA_RAW scat_data[i_ss][i_se].pha_mat_data  //CPD: changed from pha_mat_data
#define EXT_MAT_DATA_RAW scat_data[i_ss][i_se].ext_mat_data  //which wouldn't let me play with
#define ABS_VEC_DATA_RAW scat_data[i_ss][i_se].abs_vec_data  //scat_data_mono.



//! Extraction of scat_data_mono from pre-interpolated scat_data
/*! 
  \param scat_data_mono  as the WSV
  \param scat_data       as the WSV
  \param f_index         as the WSV
     
  \author Jana Mendrok
  \date   2017-07-11 
*/
void scat_data_monoExtract(//Output
                           ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                           //Input
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Index& f_index )
{
  //Initialise scat_data_mono
  scat_data_mono.resize(scat_data.nelem());

  Index this_f_index;

  // Loop over the included scattering species
  for (Index i_ss = 0; i_ss<scat_data.nelem(); i_ss++)
  {
      const Index N_se = scat_data[i_ss].nelem();

      //Initialise scat_data_mono
      scat_data_mono[i_ss].resize(N_se);

      // Loop over the included scattering elements
      for (Index i_se = 0; i_se < N_se; i_se++)
      {
          scat_data_mono[i_ss][i_se].ptype=PART_TYPE;
          scat_data_mono[i_ss][i_se].f_grid=F_DATAGRID;
          scat_data_mono[i_ss][i_se].T_grid=T_DATAGRID;
          scat_data_mono[i_ss][i_se].za_grid=ZA_DATAGRID;
          scat_data_mono[i_ss][i_se].aa_grid=AA_DATAGRID;

          //Phase matrix data
          if( PHA_MAT_DATA_RAW.nlibraries()==1 )
            this_f_index = 0;
          else
            this_f_index = f_index;
          scat_data_mono[i_ss][i_se].pha_mat_data.resize(1,
                                                  PHA_MAT_DATA_RAW.nvitrines(),
                                                  PHA_MAT_DATA_RAW.nshelves(),
                                                  PHA_MAT_DATA_RAW.nbooks(),
                                                  PHA_MAT_DATA_RAW.npages(),
                                                  PHA_MAT_DATA_RAW.nrows(),
                                                  PHA_MAT_DATA_RAW.ncols());
          scat_data_mono[i_ss][i_se].pha_mat_data(0, joker, joker, joker,
                                                     joker, joker, joker) =
            PHA_MAT_DATA_RAW(this_f_index, joker, joker, joker,
                             joker, joker, joker );

          //Extinction matrix data
          if( EXT_MAT_DATA_RAW.nshelves()==1 )
            this_f_index = 0;
          else
            this_f_index = f_index;
          scat_data_mono[i_ss][i_se].ext_mat_data.resize(1,
                                                  EXT_MAT_DATA_RAW.nbooks(),
                                                  EXT_MAT_DATA_RAW.npages(),
                                                  EXT_MAT_DATA_RAW.nrows(),
                                                  EXT_MAT_DATA_RAW.ncols());
          scat_data_mono[i_ss][i_se].ext_mat_data(0, joker,
                                                  joker, joker, joker) =
            EXT_MAT_DATA_RAW(this_f_index, joker, joker, joker, joker );

          //Absorption vector data
          if(  ABS_VEC_DATA_RAW.nshelves()==1 )
            this_f_index = 0;
          else
            this_f_index = f_index;
          scat_data_mono[i_ss][i_se].abs_vec_data.resize(1,
                                                  ABS_VEC_DATA_RAW.nbooks(),
                                                  ABS_VEC_DATA_RAW.npages(),
                                                  ABS_VEC_DATA_RAW.nrows(),
                                                  ABS_VEC_DATA_RAW.ncols());
          scat_data_mono[i_ss][i_se].abs_vec_data(0, joker,
                                                  joker, joker, joker) =
            ABS_VEC_DATA_RAW(this_f_index, joker, joker, joker, joker );
      }
  }
}


//! Transformation of absorption vector.
/*! 
  In the single scattering database the data of the absorption vector is 
  stored in different coordinate systems, depending on the type (ptype) of 
  the scattering particle (scattering element).

  See AUG for information about the different classifications of the scattering
  elements.

  Output and Input:
  \param abs_vec_lab Absorption vector in Laboratory frame.
  Input:
  \param abs_vec_data Absorption vector in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Type of scattering element.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.
     
  \author Claudia Emde
  \date   2003-05-24 
*/
void abs_vecTransform(//Output and Input
                      VectorView abs_vec_lab,
                      //Input
                      ConstTensor3View abs_vec_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid _U_,
                      const PType& ptype,
                      const Numeric& za_sca _U_,
                      const Numeric& aa_sca _U_,
                      const Verbosity& verbosity)
{
  const Index stokes_dim = abs_vec_lab.nelem();
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  switch (ptype){
      
    case PTYPE_GENERAL:
    {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         issue as long as only PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND are used,
         but will be a problem for PTYPE_GENERAL, ie needs to be fixed BEFORE
         adding PTYPE_GENERAL support (see AUG appendix for more info).
      */

      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
      break;
    }
    case PTYPE_TOTAL_RND:
    {
      // The first element of the vector corresponds to the absorption 
      // coefficient which is stored in the database, the others are 0.
      
      abs_vec_lab = 0;
      
      abs_vec_lab[0] = abs_vec_data(0,0,0);
      break;
    }
      
    case PTYPE_AZIMUTH_RND://Added by Cory Davis 9/12/03
    {
      assert (abs_vec_data.ncols() == 2);
      
      // In the case of azimuthally randomly oriented particles, only the first
      // two elements of the absorption coefficient vector are non-zero.
      // These values are dependent on the zenith angle of propagation.
      
      // 1st interpolate data by za_sca
      GridPos gp;
      Vector itw(2);
      
      gridpos(gp,za_datagrid,za_sca);
      interpweights(itw,gp);
      abs_vec_lab = 0;
      abs_vec_lab[0] = interp(itw,abs_vec_data(Range(joker),0,0),gp);
      
      if( stokes_dim == 1 ){
        break;
      }
      abs_vec_lab[1] = interp(itw,abs_vec_data(Range(joker),0,1),gp);
      break;
    }
    default:
    {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
  }  
}


//! Transformation of extinction matrix.
/*! 
  In the single scattering database the data of the extinction matrix is 
  stored in different coordinate systems, depending on the type (ptype) of 
  the scattering particle (scattering element).

  See AUG for information about the different classifications of the scattering
  elements.

  Output and Input:
  \param ext_mat_lab Extinction matrix in Laboratory frame.
  Input:
  \param ext_mat_data Extinction matrix in database.
  \param za_datagrid Zenith angle grid in the database.
  \param aa_datagrid Zenith angle grid in the database.
  \param ptype Type of scattering element.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.

  \author Claudia Emde
  \date   2003-05-24 
*/
void ext_matTransform(//Output and Input
                      MatrixView ext_mat_lab,
                      //Input
                      ConstTensor3View ext_mat_data,
                      ConstVectorView za_datagrid,
                      ConstVectorView aa_datagrid _U_,
                      const PType& ptype,
                      const Numeric& za_sca,
                      const Numeric& aa_sca _U_,
                      const Verbosity& verbosity)
{
  const Index stokes_dim = ext_mat_lab.ncols();

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  switch (ptype){
      
    case PTYPE_GENERAL:
    {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         issue as long as only PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND are used,
         but will be a problem for PTYPE_GENERAL, ie needs to be fixed BEFORE
         adding PTYPE_GENERAL support (see AUG appendix for more info).
      */

      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
      break;
    }
    case PTYPE_TOTAL_RND:
    {
      assert (ext_mat_data.ncols() == 1);
      
      // In the case of randomly oriented particles the extinction matrix is 
      // diagonal. The value of each element of the diagonal is the
      // extinction cross section, which is stored in the database.
      
      ext_mat_lab = 0.;
      
      ext_mat_lab(0,0) = ext_mat_data(0,0,0);
      
      
      if( stokes_dim == 1 ){
        break;
      }
      
      ext_mat_lab(1,1) = ext_mat_data(0,0,0);
      
      if( stokes_dim == 2 ){
        break;
      }
      
      ext_mat_lab(2,2) = ext_mat_data(0,0,0);
      
      if( stokes_dim == 3 ){
        break;
      }
      
      ext_mat_lab(3,3) = ext_mat_data(0,0,0);
      break;
    }
      
    case PTYPE_AZIMUTH_RND://Added by Cory Davis 9/12/03
    {
      assert (ext_mat_data.ncols() == 3);
      
      // In the case of azimuthally randomly oriented particles, the extinction
      // matrix has only 3 independent non-zero elements Kjj, K12=K21, and K34=-K43.
      // These values are dependent on the zenith angle of propagation.
      
      // 1st interpolate data by za_sca
      GridPos gp;
      Vector itw(2);
      Numeric Kjj;
      Numeric K12;
      Numeric K34;
      
      gridpos(gp,za_datagrid,za_sca);
      interpweights(itw,gp);
      
      ext_mat_lab=0.0;
      Kjj=interp(itw,ext_mat_data(Range(joker),0,0),gp);
      ext_mat_lab(0,0)=Kjj;
      
      if( stokes_dim == 1 ){
        break;
      }
      
      K12=interp(itw,ext_mat_data(Range(joker),0,1),gp);
      ext_mat_lab(1,1)=Kjj;
      ext_mat_lab(0,1)=K12;
      ext_mat_lab(1,0)=K12;
      
      if( stokes_dim == 2 ){
        break;
      }
      
      ext_mat_lab(2,2)=Kjj;
      
      if( stokes_dim == 3 ){
        break;
      }
      
      K34=interp(itw,ext_mat_data(Range(joker),0,2),gp);
      ext_mat_lab(2,3)=K34;
      ext_mat_lab(3,2)=-K34;
      ext_mat_lab(3,3)=Kjj;
      break;
      
    }
    default:
    {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
  }
}  
 

//! Transformation of phase matrix.
/*! 
  In the single scattering database the data of the phase matrix is 
  stored in different coordinate systems, depending on the type (ptype) of 
  the scattering particle (scattering element).

  See AUG for information about the different classifications of the scattering
  elements.

  \param[out,in] pha_mat_lab   Phase matrix in Laboratory frame.
  \param[in]     pha_mat_data  Phase matrix in database.
  \param[in]     za_datagrid   Zenith angle grid in the database.
  \param[in]     aa_datagrid   Zenith angle grid in the database.
  \param[in]     ptype Type of scattering element.
  \param[in]     za_sca_idx    Index of zenith angle of scattered direction
                                 within scat_za_grid.
  \param[in]     aa_sca_idx    Index of azimuth angle of scattered direction
                                 within scat_aa_grid.
  \param[in]     za_inc_idx    Index of zenith angle of incoming direction
                                 within scat_za_grid.
  \param[in]     aa_inc_idx    Index of azimuth angle of incoming direction
                                 within scat_aa_grid.
  \param[in]     scat_za_grid  Grid of zenith angles to extract pha_mat for.
  \param[in]     scat_aa_grid  Grid of azimuth angles to extract pha_mat for.
  
  \author Claudia Emde
  \date   2003-08-19
*/
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
                      const Verbosity& verbosity)
{
  const Index stokes_dim = pha_mat_lab.ncols();
  
  Numeric za_sca = scat_za_grid[za_sca_idx]; 
  Numeric aa_sca = scat_aa_grid[aa_sca_idx];
  Numeric za_inc = scat_za_grid[za_inc_idx]; 
  Numeric aa_inc = scat_aa_grid[aa_inc_idx];
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  switch (ptype){
      
    case PTYPE_GENERAL:
    {
      /*
         TO ANY DEVELOPER:
         current usage of coordinate systems in scattering solvers (RT and SSD
         extraction) and general radiative transfer is not consistent. Not an
         issue as long as only PTYPE_TOTAL_RND and PTYPE_AZIMUTH_RND are used,
         but will be a problem for PTYPE_GENERAL, ie needs to be fixed BEFORE
         adding PTYPE_GENERAL support (see AUG appendix for more info).
      */

      CREATE_OUT0;
      out0 << "Case PTYPE_GENERAL not yet implemented. \n"; 
      break;
    }
    case PTYPE_TOTAL_RND:
    {
      // Calculate the scattering and interpolate the data on the scattering
      // angle:
      
      Vector pha_mat_int(6);
      Numeric theta_rad;
      
      // Interpolation of the data on the scattering angle:
      interpolate_scat_angle(pha_mat_int, theta_rad, pha_mat_data,
                             za_datagrid, za_sca, aa_sca,
                             za_inc, aa_inc);
      
      // Calculate the phase matrix in the laboratory frame:
      pha_mat_labCalc(pha_mat_lab, pha_mat_int, za_sca, aa_sca, za_inc, 
                      aa_inc, theta_rad);
      
      break;
    }
      
    case PTYPE_AZIMUTH_RND://Added by Cory Davis
                           //Data is already stored in the laboratory frame,
                           //but it is compressed a little.  Details elsewhere.
    {
      assert (pha_mat_data.ncols()==16);
      assert (pha_mat_data.npages()==za_datagrid.nelem());
      Numeric delta_aa=aa_sca-aa_inc+(aa_sca-aa_inc<-180)*360-
      (aa_sca-aa_inc>180)*360;//delta_aa corresponds to the "books" 
                              //dimension of pha_mat_data
      GridPos za_sca_gp;
      GridPos delta_aa_gp;
      GridPos za_inc_gp;
      Vector itw(8);
      
      gridpos(delta_aa_gp,aa_datagrid,abs(delta_aa));
      gridpos(za_inc_gp,za_datagrid,za_inc);
      gridpos(za_sca_gp,za_datagrid,za_sca);
      
      interpweights(itw,za_sca_gp,delta_aa_gp,za_inc_gp);
      
      pha_mat_lab(0,0)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,0),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      if( stokes_dim == 1 ){
        break;
      }
      pha_mat_lab(0,1)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,1),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      pha_mat_lab(1,0)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,4),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      pha_mat_lab(1,1)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,5),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      if( stokes_dim == 2 ){
        break;
      }
      if (delta_aa>=0)
      {
        pha_mat_lab(0,2)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,2),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(1,2)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,6),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(2,0)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,8),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(2,1)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,9),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
      }
      else
      {
        pha_mat_lab(0,2)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,2),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(1,2)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,6),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(2,0)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,8),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(2,1)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,9),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
      }                             
      pha_mat_lab(2,2)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,10),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      if( stokes_dim == 3 ){
        break;
      }
      if (delta_aa>=0)
      {
        pha_mat_lab(0,3)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,3),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(1,3)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,7),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(3,0)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,12),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(3,1)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                 Range(joker),0,13),
                                za_sca_gp,delta_aa_gp,za_inc_gp);
      }
      else
      {
        pha_mat_lab(0,3)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,3),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(1,3)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,7),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(3,0)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,12),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
        pha_mat_lab(3,1)=-interp(itw,pha_mat_data(Range(joker),Range(joker),
                                                  Range(joker),0,13),
                                 za_sca_gp,delta_aa_gp,za_inc_gp);
      }
      pha_mat_lab(2,3)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,11),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      pha_mat_lab(3,2)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,14),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      pha_mat_lab(3,3)=interp(itw,pha_mat_data(Range(joker),Range(joker),
                                               Range(joker),0,15),
                              za_sca_gp,delta_aa_gp,za_inc_gp);
      break;  
    }

    default:
    {
      CREATE_OUT0;
      out0 << "Not all ptype cases are implemented\n";
    }
  }
    
}



//! Derive extinction matrix from absorption vector.
/*! 
  In case, when only absorption of a scattering element shall be considered, and
  the scattering is negelected, the extinction matrix is set from the absorption
  vector only.

  Extinction matrix is set the following way:
  
  K11 = K22 = K33 = K44 = a1
  K12 = K21 = a2
  K13 = K31 = a3
  K14 = K41 = a4

  Other elements remain 0.
  However, note that the other elements might be supposed to contain non-zero
  values as well. We couldn't find an appropriate solution yet, though, and
  these elements are expected to be rather small. Also, using the absorption
  part only is anyway only a rough approximation. Hence, we deem the assumption
  of setting the remaining elements to zero reasonable.

  Output and Input:
  \param ext_mat     Extinction matrix.
  Input: 
  \param abs_vec     Absorption vector.
  \param stokes_dim  as the WSV.
  
  \author Jana Mendrok
  \date   2013-04-30 
*/
void ext_matFromabs_vec(//Output:
                        MatrixView ext_mat,
                        //Input:
                        ConstVectorView abs_vec,
                        const Index& stokes_dim)
{
  assert( stokes_dim>=1  &&  stokes_dim<=4 );
  assert( ext_mat.nrows() == stokes_dim );
  assert( ext_mat.ncols() == stokes_dim );
  assert( abs_vec.nelem() == stokes_dim );

  // first: diagonal elements
  for (Index is=0; is<stokes_dim; is++)
    {
      ext_mat(is,is) += abs_vec[0];
    }
  // second: off-diagonal elements, namely first row and column
  for (Index is=1; is<stokes_dim; is++)
    {
      ext_mat(0,is) += abs_vec[is];
      ext_mat(is,0) += abs_vec[is];
    }
}



//! Interpolate data on the scattering angle.
/*! 
  This function is used for the transformation of the phase matrix 
  from scattering frame to the laboratory frame for randomly oriented
  scattering media (case PTYPE_MACRO_ISO).

  The scattering angle is calculated from the angles defining
  the directions of the incoming and scattered radiation.
  After that the data (which is stored in the data files as a function
  of the scattering angle) is interpolated on the calculated 
  scattering angle.

  \param[out] pha_mat_int     Interpolated phase matrix.
  \param[in]  pha_mat_data    Phase matrix in database.
  \param[in]  za_sca_idx      Index of zenith angle of scattered direction.
  \param[in]  aa_sca_idx      Index of azimuth angle of scattered direction.
  \param[in]  za_inc_idx      Zenith angle of incoming direction.
  \param[in]  aa_inc_idx      Azimuth angle of incoming direction.
  \param[in]  scat_theta_gps  Array of gridposizions for scattering angle.
  \param[in]  scat_theta_itws Interpolation weights belonging to the scattering
                              angles.
     
  \author Claudia Emde
  \date   2003-05-13 
*/
void interpolate_scat_angleDOIT(//Output:
                            VectorView pha_mat_int,
                            //Input:
                            ConstTensor5View pha_mat_data,
                            const Index& za_sca_idx,
                            const Index& aa_sca_idx,
                            const Index& za_inc_idx,
                            const Index& aa_inc_idx,
                            const ArrayOfArrayOfArrayOfArrayOfGridPos&
                                scat_theta_gps,
                            ConstTensor5View scat_theta_itws
                            )
{
  
  ConstVectorView itw = scat_theta_itws(za_sca_idx, aa_sca_idx, za_inc_idx, aa_inc_idx, joker);
  
  for (Index i = 0; i < 6; i++)
    {
      pha_mat_int[i] = interp(itw, pha_mat_data(joker, 0, 0, 0, i), 
                             scat_theta_gps[za_sca_idx][aa_sca_idx][za_inc_idx][aa_inc_idx]);
    }
  
} 


//! Interpolate data on the scattering angle.
/*! 
  This function is used for the transformation of the phase matrix 
  from scattering frame to the laboratory frame for randomly oriented
  scattering media (case PTYPE_MACRO_ISO).

  The scattering angle is calculated from the angles defining
  the directions of the incoming and scattered radiation.
  After that the data (which is stored in the data files as a function
  of the scattering angle) is interpolated on the calculated 
  scattering angle.

  \param[out] pha_mat_int  Interpolated phase matrix.
  \param[out] theta_rad    Scattering angle [rad].
  \param[in]  pha_mat_data Phase matrix in database.
  \param[in]  za_datagrid  Zenith angle grid in the database.
  \param[in]  za_sca       Zenith angle of scattered direction [rad].
  \param[in]  aa_sca       Azimuth angle of scattered direction [rad].
  \param[in]  za_inc       Zenith angle of incoming direction [rad].
  \param[in]  aa_inc       Azimuth angle of incoming direction [rad].
     
  \author Claudia Emde
  \date   2003-08-19
*/
void interpolate_scat_angle(//Output:
                            VectorView pha_mat_int,
                            Numeric& theta_rad,
                            //Input:
                            ConstTensor5View pha_mat_data,
                            ConstVectorView za_datagrid,
                            const Numeric& za_sca,
                            const Numeric& aa_sca,
                            const Numeric& za_inc,
                            const Numeric& aa_inc
                            )
{
  Numeric ANG_TOL=1e-7;

  //Calculate scattering angle from incident and scattered directions.
  //The two special cases are implemented here to avoid NaNs that can 
  //sometimes occur in in the acos... formula in forward and backscatter
  //cases. CPD 5/10/03.
  //
  // Consider not only aa_sca-aa_inc ~= 0, but also aa_sca-aa_inc ~= 360.
  // GH 2011-05-31
  
  if( (abs(aa_sca-aa_inc)<ANG_TOL) || (abs(abs(aa_sca-aa_inc)-360) < ANG_TOL) )
    {
      theta_rad=DEG2RAD*abs(za_sca-za_inc);
    }
  else if (abs(abs(aa_sca-aa_inc)-180)<ANG_TOL)
    {
      theta_rad=DEG2RAD*(za_sca+za_inc);
      if (theta_rad>PI){theta_rad=2*PI-theta_rad;}
    }
  else
    {
      const Numeric za_sca_rad = za_sca * DEG2RAD;
      const Numeric za_inc_rad = za_inc * DEG2RAD;
      const Numeric aa_sca_rad = aa_sca * DEG2RAD;
      const Numeric aa_inc_rad = aa_inc * DEG2RAD;
      
      // cout << "Interpolation on scattering angle" << endl;
      assert (pha_mat_data.ncols() == 6);
      // Calculation of the scattering angle:
      theta_rad = acos(cos(za_sca_rad) * cos(za_inc_rad) + 
                       sin(za_sca_rad) * sin(za_inc_rad) * 
                       cos(aa_sca_rad - aa_inc_rad));
   }
      const Numeric theta = RAD2DEG * theta_rad;
      
  // Interpolation of the data on the scattering angle:
 
      GridPos thet_gp;
      gridpos(thet_gp, za_datagrid, theta);
      
      Vector itw(2);
      interpweights(itw, thet_gp);
      
      for (Index i = 0; i < 6; i++)
        {
          pha_mat_int[i] = interp(itw, pha_mat_data(joker, 0, 0, 0, i), 
                              thet_gp);
    }
} 




//! Calculate phase matrix in laboratory coordinate system.
/*! 
  Transformation function for the phase matrix for the case of
  randomly oriented particles (case PTYPE_MACRO_ISO).
  
  Some of the formulas can be found in 

  Mishchenkho: "Scattering, Absorption and Emission of Light 
  by Small Particles", Cambridge University Press, 2002
  Capter 4

  The full set of formulas will be documented in AUG.

  Output and Input:
  \param pha_mat_lab Phase matrix in laboratory frame.
  Input: 
  \param pha_mat_int Interpolated phase matrix.
  \param za_sca Zenith angle of scattered direction.
  \param aa_sca Azimuth angle of scattered direction.
  \param za_inc Zenith angle of incoming direction.
  \param aa_inc Azimuth angle of incoming direction.
  \param theta_rad Scattering angle [rad].
  
  \author Claudia Emde
  \date   2003-05-13 
*/
void pha_mat_labCalc(//Output:
                      MatrixView pha_mat_lab,
                      //Input:
                      ConstVectorView pha_mat_int,
                      const Numeric& za_sca,
                      const Numeric& aa_sca,
                      const Numeric& za_inc,
                      const Numeric& aa_inc,
                      const Numeric& theta_rad)
{
  Numeric za_sca_rad = za_sca * DEG2RAD;
  Numeric za_inc_rad = za_inc * DEG2RAD;
  Numeric aa_sca_rad = aa_sca * DEG2RAD;
  Numeric aa_inc_rad = aa_inc * DEG2RAD;

  const Index stokes_dim = pha_mat_lab.ncols();

  // For stokes_dim = 1, we only need Z11=F11:
  pha_mat_lab(0,0) = F11;

  if( isnan(pha_mat_lab(0,0)) )
    {
      throw runtime_error(
        "NaN value(s) detected in *pha_mat_labCalc* (0,0). Could the "
        "input data contain NaNs? Please check with *scat_dataCheck*. If "
        "input data are OK and you critically need the ongoing calculations, "
        "try to change the observation LOS slightly. If you can reproduce "
        "this error, please contact Patrick in order to help tracking down "
        "the reason to this problem. If you see this message occasionally "
        "when doing MC calculations, it should not be critical. This path "
        "sampling will be rejected and replaced with a new one." );
    }
  
  if( stokes_dim > 1 ){

    const Numeric ANGTOL_RAD = 1e-6; //CPD: this constant is used to adjust
                                     //zenith angles close to 0 and PI.  This is
                                     //also used to avoid float == float statements.

    //
    // Several cases have to be considered:
    //

    if(    ( abs(theta_rad) < ANGTOL_RAD )            // forward scattering
        || ( abs(theta_rad-PI) < ANGTOL_RAD )         // backward scattering
        || ( abs(aa_inc_rad-aa_sca_rad) < ANGTOL_RAD ) // inc and sca on meridian
        || ( abs(abs(aa_inc_rad-aa_sca_rad)-360.) < ANGTOL_RAD ) //   "
        || ( abs(abs(aa_inc_rad-aa_sca_rad)-180.) < ANGTOL_RAD ) //   "
      )
      {
        pha_mat_lab(0,1) = F12;
        pha_mat_lab(1,0) = F12;
        pha_mat_lab(1,1) = F22;
        
        if( stokes_dim > 2 ){
          pha_mat_lab(0,2) = 0;
          pha_mat_lab(1,2) = 0;
          pha_mat_lab(2,0) = 0;
          pha_mat_lab(2,1) = 0;
          pha_mat_lab(2,2) = F33;
          
          if( stokes_dim > 3 ){
            pha_mat_lab(0,3) = 0;
            pha_mat_lab(1,3) = 0;
            pha_mat_lab(2,3) = F34;
            pha_mat_lab(3,0) = 0;
            pha_mat_lab(3,1) = 0;
            pha_mat_lab(3,2) = -F34;
            pha_mat_lab(3,3) = F44;
          }
        }
      }
   
   else 
     {
       Numeric sigma1;
       Numeric sigma2;

       Numeric s1, s2;

       // In these cases we have to take limiting values.
 
       if (za_inc_rad < ANGTOL_RAD)
         {
           sigma1 = PI + aa_sca_rad - aa_inc_rad;
           sigma2 = 0;
         }
       else if (za_inc_rad > PI-ANGTOL_RAD)
         {
           sigma1 = aa_sca_rad - aa_inc_rad;
           sigma2 = PI; 
         }
       else if (za_sca_rad < ANGTOL_RAD)
         {
           sigma1 = 0;
           sigma2 = PI + aa_sca_rad - aa_inc_rad;
         }
       else if (za_sca_rad > PI - ANGTOL_RAD)
         {
           sigma1 = PI;
           sigma2 = aa_sca_rad - aa_inc_rad; 
         }
       else
         {
           s1 = (cos(za_sca_rad) - cos(za_inc_rad) * cos(theta_rad))
              /(sin(za_inc_rad)*sin(theta_rad));
           s2 = (cos(za_inc_rad) - cos(za_sca_rad) * cos (theta_rad))/
             (sin(za_sca_rad)*sin(theta_rad)); 
       
           sigma1 =  acos(s1);
           sigma2 =  acos(s2);
           
           // Arccos is only defined in the range from -1 ... 1
           // Numerical problems can appear for values close to 1 or -1
           // this (also) catches the case when inc and sca are on one meridian
           if ( isnan(sigma1) || isnan(sigma2) )
             {
               if ( abs(s1 - 1) < ANGTOL_RAD)
                 sigma1 = 0;
               if ( abs(s1 + 1) < ANGTOL_RAD)
                 sigma1 = PI;
               if ( abs(s2 - 1) < ANGTOL_RAD)
                 sigma2 = 0;
               if ( abs(s2 + 1) < ANGTOL_RAD)
                 sigma2 = PI;
             }
         }
      
       const Numeric C1 = cos(2*sigma1);
       const Numeric C2 = cos(2*sigma2);
        
       const Numeric S1 = sin(2*sigma1); 
       const Numeric S2 = sin(2*sigma2);
        
        pha_mat_lab(0,1) = C1 * F12;
        pha_mat_lab(1,0) = C2 * F12;
        pha_mat_lab(1,1) = C1 * C2 * F22 - S1 * S2 * F33;

        //assert(!isnan(pha_mat_lab(0,1)));        
        //assert(!isnan(pha_mat_lab(1,0)));
        //assert(!isnan(pha_mat_lab(1,1)));
        if( isnan(pha_mat_lab(0,1))  ||  isnan(pha_mat_lab(1,0)) ||
            isnan(pha_mat_lab(1,1)) )
          {
            throw runtime_error(
            "NaN value(s) detected in *pha_mat_labCalc* (0/1,1). Could the "
            "input data contain NaNs? Please check with *scat_dataChack*. If "
            "input data are OK  and you critically need the ongoing calculations, "
            "try to change the observation LOS slightly. If you can reproduce "
            "this error, please contact Patrick in order to help tracking down "
            "the reason to this problem. If you see this message occasionally "
            "when doing MC calculations, it should not be critical. This path "
            "sampling will be rejected and replaced with a new one." );
          }

        if( stokes_dim > 2 ){
          /*CPD: For skokes_dim > 2 some of the transformation formula 
            for each element have a different sign depending on whether or
            not 0<aa_scat-aa_inc<180.  For details see pages 94 and 95 of 
            Mishchenkos chapter in : 
            Mishchenko, M. I., and L. D. Travis, 2003: Electromagnetic 
            scattering by nonspherical particles. In Exploring the Atmosphere 
            by Remote Sensing Techniques (R. Guzzi, Ed.), Springer-Verlag, 
            Berlin, pp. 77-127. 
            This is available at http://www.giss.nasa.gov/~crmim/publications/ */
          Numeric delta_aa=aa_sca-aa_inc+(aa_sca-aa_inc<-180)*360-
            (aa_sca-aa_inc>180)*360;
          if(delta_aa>=0)
            {
              pha_mat_lab(0,2) = S1 * F12;
              pha_mat_lab(1,2) = S1 * C2 * F22 + C1 * S2 * F33;
              pha_mat_lab(2,0) = -S2 * F12;
              pha_mat_lab(2,1) = -C1 * S2 * F22 - S1 * C2 * F33;
            }
          else
            {
              pha_mat_lab(0,2) = -S1 * F12;
              pha_mat_lab(1,2) = -S1 * C2 * F22 - C1 * S2 * F33;
              pha_mat_lab(2,0) = S2 * F12;
              pha_mat_lab(2,1) = C1 * S2 * F22 + S1 * C2 * F33;
            }
          pha_mat_lab(2,2) = -S1 * S2 * F22 + C1 * C2 * F33;

          if( stokes_dim > 3 ){
            if(delta_aa>=0)
              {
                pha_mat_lab(1,3) = S2 * F34;
                pha_mat_lab(3,1) = S1 * F34;
              }
            else
              {
                pha_mat_lab(1,3) = -S2 * F34;
                pha_mat_lab(3,1) = -S1 * F34;
              }
            pha_mat_lab(0,3) = 0;
            pha_mat_lab(2,3) = C2 * F34;
            pha_mat_lab(3,0) = 0;
            pha_mat_lab(3,2) = -C1 * F34;
            pha_mat_lab(3,3) = F44;
          }
        }     
     }
   }
}
     

ostream& operator<< (ostream& os, const SingleScatteringData& /*ssd*/)
{
  os << "SingleScatteringData: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream& os, const ArrayOfSingleScatteringData& /*assd*/)
{
  os << "ArrayOfSingleScatteringData: Output operator not implemented";
  return os;
}

ostream& operator<< (ostream& os, const ScatteringMetaData& /*ssd*/)
{
  os << "ScatteringMetaData: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream& os, const ArrayOfScatteringMetaData& /*assd*/)
{
  os << "ArrayOfScatteringMetaData: Output operator not implemented";
  return os;
}


//! Get optical properties from propmat_clearsky
/*!
  This turns propmat_clearsky into the extinction matrix
  and absorption vector for use when these are important.

  Internal function to replace the old opt_prop_gas_agenda.
  
  Output and Input:
  \param ext_mat Extinction matrix.
  \param abs_vec Absorption vector.
  Input:
  \param propmat_clearsky as the WSV.

  \author Richard Larsson
  \date   2012-07-24
*/
void opt_prop_sum_propmat_clearsky(//Output:
                                      Tensor3&         ext_mat,
                                      Matrix&          abs_vec,
                                      //Input:
                                      const Tensor4    propmat_clearsky)
{

    Index stokes_dim = propmat_clearsky.ncols();
    
    Index freq_dim = propmat_clearsky.npages();

    // old abs_vecInit
    abs_vec.resize( freq_dim, stokes_dim );
    abs_vec = 0;                  // Initialize to zero!

    // old ext_matInit
    ext_mat.resize( freq_dim, stokes_dim, stokes_dim );
    ext_mat = 0;                  // Initialize to zero!

   // old ext_matAddGas and abs_vecAddGas for 0 vector and matrix
    for ( Index iv=0; iv<freq_dim; ++iv )
        for ( Index is1=0; is1<stokes_dim; ++is1 )
        {
            abs_vec(iv,is1) += propmat_clearsky(joker,iv,is1,0).sum();
            for ( Index is2=0; is2<stokes_dim; ++is2 )
                ext_mat(iv,is1,is2) += propmat_clearsky(joker,iv,is1,is2).sum();
        }
}


//! Convert ptype name to enum value
/*!
 Returns the PType enum value for the given String.

 This is the conversion for SingleScatteringData version 2.

 \param[in]  ptype_string  Particle type name
 \return     PType enum value

 \author Oliver Lemke
 */
PType PTypeFromString(const String& ptype_string)
{
    PType ptype;
    if (ptype_string == "general")
        ptype = PTYPE_GENERAL;
    else if (ptype_string == "totally_random")
        ptype = PTYPE_TOTAL_RND;
    else if (ptype_string == "azimuthally_random")
        ptype = PTYPE_AZIMUTH_RND;
    else
    {
        ostringstream os;
        os << "Unknown ptype: " << ptype_string << endl
           << "Valid types are: general, totally_random and "
           << "azimuthally_random.";
        throw std::runtime_error(os.str());
    }

    return ptype;
}


//! Convert ptype name to enum value
/*!
 Returns the PType enum value for the given String.
 
 This is the conversion for SingleScatteringData version 2.

 \param[in]  ptype_string  Particle type name
 \return     PType enum value

 \author Oliver Lemke
 */
PType PType2FromString(const String& ptype_string)
{
    PType ptype;
    if (ptype_string == "general")
        ptype = PTYPE_GENERAL;
    else if (ptype_string == "macroscopically_isotropic")
        ptype = PTYPE_TOTAL_RND;
    else if (ptype_string == "horizontally_aligned")
        ptype = PTYPE_AZIMUTH_RND;
    else
    {
        ostringstream os;
        os << "Unknown ptype: " << ptype_string << endl
           << "Valid types are: general, macroscopically_isotropic and "
           << "horizontally_aligned.";
        throw std::runtime_error(os.str());
    }

    return ptype;
}


//! Convert particle type enum value to String.
/*!
 Returns the PType enum value for the given String.

 \param[in]  ptype  Particle type
 \return     String representation of PType

 \author Oliver Lemke
 */
String PTypeToString(const PType& ptype)
{
    String ptype_string;

    switch (ptype) {
        case PTYPE_GENERAL:
            ptype_string = "general";
            break;
        case PTYPE_TOTAL_RND:
            ptype_string = "totally_random";
            break;
        case PTYPE_AZIMUTH_RND:
            ptype_string = "azimuthally_random";
            break;
        default:
            ostringstream os;
            os << "Internal error: Cannot map PType enum value "
            << ptype << " to String.";
            throw std::runtime_error(os.str());
            break;
    }

    return ptype_string;
}


//! Convert azimuthally-random oriented SingleScatteringData to latest version.
/*!
 Converts SingleScatteringData to version 3.

 \param[in,out]  ssd  SingleScatteringData

 \author Jana Mendrok
*/
void ConvertAzimuthallyRandomSingleScatteringData(SingleScatteringData& ssd)
{
  // First check that input fulfills requirements on older data formats:
  // 1) Is za_grid symmetric and includes 90deg?
  Index nza = ssd.za_grid.nelem();
  for (Index i=0; i<nza/2; i++)
    {

      if (!is_same_within_epsilon(180.-ssd.za_grid[nza-1-i],
                                  ssd.za_grid[i],
                                  2*DBL_EPSILON))
        {
            ostringstream os;
            os << "Zenith grid of azimuthally_random single scattering data\n"
               << "is not symmetric with respect to 90degree.";
            throw std::runtime_error(os.str());
        }
    }
    if (!is_same_within_epsilon(ssd.za_grid[nza/2], 90., 2*DBL_EPSILON))
        {
            ostringstream os;
            os << "Zenith grid of azimuthally_random single scattering data\n"
               << "does not contain 90 degree grid point.";
            throw std::runtime_error(os.str());
        }

  // 2) Are data sizes correct?
  ostringstream os_pha_mat;
  os_pha_mat << "pha_mat ";
  ostringstream os_ext_mat;
  os_ext_mat << "ext_mat ";
  ostringstream os_abs_vec;
  os_abs_vec << "abs_vec ";
  chk_size(os_pha_mat.str(), ssd.pha_mat_data,
           ssd.f_grid.nelem(), ssd.T_grid.nelem(),
           ssd.za_grid.nelem(), ssd.aa_grid.nelem(),
           ssd.za_grid.nelem()/2+1, 1, 16); 

  chk_size(os_ext_mat.str(), ssd.ext_mat_data,
           ssd.f_grid.nelem(), ssd.T_grid.nelem(),
           ssd.za_grid.nelem()/2+1, 1, 3);
    
  chk_size(os_abs_vec.str(), ssd.abs_vec_data,
           ssd.f_grid.nelem(), ssd.T_grid.nelem(),
           ssd.za_grid.nelem()/2+1, 1, 2);

  // Now that we are sure that za_grid is properly symmetric, we just need to
  // copy over the data (ie no interpolation).
  Tensor5 tmpT5 = ssd.abs_vec_data;
  ssd.abs_vec_data.resize(tmpT5.nshelves(),tmpT5.nbooks(),
                          ssd.za_grid.nelem(),
                          tmpT5.nrows(),tmpT5.ncols());
  ssd.abs_vec_data(joker,joker,Range(0,nza/2+1),joker,joker) = tmpT5;
  for (Index i=0; i<nza/2; i++)
    {
      ssd.abs_vec_data(joker,joker,nza-1-i,joker,joker) =
        tmpT5(joker,joker,i,joker,joker);
    }

  tmpT5 = ssd.ext_mat_data;
  ssd.ext_mat_data.resize(tmpT5.nshelves(),tmpT5.nbooks(),
                          ssd.za_grid.nelem(),
                          tmpT5.nrows(),tmpT5.ncols());
  ssd.ext_mat_data(joker,joker,Range(0,nza/2+1),joker,joker) = tmpT5;
  for (Index i=0; i<nza/2; i++)
    {
      ssd.ext_mat_data(joker,joker,nza-1-i,joker,joker) =
        tmpT5(joker,joker,i,joker,joker);
    }

  Tensor7 tmpT7 = ssd.pha_mat_data;
  ssd.pha_mat_data.resize(tmpT7.nlibraries(),tmpT7.nvitrines(),
                          tmpT7.nshelves(),tmpT7.nbooks(),
                          ssd.za_grid.nelem(),
                          tmpT7.nrows(),tmpT7.ncols());
  ssd.pha_mat_data(joker,joker,joker,joker,Range(0,nza/2+1),joker,joker) = tmpT7;

  // scatt. matrix elements 13,23,31,32 and 14,24,41,42 (=elements 2,6,8,9 and
  // 3,7,12,13 in ARTS' flattened format, respectively) change sign.
  tmpT7(joker,joker,joker,joker,joker,joker,Range(2,2)) *= -1.;
  tmpT7(joker,joker,joker,joker,joker,joker,Range(6,4)) *= -1.;
  tmpT7(joker,joker,joker,joker,joker,joker,Range(12,2)) *= -1.;

  // For second half of incident polar angles (>90deg), we need to mirror the
  // original data in both incident and scattered polar angle around 90deg "planes".
  for (Index i=0; i<nza/2; i++)
    for (Index j=0; j<nza; j++)
      ssd.pha_mat_data(joker,joker,nza-1-j,joker,nza-1-i,joker,joker) =
        tmpT7(joker,joker,j,joker,i,joker,joker);

}


//! Convert particle ssd method name to enum value
/*!
 Returns the ParticleSSDMethod enum value for the given String.

 \param[in]  particle_ssdmethod_string  Particle SSD method name
 \return     ParticleSSDMethod enum value

 \author Oliver Lemke
 */
ParticleSSDMethod ParticleSSDMethodFromString(const String& particle_ssdmethod_string)
{
    ParticleSSDMethod particle_ssdmethod;
    if (particle_ssdmethod_string == "tmatrix")
        particle_ssdmethod = PARTICLE_SSDMETHOD_TMATRIX;
    else
    {
        ostringstream os;
        os << "Unknown particle SSD method: " << particle_ssdmethod_string << endl
        << "Valid methods: tmatrix";
        throw std::runtime_error(os.str());
    }

    return particle_ssdmethod;
}


//! Convert particle type enum value to String.
/*!
 Returns the PType enum value for the given String.

 \param[in]  ptype  Particle type
 \return     String representation of ParticleSSDMethod

 \author Oliver Lemke
 */
String PTypeToString(const ParticleSSDMethod& particle_ssdmethod)
{
    String particle_ssdmethod_string;

    switch (particle_ssdmethod) {
        case PARTICLE_SSDMETHOD_TMATRIX:
            particle_ssdmethod_string = "tmatrix";
            break;
        default:
            ostringstream os;
            os << "Internal error: Cannot map ParticleSSDMethod enum value "
            << particle_ssdmethod << " to String.";
            throw std::runtime_error(os.str());
            break;
    }

    return particle_ssdmethod_string;
}

