/* Copyright (C) 2002-2012
   Sreerekha T.R. <rekha@uni-bremen.de>
   Claudia Emde <claudia.emde@dlr.de>
   Cory Davies <cory@met.ed.ac.uk>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   m_optproperties.cc
  \author Sreerekha T.R. <rekha@uni-bremen.de>, 
          Claudia Emde <claudia.emde@dlr.de>
          Cory Davies <cory@met.ed.ac.uk>
  \date   Mon Jun 10 11:19:11 2002 
  \brief  This filecontains workspace methods for calculating the optical 
  properties for the radiative transfer. 

  Optical properties are the extinction matrix, absorption vector and
  scattering vector.  The optical properties for the gases can be
  calculated with or without Zeeman effect.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts.h"
#include "exceptions.h"
#include "array.h"
#include "matpackIII.h"
#include "matpackVII.h"
#include "logic.h"
#include "interpolation.h"
#include "messages.h"
#include "xml_io.h"
#include "optproperties.h"
#include "math_funcs.h"
#include "sorting.h"
#include "check_input.h"
#include "auto_md.h" 

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;

#define PART_TYPE scat_data_array[i_pt].particle_type
#define F_DATAGRID scat_data_array[i_pt].f_grid
#define T_DATAGRID scat_data_array[i_pt].T_grid
#define ZA_DATAGRID scat_data_array[i_pt].za_grid
#define AA_DATAGRID scat_data_array[i_pt].aa_grid
#define PHA_MAT_DATA_RAW scat_data_array[i_pt].pha_mat_data  //CPD: changed from pha_mat_data
#define EXT_MAT_DATA_RAW scat_data_array[i_pt].ext_mat_data  //which wouldn't let me play with
#define ABS_VEC_DATA_RAW scat_data_array[i_pt].abs_vec_data  //scat_data_array_mono.
#define PND_LIMIT 1e-12 // If particle number density is below this value, 
                        // no transformations will be performed


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromData( // Output:
                         Tensor5& pha_mat_spt,
                         // Input:
                         const ArrayOfSingleScatteringData& scat_data_array,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index,
                         const Index& f_index,
                         const Vector& f_grid,
                         const Numeric& rtp_temperature,
                         const Tensor4& pnd_field, 
                         const Index& scat_p_index,
                         const Index& scat_lat_index,
                         const Index& scat_lon_index,
                         const Verbosity& verbosity
                         )
{
  CREATE_OUT3;
  
  out3 << "Calculate *pha_mat_spt* from database\n";

  const Index N_pt = scat_data_array.nelem();
  const Index stokes_dim = pha_mat_spt.ncols();

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  assert( pha_mat_spt.nshelves() == N_pt );
  
  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor5 pha_mat_data_int;
  

  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // If the particle number density at a specific point in the atmosphere for 
      // the i_pt particle type is zero, we don't need to do the transfromation!
      if (pnd_field(i_pt, scat_p_index, scat_lat_index, scat_lon_index) > PND_LIMIT)
        {

          // First we have to transform the data from the coordinate system 
          // used in the database (depending on the kind of particle type 
          // specified by *particle_type*) to the laboratory coordinate sytem.
      
          // Frequency interpolation:
     
          // The data is interpolated on one frequency. 
          pha_mat_data_int.resize(PHA_MAT_DATA_RAW.nshelves(), PHA_MAT_DATA_RAW.nbooks(),
                                  PHA_MAT_DATA_RAW.npages(), PHA_MAT_DATA_RAW.nrows(), 
                                  PHA_MAT_DATA_RAW.ncols());

      
          // Gridpositions:
          GridPos freq_gp;
          gridpos(freq_gp, F_DATAGRID, f_grid[f_index]);

          GridPos t_gp;
          gridpos(t_gp, T_DATAGRID, rtp_temperature);

          // Interpolationweights:
          Vector itw(4);
          interpweights(itw, freq_gp, t_gp);
     
          for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA_RAW.nshelves() ; i_za_sca++)
            {
              for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA_RAW.nbooks(); i_aa_sca++)
                {
                  for (Index i_za_inc = 0; i_za_inc < PHA_MAT_DATA_RAW.npages(); 
                       i_za_inc++)
                    {
                      for (Index i_aa_inc = 0; i_aa_inc < PHA_MAT_DATA_RAW.nrows(); 
                           i_aa_inc++)
                        {  
                          for (Index i = 0; i < PHA_MAT_DATA_RAW.ncols(); i++)
                            {
                              pha_mat_data_int(i_za_sca, 
                                               i_aa_sca, i_za_inc, 
                                               i_aa_inc, i) =
                                interp(itw,
                                       PHA_MAT_DATA_RAW(joker, joker, i_za_sca, 
                                                        i_aa_sca, i_za_inc, 
                                                        i_aa_inc, i),
                                       freq_gp, t_gp);
                            }
                        }
                    }
                }
            }
      
          // Do the transformation into the laboratory coordinate system.
          for (Index za_inc_idx = 0; za_inc_idx < scat_za_grid.nelem(); 
               za_inc_idx ++)
            {
              for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem(); 
                   aa_inc_idx ++) 
                {
                  pha_matTransform(pha_mat_spt
                                   (i_pt, za_inc_idx, aa_inc_idx, joker, joker),
                                   pha_mat_data_int, 
                                   ZA_DATAGRID, AA_DATAGRID,
                                   PART_TYPE, scat_za_index, scat_aa_index,
                                   za_inc_idx, 
                                   aa_inc_idx, scat_za_grid, scat_aa_grid,
                                   verbosity);
                }
            }
        }
    }
}
  

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromDataDOITOpt(// Output:
                                Tensor5& pha_mat_spt,
                                // Input:
                                const ArrayOfTensor7& pha_mat_sptDOITOpt,
                                const ArrayOfSingleScatteringData& scat_data_array_mono,
                                const Index& doit_za_grid_size,
                                const Vector& scat_aa_grid,
                                const Index& scat_za_index, // propagation directions
                                const Index& scat_aa_index,
                                const Numeric& rtp_temperature,
                                const Tensor4&  pnd_field, 
                                const Index& scat_p_index,
                                const Index&  scat_lat_index,
                                const Index& scat_lon_index,
                                const Verbosity&)
{
  // atmosphere_dim = 3
  if (pnd_field.ncols() > 1)
    {
      assert(pha_mat_sptDOITOpt.nelem() == scat_data_array_mono.nelem());
      // I assume that if the size is o.k. for one particle type is will 
      // also be o.k. for more particle types. 
      assert(pha_mat_sptDOITOpt[0].nlibraries() == scat_data_array_mono[0].T_grid.nelem());
      assert(pha_mat_sptDOITOpt[0].nvitrines() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].nshelves() == scat_aa_grid.nelem() );
      assert(pha_mat_sptDOITOpt[0].nbooks() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].npages() == scat_aa_grid.nelem()); 
    }
  
  // atmosphere_dim = 1, only zenith angle required for scattered directions. 
  else if ( pnd_field.ncols() == 1 )
    {
      //assert(is_size(scat_theta, doit_za_grid_size, 1,
      //                doit_za_grid_size, scat_aa_grid.nelem()));
      
      assert(pha_mat_sptDOITOpt.nelem() == scat_data_array_mono.nelem());
      // I assume that if the size is o.k. for one particle type is will 
      // also be o.k. for more particle types. 
      assert(pha_mat_sptDOITOpt[0].nlibraries() == scat_data_array_mono[0].T_grid.nelem());
      assert(pha_mat_sptDOITOpt[0].nvitrines() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].nshelves() == 1);
      assert(pha_mat_sptDOITOpt[0].nbooks() == doit_za_grid_size);
      assert(pha_mat_sptDOITOpt[0].npages() == scat_aa_grid.nelem()); 
    }
  
  assert(doit_za_grid_size > 0);
  
  // Create equidistant zenith angle grid
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);  
  
  const Index N_pt = scat_data_array_mono.nelem();
  const Index stokes_dim = pha_mat_spt.ncols();
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  assert( pha_mat_spt.nshelves() == N_pt );

  GridPos T_gp;
  Vector itw(2);
  
  // Initialisation
  pha_mat_spt = 0.;

  // Do the transformation into the laboratory coordinate system.
  for (Index i_pt = 0; i_pt < N_pt; i_pt ++)
    {
      // If the particle number density at a specific point in the atmosphere  
      // for the i_pt particle type is zero, we don't need to do the 
      // transfromation!
      if (pnd_field(i_pt, scat_p_index, scat_lat_index, scat_lon_index) > PND_LIMIT) //TRS
        {
          if( scat_data_array_mono[i_pt].T_grid.nelem() > 1)
            {
              ostringstream os;
              os << "The temperature grid of the scattering data does not cover the \n"
                  "atmospheric temperature at cloud location. The data should \n"
                  "include the value T="<< rtp_temperature << " K. \n";
              chk_interpolation_grids(os.str(),	scat_data_array_mono[i_pt].T_grid, rtp_temperature);
              
              // Gridpositions:
              gridpos(T_gp, scat_data_array_mono[i_pt].T_grid, rtp_temperature); 
              // Interpolationweights:
              interpweights(itw, T_gp);
            }
          
          
          
          for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
               za_inc_idx ++)
            {
              for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
                   aa_inc_idx ++) 
                {
                  if( scat_data_array_mono[i_pt].T_grid.nelem() == 1)
                    {
                      pha_mat_spt(i_pt, za_inc_idx, aa_inc_idx, joker, joker) =
                        pha_mat_sptDOITOpt[i_pt](0, scat_za_index,
                                                 scat_aa_index, za_inc_idx, 
                                                 aa_inc_idx, joker, joker);
                    }
                  
                  // Temperature interpolation
                  else
                    {
                      for (Index i = 0; i< stokes_dim; i++)
                        {
                          for (Index j = 0; j< stokes_dim; j++)
                            {
                              pha_mat_spt(i_pt, za_inc_idx, aa_inc_idx, i, j)=
                                interp(itw,pha_mat_sptDOITOpt[i_pt]
                                       (joker, scat_za_index,
                                        scat_aa_index, za_inc_idx, 
                                        aa_inc_idx, i, j) , T_gp);
                            }
                        }
                    }
                }
            }
        }// TRS
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromData(// Output and Input:
                          Tensor3& ext_mat_spt,
                          Matrix& abs_vec_spt,
                          // Input:
                          const ArrayOfSingleScatteringData& scat_data_array,
                          const Vector& scat_za_grid,
                          const Vector& scat_aa_grid,
                          const Index& scat_za_index, // propagation directions
                          const Index& scat_aa_index,
                          const Index& f_index,
                          const Vector& f_grid,
                          const Numeric& rtp_temperature, 
                          const Tensor4& pnd_field, 
                          const Index& scat_p_index,
                          const Index& scat_lat_index,
                          const Index& scat_lon_index,
                          const Verbosity& verbosity)
{
  
  const Index N_pt = scat_data_array.nelem();
  const Index stokes_dim = ext_mat_spt.ncols();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  assert( ext_mat_spt.npages() == N_pt );
  assert( abs_vec_spt.nrows() == N_pt );

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor3 ext_mat_data_int;
  Tensor3 abs_vec_data_int;
  
  // Initialisation
  ext_mat_spt = 0.;
  abs_vec_spt = 0.;


  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // If the particle number density at a specific point in the atmosphere for 
      // the i_pt particle type is zero, we don't need to do the transfromation

      if (pnd_field(i_pt, scat_p_index, scat_lat_index, scat_lon_index) > PND_LIMIT)
        {
          // First we have to transform the data from the coordinate system 
          // used in the database (depending on the kind of particle type 
          // specified by *particle_type*) to the laboratory coordinate sytem.
      
          // Frequency interpolation:
     
          // The data is interpolated on one frequency. 
          //
          // Resize the variables for the interpolated data:
          //
          ext_mat_data_int.resize(EXT_MAT_DATA_RAW.npages(),
                                  EXT_MAT_DATA_RAW.nrows(), 
                                  EXT_MAT_DATA_RAW.ncols());
          //
          abs_vec_data_int.resize(ABS_VEC_DATA_RAW.npages(),
                                  ABS_VEC_DATA_RAW.nrows(), 
                                  ABS_VEC_DATA_RAW.ncols());
      
      
          // Gridpositions:
          GridPos freq_gp;
          gridpos(freq_gp, F_DATAGRID, f_grid[f_index]); 
          GridPos t_gp;
          Vector itw;
          
          if ( T_DATAGRID.nelem() > 1)
            {
              ostringstream os;
              os << "The temperature grid of the scattering data does not cover the \n"
                    "atmospheric temperature at cloud location. The data should \n"
                    "include the value T="<< rtp_temperature << " K. \n";
              chk_interpolation_grids(os.str(),	T_DATAGRID, rtp_temperature);
              
              gridpos(t_gp, T_DATAGRID, rtp_temperature);
          
              // Interpolationweights:
              itw.resize(4);
              interpweights(itw, freq_gp, t_gp);
          
              for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA_RAW.npages();
                   i_za_sca++)
                {
                  for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA_RAW.nrows(); 
                      i_aa_sca++)
                    {
                      //
                      // Interpolation of extinction matrix:
                      //
                      for (Index i = 0; i < EXT_MAT_DATA_RAW.ncols(); i++)
                        {
                          ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                            interp(itw, EXT_MAT_DATA_RAW(joker, joker, 
                                                         i_za_sca, i_aa_sca, i),
                                   freq_gp, t_gp);
                        }
                    }
                }

              for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA_RAW.npages();
                   i_za_sca++)
                {
                  for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA_RAW.nrows(); 
                      i_aa_sca++)
                    {
                      //
                      // Interpolation of absorption vector:
                      //
                      for (Index i = 0; i < ABS_VEC_DATA_RAW.ncols(); i++)
                        {
                          abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                            interp(itw, ABS_VEC_DATA_RAW(joker, joker, i_za_sca, 
                                                         i_aa_sca, i),
                                   freq_gp, t_gp);
                        }
                    }
                }
            }
          else 
            {
              // Interpolationweights:
              itw.resize(2);
              interpweights(itw, freq_gp);
          
              for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA_RAW.npages();
                   i_za_sca++)
                {
                  for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA_RAW.nrows(); 
                      i_aa_sca++)
                    {
                      //
                      // Interpolation of extinction matrix:
                      //
                      for (Index i = 0; i < EXT_MAT_DATA_RAW.ncols(); i++)
                        {
                          ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                            interp(itw, EXT_MAT_DATA_RAW(joker, 0, 
                                                         i_za_sca, i_aa_sca, i),
                                   freq_gp);
                        }
                    }
                }

              for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA_RAW.npages();
                   i_za_sca++)
                {
                  for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA_RAW.nrows(); 
                      i_aa_sca++)
                    {
                      //
                      // Interpolation of absorption vector:
                      //
                      for (Index i = 0; i < ABS_VEC_DATA_RAW.ncols(); i++)
                        {
                          abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                            interp(itw, ABS_VEC_DATA_RAW(joker, 0, i_za_sca, 
                                                         i_aa_sca, i),
                                   freq_gp);
                        }
                    }
                }
            } 
      

          //
          // Do the transformation into the laboratory coordinate system.
          //
          // Extinction matrix:
          //
     
  
          ext_matTransform(ext_mat_spt(i_pt, joker, joker),
                           ext_mat_data_int,
                           ZA_DATAGRID, AA_DATAGRID, PART_TYPE,
                           za_sca, aa_sca,
                           verbosity);
          // 
          // Absorption vector:
          //
          abs_vecTransform(abs_vec_spt(i_pt, joker),
                           abs_vec_data_int,
                           ZA_DATAGRID, AA_DATAGRID, PART_TYPE,
                           za_sca, aa_sca, verbosity);                
        }

    }
}
                          

/* Workspace method: Doxygen documentation will be auto-generated */
void ext_matAddPart(Tensor3& ext_mat,
                    const Tensor3& ext_mat_spt,
                    const Tensor4& pnd_field,
                    const Index& atmosphere_dim,
                    const Index& scat_p_index,
                    const Index& scat_lat_index,
                    const Index& scat_lon_index,
                    const Verbosity&)
                     
{
  Index N_pt = ext_mat_spt.npages();
  Index stokes_dim = ext_mat_spt.nrows();
  
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);

  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error(
                        "The dimension of stokes vector can be "
                        "only 1,2,3, or 4");
  }
  if ( ext_mat_spt.ncols() != stokes_dim){
    
    throw runtime_error(" The columns of ext_mat_spt should "
                        "agree to stokes_dim");
  }

  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
        { 
          
          // now the last two loops over the stokes dimension.
          for (Index m = 0; m < stokes_dim; m++)
            {
              for (Index n = 0; n < stokes_dim; n++)
                //summation of the product of pnd_field and 
                //ext_mat_spt.
                ext_mat_part(m, n) += 
                  (ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
            }
        }

      //Add particle extinction matrix to *ext_mat*.
      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;
    }
 
  if (atmosphere_dim == 3)
    {
      
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
        { 
          
          // now the last two loops over the stokes dimension.
          for (Index m = 0; m < stokes_dim; m++)
            {
              for (Index n = 0; n < stokes_dim; n++)
                //summation of the product of pnd_field and 
                //ext_mat_spt.
                ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
                                        pnd_field(l, scat_p_index, 
                                                  scat_lat_index, 
                                                  scat_lon_index));
              
            } 
        }

      //Add particle extinction matrix to *ext_mat*.
      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;

    }

} 


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_vecAddPart(Matrix& abs_vec,
                    const Matrix& abs_vec_spt,
                    const Tensor4& pnd_field,
                    const Index& atmosphere_dim,
                    const Index& scat_p_index,
                    const Index& scat_lat_index,
                    const Index& scat_lon_index,
                    const Verbosity&)
                    
{
  Index N_pt = abs_vec_spt.nrows();
  Index stokes_dim = abs_vec_spt.ncols();

  Vector abs_vec_part(stokes_dim, 0.0);

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  } 
 
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
        {
          // now the loop over the stokes dimension.
          //(CE:) in the middle was l instead of m
          for (Index m = 0; m < stokes_dim; ++m)
             //summation of the product of pnd_field and 
            //abs_vec_spt.
            abs_vec_part[m] += 
              (abs_vec_spt(l, m) * pnd_field(l, scat_p_index, 0, 0));
          
        }
      //Add the particle absorption
      abs_vec(0, Range(joker)) += abs_vec_part;
    }
  
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
        {
          
          // now the loop over the stokes dimension.
          for (Index m = 0; m < stokes_dim; ++m)
             //summation of the product of pnd_field and 
            //abs_vec_spt.
            abs_vec_part[m] += (abs_vec_spt(l, m) *
                                pnd_field(l, scat_p_index,
                                          scat_lat_index, 
                                          scat_lon_index));
          
        }
      //Add the particle absorption
      abs_vec(0,Range(joker)) += abs_vec_part;
    }
} 


/* Workspace method: Doxygen documentation will be auto-generated */
void ext_matInit(Tensor3&         ext_mat,
                 const Vector&    f_grid,
                 const Index&     stokes_dim,
                 const Index&     f_index,
                 const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  Index freq_dim;

  if( f_index < 0 )
    freq_dim = f_grid.nelem();
  else
    freq_dim = 1;
 
  ext_mat.resize( freq_dim,
                  stokes_dim,
                  stokes_dim );
  ext_mat = 0;                  // Initialize to zero!

  out2 << "Set dimensions of ext_mat as ["
       << freq_dim << ","
       << stokes_dim << ","
       << stokes_dim << "] and initialized to 0.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ext_matAddGas(Tensor3&      ext_mat,
                   const Tensor4& propmat_clearsky,
                   const Verbosity&)
{
  // Number of Stokes parameters:
  const Index stokes_dim = ext_mat.ncols();

  // The second dimension of ext_mat must also match the number of
  // Stokes parameters:
  if ( stokes_dim != ext_mat.nrows() )
    throw runtime_error("Row dimension of ext_mat inconsistent with "
                        "column dimension.");
  if ( stokes_dim != propmat_clearsky.ncols() )
    throw runtime_error("Col dimension of propmat_clearsky "
                        "inconsistent with col dimension in ext_mat.");

  // Number of frequencies:
  const Index f_dim = ext_mat.npages();

  // This must be consistent with the second dimension of
  // propmat_clearsky. Check this:
  if ( f_dim != propmat_clearsky.npages() )
    throw runtime_error("Frequency dimension of ext_mat and propmat_clearsky\n"
                        "are inconsistent in ext_matAddGas.");

  // Sum up absorption over all species.
  // This gives us an absorption vector for all frequencies. Of course
  // this includes the special case that there is only one frequency.
  Tensor3 abs_total(f_dim,stokes_dim,stokes_dim);
  abs_total = 0;
  
//   for ( Index i=0; i<f_dim; ++i )
//     abs_total[i] = abs_scalar_gas(i,joker).sum();
  for ( Index iv=0; iv<f_dim; ++iv )
        for ( Index is1=0; is1<stokes_dim; ++is1 )
              for ( Index is2=0; is2<stokes_dim; ++is2 )
                    abs_total(iv,is1,is2) += propmat_clearsky(joker,iv,is1,is2).sum();
  
    // Add the absorption value to all the elements:
      ext_mat += abs_total;
      
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_vecInit(Matrix&       abs_vec,
                 const Vector& f_grid,
                 const Index&  stokes_dim,
                 const Index&  f_index,
                 const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  Index freq_dim;

  if( f_index < 0 )
    freq_dim = f_grid.nelem();
  else
    freq_dim = 1;
 
  abs_vec.resize( freq_dim,
                  stokes_dim );
  abs_vec = 0;                  // Initialize to zero!

  out2 << "Set dimensions of abs_vec as ["
       << freq_dim << ","
       << stokes_dim << "] and initialized to 0.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_vecAddGas(Matrix&       abs_vec,
                   const Tensor4& propmat_clearsky,
                   const Verbosity&)
{
  // Number of frequencies:
  const Index f_dim = abs_vec.nrows();
  const Index stokes_dim = abs_vec.ncols();
  
  // This must be consistent with the second dimension of
  // propmat_clearsky. Check this:
  if ( f_dim != propmat_clearsky.npages() )
    throw runtime_error("Frequency dimension of abs_vec and propmat_clearsky\n"
                        "are inconsistent in abs_vecAddGas.");
  if ( stokes_dim != propmat_clearsky.ncols() )
    throw runtime_error("Stokes dimension of abs_vec and propmat_clearsky\n"
                        "are inconsistent in abs_vecAddGas.");
    
  // Loop all frequencies. Of course this includes the special case
  // that there is only one frequency.
  for ( Index i=0; i<f_dim; ++i )
    {
      // Sum up the columns of propmat_clearsky and add to the first
      // element of abs_vec.
      for(Index is = 0; is < stokes_dim;is++)
        abs_vec(i,is) += propmat_clearsky(joker,i,is,0).sum();
    }

  // We don't have to do anything about higher elements of abs_vec,
  // since scalar gas absorption only influences the first element.
}


/* Workspace method: Doxygen documentation will be auto-generated */
/*
void ext_matAddGasZeeman( Tensor3&      ext_mat,
                          const Tensor3&  ext_mat_zee,
                          const Verbosity&)
{
  // Number of Stokes parameters:
  const Index stokes_dim = ext_mat.ncols();

  // The second dimension of ext_mat must also match the number of
  // Stokes parameters:
  if ( stokes_dim != ext_mat.nrows() )
    throw runtime_error("Row dimension of ext_mat inconsistent with "
                        "column dimension."); 

  for ( Index i=0; i<stokes_dim; ++i )
    {
      for ( Index j=0; j<stokes_dim; ++j )
        {
          // Add the zeeman extinction to extinction matrix.
          ext_mat(joker,i,j) += ext_mat_zee(joker, i, j);
        }
      
    }
}
*/


/* Workspace method: Doxygen documentation will be auto-generated */
/*
void abs_vecAddGasZeeman( Matrix&      abs_vec,
                          const Matrix& abs_vec_zee,
                          const Verbosity&)
{
  // Number of Stokes parameters:
  const Index stokes_dim = abs_vec_zee.ncols();
  // that there is only one frequency.
  for ( Index j=0; j<stokes_dim; ++j )
    {
      abs_vec(joker,j) += abs_vec_zee(joker,j);
    }
}
*/


/* Workspace method: Doxygen documentation will be auto-generated */
void pha_matCalc(Tensor4& pha_mat,
                 const Tensor5& pha_mat_spt,
                 const Tensor4& pnd_field,
                 const Index& atmosphere_dim,
                 const Index& scat_p_index,
                 const Index& scat_lat_index,
                 const Index& scat_lon_index,
                 const Verbosity&)
{

  Index N_pt = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();
  Index stokes_dim = pha_mat_spt.nrows();
 
  pha_mat.resize(Nza, Naa, stokes_dim, stokes_dim);

  // Initialisation
  pha_mat = 0.0;
          
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index pt_index = 0; pt_index < N_pt; ++ pt_index)
        {
          // these are loops over zenith angle and azimuth angle
          for (Index za_index = 0; za_index < Nza; ++ za_index)
            {
              for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
                {
                  
                  // now the last two loops over the stokes dimension.
                  for (Index stokes_index_1 = 0; stokes_index_1 < stokes_dim; 
                       ++  stokes_index_1)
                    {
                      for (Index stokes_index_2 = 0; stokes_index_2 < stokes_dim;
                           ++ stokes_index_2)
                         //summation of the product of pnd_field and 
                          //pha_mat_spt.
                        pha_mat(za_index, aa_index,  
                                     stokes_index_1, stokes_index_2) += 
                          
                          (pha_mat_spt(pt_index, za_index, aa_index,  
                                       stokes_index_1, stokes_index_2) * 
                           pnd_field(pt_index,scat_p_index, 0, 0));
                    }
                }
            }
        }
    }
          
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index pt_index = 0; pt_index < N_pt; ++ pt_index)
        {
          
          // these are loops over zenith angle and azimuth angle
          for (Index za_index = 0; za_index < Nza; ++ za_index)
            {
              for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
                {
                  
                  // now the last two loops over the stokes dimension.
                  for (Index i = 0;  i < stokes_dim; ++  i)
                    {
                      for (Index j = 0; j < stokes_dim; ++ j)
                        {
                          //summation of the product of pnd_field and 
                          //pha_mat_spt.
                          pha_mat(za_index, aa_index, i,j ) += 
                            (pha_mat_spt(pt_index, za_index, aa_index, i, j) * 
                             pnd_field(pt_index, scat_p_index,  
                                       scat_lat_index, scat_lon_index));
                          
                          
                        } 
                    }	
                }
            }		
        }	
    }		
}



/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_arrayCheck(//Input:
                        const ArrayOfSingleScatteringData& scat_data_array,
                        const Numeric& threshold,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;

/* JM121024: we do not really need to write the scatt_data to file again. we
             usually have just read them in a couple of commands before!?
             if wanted/needed for debug cases, just uncomment the 2 lines below */
//  xml_write_to_file("SingleScatteringData", scat_data_array, FILE_TYPE_ASCII,
//                    verbosity);
  
  const Index N_pt = scat_data_array.nelem();
  
  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      out2 << "  particle " << i_pt << "\n";

      switch (PART_TYPE){

        case PARTICLE_TYPE_MACROS_ISO:
          {
            for (Index f = 0; f < F_DATAGRID.nelem(); f++)
              {
                out2 << "frequency " << F_DATAGRID[f] << "Hz\n";
                for (Index t = 0; t < T_DATAGRID.nelem(); t++)
                  {
                    out2 << "Temperature " << T_DATAGRID[t] << "K\n";

                    Numeric Csca = AngIntegrate_trapezoid
                      (PHA_MAT_DATA_RAW(f, t, joker, 0, 0, 0, 0), ZA_DATAGRID);

                    Numeric Cext_data = EXT_MAT_DATA_RAW(f,t,0,0,0);
 
                    Numeric Cabs = Cext_data - Csca;

                    Numeric Cabs_data = ABS_VEC_DATA_RAW(f,t,0,0,0);

                    Numeric Csca_data = Cext_data - Cabs_data;


                    out2 << "  Coefficients in database: "
                         << "Cext: " << Cext_data << " Cabs: " << Cabs_data
                         << " Csca: " << Csca_data << "\n"
                         << "  Calculated coefficients: "
                         << "Cabs calc: " << Cabs   
                         << " Csca calc: " << Csca << "\n"
                         << "  Deviations "
                         << "Cabs: " << 1e2*Cabs/Cabs_data-1e2
                         << "% Csca: " << 1e2*Csca/Csca_data-1e2
                         << "% Alb: " << (Csca-Csca_data)/Cext_data << "\n";


//                    if (abs(Csca/Csca_data-1.)*Csca_data/Cext_data > threshold)
//                  equivalent to the above (it's actually the (absolute) albedo
//                  deviation!)
                    if (abs(Csca-Csca_data)/Cext_data > threshold)
                      {
                        ostringstream os;
                        os << "  Deviations in scat_data_array too large:\n"
                           << "  scat dev [%] " << 1e2*Csca/Csca_data-1e2
                           << " at albedo of " << Csca_data/Cext_data << "\n"
                           << "  Check entry for particle " << i_pt << " at "
                           << f << ".frequency and " << t << ".temperature!\n";
                        throw runtime_error( os.str() );
                      }
                  }
               }
            break;
          }

        default:
          {
            CREATE_OUT0;
            out0 << "  WARNING:\n"
                 << "  scat_data_array consistency check not implemented (yet?!) for\n"
                 << "  particle type " << PART_TYPE << "!\n";
          }
      }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void DoitScatteringDataPrepare(//Output:
                               ArrayOfTensor7& pha_mat_sptDOITOpt,
                               ArrayOfSingleScatteringData& scat_data_array_mono,
                               //Input:
                               const Index& doit_za_grid_size,
                               const Vector& scat_aa_grid,
                               const ArrayOfSingleScatteringData&
                               scat_data_array,
                               const Vector& f_grid,
                               const Index& f_index,
                               const Index& atmosphere_dim,
                               const Index& stokes_dim,
                               const Verbosity& verbosity)
{
  
  // Interpolate all the data in frequency
  scat_data_array_monoCalc(scat_data_array_mono, scat_data_array, f_grid, f_index, verbosity);
  
  // For 1D calculation the scat_aa dimension is not required:
  Index N_aa_sca;
  if(atmosphere_dim == 1)
    N_aa_sca = 1;
  else
    N_aa_sca = scat_aa_grid.nelem();
  
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);

  assert( scat_data_array.nelem() == scat_data_array_mono.nelem() );
  
  Index N_pt = scat_data_array.nelem();  
  
  pha_mat_sptDOITOpt.resize(N_pt);

  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      Index N_T = scat_data_array_mono[i_pt].T_grid.nelem();
      pha_mat_sptDOITOpt[i_pt].resize(N_T, doit_za_grid_size, N_aa_sca, 
                                      doit_za_grid_size, scat_aa_grid.nelem(), 
                                      stokes_dim, stokes_dim);
      
      //    Initialize:
      pha_mat_sptDOITOpt[i_pt]= 0.;
   
      // Calculate all scattering angles for all combinations of incoming 
      // and scattered directions and interpolation.
      for (Index t_idx = 0; t_idx < N_T; t_idx ++)
        {
          // These are the scattered directions as called in *scat_field_calc*
          for (Index za_sca_idx = 0; za_sca_idx < doit_za_grid_size; za_sca_idx ++)
            {
              for (Index aa_sca_idx = 0; aa_sca_idx < N_aa_sca; aa_sca_idx ++)
                {
                  // Integration is performed over all incoming directions
                  for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
                       za_inc_idx ++)
                    {
                      for (Index aa_inc_idx = 0; aa_inc_idx <
                             scat_aa_grid.nelem();
                           aa_inc_idx ++)
                        {
                          pha_matTransform(pha_mat_sptDOITOpt[i_pt]
                                           (t_idx, za_sca_idx,            
                                            aa_sca_idx, za_inc_idx, aa_inc_idx,
                                            joker, joker),
                                           scat_data_array_mono[i_pt].
                                           pha_mat_data
                                           (0,t_idx,joker,joker,joker,
                                            joker,joker),
                                           scat_data_array_mono[i_pt].za_grid,
                                           scat_data_array_mono[i_pt].aa_grid,
                                           scat_data_array_mono[i_pt].particle_type,
                                           za_sca_idx,
                                           aa_sca_idx,
                                           za_inc_idx,
                                           aa_inc_idx,
                                           za_grid,
                                           scat_aa_grid,
                                           verbosity);
                        }
                    }
                }
            }
        }
    }
 } 


/* Workspace method: Doxygen documentation will be auto-generated */
void scat_data_array_monoCalc(ArrayOfSingleScatteringData& scat_data_array_mono,
                        const ArrayOfSingleScatteringData& scat_data_array,
                        const Vector& f_grid,
                        const Index& f_index,
                        const Verbosity&)
{
  //Extrapolation factor:
  //const Numeric extpolfac = 0.5;
  
  // Check, whether single scattering data contains the right frequencies:
  // The check was changed to allow extrapolation at the boundaries of the 
  // frequency grid.
  for (Index i = 0; i<scat_data_array.nelem(); i++)
    {
      // check with extrapolation
      chk_interpolation_grids("scat_data_array.f_grid to f_grid",
			      scat_data_array[i].f_grid,
			      f_grid[f_index]);
      
      // old check without extrapolation
      /*if (scat_data_array[i].f_grid[0] > f_grid[f_index] || 
          scat_data_array[i].f_grid[scat_data_array[i].f_grid.nelem()-1] < f_grid[f_index])
        {
          ostringstream os;
          os << "Frequency of the scattering calculation " << f_grid[f_index] 
             << " GHz is not contained \nin the frequency grid of the " << i+1 
             << "the single scattering data file \n(*ParticleTypeAdd*). "
             << "Range:"  << scat_data_array[i].f_grid[0]/1e9 <<" - "
             << scat_data_array[i].f_grid[scat_data_array[i].f_grid.nelem()-1]/1e9
             <<" GHz \n";
          throw runtime_error( os.str() );
        }*/
    }

  const Index N_pt = scat_data_array.nelem();
  
  //Initialise scat_data_array_mono
  scat_data_array_mono.resize(N_pt);
  
  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // Gridpositions:
      GridPos freq_gp;
      gridpos(freq_gp, F_DATAGRID, f_grid[f_index]); 
      
      // Interpolationweights:
      Vector itw(2);
      interpweights(itw, freq_gp);

      //Stuff that doesn't need interpolating
      scat_data_array_mono[i_pt].particle_type=PART_TYPE;
      scat_data_array_mono[i_pt].f_grid.resize(1);
      scat_data_array_mono[i_pt].f_grid=f_grid[f_index];
      scat_data_array_mono[i_pt].T_grid=scat_data_array[i_pt].T_grid;
      scat_data_array_mono[i_pt].za_grid=ZA_DATAGRID;
      scat_data_array_mono[i_pt].aa_grid=AA_DATAGRID;
          
      //Phase matrix data
      scat_data_array_mono[i_pt].pha_mat_data.resize(1,
                                               PHA_MAT_DATA_RAW.nvitrines(),
                                               PHA_MAT_DATA_RAW.nshelves(),
                                               PHA_MAT_DATA_RAW.nbooks(),
                                               PHA_MAT_DATA_RAW.npages(),
                                               PHA_MAT_DATA_RAW.nrows(),
                                               PHA_MAT_DATA_RAW.ncols());
      
      for (Index t_index = 0; t_index < PHA_MAT_DATA_RAW.nvitrines(); t_index ++)
        {
          for (Index i_za_sca = 0; i_za_sca < PHA_MAT_DATA_RAW.nshelves();
               i_za_sca++)
            {
              for (Index i_aa_sca = 0; i_aa_sca < PHA_MAT_DATA_RAW.nbooks();
                   i_aa_sca++)
                {
                  for (Index i_za_inc = 0; i_za_inc < 
                         PHA_MAT_DATA_RAW.npages(); 
                       i_za_inc++)
                    {
                      for (Index i_aa_inc = 0; 
                           i_aa_inc < PHA_MAT_DATA_RAW.nrows(); 
                           i_aa_inc++)
                        {  
                          for (Index i = 0; i < PHA_MAT_DATA_RAW.ncols(); i++)
                            {
                              scat_data_array_mono[i_pt].pha_mat_data(0, t_index, 
                                                                i_za_sca, 
                                                                i_aa_sca,
                                                                i_za_inc, 
                                                                i_aa_inc, i) =
                                interp(itw,
                                       PHA_MAT_DATA_RAW(joker, t_index,
                                                        i_za_sca, 
                                                        i_aa_sca, i_za_inc, 
                                                        i_aa_inc, i),
                                       freq_gp);
                            }
                        }
                    }
                }
            }
          //Extinction matrix data
          scat_data_array_mono[i_pt].ext_mat_data.resize(1, T_DATAGRID.nelem(), 
                                                   EXT_MAT_DATA_RAW.npages(),
                                                   EXT_MAT_DATA_RAW.nrows(),
                                                   EXT_MAT_DATA_RAW.ncols());
          for (Index i_za_sca = 0; i_za_sca < EXT_MAT_DATA_RAW.npages();
               i_za_sca++)
            {
              for(Index i_aa_sca = 0; i_aa_sca < EXT_MAT_DATA_RAW.nrows();
                  i_aa_sca++)
                {
                  //
                  // Interpolation of extinction matrix:
                  //
                  for (Index i = 0; i < EXT_MAT_DATA_RAW.ncols(); i++)
                    {
                      scat_data_array_mono[i_pt].ext_mat_data(0, t_index, 
                                                        i_za_sca, i_aa_sca, i)
                        = interp(itw, EXT_MAT_DATA_RAW(joker, t_index, i_za_sca,
                                                       i_aa_sca, i),
                                 freq_gp);
                    }
                }
            }
          //Absorption vector data
          scat_data_array_mono[i_pt].abs_vec_data.resize(1, T_DATAGRID.nelem(),
                                                   ABS_VEC_DATA_RAW.npages(),
                                                   ABS_VEC_DATA_RAW.nrows(),
                                                   ABS_VEC_DATA_RAW.ncols());
          for (Index i_za_sca = 0; i_za_sca < ABS_VEC_DATA_RAW.npages() ;
               i_za_sca++)
            {
              for(Index i_aa_sca = 0; i_aa_sca < ABS_VEC_DATA_RAW.nrows();
                  i_aa_sca++)
                {
                  //
                  // Interpolation of absorption vector:
                  //
                  for (Index i = 0; i < ABS_VEC_DATA_RAW.ncols(); i++)
                    {
                      scat_data_array_mono[i_pt].abs_vec_data(0, t_index, i_za_sca,
                                                        i_aa_sca, i) =
                        interp(itw, ABS_VEC_DATA_RAW(joker, t_index, i_za_sca,
                                                     i_aa_sca, i),
                               freq_gp);
                    }
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void opt_prop_sptFromMonoData(// Output and Input:
                              Tensor3& ext_mat_spt,
                              Matrix& abs_vec_spt,
                              // Input:
                              const ArrayOfSingleScatteringData& scat_data_array_mono,
                              const Vector& scat_za_grid,
                              const Vector& scat_aa_grid,
                              const Index& scat_za_index, // propagation directions
                              const Index& scat_aa_index,
                              const Numeric& rtp_temperature,
                              const Tensor4& pnd_field, 
                              const Index& scat_p_index,
                              const Index& scat_lat_index,
                              const Index& scat_lon_index,
                              const Verbosity& verbosity)
{
  const Index N_pt = scat_data_array_mono.nelem();
  const Index stokes_dim = ext_mat_spt.ncols();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];
  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }
  
  assert( ext_mat_spt.npages() == N_pt );
  assert( abs_vec_spt.nrows() == N_pt );

  // Initialisation
  ext_mat_spt = 0.;
  abs_vec_spt = 0.;

  GridPos t_gp;
  
  Vector itw(2);
  
  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // If the particle number density at a specific point in the atmosphere for 
      // the i_pt particle type is zero, we don't need to do the transfromation!
      if (pnd_field(i_pt, scat_p_index, scat_lat_index, scat_lon_index) > PND_LIMIT)
        {
 
          // First we have to transform the data from the coordinate system 
          // used in the database (depending on the kind of particle type 
          // specified by *particle_type*) to the laboratory coordinate sytem. 
          
          //
          // Do the transformation into the laboratory coordinate system.
          //
          // Extinction matrix:
          //
          Index ext_npages = scat_data_array_mono[i_pt].ext_mat_data.npages();  
          Index ext_nrows = scat_data_array_mono[i_pt].ext_mat_data.nrows();  
          Index ext_ncols = scat_data_array_mono[i_pt].ext_mat_data.ncols();  
          Index abs_npages = scat_data_array_mono[i_pt].abs_vec_data.npages();  
          Index abs_nrows = scat_data_array_mono[i_pt].abs_vec_data.nrows();  
          Index abs_ncols = scat_data_array_mono[i_pt].abs_vec_data.ncols();  
          Tensor3 ext_mat_data1temp(ext_npages,ext_nrows,ext_ncols,0.0);
          Tensor3 abs_vec_data1temp(abs_npages,abs_nrows,abs_ncols,0.0);

          //Check that scattering data temperature range covers required temperature
          ConstVectorView t_grid = scat_data_array_mono[i_pt].T_grid;
      
          if (t_grid.nelem() > 1)
            {
              //   if ((rtp_temperature<t_grid[0])||(rtp_temperature>t_grid[t_grid.nelem()-1]))
              //             {
              //               throw runtime_error("rtp_temperature outside scattering data temperature range");
              //             }
          
              //interpolate over temperature
              gridpos(t_gp, scat_data_array_mono[i_pt].T_grid, rtp_temperature);
              interpweights(itw, t_gp);
              for (Index i_p = 0; i_p < ext_npages ; i_p++)
                {
                  for (Index i_r = 0; i_r < ext_nrows ; i_r++)
                    {
                      for (Index i_c = 0; i_c < ext_ncols ; i_c++)
                        {
                          ext_mat_data1temp(i_p,i_r,i_c)=interp(itw,
                                                                scat_data_array_mono[i_pt].ext_mat_data(0,joker,i_p,i_r,i_c),t_gp);
                        }
                    }
                }
            }
          else 
            {
              ext_mat_data1temp = 
                scat_data_array_mono[i_pt].ext_mat_data(0,0,joker,joker,joker);
            }
      
          ext_matTransform(ext_mat_spt(i_pt, joker, joker),
                           ext_mat_data1temp,
                           scat_data_array_mono[i_pt].za_grid, 
                           scat_data_array_mono[i_pt].aa_grid, 
                           scat_data_array_mono[i_pt].particle_type,
                           za_sca, aa_sca,
                           verbosity);
          // 
          // Absorption vector:
          //
     
          if (t_grid.nelem() > 1)
            {
              //interpolate over temperature
              for (Index i_p = 0; i_p < abs_npages ; i_p++)
                {
                  for (Index i_r = 0; i_r < abs_nrows ; i_r++)
                    {
                      for (Index i_c = 0; i_c < abs_ncols ; i_c++)
                        {
                          abs_vec_data1temp(i_p,i_r,i_c)=interp(itw,
                                                                scat_data_array_mono[i_pt].abs_vec_data(0,joker,i_p,i_r,i_c),t_gp);
                        }
                    }
                }
            }
          else
            {
              abs_vec_data1temp = 
                scat_data_array_mono[i_pt].abs_vec_data(0,0,joker,joker,joker);
            }
      
          abs_vecTransform(abs_vec_spt(i_pt, joker),
                           abs_vec_data1temp,
                           scat_data_array_mono[i_pt].za_grid, 
                           scat_data_array_mono[i_pt].aa_grid, 
                           scat_data_array_mono[i_pt].particle_type,
                           za_sca, aa_sca,
                           verbosity);                
        }

    }
}
 

/* Workspace method: Doxygen documentation will be auto-generated */
void pha_mat_sptFromMonoData(// Output:
                             Tensor5& pha_mat_spt,
                             // Input:
                             const ArrayOfSingleScatteringData& scat_data_array_mono,
                             const Index& doit_za_grid_size,
                             const Vector& scat_aa_grid,
                             const Index& scat_za_index, // propagation directions
                             const Index& scat_aa_index,
                             const Numeric& rtp_temperature,
                             const Tensor4& pnd_field, 
                             const Index& scat_p_index,
                             const Index& scat_lat_index,
                             const Index& scat_lon_index,
                             const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  out3 << "Calculate *pha_mat_spt* from scat_data_array_mono. \n";
  
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size); 

  const Index N_pt = scat_data_array_mono.nelem();
  const Index stokes_dim = pha_mat_spt.ncols();

 

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                        "must be 1,2,3 or 4");
  }
  
  assert( pha_mat_spt.nshelves() == N_pt );
  
  GridPos T_gp;
  Vector itw(2);

  // Initialisation
  pha_mat_spt = 0.;
  
  for (Index i_pt = 0; i_pt < N_pt; i_pt ++)
    { 
      // If the particle number density at a specific point in the atmosphere 
      // for the i_pt particle type is zero, we don't need to do the 
      // transfromation!
      if (pnd_field(i_pt, scat_p_index, scat_lat_index, scat_lon_index) >
          PND_LIMIT)
        { 
          // Temporary phase matrix wich icludes the all temperatures.
          Tensor3 pha_mat_spt_tmp(scat_data_array_mono[i_pt].T_grid.nelem(), 
                                  pha_mat_spt.nrows(), pha_mat_spt.ncols());
  
          pha_mat_spt_tmp = 0.; 
      
          if( scat_data_array_mono[i_pt].T_grid.nelem() > 1)
            {
              ostringstream os;
              os << "The temperature grid of the scattering data does not cover the \n"
                    "atmospheric temperature at cloud location. The data should \n"
                    "include the value T="<< rtp_temperature << " K. \n";
              chk_interpolation_grids(os.str(),	scat_data_array_mono[i_pt].T_grid, rtp_temperature);
              
              // Gridpositions:
              gridpos(T_gp, scat_data_array_mono[i_pt].T_grid, rtp_temperature); 
              // Interpolationweights:
              interpweights(itw, T_gp);
            }
      
          // Do the transformation into the laboratory coordinate system.
          for (Index za_inc_idx = 0; za_inc_idx < doit_za_grid_size;
               za_inc_idx ++)
            {
              for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
                   aa_inc_idx ++) 
                {
                  for (Index t_idx = 0; t_idx < 
                         scat_data_array_mono[i_pt].T_grid.nelem();
                       t_idx ++)
                    {
                      pha_matTransform( pha_mat_spt_tmp(t_idx, joker, joker),
                                        scat_data_array_mono[i_pt].
                                        pha_mat_data
                                        (0,0,joker,joker,joker,
                                         joker,joker),
                                        scat_data_array_mono[i_pt].za_grid, 
                                        scat_data_array_mono[i_pt].aa_grid,
                                        scat_data_array_mono[i_pt].particle_type,
                                        scat_za_index, scat_aa_index, 
                                        za_inc_idx, 
                                        aa_inc_idx, za_grid, scat_aa_grid,
                                        verbosity );
                    }
                  // Temperature interpolation
                  if( scat_data_array_mono[i_pt].T_grid.nelem() > 1)
                    {
                      for (Index i = 0; i< stokes_dim; i++)
                        {
                          for (Index j = 0; j< stokes_dim; j++)
                            {
                              pha_mat_spt(i_pt, za_inc_idx, aa_inc_idx, i, j)=
                                interp(itw, pha_mat_spt_tmp(joker, i, j), T_gp);
                            }
                        }
                    }
                  else // no temperatue interpolation required
                    {
                      pha_mat_spt(i_pt, za_inc_idx, aa_inc_idx, joker, joker) =
                        pha_mat_spt_tmp(0, joker, joker);
                    }
                }
            }
        }
    }
}

