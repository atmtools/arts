/* Copyright (C) 2003-2012 Mattias Ekström <ekstrom@rss.chalmers.se>

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
  \file   sensor.h
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2003-02-28

  \brief  Sensor modelling functions.

   This file contains the definition of the functions in sensor.cc
   that are of interest elsewhere.
*/


#ifndef sensor_h
#define sensor_h

#include "arts.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "math_funcs.h"
#include "matpackI.h"
#include "matpackII.h"


/*===========================================================================
  === Functions from sensor.cc
  ===========================================================================*/

void antenna1d_matrix(      
           Sparse&   H,
      const Index&   antenna_dim,
   ConstVectorView   antenna_dza,
    const GriddedField4&   antenna_response,
   ConstVectorView   za_grid,
   ConstVectorView   f_grid,
       const Index   n_pol,
       const Index   do_norm );

void gaussian_response_autogrid(
           Vector&   x,
           Vector&   y,
    const Numeric&   x0,
    const Numeric&   fwhm,
    const Numeric&   xwidth_si,
    const Numeric&   dx_si );

void gaussian_response(
          Vector&    y,
    const Vector&    x,
    const Numeric&   x0,
    const Numeric&   fwhm );

void mixer_matrix(
           Sparse&   H,
           Vector&   f_mixer,
    const Numeric&   lo,
    const GriddedField1&   filter,
   ConstVectorView   f_grid,
      const Index&   n_pol,
      const Index&   n_sp,
      const Index&   do_norm );

void mueller_rotation(
          Sparse&    H,
    const Index&     stokes_dim,
    const Numeric&   rotangle );

void met_mm_polarisation_hmatrix(Sparse& H,
                                 const ArrayOfString& pol,
                                 const Numeric dza,
                                 const Index stokes_dim,
                                 const String& iy_unit);

void sensor_aux_vectors(
               Vector&   sensor_response_f,
         ArrayOfIndex&   sensor_response_pol,
               Matrix&   sensor_response_dlos,
       ConstVectorView   sensor_response_f_grid,
   const ArrayOfIndex&   sensor_response_pol_grid,
       ConstMatrixView   sensor_response_dlos_grid );

void spectrometer_matrix( 
           Sparse&         H,
   ConstVectorView         ch_f,
   const ArrayOfGriddedField1&   ch_response,
   ConstVectorView         sensor_f,
      const Index&         n_pol,
      const Index&         n_sp,
      const Index&         do_norm );

void stokes2pol( 
            ArrayOfVector&  s2p,
      const Numeric&        w );

void stokes2pol(
        Vector&   w,
  const Index&    stokes_dim,
  const Index&    ipol_1based,
  const Numeric   nv = 1 );

void find_effective_channel_boundaries(// Output:
                                       Vector& fmin,
                                       Vector& fmax,
                                       // Input:
                                       const Vector& f_backend,
                                       const ArrayOfGriddedField1& backend_channel_response,
                                       const Numeric& delta,
                                       const Verbosity& verbosity);



void integration_func_by_vecmult(
        VectorView   h,
   ConstVectorView   f,
   ConstVectorView   x_f_in,
   ConstVectorView   x_g_in );

void integration_bin_by_vecmult(
        VectorView   h,
   ConstVectorView   x_g_in,
   const Numeric&    limit1, 
   const Numeric&    limit2 );

void summation_by_vecmult(
        VectorView   h,
   ConstVectorView   f,
   ConstVectorView   x_f,
   ConstVectorView   x_g,
     const Numeric   x1,
     const Numeric   x2 );



#endif  // sensor_h
