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
#include "matpack_data.h"
#include "matpack_sparse.h"
#include "messages.h"

/*===========================================================================
  === Functions from sensor.cc
  ===========================================================================*/

//! antenna1d_matrix
/*!
  Core function for setting up the response matrix for 1D antenna cases.
  
  Main task is to extract correct antenna pattern, including frequency 
  interpolation. Actual weights are calculated in *integration_func_by_vectmult*.


   \param   H            The antenna transfer matrix
   \param   antenna_dim  As the WSV with the same name
   \param   antenna_dza  The zenith angle column of *antenna_dlos*.
   \param   antenna_response  As the WSV with the same name
   \param   za_grid      Zenith angle grid for pencil beam calculations
   \param   f_grid       Frequency grid for monochromatic calculations
   \param   n_pol        Number of polarisation states
   \param   do_norm      Flag whether response should be normalised

   \author Mattias Ekstr�m / Patrick Eriksson
   \date   2003-05-27 / 2008-06-17
*/
void antenna1d_matrix(Sparse& H,
                      const Index& antenna_dim,
                      ConstVectorView antenna_dza,
                      const GriddedField4& antenna_response,
                      ConstVectorView za_grid,
                      ConstVectorView f_grid,
                      const Index n_pol,
                      const Index do_norm);



//! antenna2d_interp_gridded_dlos
/*!
  The radiances are treated as a bi-linear function, but the antenna response
  is treated as step-wise constant function (in contrast to 1D). See also
  built-in doc.

   \param   H                 The antenna transfer matrix
   \param   antenna_dim       As the WSV with the same name
   \param   antenna_dlos      As the WSV with the same name
   \param   antenna_response  As the WSV with the same name
   \param   mblock_dlos       As the WSV with the same name
   \param   f_grid            Frequency grid for monochromatic calculations
   \param   n_pol             Number of polarisation states

   \author  Patrick Eriksson
   \date   2020-09-01
*/
void antenna2d_gridded_dlos(Sparse& H,
                            const Index& antenna_dim,
                            ConstMatrixView antenna_dlos,
                            const GriddedField4& antenna_response,
                            ConstMatrixView mblock_dlos,
                            ConstVectorView f_grid,
                            const Index n_pol);



//! antenna2d_interp_response
/*!
  The antenna pattern is interpolated to the dlos directions and solid
  beam angles are applied. See also built-in doc.

   \param   H                 The antenna transfer matrix
   \param   antenna_dim       As the WSV with the same name
   \param   antenna_dlos      As the WSV with the same name
   \param   antenna_response  As the WSV with the same name
   \param   mblock_dlos       As the WSV with the same name
   \param   solid_angles      The solid angle of each dlos
   \param   f_grid            Frequency grid for monochromatic calculations
   \param   n_pol             Number of polarisation states

   \author  Patrick Eriksson
   \date   2018-09-12
*/
void antenna2d_interp_response(Sparse& H,
                               const Index& antenna_dim,
                               ConstMatrixView antenna_dlos,
                               const GriddedField4& antenna_response,
                               ConstMatrixView mblock_dlos,
                               ConstVectorView solid_angles,
                               ConstVectorView f_grid,
                               const Index n_pol);


//! mixer_matrix
/*!
   Sets up the sparse matrix that models the response from sideband filtering
   and the mixer.

   The size of the transfer matrix is changed in the function
   as follows:
     nrows = f_mixer.nelem()
     ncols = f_grid.nelem()

   The returned frequencies are given in IF, so both primary and mirror band
   is converted down.

   \param   H         The mixer/sideband filter transfer matrix
   \param   f_mixer   The frequency grid of the mixer
   \param   lo        The local oscillator frequency
   \param   filter    The sideband filter data. See *sideband_response*
                      for format and constraints.
   \param   f_grid    The original frequency grid of the spectrum
   \param   n_pol     The number of polarisations to consider
   \param   n_sp      The number of spectra (viewing directions)
   \param   do_norm   Flag whether rows should be normalised

   \author Mattias Ekstr�m / Patrick Eriksson
   \date   2003-05-27 / 2008-06-17
*/
void mixer_matrix(Sparse& H,
                  Vector& f_mixer,
                  const Numeric& lo,
                  const GriddedField1& filter,
                  ConstVectorView f_grid,
                  const Index& n_pol,
                  const Index& n_sp,
                  const Index& do_norm);



/** Calculate polarisation H-matrix
 
 Takes into account instrument channel polarisation and zenith angle.

 \param[out] H          Polarisation matrix
 \param[in]  mm_pol     Instrument channel polarisations
 \param[in]  dza        Zenith angle, from reference direction
 \param[in]  stokes_dim Workspace variable
 \param[in]  iy_unit    Workspace variable
 */
void met_mm_polarisation_hmatrix(Sparse& H,
                                 const ArrayOfString& pol,
                                 const Numeric dza,
                                 const Index stokes_dim,
                                 const String& iy_unit);



//! sensor_aux_vectors
/*!
   Sets up the the auxiliary vectors for sensor_response.

   The function assumes that all grids are common, and the full 
   vectors are just the grids repeated.

   \param   sensor_response_f          As the WSV with same name
   \param   sensor_response_pol        As the WSV with same name
   \param   sensor_response_za         As the WSV with same name
   \param   sensor_response_aa         As the WSV with same name
   \param   sensor_response_f_grid     As the WSV with same name
   \param   sensor_response_pol_grid   As the WSV with same name
   \param   sensor_response_dlos_grid  As the WSV with same name

   \author Patrick Eriksson
   \date   2008-06-09
*/
void sensor_aux_vectors(Vector& sensor_response_f,
                        ArrayOfIndex& sensor_response_pol,
                        Matrix& sensor_response_dlos,
                        ConstVectorView sensor_response_f_grid,
                        const ArrayOfIndex& sensor_response_pol_grid,
                        ConstMatrixView sensor_response_dlos_grid);



//! spectrometer_matrix
/*!
   Constructs the sparse matrix that multiplied with the spectral values
   gives the spectra from the spectrometer.

   The input to the function corresponds mainly to WSVs. See f_backend and
   backend_channel_response for how the backend response is specified.

   \param   H             The response matrix.
   \param   ch_f          Corresponds directly to WSV f_backend.
   \param   ch_response   Corresponds directly to WSV backend_channel_response.
   \param   sensor_f      Corresponds directly to WSV sensor_response_f_grid.
   \param   n_pol         The number of polarisations.
   \param   n_sp          The number of spectra (viewing directions).
   \param   do_norm       Corresponds directly to WSV sensor_norm.

   \author Mattias Ekstr�m and Patrick Eriksson
   \date   2003-08-26 / 2008-06-10
*/
void spectrometer_matrix(Sparse& H,
                         ConstVectorView ch_f,
                         const ArrayOfGriddedField1& ch_response,
                         ConstVectorView sensor_f,
                         const Index& n_pol,
                         const Index& n_sp,
                         const Index& do_norm);



//! stokes2pol
/*!
   Sets up a vector to convert the Stokes vector to different polarsiations.

   The measured value is the sum of the element product of the conversion
   vector and the Stokes vector. Schematically:

   y[iout] = sum( w.*iy(iin,joker)

   Vectors for I, Q, U and V are always normalised to have unit length (one
   value is 1, the remaining ones zero). The first element of remaining vectors
   is set to nv (and other values normalised accordingly), to allow that
   calibartion and other normalisation effects can be incorporated.

   \param   s2p           Array of conversion vectors.
   \param   nv            Norm value for polarisations beside I, Q, U and V.

   \author Patrick Eriksson
   \date   2011-11-01 and 2018-03-16
*/
void stokes2pol(VectorView w,
                const Index& stokes_dim,
                const Index& ipol_1based,
                const Numeric nv = 1);



//! Calculate channel boundaries from instrument response functions.
/*!
  This function finds out the unique channel boundaries from
  f_backend and backend_channel_response. This is not a trivial task,
  since channels may overlap, or may be sorted in a strange way. The
  function tries to take care of all that. If channels overlap, they
  are combined to one continuous frequency region. therefore the
  number of elements in the output vectors fmin and fmax can be lower
  than the number of elements in f_backend and
  backend_channel_response. 

  The function also does consistency checking on the two input
  variables.
 
  The output vectors fmin and fmax will be monotonically increasing.

  \author Stefan Buehler
  
  \param[out] fmin                      Vector of lower boundaries of instrument channels.
  \param[out] fmax                      Vector of upper boundaries of instrument channels.
  \param[in]  f_backend                 Nominal backend frequencies.
  \param[in]  backend_channel_response  Channel response, relative to nominal frequencies.
  \param[in]  delta                     Extra margin on both sides of each band. Has a 
                                        default value of 0.
*/
void find_effective_channel_boundaries(  // Output:
    Vector& fmin,
    Vector& fmax,
    // Input:
    const Vector& f_backend,
    const ArrayOfGriddedField1& backend_channel_response,
    const Numeric& delta,
    const Verbosity& verbosity);



//! integration_func_by_vecmult
/*!
   Calculates the (row) vector that multiplied with an unknown (column) vector
   approximates the integral of the product between the functions represented
   by the two vectors: h*g = integral( f(x)*g(x) dx )

   Basic principle follows Eriksson et al., Efficient forward modelling by
   matrix representation of sensor responses, Int. J. Remote Sensing, 27,
   1793-1808, 2006. However, while in Eriksson et al. the product between f*g
   is assumed to vary linearly between the grid point, the expressions applied
   here are more advanced and are completly exact as long as f and g are
   piece-wise linear functions. The product f*g is then a quadratic funtion
   between the grid points.

   \param   h       The multiplication (row) vector.
   \param   f       The values of function f(x).
   \param   x_f_in  The grid points of function f(x). Must be increasing.
   \param   x_g_in  The grid points of function g(x). Can be increasing or 
                    decreasing. Must cover a wider range than x_f (in
                    both ends).

   \author Mattias Ekstr�m and Patrick Eriksson
   \date   2003-02-13 / 2008-06-12
*/
void integration_func_by_vecmult(VectorView h,
                                 ConstVectorView f,
                                 ConstVectorView x_f_in,
                                 ConstVectorView x_g_in);



//! integration_bin_by_vecmult
/*!
   Calculates the (row) vector that multiplied with an unknown (column) vector,
   g, approximates the integral between limit1 and limit2, where limit1 >=
   limit2.

   This can be seen as a special case of what is handled by 
   *integration_func_by_vecmult*, where the function f is a boxcar function. Or
   expressed differently, the function g is "binned" between limit1 and limit2.

   The limits must be inside the range the x_g.

   \param   h       The multiplication (row) vector.
   \param   x_g_in  The grid points of function g(x). Can be increasing or 
                    decreasing. Must cover a wider range than the boxcar
                    function (in both ends).
   \param   limit1  The lower integration limit.
   \param   limit2  The upper integration limit.

   \author Patrick Eriksson
   \date   2017-06-02
*/
void integration_bin_by_vecmult(VectorView h,
                                ConstVectorView x_g_in,
                                const Numeric& limit1,
                                const Numeric& limit2);


//! summation_by_vecmult
/*!
   Calculates the (row) vector that multiplied with an unknown
   (column) vector approximates the sum of the product 
   between the functions at two points.

   E.g. h*g = f(x1)*g(x1) + f(x2)*g(x2)

   The typical application is to set up the combined response matrix
   for mixer and sideband filter.

   See Eriksson et al., Efficient forward modelling by matrix
   representation of sensor responses, Int. J. Remote Sensing, 27,
   1793-1808, 2006, for details.

   No normalisation of the response is made.

   \param   h     The summation (row) vector.
   \param   f     Sideband response.
   \param   x_f   The grid points of function f(x).
   \param   x_g   The grid for spectral values (normally equal to f_grid) 
   \param   x1    Point 1
   \param   x2    Point 2

   \author Mattias Ekstr�m / Patrick Eriksson
   \date   2003-05-26 / 2008-06-17
*/
void summation_by_vecmult(VectorView h,
                          ConstVectorView f,
                          ConstVectorView x_f,
                          ConstVectorView x_g,
                          const Numeric x1,
                          const Numeric x2);

#endif  // sensor_h
