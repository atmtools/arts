/* Copyright (C) 2003-2008 Mattias Ekström <ekstrom@rss.chalmers.se>

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

/*!a
  \file   sensor.cc
  \author Mattias Ekström <ekstrom@rss.chalmers.se>
  \date   2003-02-27

  \brief  Functions related to sensor modelling.

  Functions to model sensor behaviour and integration calculated as vector
  multiplication.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <list>
#include "arts.h"
#include "logic.h"
#include "matpackI.h"
#include "matpackII.h"
#include "messages.h"
#include "sensor.h"

extern const Numeric PI;
extern const Index GFIELD1_F_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_F_GRID;
extern const Index GFIELD4_ZA_GRID;
extern const Index GFIELD4_AA_GRID;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

void antenna1d_matrix(      
           Sparse&   H,
#ifndef NDEBUG
      const Index&   antenna_dim,
#else
      const Index&   antenna_dim _U_,
#endif
   ConstMatrixView   antenna_los,
    const GField4&   antenna_response,
   ConstVectorView   za_grid,
   ConstVectorView   f_grid,
       const Index   n_pol,
       const Index   do_norm )
{
  // Number of input za and frequency angles
  const Index n_za = za_grid.nelem();
  const Index n_f  = f_grid.nelem();

  // Calculate number of antenna beams
  const Index n_ant = antenna_los.nrows();

  // Asserts for variables beside antenna_response
  assert( antenna_dim == 1 );
  assert( antenna_los.ncols() == antenna_dim );
  assert( n_za >= 2 );
  assert( n_pol >= 1 );
  assert( do_norm >= 0  &&  do_norm <= 1 );
  
  // Extract antenna_response grids
  const Index n_ar_pol = 
                  antenna_response.get_string_grid(GFIELD4_FIELD_NAMES).nelem();
  ConstVectorView aresponse_f_grid = 
                  antenna_response.get_numeric_grid(GFIELD4_F_GRID);
  ConstVectorView aresponse_za_grid = 
                  antenna_response.get_numeric_grid(GFIELD4_ZA_GRID);
  DEBUG_ONLY( const Index n_ar_aa = 
                  antenna_response.get_numeric_grid(GFIELD4_AA_GRID).nelem(); )

  //
  const Index n_ar_f  = aresponse_f_grid.nelem();
  const Index n_ar_za = aresponse_za_grid.nelem();
  const Index pol_step = n_ar_pol > 1;
  
  // Asserts for antenna_response
  assert( n_ar_pol == 1  ||  n_ar_pol == n_pol );
  assert( n_ar_f );
  assert( n_ar_za > 1 );
  assert( n_ar_aa == 1 );

  // If response data extend outside za_grid is checked in 
  // sensor_integration_vector
  

  // Some size(s)
  const Index nfpol = n_f * n_pol;  

  // Resize H
  H.resize( n_ant*nfpol, n_za*nfpol );

  // Storage vectors for response weights
  Vector hrow( H.ncols(), 0.0 );
  Vector hza( n_za, 0.0 );

  // Antenna response to apply (possibly obtained by frequency interpolation)
  Vector aresponse( n_ar_za, 0.0 );


  // Antenna beam loop
  for( Index ia=0; ia<n_ant; ia++ )
    {
      Vector shifted_aresponse_za_grid  = aresponse_za_grid;
             shifted_aresponse_za_grid += antenna_los(ia,0);


      // Order of loops assumes that the antenna response more often
      // changes with frequency than for polarisation

      // Frequency loop
      for( Index f=0; f<n_f; f++ )
        {

          // Polarisation loop
          for( Index ip=0; ip<n_pol; ip++ )
            {
              // Determine antenna pattern to apply
              //
              // Interpolation needed only if response has a frequency grid
              // New antenna for each loop of response changes with polarisation
              //
              Index new_antenna = 1; 
              //
              if( n_ar_f > 1 )
                {
                  // Interpolation (do this in "green way")
                  ArrayOfGridPos gp_f( 1 ), gp_za(n_za);
                  gridpos( gp_f, aresponse_f_grid, Vector(1,f_grid[f]) );
                  gridpos( gp_za, aresponse_za_grid, aresponse_za_grid );
                  Tensor3 itw( 1, n_za, 4 );
                  interpweights( itw, gp_f, gp_za );
                  Matrix aresponse_matrix(1,n_za);
                  interp( aresponse_matrix, itw, 
                          antenna_response(ip,joker,joker,0), gp_f, gp_za );
                  aresponse = aresponse_matrix(0,joker);
                }
              else if( pol_step )   // Response changes with polarisation
                {
                  aresponse = antenna_response(ip,0,joker,0);
                }
              else if( f == 0 )  // Same response for all f and polarisations
                {
                  aresponse = antenna_response(0,0,joker,0);
                }
              else
                {
                  new_antenna = 0;
                }

              // Calculate response weights
              if( new_antenna )
                {
                  sensor_integration_vector( hza, aresponse,
                                             shifted_aresponse_za_grid,
                                             za_grid );
                  // Normalisation?
                  if( do_norm )
                    { hza /= hza.sum(); }
                }

              // Put weights into H
              //
              const Index ii = f*n_pol + ip;
              //
              hrow[ Range(ii,n_za,nfpol) ] = hza;
              //
              H.insert_row( ia*nfpol+ii, hrow );
              //
              hrow = 0;
            }
        }
    }
}



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

   \author Mattias Ekström / Patrick Eriksson
   \date   2003-05-27 / 2008-06-17
*/
void mixer_matrix(
           Sparse&   H,
           Vector&   f_mixer,
    const Numeric&   lo,
    const GField1&   filter,
   ConstVectorView   f_grid,
      const Index&   n_pol,
      const Index&   n_sp,
      const Index&   do_norm )
{
  // Frequency grid of for sideband response specification
  ConstVectorView filter_grid = filter.get_numeric_grid(GFIELD1_F_GRID);

  DEBUG_ONLY( const Index nrp = filter.nelem(); )

  // Asserts
  assert( lo > f_grid[0] );
  assert( lo < last(f_grid) );
  assert( filter_grid.nelem() == nrp );
  assert( fabs(last(filter_grid)+filter_grid[0]) < 1e3 );
  // If response data extend outside f_grid is checked in sensor_summation_vector

  // Find indices in f_grid where f_grid is just below and above the
  // lo frequency.
  Index i_low = 0, i_high = f_grid.nelem()-1, i_mean;
  while( i_high-i_low > 1 )
    {
      i_mean = (Index) (i_high+i_low)/2;
      if (f_grid[i_mean]<lo)
        { 
          i_low = i_mean; 
        }
      else
        {
          i_high = i_mean;
        }
    }
  if (f_grid[i_high]==lo)
    {
      i_high++;
    }

  // Determine IF limits for new frequency grid
  const Numeric lim_low  = max( lo-f_grid[i_low], f_grid[i_high]-lo );
  const Numeric lim_high = -filter_grid[0];

  // Convert sidebands to IF and use list to make a unique sorted
  // vector, this sorted vector is f_mixer.
  list<Numeric> l_mixer;
  for( Index i=0; i<f_grid.nelem(); i++ )
    {
      if( fabs(f_grid[i]-lo)>=lim_low && fabs(f_grid[i]-lo)<=lim_high )
        {
          l_mixer.push_back(fabs(f_grid[i]-lo));
        }
    }
  l_mixer.push_back(lim_high);   // Not necessarily a point in f_grid
  l_mixer.sort();
  l_mixer.unique();
  f_mixer.resize((Index) l_mixer.size());
  Index e=0;
  for (list<Numeric>::iterator li=l_mixer.begin(); li != l_mixer.end(); li++)
    {
      f_mixer[e] = *li;
      e++;
    }

  // Resize H
  H.resize( f_mixer.nelem()*n_pol*n_sp, f_grid.nelem()*n_pol*n_sp );

  // Calculate the sensor summation vector and insert the values in the
  // final matrix taking number of polarisations and zenith angles into
  // account.
  Vector row_temp( f_grid.nelem() );
  Vector row_final( f_grid.nelem()*n_pol*n_sp );
  //
  Vector if_grid  = f_grid;
         if_grid -= lo;
  //
  for( Index i=0; i<f_mixer.nelem(); i++ ) 
    {
      sensor_summation_vector( row_temp, filter, filter_grid, 
                               if_grid, f_mixer[i], -f_mixer[i] );

      // Normalise if flag is set
      if (do_norm)
        row_temp /= row_temp.sum();

      // Loop over number of polarisations
      for (Index p=0; p<n_pol; p++)
        {
          // Loop over number of zenith angles/antennas
          for (Index a=0; a<n_sp; a++)
            {
              // Distribute elements of row_temp to row_final.
              row_final = 0.0;
              row_final[Range(a*f_grid.nelem()*n_pol+p,f_grid.nelem(),n_pol)]
                                                                     = row_temp;
              H.insert_row(a*f_mixer.nelem()*n_pol+p+i*n_pol,row_final);
            }
        }
    }
}



//! sensor_aux_vectors
/*!
   Sets up the the auxiliary vectors for sensor_response.

   The function assumes that all grids are common, and the aux vectors
   are just the grids repeated

   \param   sensor_response_f          As the WSV with same name
   \param   sensor_response_pol        As the WSV with same name
   \param   sensor_response_za         As the WSV with same name
   \param   sensor_response_aa         As the WSV with same name
   \param   sensor_response_f_grid     As the WSV with same name
   \param   sensor_response_pol_grid   As the WSV with same name
   \param   sensor_response_za_grid    As the WSV with same name
   \param   sensor_response_aa_grid    As the WSV with same name

   \author Patrick Eriksson
   \date   2008-06-09
*/
void sensor_aux_vectors(
               Vector&   sensor_response_f,
         ArrayOfIndex&   sensor_response_pol,
               Vector&   sensor_response_za,
               Vector&   sensor_response_aa,
       ConstVectorView   sensor_response_f_grid,
   const ArrayOfIndex&   sensor_response_pol_grid,
       ConstVectorView   sensor_response_za_grid,
       ConstVectorView   sensor_response_aa_grid )
{
  // Sizes
  const Index nf       = sensor_response_f_grid.nelem();
  const Index npol     = sensor_response_pol_grid.nelem();
  const Index nza      = sensor_response_za_grid.nelem();
        Index naa      = sensor_response_aa_grid.nelem();
        Index empty_aa = 0;
  //
  if( naa == 0 )
    {
      empty_aa = 1;
      naa      = 1; 
    }
  //
  const Index n = nf * npol * nza * naa;

  // Allocate
  sensor_response_f.resize( n );
  sensor_response_pol.resize( n );
  sensor_response_za.resize( n );
  if( empty_aa )
    { sensor_response_aa.resize( 0 ); }
  else
    { sensor_response_aa.resize( n ); }
  
  // Fill
  for( Index iaa=0; iaa<naa; iaa++ )
    {
      const Index i1 = iaa * nza * nf * npol;
      //
      for( Index iza=0; iza<nza; iza++ )
        {
          const Index i2 = i1 + iza * nf * npol;
          //
          for( Index ifr=0; ifr<nf; ifr++ ) 
            {
              const Index i3 = i2 + ifr * npol;
              //
              for( Index ip=0; ip<npol; ip++ )
                {
                  const Index i = i3 + ip;
                  //
                  sensor_response_f[i]   = sensor_response_f_grid[ifr];
                  sensor_response_pol[i] = sensor_response_pol_grid[ip];
                  sensor_response_za[i]  = sensor_response_za_grid[iza];
                  if( !empty_aa )
                    { sensor_response_aa[i] = sensor_response_aa_grid[iaa]; }
                }
            }
        }
    }
}



//! sensor_integration_vector
/*!
   Calculates the (row) vector that multiplied with an unknown
   (column) vector approximates the integral of the product
   between the functions represented by the two vectors.

   E.g. h*g = integral( f(x)*g(x) dx )

   See Eriksson et al., Efficient forward modelling by matrix
   representation of sensor responses, Int. J. Remote Sensing, 27,
   1793-1808, 2006, for details.

   The grids are internally normalised to cover the range [0,1] for
   increased numerical stability.

   \param   h       The multiplication (row) vector.
   \param   f       The values of function f(x).
   \param   x_f_in  The grid points of function f(x). Must be increasing.
   \param   x_g_in  The grid points of function g(x). Can be increasing or 
                    decreasing. Must cover a wider range than x_ft (in
                    both ends).

   \author Mattias Ekström and Patrick Eriksson
   \date   2003-02-13 / 2008-06-12
*/
void sensor_integration_vector(
        VectorView   h,
   ConstVectorView   f,
   ConstVectorView   x_f_in,
   ConstVectorView   x_g_in )
{
  // Basic sizes 
  const Index nf = x_f_in.nelem();
  const Index ng = x_g_in.nelem();

  // Asserts
  assert( h.nelem() == ng );
  assert( f.nelem() == nf );
  assert( is_increasing( x_f_in ) );
  assert( is_increasing( x_g_in ) || is_decreasing( x_g_in ) );
  // More asserts below

  // Copy grids, handle reversed x_g and normalise to cover the range
  // [0 1]. This is necessary to avoid numerical problems for
  // frequency grids (e.g. experienced for a case with frequencies
  // around 501 GHz).
  //
  Vector x_g         = x_g_in;
  Vector x_f         = x_f_in;
  Index  xg_reversed = 0;
  //
  if( is_decreasing( x_g ) )
    {
      xg_reversed = 1;
      Vector tmp  = x_g[Range(ng-1,ng,-1)];   // Flip order
      x_g         = tmp;
    }
  //
  assert( x_g[0]    <= x_f[0] );
  assert( x_g[ng-1] >= x_f[nf-1] );
  //
  const Numeric xmin = x_g[0];
  const Numeric xmax = x_g[ng-1];
  //
  x_f -= xmin;
  x_g -= xmin;
  x_f /= xmax - xmin;
  x_g /= xmax - xmin;

  //Create a reference grid vector, x_ref that containing the values
  //of x_f and x_g strictly sorted. Only g points inside the f range
  //are of concern.
  list<Numeric> l_x;
  for( Index it=0; it<nf; it++ )
    l_x.push_back(x_f[it]);
  for (Index it=0; it<ng; it++) 
    {
      if( x_g[it]>x_f[0] && x_g[it]<x_f[x_f.nelem()-1] )
        l_x.push_back(x_g[it]);
    }

  l_x.sort();
  l_x.unique();

  Vector x_ref(l_x.size());
  Index e=0;
  for (list<Numeric>::iterator li=l_x.begin(); li != l_x.end(); li++) {
    x_ref[e] = *li;
    e++;
  }

  //Initiate output vector, with equal size as x_g, with zeros.
  //Start calculations
  h = 0.0;
  Index i_f = 0;
  Index i_g = 0;
  //i = 0;
  Numeric dx,a0,b0,c0,a1,b1,c1,x3,x2,x1;
  //while( i_g < ng && i_f < x_f.nelem() ) {
  for( Index i=0; i<x_ref.nelem()-1; i++ ) {
    //Find for what index in x_g (which is the same as for h) and f
    //calculation corresponds to
    while( x_g[i_g+1] <= x_ref[i] ) {
      i_g++;
    }
    while( x_f[i_f+1] <= x_ref[i] ) {
     i_f++;
    }

    //If x_ref[i] is out of x_f's range then that part of the integral
    //is set to 0, so no calculations will be done
    if( x_ref[i] >= x_f[0] && x_ref[i] < x_f[x_f.nelem()-1] ) {
      //Product of steps in x_f and x_g
      dx = (x_f[i_f+1] - x_f[i_f]) * (x_g[i_g+1] - x_g[i_g]);

      //Calculate a, b and c coefficients; h[i]=ax^3+bx^2+cx
      a0 = (f[i_f] - f[i_f+1]) / 3;
      b0 = (-f[i_f]*(x_g[i_g+1]+x_f[i_f+1])+f[i_f+1]*(x_g[i_g+1]+x_f[i_f]))
           /2;
      c0 = f[i_f]*x_f[i_f+1]*x_g[i_g+1]-f[i_f+1]*x_f[i_f]*x_g[i_g+1];

      a1 = -a0;
      b1 = (f[i_f]*(x_g[i_g]+x_f[i_f+1])-f[i_f+1]*(x_g[i_g]+x_f[i_f]))/2;
      c1 = -f[i_f]*x_f[i_f+1]*x_g[i_g]+f[i_f+1]*x_f[i_f]*x_g[i_g];

      x3 = pow(x_ref[i+1],3) - pow(x_ref[i],3);
      x2 = pow(x_ref[i+1],2) - pow(x_ref[i],2);
      x1 = x_ref[i+1]-x_ref[i];

      //Calculate h[i] and h[i+1] increment
      h[i_g] += (a0*x3+b0*x2+c0*x1) / dx;
      h[i_g+1] += (a1*x3+b1*x2+c1*x1) / dx;

    }
  }

  // Flip back if x_g was decreasing
  if( xg_reversed )
    {
      Vector tmp = h[Range(ng-1,ng,-1)];   // Flip order
      h = tmp;
    }
}



//! sensor_summation_vector
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

   \author Mattias Ekström / Patrick Eriksson
   \date   2003-05-26 / 2008-06-17
*/
void sensor_summation_vector(
        VectorView   h,
   ConstVectorView   f,
   ConstVectorView   x_f,
   ConstVectorView   x_g,
     const Numeric   x1,
     const Numeric   x2 )
{
  // Asserts
  assert( h.nelem() == x_g.nelem() );
  assert( f.nelem() == x_f.nelem() );
  assert( x_g[0]    <= x_f[0] );
  assert( last(x_g) >= last(x_f) );
  assert( x1        >= x_f[0] );
  assert( x2        >= x_f[0] );
  assert( x1        <= last(x_f) );
  assert( x2        <= last(x_f) );

  // Determine grid positions for point 1 (both with respect to f and g grids)
  // and interpolate response function.
  ArrayOfGridPos gp1g(1), gp1f(1);
  gridpos( gp1g, x_g, x1 );
  gridpos( gp1f, x_f, x1 );
  Matrix itw1(1,2);
  interpweights( itw1, gp1f );
  Numeric f1;
  interp( f1, itw1, f, gp1f );

  // Same for point 2
  ArrayOfGridPos gp2g(1), gp2f(1);
  gridpos( gp2g, x_g, x2 );
  gridpos( gp2f, x_f, x2 );
  Matrix itw2(1,2);
  interpweights( itw2, gp2f );
  Numeric f2;
  interp( f2, itw2, f, gp2f );

  //Initialise h at zero and store calculated weighting components
  h = 0.0;
  h[gp1g[0].idx]   += f1 * gp1g[0].fd[1];
  h[gp1g[0].idx+1] += f1 * gp1g[0].fd[0];
  h[gp2g[0].idx]   += f2 * gp2g[0].fd[1];
  h[gp2g[0].idx+1] += f2 * gp2g[0].fd[0];
}



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

   \author Mattias Ekström and Patrick Eriksson
   \date   2003-08-26 / 2008-06-10
*/
void spectrometer_matrix( 
           Sparse&         H,
   ConstVectorView         ch_f,
   const ArrayOfGField1&   ch_response,
   ConstVectorView         sensor_f,
      const Index&         n_pol,
      const Index&         n_sp,
      const Index&         do_norm )
{
  // Check if matrix has one frequency column or one for every channel
  // frequency
  //
  assert( ch_response.nelem()==1 || ch_response.nelem()==ch_f.nelem() );
  //
  Index freq_full = ch_response.nelem() > 1;

  // If response data extend outside sensor_f is checked in 
  // sensor_integration_vector

  // Reisze H
  //
  const Index   nin_f  = sensor_f.nelem();
  const Index   nout_f = ch_f.nelem();
  const Index   nin    = n_sp * nin_f  * n_pol;
  const Index   nout   = n_sp * nout_f * n_pol;
  //
  H.resize( nout, nin );

  // Calculate the sensor integration vector and put values in the temporary
  // vector, then copy vector to the transfer matrix
  //
  Vector ch_response_f;
  Vector weights( nin_f );
  Vector weights_long( nin, 0.0 );
  //
  for( Index ifr=0; ifr<nout_f; ifr++ ) 
    {
      const Index irp = ifr * freq_full;

      //The spectrometer response is shifted for each centre frequency step
      ch_response_f  = ch_response[irp].get_numeric_grid(GFIELD1_F_GRID);
      ch_response_f += ch_f[ifr];

      // Call sensor_integration_vector and store it in the temp vector
      sensor_integration_vector( weights, ch_response[irp],
                                 ch_response_f, sensor_f );

      // Normalise if flag is set
      if( do_norm )
        weights /= weights.sum();

      // Loop over polarisation and spectra (viewing directions)
      // Weights change only with frequency
      for( Index sp=0; sp<n_sp; sp++ ) 
        {
          for( Index pol=0; pol<n_pol; pol++ ) 
            {
              // Distribute the compact weight vector into weight_long
              weights_long[Range(sp*nin_f*n_pol+pol,nin_f,n_pol)] = weights;

              // Insert temp_long into H at the correct row
              H.insert_row( sp*nout_f*n_pol + ifr*n_pol + pol, weights_long );

              // Reset weight_long to zero.
              weights_long = 0.0;
            }
        }
    }
}






//--- Stuff not yet updated --------------------------------------------------







//! polarisation_matrix
/*!
   Sets up the polarisation transfer matrix from stokes vectors describing
   the sensor polarisation.

   The sensor polarisation matrix is here multiplied 0.5 to get intensities.

   The size of the transfer matrix has to be set up before calling the function
   as follows:
     nrows = number of polarisations times frequencies and angles
     ncols = stokes dimension times frequencies and angles.

   \param   H         The polarisation transfer matrix
   \param   pol       The polarisation matrix
   \param   n_f       The number of frequencies
   \param   n_za      The number of zenith angles/antennas
   \param   dim       The stokes dimension

   \author Mattias Ekström
   \date   2004-06-02

void polarisation_matrix(
              Sparse&   H,
      ConstMatrixView   pol,
          const Index   n_f,
          const Index   n_za,
          const Index   dim )
{
  // Assert size of H and pol
  assert( H.nrows()==pol.nrows()*n_f*n_za );
  assert( H.ncols()==dim*n_f*n_za );
  assert( pol.ncols()==dim );

  Index n_pol = pol.nrows();
  Matrix pol_half = pol;
  pol_half *= 0.5;

  // Loop over angles
  for (Index za=0; za<n_za; za++) {
    //FIXME: Add rotation here?

    // Loop over frequencies
    for (Index f=0; f<n_f; f++) {

      // Loop over stokes dimensions
      for (Index d=0; d<dim; d++) {

        // Loop over polarisations
        for (Index p=0; p<n_pol; p++) {

          if ( pol(p,d)!=0.0 )
            H.rw(za*n_f*n_pol+f*n_pol+p,za*n_f*dim+f*dim+d)=pol_half(p,d);
        }
      }
    }
  }
}
*/


//! rotation_matrix
/*!
   Sets up the rotation transfer matrix from the sensor rotation vector.

   The sensor rotation vector should contain the rotation for each
   direction. It is coupled with the antenna line-of-sight and has to have the
   same number of elements/rows.

   The size of the transfer matrix has to be set up before calling the function
   and it is a quadratic matrix with sizes equal the product of stokes_dim and
   number of antenna line-of-sight (number of rotations).

   \param   H         The polarisation transfer matrix
   \param   rot       The polarisation matrix
   \param   n_f       The number of frequencies
   \param   dim       The stokes dimension

   \author Mattias Ekström
   \date   2004-06-02

void rotation_matrix(
              Sparse&   H,
      ConstVectorView   rot,
          const Index   n_f,
          const Index   dim )
{
  // Assert that the matrix has the right size
  assert( H.nrows()==H.ncols() );
  assert( H.nrows()==dim*n_f*rot.nelem() );

  // Setup the L matrix for each rotation and distribute the elements for
  // all frequencies in the rotation matrix.
  Matrix L(dim,dim,0.0);
  if( dim==4 )
    L(3,3) = 1.0;
  L(0,0) = 1.0;
  Index tmp;

  for( Index rit=0; rit<rot.nelem(); rit++ ) {
    L(1,1) = cos(2*rot[rit]*PI/180.);
    L(2,2) = L(1,1);
    L(1,2) = sin(2*rot[rit]*PI/180.);
    L(2,1) = -L(1,2);

    for( Index fit=0; fit<n_f; fit++ ) {
      for( Index Lcit=0; Lcit<dim; Lcit++ ) {
        for( Index Lrit=0; Lrit<dim; Lrit++ ) {
          tmp = (rit*n_f+fit)*dim;
          H.rw(tmp+Lrit,tmp+Lcit)=L(Lrit,Lcit);
        }
      }
    }
  }
}
*/



