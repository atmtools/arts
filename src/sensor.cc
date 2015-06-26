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

/*!
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
#include <stdexcept>
#include "arts.h"
#include "logic.h"
#include "make_vector.h"
#include "matpackI.h"
#include "matpackII.h"
#include "messages.h"
#include "sorting.h"
#include "sensor.h"

extern const Numeric PI;
extern const Numeric NAT_LOG_2;
extern const Numeric DEG2RAD;
extern const Index GFIELD1_F_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_F_GRID;
extern const Index GFIELD4_ZA_GRID;
extern const Index GFIELD4_AA_GRID;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! antenna1d_matrix
/*!
  Core function for setting up the response matrix for 1D antenna cases.
  
  Main task is to extract correct antenna pattern, including frequency 
  interpolation. Actual weights are calculated in *sensor_integration_vector*.


   \param   H            The antenna transfer matrix
   \param   antenna_dim  As the WSV with the same name
   \param   antenna_dza  The zenith angle column of *antenna_dlos*.
   \param   antenna_response  As the WSV with the same name
   \param   za_grid      Zenith angle grid for pencil beam calculations
   \param   f_grid       Frequency grid for monochromatic calculations
   \param   n_pol        Number of polarisation states
   \param   do_norm      Flag whether response should be normalised

   \author Mattias Ekström / Patrick Eriksson
   \date   2003-05-27 / 2008-06-17
*/
void antenna1d_matrix(      
           Sparse&   H,
#ifndef NDEBUG
      const Index&   antenna_dim,
#else
      const Index&   antenna_dim _U_,
#endif
   ConstVectorView   antenna_dza,
    const GriddedField4&   antenna_response,
   ConstVectorView   za_grid,
   ConstVectorView   f_grid,
       const Index   n_pol,
       const Index   do_norm )
{
  // Number of input pencil beam directions and frequency
  const Index n_za = za_grid.nelem();
  const Index n_f  = f_grid.nelem();

  // Calculate number of antenna beams
  const Index n_ant = antenna_dza.nelem();

  // Asserts for variables beside antenna_response
  assert( antenna_dim == 1 );
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
             shifted_aresponse_za_grid += antenna_dza[ia];

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
              if( n_ar_f > 1 )  // If multiple freqs in response, new interpolation
                {               // for each freq and pol
                  // Interpolation (do this in "green way")
                  ArrayOfGridPos gp_f( 1 ), gp_za(n_ar_za);
                  gridpos( gp_f, aresponse_f_grid, Vector(1,f_grid[f]) );
                  gridpos( gp_za, aresponse_za_grid, aresponse_za_grid );
                  Tensor3 itw( 1, n_ar_za, 4 );
                  interpweights( itw, gp_f, gp_za );
                  Matrix aresponse_matrix(1,n_ar_za);
                  interp( aresponse_matrix, itw, 
                         antenna_response.data(ip,joker,joker,0), gp_f, gp_za );
                  aresponse = aresponse_matrix(0,joker);
                }
              else if( pol_step )   // Response changes with polarisation,
                {                   // but not frequency
                  aresponse = antenna_response.data(ip,0,joker,0);
                }
              else if( f == 0 )  // Same response for all f and polarisations
                {
                  aresponse = antenna_response.data(0,0,joker,0);
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



// An old version of the function above

//! antenna2d_simplified
/*!
  A first function for setting up the response matrix for 2D antenna cases.
  
  Following the ARTS definitions, a bi-linear variation (in za and aa
  dimensions) for both antenna pattern and pencil beam radiances
  should be assumed here. This function does not handle this. It
  performs instead a series of 1D antenna calculations and "sums up"
  the results. In this summation, both antenna pattern and radiances
  are assumed to constant in the azimuthal direction around each point
  in aa_grid (corresponding to mblock_aa_grid). That is, for azimuth,
  a step-wise function is used instead of a piecewise linear one.

   \param   H            The antenna transfer matrix
   \param   antenna_dim  As the WSV with the same name
   \param   antenna_los  As the WSV with the same name
   \param   antenna_response  As the WSV with the same name
   \param   za_grid      Zenith angle grid for pencil beam calculations
   \param   aa_grid      Azimuth angle grid for pencil beam calculations
   \param   f_grid       Frequency grid for monochromatic calculations
   \param   n_pol        Number of polarisation states
   \param   do_norm      Flag whether response should be normalised

   \author Patrick Eriksson
   \date   2009-09-16
*/
/*
void antenna2d_simplified(      
           Sparse&   H,
#ifndef NDEBUG
      const Index&   antenna_dim,
#else
      const Index&   antenna_dim _U_,
#endif
   ConstMatrixView   antenna_los,
    const GriddedField4&   antenna_response,
   ConstVectorView   za_grid,
   ConstVectorView   aa_grid,
   ConstVectorView   f_grid,
       const Index   n_pol,
       const Index   do_norm )
{
  // Sizes
  const Index n_f      = f_grid.nelem();
  const Index nfpol    = n_f * n_pol;  
  const Index n_aa     = aa_grid.nelem();
  const Index n_za     = za_grid.nelem();
  const Index n_ant    = antenna_los.nrows();
  const Index n_ar_pol = antenna_response.data.nbooks();
  const Index n_ar_f   = antenna_response.data.npages();
  const Index n_ar_za  = antenna_response.data.nrows();

  // Asserts for variables beside antenna_response (not done in antenna1d_matrix)
  assert( antenna_dim == 2 );
  assert( n_aa >= 2 );
  assert( do_norm >= 0  &&  do_norm <= 1 );

  // Make copy of antenna response suitable as input to antenna1d_matrix
  //
  GriddedField4 aresponse = antenna_response;
  //
  ConstVectorView response_aa_grid = 
                  antenna_response.get_numeric_grid(GFIELD4_AA_GRID);
  //
  aresponse.resize( n_ar_pol, n_ar_f, n_ar_za, 1 );
  aresponse.set_grid( GFIELD4_AA_GRID, Vector(1,0) );

  // Resize H
  H.resize( n_ant*nfpol, n_aa*n_za*nfpol );

  // Loop antenna_los
  for( Index il=0; il<n_ant; il++ )
    {

      // Set up an ArrayOfVector that can hold all data for one antenna_los
      ArrayOfVector hrows(nfpol);
      for( Index row=0; row<nfpol; row++ )
        {
          hrows[row].resize(n_aa*n_za*nfpol);
          hrows[row] = 0;     // To get correct value for aa_grid 
        }                     // points outside response aa grid
 
      // Loop azimuth angles
      for( Index ia=0; ia<n_aa; ia++ )
        {
          const Numeric aa_point = aa_grid[ia] - antenna_los(il,1);

          if( aa_point >= response_aa_grid[0]  &&  
                                            aa_point <= last(response_aa_grid) )
            {
              // Interpolate antenna patterns to aa_grid[ia] 
              // Use grid position function to find weights 
              //
              ArrayOfGridPos gp( 1 );
              gridpos( gp, response_aa_grid, Vector(1,aa_point) );
              //
              for( Index i4=0; i4<n_ar_pol; i4++ )
                {
                  for( Index i3=0; i3<n_ar_f; i3++ )
                    {
                      for( Index i2=0; i2<n_ar_za; i2++ )
                        {
                          aresponse.data(i4,i3,i2,0) = 
                            gp[0].fd[1] * antenna_response.data(i4,i3,i2,gp[0].idx) +
                            gp[0].fd[0] * antenna_response.data(i4,i3,i2,gp[0].idx+1);
                        }  
                    }   
                }
 
              // Find the aa width for present angle 
              //
              // Lower and upper end of "coverage" for present aa angle
              Numeric aa_low = response_aa_grid[0];
              if( ia > 0 ) 
                { 
                  const Numeric aam = antenna_los(il,1) + 
                                          ( aa_grid[ia] + aa_grid[ia-1] ) / 2.0;
                  if( aam > aa_low )
                    { aa_low = aam; };
                }
              Numeric aa_high = last(response_aa_grid);
              if( ia < n_aa-1 )
                { 
                  const Numeric aam = antenna_los(il,1) + 
                                          ( aa_grid[ia+1] + aa_grid[ia] ) / 2.0;
                  if( aam < aa_high )
                    { aa_high = aam; };
                }
              //
              const Numeric aa_width = aa_high - aa_low;

              // Do 1D calculations
              //
              Sparse Hza;
              //
              antenna1d_matrix( Hza, 1, antenna_los(il,0), 
                                         aresponse, za_grid, f_grid, n_pol, 0 );

              for( Index row=0; row<nfpol; row++ )
                { 
                  for( Index iz=0; iz<n_za; iz++ )
                    {
                      for( Index i=0; i<nfpol; i++ )
                        {
                          hrows[row][(iz*n_aa+ia)*nfpol+i] = 
                                                 aa_width * Hza(row,iz*nfpol+i);
                        }
                    }
                }   
            }  // if-statement
        }  // aa loop

      // Move results to H
      for( Index row=0; row<nfpol; row++ )
        {
          if( do_norm )
            { 
              hrows[row] /= hrows[row].sum(); 
            }
          H.insert_row( il*nfpol+row, hrows[row] ); 
        }

    } // antenna_los loop
}
*/


//! gaussian_response_autogrid
/*!
   Returns a 1D gaussian response with a suitable grid

   First a grid is generated. The grid is si*[-xwidth_si:dx:xwidth_si],
   where si is the "standard deviation" corresponding to the FWHM, and
   dx is biggest possible value < dx_si, to enusre an symmetric grid wth end
   points exactly at xwidth_si.
   That is, width and spacing of the grid is specified in terms of number of 
   standard deviations. If xwidth_si is set to 2, the response will cover
   about 95% the complete response. For xwidth_si=3, about 99% is covered.

   y is the gaussian function on grid x, with max at x0 and width following
   fwhm.

   \param   x           Grid generated.
   \param   y           Calculated response.
   \param   x0          The x-position of response centre/max.
   \param   fwhm        The full width at half-maximum of the response
   \param   xwidth_si   The one-sided width of x. See above.
   \param   dx_si       The grid step size of x. See above.

   \author Patrick Eriksson
   \date   2009-09-20
*/
void gaussian_response_autogrid(
           Vector&   x,
           Vector&   y,
    const Numeric&   x0,
    const Numeric&   fwhm,
    const Numeric&   xwidth_si,
    const Numeric&   dx_si )
{
  assert( dx_si <= xwidth_si );

  const Numeric si = fwhm / ( 2 * sqrt( 2 * NAT_LOG_2 ) );

  // Number of points needed to enure that spacing is max dx_si
  const Index n = (Index) floor(2*xwidth_si/dx_si) + 1;

  // Grid for response
  const Numeric dd = si * xwidth_si;
  nlinspace( x, -dd, dd, n );

  // Calculate response
  gaussian_response( y, x, x0, fwhm );
}



//! gaussian_response
/*!
   Returns a 1D gaussian response

   y is the gaussian function on grid x, with max at x0 and width following
   fwhm.

   \param   y           Calculated response.
   \param   x           Grid.
   \param   x0          The x-position of response centre/max.
   \param   fwhm        The full width at half-maximum of the response

   \author Patrick Eriksson
   \date   2009-09-20
*/
void gaussian_response(
          Vector&    y,
    const Vector&    x,
    const Numeric&   x0,
    const Numeric&   fwhm )
{
  const Numeric si = fwhm / ( 2 * sqrt( 2 * NAT_LOG_2 ) );
  const Numeric a = 1 / ( si * sqrt( 2 * PI ) );
  const Index   n = x.nelem();

  y.resize( n );
  //
  for( Index i=0; i<n; i++ )
    y[i] = a * exp( -0.5 * pow((x[i]-x0)/si,2.0) );
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
    const GriddedField1&   filter,
   ConstVectorView   f_grid,
      const Index&   n_pol,
      const Index&   n_sp,
      const Index&   do_norm )
{
  // Frequency grid of for sideband response specification
  ConstVectorView filter_grid = filter.get_numeric_grid(GFIELD1_F_GRID);

  DEBUG_ONLY( const Index nrp = filter.data.nelem(); )

  // Asserts
  assert( lo > f_grid[0] );
  assert( lo < last(f_grid) );
  assert( filter_grid.nelem() == nrp );
  assert( fabs(last(filter_grid)+filter_grid[0]) < 1e3 );
  // If response data extend outside f_grid is checked in sensor_summation_vector

  // Find indices in f_grid where f_grid is just below and above the
  // lo frequency.
  /*
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
  const Numeric lim_low  = max( lo-f_grid[i_low], f_grid[i_high]-lo );
  */

  // Determine IF limits for new frequency grid
  const Numeric lim_low  = 0;
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
      sensor_summation_vector( row_temp, filter.data, filter_grid, 
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



//! mueller_rotation
/*!
   Returns the Mueller matrix for a rotation of the coordinate system defining
   H and V directions.

   The function follows Eq 9 in the sensor response article (Eriksson et al,
   Efficient forward modelling by matrix representation of sensor responses,
   IJRS, 2006).

   The sparse matrix H is not sized by the function, in order to save time for
   repeated usage. Before first call of this function, size H as
   H.resize( stokes_dim, stokes_dim );
   The H returned of this function can be used as input for later calls. That
   is, no need to repeat the resize command above.

   \param   H           Mueller matrix
   \param   stokes_dim  Stokes dimensionality
   \param   rotangle    Rotation angle.

   \author Patrick Eriksson
   \date   2014-09-23
*/
void mueller_rotation(
          Sparse&    H,
    const Index&     stokes_dim,
    const Numeric&   rotangle )
{
  assert( stokes_dim > 1 );
  assert( H.nrows() == stokes_dim );
  assert( H.ncols() == stokes_dim );
  assert( H(0,1) == 0 );
  assert( H(1,0) == 0 );

  H.rw(0,0) = 1;
  const Numeric a = cos( 2*DEG2RAD*rotangle );
  H.rw(1,1) = a;
  if( stokes_dim > 2 )
    {
      assert( H(2,0) == 0 );
      assert( H(0,2) == 0 );

      const Numeric b = sin( 2*DEG2RAD*rotangle );
      H.rw(1,2) = b;
      H.rw(2,1) = -b;
      H.rw(2,2) = a;
      if( stokes_dim > 3 ) 
        { 
          // More values should be checked, but to save time we just assert one 
          assert( H(2,3) == 0 );
          H.rw(3,3) = 1; 
        }
    }
}

/** Calculate polarisation H-matrix
 
 Takes into account instrument channel polarisation and zenith angle.

 \param[out] H          Polarisation matrix
 \param[in]  mm_pol     Instrument channel polarisations
 \param[in]  dza        Zenith angle, from reference direction
 \param[in]  stokes_dim Workspace variable
 \param[in]  iy_unit    Workspace variable
 */
void met_mm_polarisation_hmatrix(Sparse& H,
                                 const ArrayOfString& mm_pol,
                                 const Numeric dza,
                                 const Index stokes_dim,
                                 const String& iy_unit)
{
    assert( stokes_dim > 1 );

    // Set "Stokes vector weights" according to iy_unit
    Numeric w = 0.5;
    if( iy_unit == "PlanckBT"  ||  iy_unit == "RJBT"  )
    { w = 1.0; }

    // Identify (basic) polarisation response and possible sensor specific
    // "rotation" 
    const Index nch = mm_pol.nelem();  // Number of channels
    ArrayOfString pol(nch);
    ArrayOfString rot(nch);
    for (Index i = 0; i < nch; i++)
    {
        if (mm_pol[i] == "AMSU-H")
        {
            rot[i] = "AMSU";
            pol[i] = "H";
        }
        else if (mm_pol[i] == "AMSU-V")
        {
            rot[i] = "AMSU";
            pol[i] = "V";
        }
        else if (mm_pol[i] == "ISMAR-H")
        {
            rot[i] = "ISMAR";
            pol[i] = "H";
        }
        else if (mm_pol[i] == "ISMAR-V")
        {
            rot[i] = "ISMAR";
            pol[i] = "V";
        }
        else if (mm_pol[i] == "MARSS-H")
        {
            rot[i] = "MARSS";
            pol[i] = "H";
        }
        else if (mm_pol[i] == "MARSS-V")
        {
            rot[i] = "MARSS";
            pol[i] = "V";
        }
        else if (mm_pol[i] == "H"   || mm_pol[i] == "V"   || 
                 mm_pol[i] == "LHC" || mm_pol[i] == "RHC" )
        {
            rot[i] = "none";
            pol[i] = mm_pol[i];
        }
        else
        {
            ostringstream os;
            os << "Unknown polarisation " << mm_pol[i];
            throw std::runtime_error(os.str());
        }
    }

    // Vectors representing standard cases of sensor polarisation response
    ArrayOfVector pv;
    stokes2pol( pv, w );

    // Complete H, for all channels
    H = Sparse( nch, nch*stokes_dim );

    for( Index i=0; i<nch; i++ )
    {
        // See stokes2pol for index order used in pv
        Index ipv = -1;
        if( pol[i] == "V" )
        { ipv = 4; }
        else if( pol[i] == "H" )
        { ipv = 5; }
        else if( pol[i] == "LHC" )  // Left hand circular
        { ipv = 8; }
        else if( pol[i] == "RHC" )  // Right hand circular
        { ipv = 9; }
        else
        { assert( 0 ); }

        // Maybe this error messages should be mofified
        if( pv[ipv].nelem() > stokes_dim )
        {
            ostringstream os;
            os << "You have selected a channel polarisation that is not "
            << "covered by present value of *stokes_dim* (the later has to be "
            << "increased).";
            throw runtime_error(os.str());
        }

        // No rotation, just plane polarisation response
        if( rot[i] == "none" )
        {
          // Here we just need to fill the row H
          Vector hrow( nch*stokes_dim, 0.0 );
          hrow[Range(i*stokes_dim,pv[ipv].nelem())] = pv[ipv];
          H.insert_row( i, hrow );      
        }
        
        // Rotation + pol-response
        else
        {
          // Rotation part
          Sparse Hrot( stokes_dim, stokes_dim );
          if (rot[i] == "AMSU")
          {
            // No idea about the sign. Not important if U=0,
            // but matter for U != 0.
            mueller_rotation( Hrot, stokes_dim, abs(dza) );
          }
          else if (rot[i] == "ISMAR")
          {
            // No rotation at -53 (= forward direction). But no idea about the
            // sign, as for AMSU above 
            mueller_rotation( Hrot, stokes_dim, dza+50 );
          }
          else if (rot[i] == "MARSS")
          {
            // MARSS special, as 48 deg between different polarisation (not 90,
            // as for ISMAR. This is best interpretation of information
            // from Stuart. Should be double-checked with him at some point. 
            if( pol[i] == "H" )
              { mueller_rotation( Hrot, stokes_dim, dza+42 ); }
            else
              { mueller_rotation( Hrot, stokes_dim, dza ); }
          }
          else
          {
            assert(0);
          }
        
          // H-matrix matching pv[ipv] (can this made in more compact way?)
          Sparse Hpol( 1, stokes_dim );
          {
            Vector hrow( stokes_dim, 0.0 );
            hrow[Range(0,pv[ipv].nelem())] = pv[ipv];
            Hpol.insert_row( 0, hrow );
          }

          // H for the individual channel
          Sparse Hc( 1, stokes_dim );
          mult( Hc, Hpol, Hrot );

          // Put Hc into H
          Vector hrow( nch*stokes_dim, 0.0 );
          const Index i0=i*stokes_dim;
          for( Index s=0; s<stokes_dim; s++ )
          {  hrow[i0+s] = Hc(0,s); }
          H.insert_row( i, hrow );      
        }
    }
}



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
void sensor_aux_vectors(
               Vector&   sensor_response_f,
         ArrayOfIndex&   sensor_response_pol,
               Matrix&   sensor_response_dlos,
       ConstVectorView   sensor_response_f_grid,
   const ArrayOfIndex&   sensor_response_pol_grid,
       ConstMatrixView   sensor_response_dlos_grid )
{
  // Sizes
  const Index nf    = sensor_response_f_grid.nelem();
  const Index npol  = sensor_response_pol_grid.nelem();
  const Index nlos  = sensor_response_dlos_grid.nrows();
  const Index n     = nf * npol * nlos;

  // Allocate
  sensor_response_f.resize( n );
  sensor_response_pol.resize( n );
  sensor_response_dlos.resize( n, sensor_response_dlos_grid.ncols() );
  
  // Fill
  for( Index ilos=0; ilos<nlos; ilos++ )
    {
      const Index i2 = ilos * nf * npol;
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
              sensor_response_dlos(i,joker) = sensor_response_dlos_grid(ilos,joker);
            }
        }
    }
}



//! sensor_integration_vector
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

  // End points of x_f
  Numeric xfmin = x_f_in[0];
  Numeric xfmax = x_f_in[nf-1];

  // Handle possibly reversed x_g.
  Vector x_g;
  Index  xg_reversed = 0;
  //
  if( is_decreasing( x_g_in ) )
    {
      x_g = x_g_in[Range(ng-1,ng,-1)];
      xg_reversed = 1;
    }
  else
    { 
      x_g = x_g_in; 
    }
  //
  assert( x_g[0]    <= xfmin );
  assert( x_g[ng-1] >= xfmax );

  // Normalise grids so x_f covers [0,1]
  const Numeric df = xfmax - xfmin;
  Vector x_f(nf);
  //
  for( Index i=0; i<nf; i++ )
    { x_f[i] = ( x_f_in[i] - xfmin ) / df; }
  for( Index i=0; i<ng; i++ )
    { x_g[i] = ( x_g[i] - xfmin ) / df; }
  xfmin = 0;
  xfmax = 1;
  // To test without normalisation, comment out lines above and use:
  //const Numeric df  = 1;
  //const Vector  x_f = x_f_in;
  
  // Create a reference grid vector, x_ref that containing the values
  // of x_f and x_g strictly sorted. Only g points inside the f range
  // are of concern.
  list<Numeric> l_x;
  for( Index it=0; it<nf; it++ )
    { l_x.push_back(x_f[it]); }
  for (Index it=0; it<ng; it++) 
    {
      if( x_g[it] > xfmin  &&  x_g[it] < xfmax )  
        { l_x.push_back(x_g[it]); }
    }
  l_x.sort();
  l_x.unique();
  //
  Vector x_ref(l_x.size());
  Index  e = 0;
  for( list<Numeric>::iterator li=l_x.begin(); li != l_x.end(); li++ ) 
    {
      x_ref[e] = *li;
      e++;
    }

  // Initiate output vector, with equal size as x_g, with zeros.
  // Start calculations
  h = 0.0;
  Index i_f = 0;
  Index i_g = 0;
  Numeric dx,a0,b0,c0,a1,b1,c1,x3,x2,x1;
  //
  for( Index i=0; i<x_ref.nelem()-1; i++ ) 
    {
      // Find for what index in x_g (which is the same as for h) and f
      // calculation corresponds to
      while( x_g[i_g+1] <= x_ref[i] ) 
        { i_g++; }
      while( x_f[i_f+1] <= x_ref[i] ) 
        { i_f++; }

      // If x_ref[i] is out of x_f's range then that part of the integral is 0,
      // and no calculations should be done
      if( x_ref[i] >= xfmin  &&  x_ref[i] < xfmax ) 
        {
          // Product of steps in x_f and x_g
          dx = ( x_f[i_f+1] - x_f[i_f] ) * ( x_g[i_g+1] - x_g[i_g] );

          // Calculate a, b and c coefficients; h[i]=ax^3+bx^2+cx
          a0 = ( f[i_f] - f[i_f+1] ) / 3.0;
          b0 = (-f[i_f]   * ( x_g[i_g+1] + x_f[i_f+1] ) + 
                 f[i_f+1] * (x_g[i_g+1]+x_f[i_f]) ) / 2.0;
          c0 = x_g[i_g+1] * ( f[i_f] * x_f[i_f+1] - f[i_f+1] * x_f[i_f] );

          a1 = -a0;
          b1 = ( f[i_f]   * ( x_g[i_g] + x_f[i_f+1] ) - 
                 f[i_f+1] * ( x_g[i_g] + x_f[i_f]) ) / 2.0;
          c1 = x_g[i_g] * ( -f[i_f] * x_f[i_f+1] + f[i_f+1] * x_f[i_f] );

          x1 = x_ref[i+1]-x_ref[i];
          // Simple way, but sensitive to numerical problems:
          //x2 = pow(x_ref[i+1],2) - pow(x_ref[i],2);
          //x3 = pow(x_ref[i+1],3) - pow(x_ref[i],3);
          // The same but a numerically better way:
          x2 = x1 * ( 2*x_ref[i] + x1 );
          x3 = x1 * ( 3*x_ref[i]*(x_ref[i]+x1) + x1*x1 );

          // Calculate h[i] and h[i+1] increment
          // (df-factor to compensate for normalisation)
          h[i_g]   += df * ( a0*x3 + b0*x2 + c0*x1 ) / dx;
          h[i_g+1] += df * ( a1*x3 + b1*x2 + c1*x1 ) / dx;
        }
    }

  // Flip back if x_g was decreasing
  if( xg_reversed )
    {
      Vector tmp = h[Range(ng-1,ng,-1)];   // Flip order
      h = tmp;
    }

  // The expressions are sensitive to numerical issues if two points in x_ref
  // are very close compared to the values in x_ref. A test trying to warn when
  // this happens:
  if( min(f) >= 0  &&  min(h) < -1e-15 )
    throw runtime_error( "Significant negative response value obtained, "
                         "despite sensor reponse is strictly positive. This "
                         "indicates numerical problems. Is there any very "
                         "small spacing of the sensor response grid?"
                         "Please, send a report to Patrick if you see this!" );
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

   \param   h       The multiplication (row) vector.
   \param   f       The values of function f(x).
   \param   x_f_in  The grid points of function f(x). Must be increasing.
   \param   x_g_in  The grid points of function g(x). Can be increasing or 
                    decreasing. Must cover a wider range than x_ft (in
                    both ends).

   \author Mattias Ekström and Patrick Eriksson
   \date   2003-02-13 / 2008-06-12
*/
void sensor_integration_vector2(
        VectorView   h,
   ConstVectorView   f,
   ConstVectorView   x_f,
   ConstVectorView   x_g_in )
{
  // Basic sizes 
  const Index nf = x_f.nelem();
  const Index ng = x_g_in.nelem();

  // Asserts
  assert( h.nelem() == ng );
  assert( f.nelem() == nf );
  assert( is_increasing( x_f ) );
  assert( is_increasing( x_g_in ) || is_decreasing( x_g_in ) );
  // More asserts below

  // End points of x_f
  const Numeric xfmin = x_f[0];
  const Numeric xfmax = x_f[nf-1];

  // Handle possibly reversed x_g.
  Vector x_g;
  Index  xg_reversed = 0;
  //
  if( is_decreasing( x_g_in ) )
    {
      x_g = x_g_in[Range(ng-1,ng,-1)];
      xg_reversed = 1;
    }
  else
    { 
      x_g = x_g_in; 
    }
  //
  assert( x_g[0]    <= xfmin );
  assert( x_g[ng-1] >= xfmax );

  // Create a reference grid vector, x_ref that containing the values
  // of x_f and x_g strictly sorted. Only g points inside the f range
  // are of concern.
  list<Numeric> l_x;
  for( Index it=0; it<nf; it++ )
    l_x.push_back(x_f[it]);
  for( Index it=0; it<ng; it++ ) 
    {
      if( x_g[it] > xfmin  &&  x_g[it] < xfmax )
        { l_x.push_back(x_g[it]); }
    }
  l_x.sort();
  l_x.unique();
  //
  Vector x_ref(l_x.size());
  Index  e = 0;
  for( list<Numeric>::iterator li=l_x.begin(); li != l_x.end(); li++ ) 
    {
      x_ref[e] = *li;
      e++;
    }

  // Initiate output vector, with equal size as x_g, with zeros.
  h = 0.0;

  // Start calculations
  //
  Index i_f = 0;
  Index i_g = 0;
  //
  for( Index i=0; i<x_ref.nelem()-1; i++ ) 
    {
      // Find for what index in x_g (which is the same as for h) and f
      // calculation corresponds to
      while( x_g[i_g+1] <= x_ref[i] ) 
        { i_g++; }
      while( x_f[i_f+1] <= x_ref[i] ) 
        { i_f++; }

      // If x_ref[i] is out of x_f's range then that part of the integral
      // is  0, and no calculations should be done
      if( x_ref[i] >= xfmin  &&  x_ref[i] < xfmax ) 
        {
          // ( x_k+1 - x_k )/2
          const Numeric dxk = 0.5 * ( x_ref[i+1] - x_ref[i] );

          // ( x_i+1 - x_i )
          const Numeric dxi = x_g[i_g+1] - x_g[i_g];

          // The ratio of the two quantities above
          const Numeric r = dxk / dxi;
        
          // Values of f at x_ref[i] and x_ref[i+1]
          const Numeric dx = x_f[i_f+1] - x_f[i_f];
          const Numeric f0 = ( f[i_f]   * (x_f[i_f+1]-x_ref[i]) +
                               f[i_f+1] * (x_ref[i]-x_f[i_f]) ) / dx;
          const Numeric f1 = ( f[i_f]   * (x_f[i_f+1]-x_ref[i+1]) +
                               f[i_f+1] * (x_ref[i+1]-x_f[i_f]) ) / dx;

          //cout << "f0/f1 = " << f0-x_ref[i] << " / " << f1-x_ref[i+1] << endl;

          // Add increaments to h
          h[i_g]   += r * ( f0 * ( x_g[i_g+1] - x_ref[i]   ) + 
                            f1 * ( x_g[i_g+1] - x_ref[i+1] ) );
          h[i_g+1] += r * ( f0 * ( x_ref[i]   - x_g[i_g]   ) + 
                            f1 * ( x_ref[i+1] - x_g[i_g]   ) );
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
   const ArrayOfGriddedField1&   ch_response,
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
      sensor_integration_vector( weights, ch_response[irp].data,
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



//! stokes2pol
/*!
   Sets up an array of vectors to convert the Stokes vector to different
   polarsiations.

   The measured value is the sum of the element product of the conversion
   vector and the Stokes vector. Schematically:

   y[iout] = sum( s2p[i].*iy(iin,joker)

   The order of the vectors follow the coding described for sensor_pol, but
   zero-based indexing is used here. That is, the first vector (s2p[0])
   corresponds to I.

   Trailing zeros are not stored, to indicate the required value of
   *stokes_dim*.

   Vectors for I, Q, U and V are always normalised to have unit length (one
   value is 1, the remaining ones zero). The first element of remaining vectors
   is set to nv (and other values normalised accordingly), to allow that
   calibartion and other normalisation effects can be incorporated.

   \param   s2p           Array of conversion vectors.
   \param   nv            Norm value for polarisations beside I, Q, U and V.

   \author Patrick Eriksson
   \date   2011-11-01
*/
void stokes2pol( 
            ArrayOfVector&  s2p,
      const Numeric&        nv )
{
  s2p.resize(10);

  s2p[0] = MakeVector( 1 );               // I
  s2p[1] = MakeVector( 0, 1 );            // Q
  s2p[2] = MakeVector( 0, 0, 1 );         // U
  s2p[3] = MakeVector( 0, 0, 0, 1 );      // V
  s2p[4] = MakeVector( nv, nv );          // Iv
  s2p[5] = MakeVector( nv, -nv );         // Ih
  s2p[6] = MakeVector( nv, 0, nv );       // I+45
  s2p[7] = MakeVector( nv, 0, -nv );      // I-45
  s2p[8] = MakeVector( nv, 0, 0, nv );    // Ilhc
  s2p[9] = MakeVector( nv, 0, 0, -nv );   // Irhc
}




// Functions by Stefan, needed for HIRS:


//! Test if two instrument channels overlap, and if so, merge them. 
/*!
  The channels boundaries are specified in two separate vectors, fmin
  and fmax. These vectors are both input and output. If merging has
  happened, they will each be one element shorter. 

  The positions of the channels to compare is given by the input
  parameters i and j. It is assumed that the minimum frequency of i
  is lower than or equal to that of j.

  Furthermore, it is assumed that i itself is lower than j.

  The range of the first channel (i) will have been extended to
  accomodate the second channel (j). The second channel will have been
  removed.

  The function also handles the updating of index j: If the two
  channels do not overlap, j is increased by one.

  Function returns true if merging has happened.

  \author Stefan Buehler
  
  \return True if channels were merged, otherwise false.
  \retval fmin Lower channel boundaries.
  \retval fmax Upper channel boundaries.
  \param i Index of first channel.
  \param j Index of second channel.
*/
bool test_and_merge_two_channels(Vector& fmin,
                                 Vector& fmax,
                                 Index i,
                                 Index j)
{
  const Index  nf = fmin.nelem();
  assert(fmax.nelem()==nf);
  assert(i>=0 && i<nf);
  assert(j>=0 && j<nf);
  assert(fmin[i]<=fmin[j]);
  assert(i<j);

  // There are three cases to consider:
  // a) The two channels are separate: fmax[i] <  fmin[j]
  // b) They overlap:                  fmax[i] >= fmin[j]
  // c) j is inside i:                 fmax[i] >  fmax[j]

  // In the easiest case (a), we do not have to do anything.
  if (fmax[i] >= fmin[j])
    {
      // We are in case (b) or (c), so we know that we have to combine
      // the channels. The new minimum frequency is fmin[i]. The new
      // maximum frequency is the larger one of the two channels we
      // are combining:
      if (fmax[j] > fmax[i])
        fmax[i] = fmax[j];

      // We now have to kick out element j from both vectors.

      // Number of elements behind j:
      Index n_behind = nf-1 - j;

      Vector dummy = fmin;
      fmin.resize(nf-1);
      fmin[Range(0,j)] = dummy[Range(0,j)];
      if (n_behind > 0)
        fmin[Range(j,n_behind)] = dummy[Range(j+1,n_behind)];
       
      dummy = fmax;
      fmax.resize(nf-1);
      fmax[Range(0,j)] = dummy[Range(0,j)];
      if (n_behind > 0)
        fmax[Range(j,n_behind)] = dummy[Range(j+1,n_behind)];

      return true;
    }

  return false;
}


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
void find_effective_channel_boundaries(// Output:
                                       Vector& fmin,
                                       Vector& fmax,
                                       // Input:
                                       const Vector& f_backend,
                                       const ArrayOfGriddedField1& backend_channel_response,
                                       const Numeric& delta,
                                       const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  // How many channels in total:
  const Index n_chan = f_backend.nelem();

  // Checks on input quantities:

  // There must be at least one channel.
  if (n_chan < 1)
    {
      ostringstream os;
      os << "There must be at least one channel.\n"
         << "(The vector *f_backend* must have at least one element.)";
      throw runtime_error(os.str());
    }

  // There must be a response function for each channel.
  if (n_chan != backend_channel_response.nelem())
    {
      ostringstream os;
      os << "Variables *f_backend_multi* and *backend_channel_response_multi*\n"
         << "must have same number of bands for each LO.";
      throw runtime_error(os.str());
    }

  // Frequency grids for response functions must be strictly increasing.
  for (Index i=0; i<n_chan; ++i)
    {
      // Frequency grid for this response function:
      const Vector& backend_f_grid = backend_channel_response[i].get_numeric_grid(0);

      if ( !is_increasing(backend_f_grid) )
        {
          ostringstream os;
          os << "The frequency grid for the backend channel response of\n"
             << "channel " << i << " is not strictly increasing.\n";
          os << "It is: " << backend_f_grid << "\n";
          throw runtime_error( os.str() );
        }
    }


  // Start the actual work.

  out2 << "  Original channel characteristics:\n"
       << "  min         nominal      max (all in Hz):\n";

  // count the number of passbands as defined by segments of filtershapes, backend_channel_response.data = [0,>0,>0,0]
  // Borders between passbands are identified as [...0,0...]

  // Get a list of original channel boundaries:
  Index numPB = 0;
  for(Index idx=0;idx<n_chan;++idx)
  {
    const Vector& backend_filter = backend_channel_response[idx].data;
    if(backend_filter.nelem()>2)
    { // only run this code when there is more then two elements in the backend
      for(Index idy =1;idy < backend_filter.nelem();++idy)
      {
        if((backend_filter[idy] >0)&& (backend_filter[idy-1]==0))
        {
          numPB ++; // Two consecutive zeros gives the border between two passbands
        }
      }
    }
    else 
    {
    	numPB ++; 
    }
  }

  Vector fmin_pb(numPB);
  Vector fmax_pb(numPB);
  Index pbIdx = 0;

  for(Index idx=0;idx<n_chan;++idx) 
  {
    // Some handy shortcuts:
    //
    // We have to find the first and last frequency where the
    // response is actually different from 0. (No point in making
    // calculations for frequencies where the response is 0.)
    //       Index j=0;
    //       while (backend_response[j] <= 0) ++j;
    //       Numeric bf_min = backend_f_grid[j];

    //       j=nf-1;
    //       while (backend_response[j] <= 0) --j;
    //       Numeric bf_max = backend_f_grid[j];
    //
    // No, aparently the sensor part want values also where the
    // response is zero. So we simply take the grid boundaries here.
    //
    // We need to add a bit of extra margin at both sides,
    // otherwise there is a numerical problem in the sensor WSMs.
    //
    // PE 081003: The accuracy for me (double on 64 bit machine) appears to
    // be about 3 Hz. Select 1 MHz to have a clear margin. Hopefully OK
    // for other machines.
    //
    // SAB 2010-04-14: The approach with a constant delta does not seem to work 
    // well in practice. What I do now is that I add a delta corresponding to a 
    // fraction of the grid spacing. But that is done outside of this function. 
    // So we set delta = 0 here.
    //
    // SAB 2010-05-03: Now we pass delta as a parameter (with a default value of 0).
    //
    // Isoz 2013-05-21: Added methods to ignore areas between passbands
    //
    const Vector& backend_f_grid   = backend_channel_response[idx].get_numeric_grid(0);
    const Vector& backend_filter = backend_channel_response[idx].data;
    if(backend_filter.nelem()>=4)// Is the passband frequency response given explicitly ? e.g. [0,>0,>0,0]
    {
      for(Index idy =1;idy < backend_filter.nelem();++idy)
      {
        if(idy==1)
        {
          fmin_pb[pbIdx] = f_backend[idx]+backend_f_grid[0];
        }
        else if((backend_filter[idy] >0) && (backend_filter[idy-1]==0))
        {
          fmin_pb[pbIdx]= f_backend[idx] + backend_f_grid[idy-1]-delta;
        }
        if((backend_filter[idy] ==0) && (backend_filter[idy-1]>0))
        {
          fmax_pb[pbIdx]= f_backend[idx] + backend_f_grid[idy]+delta;
          out2 << "  " << "fmin_pb "<< fmin_pb[pbIdx] 
            << "  " << "f_backend" <<f_backend[idx] 
            << "  " << "fmax_pb "<<fmax_pb[pbIdx] 
            << "  " << "diff " << fmax_pb[pbIdx]-fmin_pb[pbIdx]
            << "\n";
          pbIdx++;
        }
      }
      fmax_pb[pbIdx-1] = f_backend[idx]+last(backend_f_grid);
    } 
    else // Or are the passbands given implicitly - such as the default for AMSUA and MHS
    {
      fmin_pb[pbIdx]= f_backend[idx] + backend_f_grid[0]-delta ; //delta;
      fmax_pb[pbIdx]= f_backend[idx] + backend_f_grid[backend_f_grid.nelem()-1]+delta; ///delta;
      out2 << "  " << "fmin_pb "<< fmin_pb[pbIdx] 
           << "  " << "f_backend" <<f_backend[pbIdx] 
           << "  " << "fmax_pb "<<fmax_pb[pbIdx] 
           << "\n";
      pbIdx++;
    }
  }

  // The problem is that channels may be overlapping. In that case, we
  // want to create a frequency grid that covers their entire range,
  // but we do not want to duplicate any frequencies.

  // To make matters worse, one or even several channels may be
  // completely inside another very broad channel.

  // Sort channels by frequency:
  // Caveat: A channel may be higher in
  // characteristic frequency f_backend, but also wider, so that it
  // has a lower minimum frequency fmin_orig. (This is the case for
  // some HIRS channels.) We sort by the minimum frequency here, not
  // by f_backend. This is necessary for function
  // test_and_merge_two_channels to work correctly. 
  ArrayOfIndex isorted;
  get_sorted_indexes(isorted,fmin_pb);

  fmin.resize(numPB);
  fmax.resize(numPB);
  out2 <<" resize numPb "<<numPB<<"\n";
  for(Index idx = 0;idx<numPB;++idx)
  {
    fmin[idx] = fmin_pb[isorted[idx]];
    fmax[idx] = fmax_pb[isorted[idx]];
  }
  // We will be testing pairs of channels, and combine them if
  // possible. We have to test always only against the direct
  // neighbour. If that has no overlap, higher channels can not have
  // any either, due to the sorting by fmin.
  //
  // Note that fmin.nelem() changes, as the loop is
  // iterated. Nevertheless this is the correct stop condition.
  for (Index i=0; i<fmin.nelem()-1; ++i)
    {
      bool continue_checking = true;
      // The "i<fmin.nelem()" condition below is necessary, since
      // fmin.nelem() can decrease while the loop is executed, due to
      // merging. 
      while (continue_checking && i<fmin.nelem()-1)
        {
          continue_checking =
            test_and_merge_two_channels(fmin, fmax, i, i+1);

          // Function returns true if merging has taken place.
          // In this case, we have to check again.
  
        }
    }

  out2 << "  New channel characteristics:\n"
       << "  min                       max (all in Hz):\n";
  for (Index i=0; i<fmin.nelem(); ++i)
    out2 << "  " << fmin[i] << "               " << fmax[i] << "\n";

}


