/* Copyright (C)  Mattias Ekström <ekstrom@rss.chalmers.se>

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
#include "matpackI.h"
#include "matpackII.h"
//#include <stdexcept>
//#include "array.h"
//#include "auto_md.h"
//#include "check_input.h"
//#include "math_funcs.h"
//#include "interpolation.h"
#include "sensor.h"
//#include "messages.h"
//#include "mystring.h"
//#include "poly_roots.h"
//#include "special_interp.h"

  extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! antenna_transfer_matrix
/*!
   Constructs the sparse matrix that multiplied with the spectral values
   for one or several line-of-sights models the antenna transfer matrix.
   The matrix it built up of spaced row vectors, to match the format
   of the spectral values.

   The size of the antenna transfer matrix has to be set by the calling
   function, and it must be set as:
    nrows = x_f.nelem()
    ncols = x_f.nelem() * m_za.nelem()

   The number of line-of-sights is determined by the length of the
   measurement block zenith angle grid. The number of sensor response
   matrix rows don't need to match the number of frequency grid points.

   The antenna diagram must either have two columns or x_f.nelem()+1
   columns. The function will use the size of the matrix to determine
   how to handle it.

   \param   H      The antenna transfer matrix.
   \param   m_za   The measurement block grid of zenith angles.
   \param   srm	   The sensor response matrix, i.e. the antenna diagram
   \param   x_f    The frequency grid points.
   \param   ant_za The antenna zenith angle grid.
   \param   n_pol  The number of polarisations to consider.

   \author Mattias Ekström
   \date   2003-04-09
*/
void antenna_transfer_matrix( Sparse&   H,
                      ConstVectorView   m_za,
          const ArrayOfArrayOfMatrix&   diag,
                      ConstVectorView   x_f,
                      ConstVectorView   ant_za,
                         const Index&   n_pol )
{
  // Calculate number of antennas/beams
  const Index n_ant = ant_za.nelem();

  // Check that the output matrix the right size
  assert(H.nrows()==x_f.nelem()*n_ant*n_pol);
  assert(H.ncols()==m_za.nelem()*x_f.nelem()*n_pol);
  
  // Check the size of the antenna diagram array of arrays and set a flag
  // if only one angle is given or if there is a complete set for each
  // angle. Initialise also flags for polarisation and frequency.
  assert(diag.nelem()==1 || diag.nelem()==n_ant);
  Index a_step = 0;
  Index p_step = 0;
  Index f_step = 0;
  if (diag.nelem()>1)
    a_step = 1;

  // Initialise variables that will store the angle, polarisation and
  // frequency grid points, so that we can check if the same antenna
  // diagram is used in succeding loops. We need the variables to
  // start with values out of the range of the different grids, therefore
  // we initialise them to their respective grid number of elements plus
  // one.
  Index a_old = n_ant;
  Index p_old = n_pol;
  Index f_old = x_f.nelem()+1;

  // Initialise temporary vectors for storing integration vector values
  // before storing them in the final sparse matrix and start looping
  // through the viewing angles.
  Vector temp(H.ncols(), 0.0);
  Vector temp_za(m_za.nelem(), 0.0);
  for (Index a=0; a<n_ant; a++) {
    Index a_this = a*a_step;

    // Check the size of this element of diag and set a flag if only one
    // polarisation is given or if there is a complete set for each
    // polarisation.
    assert(diag[a_this].nelem()==1 || diag[a_this].nelem()==n_pol);
    if (diag[a_this].nelem()>1)
      p_step = 1;

    // Loop through the polarisation antenna diagrams.
    for (Index p=0; p<n_pol; p++) {
      Index p_this = p*p_step;

      // Check the number of columns in this matrix and set flag if one
      // column is given or if there is a complete set for each frequency.
      assert((diag[a_this])[p_this].ncols()==2 ||
             (diag[a_this])[p_this].ncols()==x_f.nelem()+1);
      if ((diag[a_this])[p_this].ncols()!=2)
        f_step = 1;

      // Loop through x_f and calculate the sensor integration vector
      // for each frequency and put values in the temp vector. For this
      // we use vectorviews where the elements are separated by number
      // of frequencies in x_f.
      for (Index f=0; f<x_f.nelem(); f++) {
        Index f_this = f*f_step;

        // Check if the antenna pointer still points to the same antenna
        // diagram, if so don't recalculate the integration vector.
        // Add the angle offset of this antenna/beam.
        Vector za_rel = (diag[a_this])[p_this](joker, 0);
        za_rel += ant_za[a];
        if (a_this!=a_old || p_this!=p_old || f_this!=f_old)
          sensor_integration_vector(temp_za,
            (diag[a_this])[p_this](joker, 1+f_this),
            za_rel, m_za);

        // Now distribute the temp_za elements into temp, where they will
        // be spread with the number of frequencies. Then insert the
        // vector into the output matrix at the specific row corresponding
        // to this frequency, polarisation and viewing direction. To do
        // we first check if the same antenna diagram applies for all
        // polarisations, i.e. p_step = 0, if so insert it n_pol times.
        Index p_step_tmp = p_this;
        //Index p_tmp = p_this;
        if (p_step==0)
          p_step_tmp = n_pol-1;
        //for (p_tmp; p_tmp<=p_step_tmp; p_tmp++) {
          temp[Range(f*n_pol+p,m_za.nelem(),x_f.nelem()*n_pol)]
          //temp[Range(f*n_pol+p_tmp,m_za.nelem(),x_f.nelem()*n_pol)]
            = temp_za;
          H.insert_row(a*n_pol*x_f.nelem()+f*n_pol+p, temp);
          temp = 0.0;
        //}

        // Store antenna diagram index for this loop so that we can
        // compare it with next step.
        a_old = a;
        p_old = p;
        f_old = f;
      }
    }
  }
}

//! mixer_transfer_matrix
/*!
   Sets up the sparse matrix that models the response from sideband filtering
   and the mixer.

   The size of the transfer matrix has to be set up before calling the function
   as follows:
     nrows = f_mixer.nelem()
     ncols = f_grid.nelem()

   The returned frequencies are given in IF, so both primary and mirror band
   is converted down.

   \param   H         The mixer/sideband filter transfer matrix
   \param   f_mixer   The frequency grid of the mixer
   \param   f_grid    The original frequency grid of the spectrum
   \param   lo        The local oscillator frequency
   \param   filter    The sideband filter matrix
   \param   n_pol     The number of polarisations to consider
   \param   n_za      The number of zenith angles/antennas

   \author Mattias Ekström
   \date   2003-05-27
*/
void mixer_transfer_matrix(
              Sparse&   H,
              Vector&   f_mixer,
      ConstVectorView   f_grid,
        const Numeric   lo,
      ConstMatrixView   filter,
          const Index   n_pol,
          const Index   n_za )
{
  // Check that the sideband filter matrix at least has two columns and
  // that its frequency grid expands outside f_grid.
  assert( filter.ncols()==2 );
  assert( filter(0,0)<=f_grid[0] );
  assert( filter(filter.nrows()-1,0)>=f_grid[f_grid.nelem()-1] );

  // Check that the lo frequency is within the f_grid
  assert( lo>f_grid[0] && lo<f_grid[f_grid.nelem()-1] );

  // Find indices in f_grid where f_grid is just below and above the
  // lo frequency.
  Index i_low = 0, i_high = f_grid.nelem()-1, i_mean;
  while (i_high-i_low>1)
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
    i_high++;

  // Calculate the cut-off limits to assure that all frequencies in IF are
  // computable, i.e. possible to interpolate, in RF.
  const Numeric lim_low = max(lo-f_grid[i_low], f_grid[i_high]-lo);
  const Numeric lim_high = min(lo-f_grid[0], f_grid[f_grid.nelem()-1]-lo);

  // Convert sidebands to IF and use std::list to make a unique sorted
  // vector, this sorted vector is f_mixer.
  list<Numeric> l_mixer;
  for (Index i=0; i<f_grid.nelem(); i++)
  {
    if (fabs(f_grid[i]-lo)>=lim_low && fabs(f_grid[i]-lo)<=lim_high)
      l_mixer.push_back(fabs(f_grid[i]-lo));
  }
  l_mixer.unique();
  l_mixer.sort();
  f_mixer.resize((Index) l_mixer.size());
  Index i=0;
  for (list<Numeric>::iterator li=l_mixer.begin(); li != l_mixer.end(); li++)
  {
    f_mixer[i] = *li;
    i++;
  }

  // Now we know the final size of H, so resize it
  H.resize( f_mixer.nelem()*n_pol*n_za, f_grid.nelem()*n_pol*n_za );

  // Calculate the sensor summation vector and insert the values in the
  // final matrix taking number of polarisations and zenith angles into
  // account.
  Vector row_temp(f_grid.nelem());
  Vector row_final(f_grid.nelem()*n_pol*n_za);
  for (Index i=0; i<f_mixer.nelem(); i++) {
    sensor_summation_vector(row_temp, f_mixer[i], f_grid, lo, filter );
    // Loop over number of polarisations
    for (Index p=0; p<n_pol; p++)
    {
      // Loop over number of zenith angles/antennas
      for (Index a=0; a<n_za; a++)
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

//! scale_antenna_diagram
/*!
   Scales a Gaussian antenna diagram for a reference frequency to match
   the new frequency.

   \param   s       The scaled antenna diagram
   \param   srm     The antenna diagram matrix
   \param   f_ref   The reference frequency
   \param   f_new   The new frequency

   \author Mattias Ekström
   \date   2003-08-14
*/
void scale_antenna_diagram(
             VectorView   sc,
        ConstMatrixView   srm,
         const Numeric&   f_ref,
         const Numeric&   f_new )
{
  // Check output vector size
  assert( sc.nelem()==srm.ncols() );

  // Calculate the scale factor
  Numeric s = f_new / f_ref;

  // Perform the scaling, by scaling the gain values
  for (Index i=0; i<srm.nrows(); i++) {
    sc[i]=pow(srm(i,1), s);
  }
}

//! sensor_integration_vector
/*!
   Calculates the (row) vector that multiplied with an unknown
   (column) vector approximates the integral of the product
   between the functions represented by the two vectors.

   E.g. h*g = integral( f(x)*g(x) dx )

   \param   h      The multiplication (row) vector.
   \param   f      The values of function f(x).
   \param   x_ftot The grid points of function f(x).
   \param   x_g    The grid points of function g(x).

   \author Mattias Ekström
   \date   2003-02-13
*/
void sensor_integration_vector(
           VectorView   h,
      ConstVectorView   f,
      ConstVectorView   x_ftot,
      ConstVectorView   x_g )
{
  //Check that vectors are sorted, ascending (no descending?)

  //Assert that h has the right size
  assert( h.nelem() == x_g.nelem() );

  //Find x_f points that lies outside the scope of x_g and remove them
  Index i1_f = 0, i2_f = x_ftot.nelem()-1;
  while( x_ftot[i1_f] < x_g[0] ) {
    i1_f++;
  }
  while( x_ftot[i2_f] > x_g[x_g.nelem()-1] ) {
    i2_f--;
  }
  Vector x_f = x_ftot[Range(i1_f, i2_f-i1_f+1)];

  //Create a reference grid vector, x_ref that containing the values of
  //x_f and x_g strictly sorted.
  list<Numeric> l_x;
  for (Index i=0; i<x_f.nelem(); i++)
  	l_x.push_back(x_f[i]);
  for (Index i=0; i<x_g.nelem(); i++) {
    if( x_g[i] > l_x.front() && x_g[i]<l_x.back() )
      l_x.push_back(x_g[i]);
  }

  l_x.sort();
  l_x.unique();

  Vector x_ref(l_x.size());
  Index i=0;
  for (list<Numeric>::iterator li=l_x.begin(); li != l_x.end(); li++) {
	x_ref[i] = *li;
	i++;
  }

  //Initiate output vector, with equal size as x_g, with zeros.
  //Start calculations
  h = 0.0;
  Index i_f = 0;
  Index i_g = 0;
  //i = 0;
  Numeric dx,a0,b0,c0,a1,b1,c1,x3,x2,x1;
  //while( i_g < x_g.nelem() && i_f < x_f.nelem() ) {
  for( Index i=0; i<x_ref.nelem()-1; i++ ) {
    //Find for what index in x_g (which is the same as for h) and f
    //calculation corresponds to
    while( x_g[i_g+1] <= x_ref[i] ) {
      i_g++;
    }
    while( x_f[i_f+1] <= x_ref[i] ) {
     i_f++;
	 //cout << "i_f: " << i_f << "\n";
    }

    //If x_ref[i] is out of x_f's range then that part of the integral
    //is set to 0, so no calculations will be done
    if( x_ref[i] >= x_f[0] && x_ref[i] < x_f[x_f.nelem()-1] ) {
      //Product of steps in x_f and x_g
      dx = (x_f[i_f+1] - x_f[i_f]) * (x_g[i_g+1] - x_g[i_g]);

      //Calculate a, b and c coefficients; h[i]=ax^3+bx^2+cx
      a0 = (f[i_f] - f[i_f+1]) / 3;
      b0 = (-f[i_f]*(x_g[i_g+1]+x_f[i_f+1])+f[i_f+1]*(x_g[i_g+1]+x_f[i_f]))/2;
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
    //i++;
  }

  //Normalize h.
  h /= h.sum();
}

//! sensor_summation_vector
/*!
   Constructs the (row) vector that sums components of another (column) vector.
   These (row) vectors are a used to set up the response matrix for mixer and
   sideband filter.

   The sideband filter respone should already be normalised before calling this
   function and its relative grid should cover the whole frequency grid.

   \param   h      The summation (row) vector.
   \param   f      The frequency in the IF band.
   \param   f_grid The grid points of function f(x).
   \param   lo     The local oscillator frequency
   \param   sfrm   The sideband filter response matrix.

   \author Mattias Ekström
   \date   2003-05-26
*/
void sensor_summation_vector(
           VectorView   h,
        const Numeric   f,
      ConstVectorView   f_grid,
        const Numeric   lo,
      ConstMatrixView   sfrm )
{
  //Check that the (row) vector has the right dimensions
  assert( h.nelem() == f_grid.nelem() );

  //Check that sfrm has the right size and that it covers f_grid
  assert( sfrm.ncols()==2 );
  assert( sfrm(0,0)<=f_grid[0] );
  assert( sfrm(sfrm.nrows()-1,0)>=f_grid[f_grid.nelem()-1] );

  //Calculate the upper and lower sideband frequencies
  const Numeric f_low = lo - f;
  const Numeric f_upp = lo + f;

  //Check that the sideband frequencies lies within the frequency grid
  assert( f_low >= f_grid[0] && f_low <= f_grid[f_grid.nelem()-1] );
  assert( f_upp >= f_grid[0] && f_upp <= f_grid[f_grid.nelem()-1] );

  //Interpolate the intensity and sideband filter response for the upper
  //frequency. This should work even if the filter grid and the frequency
  //grid are the same.
  ArrayOfGridPos gp_upp(1), gp_upp_filt(1);
  gridpos( gp_upp, f_grid, f_upp);
  gridpos( gp_upp_filt, sfrm(joker,0), f_upp);

  Matrix itw_upp(1,2);
  interpweights( itw_upp, gp_upp_filt);

  Numeric filt_upp;
  interp( filt_upp, itw_upp, sfrm(joker,1), gp_upp_filt);

  //Interpolate the intensity and sideband filter response for the lower
  //frequency. Since different grids are used for the intensity and the
  //filter, different gridpos has to be set up. Also since we don't know
  //the intensity values, only gridpos are calculated.
  ArrayOfGridPos gp_low(1), gp_low_filt(1);
  gridpos( gp_low, f_grid, f_low);
  gridpos( gp_low_filt, sfrm(joker,0), f_low);

  Matrix itw_filt(1,2);
  interpweights( itw_filt, gp_low_filt);

  Numeric filt_low;
  interp( filt_low, itw_filt, sfrm(joker,1), gp_low_filt);

  //Initialise h at zero and store calculated weighting components
  h = 0.0;
  Numeric filt_sum = filt_upp + filt_low;
  h[gp_upp[0].idx] = filt_upp/filt_sum * gp_upp[0].fd[1];
  h[gp_upp[0].idx+1] = filt_upp/filt_sum * gp_upp[0].fd[0];
  h[gp_low[0].idx] = filt_low/filt_sum * gp_low[0].fd[1];
  h[gp_low[0].idx+1] = filt_low/filt_sum * gp_low[0].fd[0];
}

//! spectrometer_transfer_matrix
/*!
   Constructs the sparse matrix that multiplied with the spectral values
   gives the spectra from the spectrometer.

   The size of the spectrometer transfer matrix has to be set by the calling
   function, and the number of rows must match the number of spectrometer
   channels times the number of viewing angles and polarisations while the
   number of columns must match the number of frequency grid points times
   the number of viewing angles and polarisations.

   The spectrometer response ArrayOfMatrix decribes how the spectrometer
   behaves on a relative frequency grid for each polarisation. Each element
   corresponds to a polarisation and in the matrices the first column
   contains the relative grid points and the following holds the response
   values. Both the array and the individual matrices work on the
   single/full principle, which means that either only one element/column
   is given and used for all polarisations/channel frequencies or every
   polarisation/frequency is given its individual response.

   \param   H             The transfer matrix.
   \param   ch_response   The spectrometer response matrix.
   \param   ch_f          The spectrometer channel centre frequencies.
   \param   sensor_f      The frequency grid points.
   \param   n_za          The number of viewing angles
   \param   n_pol         The number of polarisations

   \author Mattias Ekström
   \date   2003-08-26
*/
void spectrometer_transfer_matrix( Sparse&                H,
                                   const ArrayOfMatrix&   ch_response,
                                   ConstVectorView        ch_f,
                                   ConstVectorView        sensor_f,
                                   const Index&           n_za,
                                   const Index&           n_pol)
{
  // Check that the transfer matrix has the right size
  assert (H.nrows()==ch_f.nelem()*n_za*n_pol);
  assert (H.ncols()==sensor_f.nelem()*n_za*n_pol);

  // Check that the channel response has the right number of elements
  // and if the same response will be used for all frequencies or if
  // there is one response for each frequency.
  assert (ch_response.nelem()==1 || ch_response.nelem()==n_pol);
  Index pol_single = 0;
  if (ch_response.nelem()==1)
    pol_single = 1;

  // Allocate memory for temporary vectors, temp_long is used to store the
  // expanded result from sensor_integration_vector before inserting them in
  // the transfer matrix. The second vector, temp, is used for the output
  // from sensor_integration_vector.
  Vector temp_long(sensor_f.nelem()*n_za*n_pol, 0.0);
  Vector temp(sensor_f.nelem(), 0.0);

  // Loop over elements in ch_response and calculate responses, if the same
  // response applies to all polarisations run inner loop once and
  // distribute values. Initialise a temporary vector for the frequency
  // grid shifted by channel frequency
  for (Index p=0;p<ch_response.nelem();p++) {
    Vector ch_response_f(ch_response[p].nrows());

    //Check if matrix has one frequency column or one for every channel
    // frequency
    assert (ch_response[p].ncols()==2 ||
            ch_response[p].ncols()==ch_f.nelem()+1);
    Index freq_full = 1;
    if (ch_response[p].ncols()==2)
      freq_full = 0;

    //Calculate the sensor integration vector and put values in the temporary
    //vector, then copy vector to the transfer matrix
    for (Index i=0; i<ch_f.nelem(); i++) {
      //The spectrometer response is shifted for each centre frequency step
      ch_response_f = ch_response[p](joker,0);
      ch_response_f += ch_f[i];

      // Call sensor_integration_vector and store it in the temp vector
      sensor_integration_vector(temp,ch_response[p](joker,1+i*freq_full),
        ch_response_f,sensor_f);

      // Loop over the viewing angles, here we only need on computation
      // which is then copied to all angles.
      for (Index za=0;za<n_za;za++) {
        // Here we loop if pol_single == 1, that is if the same response
        // apply to all polarisations. If pol_single == 0 the outermost
        // loop over polarisations take care of this.
        for (Index p_tmp=0;p_tmp<=(n_pol-1)*pol_single;p_tmp++) {
          // Get the current polarisation index, this works since
          // either we are looping over p_tmp or pol_single is zero
          // and we are looping over p.
          Index p_this = p_tmp+p*(1-pol_single);

          // Distribute the compact temp vector into temp_long
          temp_long[Range(n_pol*sensor_f.nelem()*za+p_this,
            sensor_f.nelem(),n_pol)] = temp;

          // Insert temp_long into H at the correct row
          H.insert_row(n_pol*ch_f.nelem()*za+i*n_pol+p_this,temp_long);

          // Reset temp_long to zero.
          temp_long = 0.0;
        }
      }
    }
  }
}
