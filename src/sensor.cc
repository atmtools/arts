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
   measurement block grid. The number of sensor response matrix rows
   don't need to match the number of frequency grid points.

   FIXME: The antenna diagram values could be set up and scaled using the
   antenna_diagram_gaussian and scale_antenna_diagram functions.

   \param   H      The antenna transfer matrix.
   \param   m_za   The measurement block grid of zenith angles.
   \param   srm	   The sensor response matrix, i.e. the antenna diagram
   \param   x_f    The frequency grid points.

   \author Mattias Ekström
   \date   2003-04-09
*/
void antenna_transfer_matrix(
           SparseView   H,
      ConstVectorView   m_za,
      ConstMatrixView   srm,
      ConstVectorView   x_f )
{
  //Assert that the transfer matrix and the sensor response matrix has the 
  //right size
  assert( H.nrows()==x_f.nelem() && H.ncols()==m_za.nelem()*x_f.nelem() );
  assert( srm.ncols()==2 );

  //FIXME: Allocate a temporary vector to keep values before sorting them into the
  //final sparse matrix, could (should) be fixed with a SparseView(range, range) operator
  Vector temp(m_za.nelem()*x_f.nelem(), 0.0);

  //Calculate the sensor integration vector and put values in the temp vector
  //FIXME: Scale the antenna diagram?? If so, do it here.
  for (Index i=0; i<x_f.nelem(); i++) {
    sensor_integration_vector( temp[Range(i, m_za.nelem(), x_f.nelem())],
      srm(Range(joker),1), srm(Range(joker),0), m_za);
  }

  //Copy values to the antenna matrix
  for (Index j=0; j<m_za.nelem(); j++) {
    for (Index i=0; i<x_f.nelem(); i++) {
      if (temp[i+j*x_f.nelem()]!=0)
        H.rw(i, i+j*x_f.nelem()) = temp[i+j*x_f.nelem()];
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

   The primary frequency band is specified to be upper or lower, that means the
   frequencies in the frequency grid above or below the LO frequency. The image
   band when mirrored has to cover the whole primary band.

   \param   H         The mixer/sideband filter transfer matrix
   \param   f_mixer   The frequency grid of the mixer
   \param   f_grid    The original frequency grid of the spectrum
   \param   is_upper  Specifies if the primary band is the upper or not.
   \param   lo		  The local oscillator frequency
   \param   sfrm      The sideband filter response matrix

   \author Mattias Ekström
   \date   2003-05-27
*/
void mixer_transfer_matrix(
               Sparse&   H,
              Vector&   f_mixer,
      ConstVectorView   f_grid,
		   const bool   is_upper,
	    const Numeric   lo,
      ConstMatrixView   sfrm )
{
  //Check that the sideband filter matrix at least has two columns
  assert( sfrm.ncols()==2 );

  //Get min and max indices for the primary band
  GridPos gp_lo;
  gridpos( gp_lo, f_grid, lo);
  Index lo_idx = gp_lo.idx;

  Index f_low, f_high;
  if( is_upper ) {
  	f_low = lo_idx+1;
	f_high = f_grid.nelem()-1;
  } else {
    f_low = 0;
	if (gp_lo.fd[0]==0) {
	  f_high = lo_idx-1;
	} else {
	  f_high = lo_idx;
	}
  }

  //Check that the LO and sidebands span are consistent?
  assert( f_high-f_low>0 );

  //Since we want to interpolate the image frequencies, we want the image band
  //to cover the primary band.
  Vector f_im(2);
  f_im[0] = 2*lo-f_grid[f_high];
  f_im[1] = 2*lo-f_grid[f_low];

  //Find the indices for f_im_min and f_im_max in f_grid
  ArrayOfGridPos gp_im(2);
  gridpos( gp_im, f_grid, f_im);
  Index f_im_low = gp_im[0].idx+1;
  Index f_im_high = gp_im[1].idx;

  //Check that image band indices are consistent
  assert( f_im_low >= 0 );		 //FIXME: unnecessary? does gridpos this?
  assert( f_im_high-f_im_low > 0 );
  assert( f_im_high < f_grid.nelem() );

  //Set up f_mixer with the correct size and with the frequencies from both
  //the primary band and the image band strictly sorted.
  Numeric f_im_im;
  list<Numeric> l_f;
  for (Index i=f_low; i<=f_high; i++)
  	l_f.push_back(f_grid[i]);
  for (Index i=f_im_low; i<=f_im_high; i++) {
    f_im_im = 2*lo-f_grid[i];
	l_f.push_back(f_im_im);
  }

  l_f.sort();
  l_f.unique();

  f_mixer.resize((Index) l_f.size());
  Index i=0;
  for (list<Numeric>::iterator li=l_f.begin(); li != l_f.end(); li++) {
	f_mixer[i] = *li;
	i++;
  }

  //Now we know the final size of H, so resize it
  H.resize( f_mixer.nelem(), f_grid.nelem() );

  //FIXME: Allocate a temporary vector to keep values before sorting them into
  //the final sparse matrix, could (should) be fixed with a
  //SparseView(range, range) operator
  Vector temp( f_mixer.nelem()*f_grid.nelem(), 0.0);

  //Calculate the sensor summation vector and put values in the temp vector
  for (Index i=0; i<f_mixer.nelem(); i++) {
	sensor_summation_vector(
	  temp[Range(i*f_grid.nelem(), f_grid.nelem())],
      f_mixer[i], f_grid, lo, sfrm );
  }

  //Copy values to the sensor matrix
  for (Index j=0; j<f_grid.nelem(); j++) {
    for (Index i=0; i<f_mixer.nelem(); i++) {
      if (temp[i+j*f_mixer.nelem()]!=0)
        H.rw(i, j) = temp[i+j*f_mixer.nelem()];
    }
  }
}

//! scale_antenna_diagram
/*!
   Scales a Gaussian antenna diagram for a reference frequency to match
   the new frequency.

   \return          The scaled antenna diagram
   \param   srm     The antenna diagram matrix
   \param   f_ref   The reference frequency
   \param   f_new   The new frequency

   \author Mattias Ekström
   \date   2003-03-11
*/
Matrix scale_antenna_diagram(
        ConstMatrixView   srm,
         const Numeric&   f_ref,
         const Numeric&   f_new )
{
  //Initialise new vector
  Matrix srm_new = srm;

  //Get scale factor
  Numeric s = f_new / f_ref;

  //Scale
  for (Index i=0; i<srm.nrows(); i++) {
    srm_new(i,1)=pow(srm(i,1), s);
  }

  return srm_new;
}

//! sensor_integration_vector
/*!
   Calculates the (row) vector that multiplied with an unknown
   (column) vector approximates the integral of the product
   between the functions represented by the two vectors.

   E.g. h*g = integral( f(x)*g(x) dx )

   \param   h      The multiplication (row) vector.
   \param   f      The values of function f(x).
   \param   x_f    The grid points of function f(x).
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
   \param   f      The primary frequency.
   \param   f_grid The grid points of function f(x).
   \param   lo	   The local oscillator frequency
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

  //Check that sfrm has the right size
  assert( sfrm.ncols()==2 );
  //FIXME:assert( sfrm.nrows()==f_grid.nelem() );

  //Check that the primary frequency lies within the frequency grid
  assert( f >= f_grid[0] );
  assert( f <= f_grid[f_grid.nelem()-1] );

  //Calculate the image frequency
  const Numeric f_im = 2*lo-f;

  //Check that the image frequency exist within the frequency span assuming
  //that f_grid is sorted and each value is unique.
  //This should already have been made sure by the calling function.
  assert( f_im >= f_grid[0] );
  assert( f_im <= f_grid[f_grid.nelem()-1] );

  //Interpolate the intensity and sideband filter response for the primary
  //frequency. This should work even if the filter grid and the frequency grid
  //are the same.
  ArrayOfGridPos gp_prim(1), gp_prim_filt(1);
  gridpos( gp_prim, f_grid, f);
  gridpos( gp_prim_filt, sfrm(joker,0), f);

  Matrix itw_prim(1,2);
  interpweights( itw_prim, gp_prim_filt);

  Numeric filt_prim;
  interp( filt_prim, itw_prim, sfrm(joker,1), gp_prim_filt);

  //Interpolate the intensity and sideband filter response for the image
  //frequency. Since different grids are used for the intensity and the
  //filter, different gridpos has to be set up. Also since we don't know
  //the intensity values, only gridpos are calculated.
  ArrayOfGridPos gp_im(1), gp_im_filt(1);
  gridpos( gp_im, f_grid, f_im);
  gridpos( gp_im_filt, sfrm(joker,0), f_im);

  Matrix itw_filt(1,2);
  interpweights( itw_filt, gp_im_filt);

  Numeric filt_im;
  interp( filt_im, itw_filt, sfrm(joker,1), gp_im_filt);

  //Initialise h at zero and store calculated weighting components
  h = 0.0;
  Numeric filt_sum = filt_prim + filt_im;
  h[gp_prim[0].idx] = filt_prim/filt_sum * gp_prim[0].fd[1];
  h[gp_prim[0].idx+1] = filt_prim/filt_sum * gp_prim[0].fd[0];
  h[gp_im[0].idx] = filt_im/filt_sum * gp_im[0].fd[1];
  h[gp_im[0].idx+1] = filt_im/filt_sum * gp_im[0].fd[0];
}

//! spectrometer_transfer_matrix
/*!
   Constructs the sparse matrix that multiplied with the spectral values.
   There are no dependency on line-of-sights, therefor the antenna transfer
   matrix has to be applied before the spectrometer.

   The size of the spectrometer transfer matrix has to be set by the calling
   function, and the number of rows must match the number of spectrometer
   channels while the number of columns must match the number of frequency
   grid points.

   The spectrometer response matrix decribes how the spectrometer behaves on a
   relative frequency grid. The first column contains the relative grid points
   and the second holds the response values.
   NB:(FIXME?) The same spectrometer response is used for all channels, it is only
   shifted to the corresponding centre frequency.

   \param   H      The transfer matrix.
   \param   srm    The spectrometer response matrix.
   \param   x_s    The spectrometer channel centre frequencies.
   \param   x_f    The frequency grid points.

   \author Mattias Ekström
   \date   2003-05-12
*/
void spectrometer_transfer_matrix(
           SparseView   H,
      ConstMatrixView   srm,
      ConstVectorView   x_s,
      ConstVectorView   x_f )
{
  //Assert that the transfer matrix and the sensor response matrix has the
  //right size
  assert( H.nrows()==x_s.nelem() && H.ncols()==x_f.nelem() );
  assert( srm.ncols()==2 );
  
  //Allocate memory for temporary vectors
  Vector temp(x_f.nelem()*x_s.nelem(), 0.0);
  Vector x_srm(srm.nrows());

  //Calculate the sensor integration vector and put values in the temporary
  //vector, then copy vector to the transfer matrix
  for (Index i=0; i<x_s.nelem(); i++) {
    //The spectrometer response is shifted for each centre frequency step
	x_srm = srm(joker,0);
	x_srm += x_s[i];
    sensor_integration_vector( temp[Range(i*x_f.nelem(), x_f.nelem())],
	                           srm(joker,1), x_srm, x_f);
	//cout << "temp:\n" << temp[Range(i*x_f.nelem(), x_f.nelem())] << "\n";
	//cout << "x_srm:\n" << x_srm << "\n";
  }

  //FIXME: Should be fixed by a SparseView(Range, Range) -> VectorView
  //This might take more time than needed
  for (Index j=0; j<x_f.nelem(); j++) {
    for (Index i=0; i<x_s.nelem(); i++) {
      if (temp[i+j*x_s.nelem()]!=0)
	    H.rw(i,j) = temp[i+j*x_s.nelem()];
	}
  }
}

