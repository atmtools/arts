/* Copyright (C)  Patrick Eriksson <patrick@rss.chalmers.se>

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
#include <stdexcept>
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "mystring.h"
//#include "logic.h"
//#include "poly_roots.h"
//#include "special_interp.h"

  extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! sensor_integration_vector
/*!
   Calculates the (row) vector that multiplied that multiplied with an
   unknown (column) vector approximates the integral of the product
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

  //Find x_f points that lies outside the scope of x_g and remove them
  Index i1_f = 0, i2_f = x_ftot.nelem()-1;
  while( x_ftot[i1_f] < x_g[0] ) {
    i1_f++;
  }
  while( x_ftot[i2_f] > x_g[x_g.nelem()-1] ) {
    i2_f--;
  }
  Vector x_f = x_ftot[Range(i1_f, i2_f-i1_f+1)];

  //Create a reference grid vector containing all x_f and x_g
  //strictly sorted by adding the in a sorted way.
  Vector x_reftot( x_f.nelem() + x_g.nelem() );

  Index i_f = 0, i_g = 0, i=0;
  while( i_f<x_f.nelem() || i_g<x_g.nelem() ) {
    if (x_f[i_f] < x_g[i_g]) {
      x_reftot[i] = x_f[i_f];
      i_f++;
    } else if (x_f[i_f] > x_g[i_g]) {
      x_reftot[i] = x_g[i_g];
      i_g++;
    } else {
      // x_f and x_g are equal
      x_reftot[i] = x_f[i_f];
      i_f++;
      i_g++;
    }
  i++;
  }
  Vector x_ref = x_reftot[Range(0,i)];

  //Initiate output vector, with equal size as x_g, with zeros.
  //Start calculations
  h = 0.0;
  i_f = 0;
  i_g = 0;
  Numeric dx,a0,b0,c0,a1,b1,c1,x3,x2,x1;
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
  }
}

//! antenna_transfer_matrix
/*!
   Constructs the matrix that multiplied with the spectral values for
   one or several line-of-sights models the antenna transfer matrix.
   The matrix it built up of spaced row vectors, to match the format
   of the spectral values.

   The size of the antenna transfer matrix has to be set by the calling
   function, and it must be set as:
    nrows = x_f.nelem()
    ncols = x_f.nelem() * m_za.nelem()

   The number of line-of-sights is determined by the length of the
   measurement block grid. The number of antenna diagram columns and
   grid points must match, but they don't need to match the number
   of frequency grid points. However, the antenna diagram row number
   must either be equal to one, then we have the same response for all
   frequencies, or it must equal the length of the frequency grid,
   where we have different responses for each frequency.

   \param   Hb     The antenna transfer matrix.
   \param   m_za   The measurement block grid of zenith angles.
   \param   a      The antenna diagram values.
   \param   x_a    The antenna diagram grid points.
   \param   x_f    The frequency grid points.

   \author Mattias Ekström
   \date   2003-03-06
*/
void antenna_transfer_matrix(
           MatrixView   Hb,
      ConstVectorView   m_za,
      ConstMatrixView   a,
      ConstVectorView   x_a,
      ConstVectorView   x_f )
{
  //FIXME: This part below could be exchanged with an interpolation to
  //cover all cases:
  //Check dimension of antenna diagram
  Index di_a;
  if ( a.nrows() == 1 ) {
    //If a only contains one row, the we use the same antenna response
    //for all frequencies
    di_a = 0;
  } else if ( a.nrows() == x_f.nelem() ) {
    //If a has equally many rows as grid points in the frequency grid
    //use the corresponding antenna diagram for each frequency.
    di_a = 1;
  } else {
    //FIXME: Should this be an assert instead?
    ostringstream os;
    os << "The number of rows in the antenna diagram matrix must equal\n"
       << "one or the length of the frequency grid.\n";
    throw runtime_error(os.str());
  }

  //Assert that the transfer matrix has the right size, and set it to zero
  assert( Hb.nrows()==x_f.nelem() && Hb.ncols()==m_za.nelem()*x_f.nelem() );
  Hb = 0.0;

  //Loop over frequencies to calculate the sensor intergration vector
  //and put values in the antenna matrix.
  Index i_a=0;
  for (Index i=0; i<x_f.nelem(); i++) {
    sensor_integration_vector( Hb(i, Range(i, m_za.nelem(), x_f.nelem())),
      a(i_a, Range(joker)), x_a, m_za);

    i_a += di_a;
  }
}

//! antenna_diagram_gaussian
/*!
   Sets up an vector containing a standardized Gaussian antenna diagram,
   described for a certain frequency. The function is called with the 
   half-power beamwidth that determines the shape of the curve for the
   reference frequency, and a grid that sets up the antenna diagram.

   \param   a       The antenna diagram vector.
   \param   a_grid  The antenna diagram grid of angles.
   \param   theta   The antenna average width

   \author Mattias Ekström
   \date   2003-03-11
*/
void antenna_diagram_gaussian(
           VectorView   a,
      ConstVectorView   a_grid,
       const Numeric&   theta )
{
  //Assert that a has the right size
  assert( a.nelem()==a_grid.nelem() );

  //Initialise variables
  Numeric ln2 = log(2.0);

  //Loop over grid points to calculate antenna diagram
  for (Index i=0; i<a_grid.nelem(); i++) {
    a[i]=exp(-4*ln2*pow(a_grid[i]*DEG2RAD/theta,2));
  }
}

//! scale_antenna_diagram
/*!
   Scales a Gaussian antenna diagram for a reference frequency to match 
   the new frequency.

   \return          The scaled antenna diagram
   \param   a       The antenna diagram vector
   \param   f_ref   The reference frequency
   \param   f_new   The new frequency
   
   \author Mattias Ekström
   \date   2003-03-11
*/
Vector scale_antenna_diagram(
        ConstVectorView   a,
         const Numeric&   f_ref,
         const Numeric&   f_new )
{
  //Initialize new vector
  Vector a_new( a.nelem() );

  //Get scale factor
  Numeric s = f_new / f_ref;

  //Scale
  for (Index i=0; i<a.nelem(); i++) {
    a_new[i]=pow(a[i], s);
  }

  return a_new;
}

//! antenna_diagram_gaussian
/*!
   Sets up an matrix containing a standardized Gaussian antenna diagram.
   Each row describes the diagram for a certain frequency at the
   measurement block grid points.

   \param   a      The antenna diagram matrix.
   \param   m_za   The measurement block grid of zenith angles.
   \param   f      The frequencys.
   \param   alpha  The antenna efficiency typically between 0.9 and 1.4
   \param   D      The antenna dimension

   \author Mattias Ekström
   \date   2003-03-07
*/
/*
void antenna_diagram_gaussian(
           MatrixView   a,
      ConstVectorView   m_za,
      ConstVectorView   f,
       const Numeric&   alpha,
       const Numeric&   D )
{
  //Assert that a has the right size
  assert( a.nrows()==f.nelem() && a.ncols()==m_za.nelem() );

  //Initialise variables
  Numeric ln2 = log(2.0);

  //Loop over frequencies and grid points to calculate antenna diagram
  for (Index i=0; i<f.nelem(); i++) {
    for (Index j=0; i<m_za.nelem(); j++) {
      NumericPrint(i, "i");
      NumericPrint(j, "j");
      a(i,j)=exp(-4*ln2*pow(m_za[j]*DEG2RAD*D*f[i]/(alpha*SPEED_OF_LIGHT),2));
    }
  }

}
*/
