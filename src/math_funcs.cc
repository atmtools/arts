/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>

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



/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
   \file   math_funcs.cc
   \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   \date   2000-09-18 

   Contains basic mathematical functions.
*/



/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include "array.h"
#include "math_funcs.h"
#include "logic.h"
extern const Numeric DEG2RAD;



/*****************************************************************************
 *** The functions (in alphabetical order)
 *****************************************************************************/

//! integer_div
/*! 
    Performs an integer division.

    The function asserts that the reminder of the division x/y is 0.

    \return      The quotient
    \param   x   Nominator
    \param   y   Denominator

    \author Patrick Eriksson 
    \date   2002-08-11
*/
Index integer_div( const Index& x, const Index& y )
{
  assert( is_multiple( x, y ) );
  return x/y;
}



//! last
/*! 
    Returns the last value of a vector.

    \return      The last value of x.
    \param   x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last( ConstVectorView x )
{
  assert( x.nelem() > 0 );
  return x[x.nelem()-1]; 
}



//! last
/*! 
    Returns the last value of an index array.

    \return      The last value of x.
    \param   x   An index array.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Index last( const ArrayOfIndex& x )
{
  assert( x.nelem() > 0 );
  return x[x.nelem()-1]; 
}



//! linspace
/*! 
    Linearly spaced vector with specified spacing. 

    The first element of x is always start. The next value is start+step etc.
    Note that the last value can deviate from stop.
    The step can be both positive and negative. 
    (in Matlab notation: start:step:stop)

    Size of result is adjusted within this function!

    \param    x       Output: linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    step    distance between values in x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void linspace(                      
              Vector&     x,           
              const Numeric     start,    
              const Numeric     stop,        
              const Numeric     step )
{
  Index n = (Index) floor( (stop-start)/step ) + 1;
  if ( n<1 )
    n=1;
  x.resize(n);
  for ( Index i=0; i<n; i++ )
    x[i] = start + i*step;
}



//! nlinspace
/*! 
    Linearly spaced vector with specified length. 

    Returns a vector equally and linearly spaced between start and stop 
    of length n. (equals the Matlab function linspace)

    The length must be > 1.

    \param    x       Output: linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void nlinspace(         
               Vector&     x, 
               const Numeric     start,     
               const Numeric     stop,        
               const Index       n )
{
  assert( 1<n );                // Number of points must be greatere 1.
  x.resize(n);
  Numeric step = (stop-start)/(n-1) ;
  for ( Index i=0; i<n; i++ )
    x[i] = start + i*step;
}



//! nlogspace
/*! 
    Logarithmically spaced vector with specified length. 

    Returns a vector logarithmically spaced vector between start and 
    stop of length n (equals the Matlab function logspace)

    The length must be > 1.

    \param    x       Output: logarithmically spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void nlogspace(         
               Vector&     x, 
               const Numeric     start,     
               const Numeric     stop,        
               const Index         n )
{
  // Number of points must be greatere 1:
  assert( 1<n );        
  // Only positive numbers are allowed for start and stop:
  assert( 0<start );
  assert( 0<stop );

  x.resize(n);
  Numeric a = log(start);
  Numeric step = (log(stop)-a)/(n-1);
  x[0] = start;
  for ( Index i=1; i<n-1; i++ )
    x[i] = exp(a + i*step);
  x[n-1] = stop;
}



//! nlogspace
/*! 
    Logarithmically spaced vector with specified length (return version). 

    Returns a vector logarithmically spaced vector between start and 
    stop of length n (equals the Matlab function logspace)

    The length must be > 1.

    \return   x       logarithmically spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
Vector nlogspace(  
                 const Numeric start, 
                 const Numeric stop,  
                 const Index     n )
{
  Vector x;
  nlogspace( x, start, stop, n );
  return x; 
}                     

//! AngIntegrate_trapezoid
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    \param Integrand The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid 
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid(MatrixView Integrand,
                               ConstVectorView za_grid,
                               ConstVectorView aa_grid)
{

  Index n = za_grid.nelem();
  Index m = aa_grid.nelem();
  Vector res1(n);
  assert (is_size(Integrand, n, m));
  
  for (Index i = 0; i < n ; ++i)
    {
      res1[i] = 0.0;
      
      for (Index j = 0; j < m - 1; ++j)
        {
          res1[i] +=  0.5 * DEG2RAD * (Integrand(i, j) + Integrand(i, j + 1)) *
            (aa_grid[j + 1] - aa_grid[j]) * sin(za_grid[i] * DEG2RAD);
        }
    }
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i)
    {
      res += 0.5 * DEG2RAD * (res1[i] + res1[i + 1]) * 
        (za_grid[i + 1] - za_grid[i]);
    }
  
  //cout<<res<<"\n";
  return res;
}



//! sign
/*! 
    Returns the sign of a numeric value.

    The function returns 1 if the value is greater than zero, 0 if it 
    equals zero and -1 if it is less than zero.

    \return      The sign of x (see above).
    \param   x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric sign( const Numeric& x )
{
  if( x < 0 )
    return -1.0;
  else if( x == 0 )
    return 0.0;
  else
    return 1.0;
}
