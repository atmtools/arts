/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2000-09-18 

   Contains basic mathematical functions.
*/



/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include <iostream>
#include <cmath>
#include <stdexcept>
#include "array.h"
#include "math_funcs.h"
#include "logic.h"
#include "mystring.h"

extern const Numeric DEG2RAD;
extern const Numeric PI;



/*****************************************************************************
 *** The functions (in alphabetical order)
 *****************************************************************************/

//! fac
/*!
    Calculates the factorial.

    The function asserts that n must be >= 0

    \return      The factorial
    \param   n   Nominator

    \author Oliver Lemke
    \date   2003-08-15
*/
Numeric fac(const Index n)
{
  Numeric sum;

  if (n == 0) return (1.0);

  sum = 1.0;
  for (Index i = 1; i <= n; i++)
    sum *= Numeric(i);

  return(sum);
}


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



//! Lagrange Interpolation (internal function).
/*! 
  This function calculates the Lagrange interpolation of four interpolation 
  points as described in 
  <a href="http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html">
  Lagrange Interpolating Polynomial</a>.<br>
  The input are the four x-axis values [x0,x1,x2,x3] and their associated 
  y-axis values [y0,y1,y2,y3]. Furthermore the x-axis point "a" at which the 
  interpolation should be calculated must be given as input. NOTE that 
  the relation x2 =< x < x3 MUST hold!

  \param x     x-vector with four elements [x0,x1,x2,x3]
  \param y     y-vector with four elements: yj = y(xj), j=0,1,2,3
  \param a     interpolation point on the x-axis with x1 =< a < x2 

  \return FIXME

  \author Thomas Kuhn
  \date   2003-11-25
*/

Numeric LagrangeInterpol4( ConstVectorView x,
                           ConstVectorView y,
                           const Numeric a)
{
  // lowermost grid spacing on x-axis
  const Numeric Dlimit = 1.00000e-15;

  // Check that dimensions of x and y vector agree
  const Index n_x = x.nelem();
  const Index n_y = y.nelem();
  if ( (n_x != 4) || (n_y != 4) )
    {
      ostringstream os;
      os << "The vectors x and y must all have the same length of 4 elements!\n"
        << "Actual lengths:\n"
        << "x:" << n_x << ", " << "y:" << n_y << ".";
      throw runtime_error(os.str());
    }

  // assure that x1 =< a < x2 holds
  if ( (a < x[1]) || (a > x[2]) )
    {
      ostringstream os;
      os << "LagrangeInterpol4: the relation x[1] =< a < x[2] is not satisfied. " 
         << "No interpolation can be calculated.\n";
      throw runtime_error(os.str());
    };

  // calculate the Lagrange polynomial coefficients for a polynomial of the order of 3
  Numeric b[4];
  for (Index i=0 ; i < 4 ; ++i)
    {
      b[i] = 1.000e0;
      for (Index k=0 ; k < 4 ; ++k)
        {
          if ( (k != i) && (fabs(x[i]-x[k]) > Dlimit) )  
            b[i] = b[i] * ( (a-x[k]) / (x[i]-x[k]) );
        };
    };

  Numeric ya = 0.000e0;
  for (Index i=0 ; i < n_x ; ++i) ya = ya + b[i]*y[i];

  return ya;
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
    x[i] = start + (double)i*step;
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
  assert( 1<n );                // Number of points must be greater 1.
  x.resize(n);
  Numeric step = (stop-start)/((double)n-1) ;
  for ( Index i=0; i<n-1; i++ )
    x[i] = start + (double)i*step;
  x[n-1] = stop;
}
void nlinspace(
               VectorView        x,
               const Numeric     start,     
               const Numeric     stop,        
               const Index       n )
{
  Numeric step = (stop-start)/((double)n-1) ;
  for ( Index i=0; i<n-1; i++ )
    x[i] = start + (double)i*step;
  x[n-1] = stop;
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
  // Number of points must be greater than 1:
  assert( 1<n );        
  // Only positive numbers are allowed for start and stop:
  assert( 0<start );
  assert( 0<stop );

  x.resize(n);
  Numeric a = log(start);
  Numeric step = (log(stop)-a)/((double)n-1);
  x[0] = start;
  for ( Index i=1; i<n-1; i++ )
    x[i] = exp(a + (double)i*step);
  x[n-1] = stop;
}


//! AngIntegrate_trapezoid
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    \param Integrand The Matrix to be integrated
    \param za_grid   The zenith angle grid 
    \param aa_grid   The azimuth angle grid 
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid(ConstMatrixView Integrand,
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
  
  return res;
}


//! AngIntegrate_trapezoid_opti
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.

    In addition to the "old fashined" integration method, it checks whether
    the stepsize is constant. If it is, it uses a faster method, if not, it
    uses the old one.

    \param Integrand Input : The Matrix to be integrated
    \param za_grid Input : The zenith angle grid 
    \param aa_grid Input : The azimuth angle grid
    \param grid_stepsize Input : stepsize of the grid
    
    \return The resulting integral

    \author Claas Teichmann <claas@sat.physik.uni-bremen.de>
    \date 2003/05/28
*/
Numeric AngIntegrate_trapezoid_opti(ConstMatrixView Integrand,
                                    ConstVectorView za_grid,
                                    ConstVectorView aa_grid,
                                    ConstVectorView grid_stepsize)
{
  Numeric res = 0;
  if ((grid_stepsize[0] > 0) && (grid_stepsize[1] > 0))
    {
      Index n = za_grid.nelem();
      Index m = aa_grid.nelem();
      Numeric stepsize_za = grid_stepsize[0];
      Numeric stepsize_aa = grid_stepsize[1];
      Vector res1(n);
      assert (is_size(Integrand, n, m));

      Numeric temp = 0.0;
      
      for (Index i = 0; i < n ; ++i)
        {
          temp = Integrand(i, 0);
          for (Index j = 1; j < m - 1; j++)
            {
              temp += Integrand(i, j) * 2;
            }
          temp += Integrand(i, m-1);
          temp *= 0.5 * DEG2RAD * stepsize_aa * sin(za_grid[i] * DEG2RAD);
          res1[i] = temp;
        }

      res = res1[0];
      for (Index i = 1; i < n - 1; i++)
        {
          res += res1[i] * 2;
        }
      res += res1[n-1];
      res *= 0.5 * DEG2RAD * stepsize_za;
    }
  else
    {
      res = AngIntegrate_trapezoid(Integrand, za_grid, aa_grid);
    }

  return res;
}


//! AngIntegrate_trapezoid
/*! 
    Performs an integration of a matrix over all directions defined in angular
    grids using the trapezoidal integration method.
    The integrand is independant of the azimuth angle. The integration over
    the azimuth angle gives a 2*PI

    \param Integrand Input : The vector to be integrated
    \param za_grid Input : The zenith angle grid 

    \author Claas Teichmann
    \date   2003-05-13
    
    \return The resulting integral
*/
Numeric AngIntegrate_trapezoid(ConstVectorView Integrand,
                               ConstVectorView za_grid)
{

  Index n = za_grid.nelem();
  assert (is_size(Integrand, n));
  
  Numeric res = 0.0;
  for (Index i = 0; i < n - 1; ++i)
    {
      // in this place 0.5 * 2 * PI is calculated:
      res += PI * DEG2RAD * (Integrand[i]* sin(za_grid[i] * DEG2RAD) 
                             + Integrand[i + 1] * sin(za_grid[i + 1] * DEG2RAD))
        * (za_grid[i + 1] - za_grid[i]);
    }
  
  return res;
}




//! sign
/*! 
    Returns the sign of a numeric value.

    The function returns 1 if the value is greater than zero, 0 if it 
    equals zero and -1 if it is less than zero.

    \return      The sign of x (see above).
    \param   x   A Numeric.

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



/*! Modified gamma distribution
 *  
 *  Uses all four free parameters (n0, mu, la, ga) to calculate
 *    psd(D) = n0 * D^mu * exp( -la * x^ga )
 *  
 *  Reference: Eq 1 of Petty & Huang, JAS, (2011).

    \param psd       Particle number density per x-interval. Sizing of vector
                     should be done before calling the function.
    \param x       Mass or size.
    \param n0      See above.
    \param mu      See above.
    \param la      See above.
    \param ga      See above.
  
  \author Jana Mendrok, Patrick Eriksson
  \date 2017-06-07

*/
void mgd(
          VectorView  psd,
    const Vector&     x,
    const Numeric&    n0,
    const Numeric&    mu,
    const Numeric&    la,
    const Numeric&    ga )
{
  const Index nx = x.nelem();

  assert( psd.nelem() == nx );

  if( ga == 1 )
    {
      if( mu == 0 )
        {
          // Exponential distribution
          for( Index ix=0; ix<nx; ix++ )
            {
              const Numeric eterm = exp( -la*x[ix] );
              psd[ix] = n0 * eterm;
            }
        }
      else
        {
          if( mu > 10 )
            {
              ostringstream os;
              os << "Given mu is " << mu << endl
                 <<"Seems unreasonable. Have you mixed up the inputs?";
              throw runtime_error(os.str());
            }
          // Gamma distribution
          for( Index ix=0; ix<nx; ix++ )
            {
              const Numeric eterm = exp( -la*x[ix] );
              const Numeric xterm = pow( x[ix], mu );
              psd[ix] = n0 * xterm * eterm;
              psd[ix] = n0 * pow( x[ix], mu ) * exp( -la*x[ix] );
            }
        }
    }
  else
    {
      // Complete MGD
      if( mu > 10 )
        {
          ostringstream os;
          os << "Given mu is " << mu << endl
             <<"Seems unreasonable. Have you mixed up the inputs?";
          throw runtime_error(os.str());
        }
      if( ga > 10 )
        {
          ostringstream os;
          os << "Given gamma is " << ga << endl
             <<"Seems unreasonable. Have you mixed up the inputs?";
          throw runtime_error(os.str());
        }
      for( Index ix=0; ix<nx; ix++ )
        {
          const Numeric pterm = pow( x[ix], ga );
          const Numeric eterm = exp( -la * pterm );
          const Numeric xterm = pow( x[ix], mu );
          psd[ix] = n0 * xterm * eterm;
        }
    }
}



/*! Modified gamma distribution, and derivatives
 *  
 *  As mgd, but this version can also return the derivate of psd with respect 
 *  to the four parameters.   

    \param psd       Particle number density per x-interval. Sizing of vector
                     should be done before calling the function.
    \param jac_data  Container for returning jacobian data. Shall be a matrix
                     with four rows, where the rows match n0, mu, la and ga.  
                     Number of columns same as length of psd.
    \param x       Mass or size.
    \param n0      See above.
    \param mu      See above.
    \param la      See above.
    \param ga      See above.
    \param do_n0_jac  Flag to actually calculate d_psd/d_n0
    \param do_mu_jac  Flag to actually calculate d_psd/d_mu
    \param do_la_jac  Flag to actually calculate d_psd/d_la
    \param do_ga_jac  Flag to actually calculate d_psd/d_ga
  
  \author Patrick Eriksson
  \date 2017-06-07

*/
void mgd_with_derivatives(
          VectorView  psd,
          MatrixView  jac_data,
    const Vector&     x,
    const Numeric&    n0,
    const Numeric&    mu,
    const Numeric&    la,
    const Numeric&    ga,
    const bool&       do_n0_jac,
    const bool&       do_mu_jac,
    const bool&       do_la_jac,
    const bool&       do_ga_jac )
{
  const Index nx = x.nelem();

  assert( psd.nelem() == nx );
  assert( jac_data.nrows() == 4 );
  assert( jac_data.ncols() == nx );

  if( ga == 1  &&  !do_ga_jac )
    {
      if( mu == 0  &&  !do_mu_jac )
        {
          // Exponential distribution
          for( Index ix=0; ix<nx; ix++ )
            {
              const Numeric eterm = exp( -la*x[ix] );
              psd[ix] = n0 * eterm;
              if( do_n0_jac )
                { jac_data(0,ix) = eterm; }
              if( do_la_jac )
                { jac_data(2,ix) = -x[ix] * psd[ix]; }
            }
        }
      else
        {
          if( mu > 10 )
            {
              ostringstream os;
              os << "Given mu is " << mu << endl
                 <<"Seems unreasonable. Have you mixed up the inputs?";
              throw runtime_error(os.str());
            }
          // Gamma distribution
          for( Index ix=0; ix<nx; ix++ )
            {
              const Numeric eterm = exp( -la*x[ix] );
              const Numeric xterm = pow( x[ix], mu );
              psd[ix] = n0 * xterm * eterm;
              if( do_n0_jac )
                { jac_data(0,ix) = xterm * eterm; }
              if( do_mu_jac )
                { jac_data(1,ix) = log(x[ix]) * psd[ix]; }
              if( do_la_jac )
                { jac_data(2,ix) = -x[ix] * psd[ix]; }
              psd[ix] = n0 * pow( x[ix], mu ) * exp( -la*x[ix] );
            }
        }
    }
  else
    {
      // Complete MGD
      if( mu > 10 )
        {
          ostringstream os;
          os << "Given mu is " << mu << endl
             <<"Seems unreasonable. Have you mixed up the inputs?";
          throw runtime_error(os.str());
        }
      if( ga > 10 )
        {
          ostringstream os;
          os << "Given gamma is " << ga << endl
             <<"Seems unreasonable. Have you mixed up the inputs?";
          throw runtime_error(os.str());
        }
      for( Index ix=0; ix<nx; ix++ )
        {
          const Numeric pterm = pow( x[ix], ga );
          const Numeric eterm = exp( -la * pterm );
          const Numeric xterm = pow( x[ix], mu );
          psd[ix] = n0 * xterm * eterm;
          if( do_n0_jac )
            { jac_data(0,ix) = xterm * eterm; }
          if( do_mu_jac )
            { jac_data(1,ix) = log(x[ix]) * psd[ix]; }
          if( do_la_jac )
            { jac_data(2,ix) = -pterm * psd[ix]; }
          if( do_ga_jac )
            { jac_data(3,ix) = -la * pterm * log(x[ix]) * psd[ix]; }
        }
    }
}



//! Generalized Modified Gamma Distribution
/*! Returns number density per unit of 'x' as function of 'x'.
 
 \return  dN Number density as function of x.
 \param   x       Numeric
 \param   N0      Numeric, Scaling parameter
 \param   Lambda  Numeric, Shape parameter
 \param   mu      Numeric, Shape parameter
 \param   gamma   Numeric, Shape parameter
 
 \author Manfred Brath
 \date   2015-01-19
 */


Numeric mod_gamma_dist(Numeric x,
                       Numeric N0,
                       Numeric Lambda,
                       Numeric mu,
                       Numeric gamma)
{
    Numeric dN;
    
    if (x > 0. && N0 > 0. && Lambda >0. && (mu+1)/gamma > 0.)
    {
        
        //Distribution function
        dN=N0*pow(x ,mu)*exp(-Lambda*pow(x,gamma));
        
        return dN;
    }
    else
    {
        ostringstream os;
        os << "At least one argument is zero or negative.\n"
        << "Modified gamma distribution can not be calculated.\n"
        << "x      = "<< x << "\n"
        << "N0     = "<< N0 << "\n"
        << "lambda = "<< Lambda << "\n"
        << "mu     = "<< mu << "\n"
        << "gamma  = "<< gamma << "\n";
        
        throw runtime_error(os.str());
    }
}

//! unitl
/*!
    Normalises a vector to have unit length.

    The standard Euclidean norm is used (2-norm).

    param    x   In/Out: A vector.

    \author Patrick Eriksson
    \date   2012-02-12
*/
void unitl( Vector& x )
{
  assert( x.nelem() > 0 );
 
  const Numeric l = sqrt(x*x);
  for(Index i=0; i<x.nelem(); i++ )
    x[i] /= l;
}

//! flat
/*!
    Flattens a matrix to a vector

    The matrix is read from front, i.e. rows are looped first. 
    In Matlab this equals x=X(:).

    \param[out] x   The vector. Should already be sized
    \param[in]  X   The matrix.

    \author Patrick Eriksson
    \date   2015-09-09
*/
void flat( VectorView x, ConstMatrixView X )
{
  assert( x.nelem() == X.nrows()*X.ncols() );

  Index i = 0; 

  for( Index c=0; c<X.ncols(); c++ )
    {
      for( Index r=0; r<X.nrows(); r++ )
        { 
          x[i] = X(r,c);
          i += 1;
        }
    }
}

//! flat
/*!
    Converts Tensor3 to a vector

    The matrix is read from front, i.e. pages are looped first, followed by rows. 
    In Matlab this equals x=X(:).

    \param[out] x   The vector. Should already be sized
    \param[in]  X   The tensor.

    \author Patrick Eriksson
    \date   2015-09-09
*/
void flat( VectorView x, ConstTensor3View X )
{
    assert( x.nelem() == X.nrows()*X.ncols()*X.npages() );

    Index i = 0;

    for( Index c=0; c<X.ncols(); c++ )
    {
        for( Index r=0; r<X.nrows(); r++ )
        {
            for( Index p=0; p<X.npages(); p++ )
            {
                x[i] = X(p,r,c);
                i += 1;
            }
        }
    }
}

//! reshape
/*!
    Converts vector to Tensor3

    The tensor is filled from front, i.e. pages are looped first, followed by rows. 
    In Matlab this equals X = reshape( x, [ X.npages(), X.nrows(), X.ncols() ]

    \param[out] X   The tensor. Should already be sized
    \param[in]  x   The vector.

    \author Patrick Eriksson
    \date   2015-09-10
*/
void reshape( Tensor3View X, ConstVectorView x )
{
    assert( x.nelem() == X.nrows()*X.ncols()*X.npages() );

    Index i = 0;

    for( Index c=0; c<X.ncols(); c++ )
    {
        for( Index r=0; r<X.nrows(); r++ )
        {
            for( Index p=0; p<X.npages(); p++ )
            {
                X(p,r,c) = x[i];
                i += 1;
            }
        }
    }
}

//! reshape
/*!
    Converts vector to Matrix

    The matrix is filled from front, i.e. rows are looped first, followed by cols. 
    In Matlab this equals X = reshape( x, [ X.nrows(), X.ncols() ]

    \param[out] X   The matrix. Should already be sized
    \param[in]  x   The vector.

    \author Patrick Eriksson
    \date   2015-09-10
*/
void reshape( MatrixView X, ConstVectorView x )
{
    assert( x.nelem() == X.nrows()*X.ncols() );

    Index i = 0;

    for( Index c=0; c<X.ncols(); c++ )
    {
        for( Index r=0; r<X.nrows(); r++ )
        {
            X(r,c) = x[i];
            i += 1;
        }
    }
}
