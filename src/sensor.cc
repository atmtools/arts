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



/*===========================================================================
  === File description
  ===========================================================================*/

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
#include "logic.h"
#include "poly_roots.h"
#include "special_interp.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! sensor_integration_vector
/*!
   Calculates the (row) vector that multiplied that multiplied with an
   unknown (column) vector approximates the integral of the product
   between the functions represented by the two vectors.

   E.g. h*g = integral( f(x)*g(x) dx )

   \return         The multiplication (row) vector.
   \param   f      The values of function f(x).
   \param   x_f    The grid points of function f(x).
   \param   x_g    The grid points of function g(x).

   \author Mattias Ekström
   \date   2003-02-13
*/
Vector sensor_integration_vector(
        const Vector&   f,
        const Vector&   x_ftot,
        const Vector&   x_g )
{
  //Check that vectors are sorted, ascending (no descending?)

  //Find x_f points that lies outside the scope of x_g and remove them
  Index i1_f = 0, i2_f = x_g.nelem()-1;
  while( x_ftot[i1_f] < x_g[0] ) {
    i1_f++;
  }
  while( x_ftot[i2_f] > x_g[x_g.nelem()-1] ) {
    i2_f--;
  }
  Vector x_f = x_ftot[Range(i1_f, i2_f)];

  //Create a reference grid vector containing all x_f and x_g
  //strictly sorted by adding the in a sorted way.
  Vector x_ref( x_f.nelem() + x_g.nelem() );

  Index i_f = 0, i_g = 0;
  for( Index i=0; i<x_ref.nelem(); i++ ) {
    if (x_f[i_f] < x_g[i_g]) {
      x_ref[i] = x_f[i_f];
      i_f++;
    } else if (x_f[i_f] > x_g[i_g]) {
      x_ref[i] = x_g[i_g];
      i_g++;
    } else {
      // x_f and x_g are equal
      x_ref[i] = x_f[i_f];
      i_f++;
      i_g++;
    }
  }

  //Initiate output vector, with equal size as x_g, with zeros.
  //Start calculations
  Vector h(x_g.nelem(), 0.0);
  i_f = 0;
  i_g = 0;
  Numeric dx,a0,b0,c0,a1,b1,c1,x3,x2,x1;
  for( Index i=0; i<x_ref.nelem(); i++ ) {
    //Find for what index in x_g (which is the same as for h) and f
    //calculation corresponds to
    while( x_g[i_g+1] < x_ref[i] ) {
      i_g++;
    }
    while( x_f[i_f+1] < x_ref[i] ) {
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

  return h;

}


