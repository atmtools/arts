/* Copyright (C) 2002 Stefan Buehler  <sbuehler@uni-bremen.de>

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
  \file   logic.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri May  3 19:02:36 2002
  
  \brief  Logical functions.
  
  All functions here have return type bool. They all check whether some
  condition is fullfilled and return true if that is the case.

  These functions are intended to be used either inside "if" statements
  or inside "assert" statements.

  The condition should have a simple and intuitive meaning!
*/

#include "logic.h"

//! Checks if a variable equals 0 or 1.
/*!
  \return       True if the variable is 0 or 1. Otherwise false.
  \param    x   A variable of type Index.
*/
bool is_bool( const Index& x )
{
  return ( x==0 || x==1 );
}

//! Verifies that the size of x is l.
/*! 
  This function is supposed to be used together with assert like this:
  assert(is_size(x,l)) 

  \param  x The Vector to check.
  \param  n The desired length.
  \return True if the size of x is l.
*/
bool is_size( ConstVectorView   x,
	      const Index&      n ) 
{
  return( n == x.nelem() );
}

//! Verifies that the size of x is r by c.
/*! 
  \param  x The Matrix to check.
  \param  r The desired number of rows.
  \param  c The desired number of columns.
  \return True if the size of x is r x c.
*/
bool is_size( ConstMatrixView   x,
	      const Index&      r,
	      const Index&      c ) 
{
  return( r == x.nrows() &&
	  c == x.ncols()     );
}

//! Checks if a vector is sorted in ascending order.
/*!
  Duplicated values are allowed.

  \param   x   A vector.
  \return      True if sorted.
*/
bool is_sorted( ConstVectorView& x )
{
  if( x.nelem() > 1 )
    {
      for( Index i=1; i<x.nelem(); i++ )
	{
	  if( x[i] < x[i-1] )
	    return false;
	}
    }
  return true;
}

//! Checks if a vector is sorted and strictly increasing. 
/*! 
    Duplicated values are not allowed.

    \return      True if strictly increasing, otherwise false.
    \param   x   A vector.
*/
bool is_increasing( ConstVectorView& x )
{
  if( x.nelem() > 1 )
    {
      for( Index i=1; i<x.nelem(); i++ )
	{
	  if( x[i] <= x[i-1] )
	    return false;
	}
    }
  return true;
}

//! Checks if a vector is sorted in reversed order and is strictly decreasing.
/*! 
    Duplicated values are not allowed.

    \return      True if strictly decreasing, otherwise false.
    \param   x   A vector.
*/
bool is_decreasing( ConstVectorView& x )
{
  if( x.nelem() > 1 )
    {
      for( Index i=1; i<x.nelem(); i++ )
	{
	  if( x[i] >= x[i-1] )
	    return false;
	}
    }
  return true;
}


