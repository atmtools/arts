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
  \file   logic.h
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri May  3 19:10:04 2002
  
  \brief  Header file for logic.cc
*/

#ifndef logic_h
#define logic_h

#include "arts.h"
#include "matpackVII.h"

bool is_bool( const Index x );

bool is_size( ConstVectorView   x,
 	      const Index       l );

bool is_sorted( ConstVectorView& x );

bool is_increasing( ConstVectorView& x );

bool is_decreasing( ConstVectorView& x );


////////////////////////////////////////////////////////////////////////////
//   Template functions (have to be here in the .h file).
////////////////////////////////////////////////////////////////////////////

//! Verifies that the size of x is l.
/*! 
  This function is supposed to be used together with assert like this:
  assert(is_size(x,l)). It works for any array type.

  \param  x The Array to check.
  \param  l The desired length.
  \return True if the size of x is l.
*/
template< class T >
bool is_size( const Array<T>& x,
	      const Index     n ) 
{
  return( n == x.nelem() );
}

#endif // logic_h
