/* Copyright (C) 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

/**
   \file   matpackI.cc

   Most of the implementation of matpack is in the header file, in
   order to allow inlining. But those methods that also require Array
   are here, in order to avoid that matpackI.h depends on array.h.

   \author Stefan Buehler
   \date   2001-09-15
*/

#include "matpackI.h"
#include "array.h"

/** Assignment operator from Array<Numeric>. This copies the data from
    an Array<Numeric> to this VectorView. Dimensions must agree! 
    Resizing would destroy the selection that we might have done in
    this VectorView by setting its range. 

    Array<Numeric> can be useful to collect things in, because there
    is a .push_back method for it. Then, after collecting we usually
    have to transfer the content to a Vector. With this assignment
    operator that's easy. */
VectorView VectorView::operator=(const Array<Numeric>& v)
{
  //  cout << "Assigning VectorView from Array<Numeric>.\n";

  // Check that sizes are compatible:
  assert(mrange.mextent==v.nelem());

  // Iterators for Array:
  Array<Numeric>::const_iterator i=v.begin();
  const Array<Numeric>::const_iterator e=v.end();
  // Iterator for Vector:
  Iterator1D target = begin();

  for ( ; i!=e ; ++i,++target )
    *target = *i;

  return *this;
}

