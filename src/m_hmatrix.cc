/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   m_hmatrix.cc

   This file contains all workspace methods (except IO) for H_total and H_data.

   \author Patrick Eriksson
   \date 2000-10-06 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "hmatrix.h"
#include "messages.h"          
#include "wsv.h"          
#include "file.h"




////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
//   Apply H on data
////////////////////////////////////////////////////////////////////////////

// The same code is found in both fuctions (vector and matrix). This to
// be more efficient. Accordingly, a change should be done for both
// functions.

//// h_apply (on vector) ///////////////////////////////////////////////
//
/** Applies a H matrix on a vector (spectra).

    \retval y2      new vector, y2 = h*y1
    \param  h       H matrix
    \param  y1      original vector

    @exception logic_error The issparse flag is not 0 or 1.

    \author Patrick Eriksson 
    \date   2000-10-22 
*/
void hApply (
              VECTOR&     y2, 
        const Hmatrix&    h,
        const VECTOR&     y1 )
{
  if ( h.issparse == 0 )
  {
    if ( h.full.dim(2) != y1.dim() )
      throw runtime_error("Size of h and length of y do not match.");
    y2 = h.full * y1;
  }

  else if ( h.issparse == 1 )
    throw runtime_error("H set to be SPARSE.");

  else
    throw logic_error("The isspare flag can only be 0 or 1.");
}



//// h_apply (on matrix) ///////////////////////////////////////////////
//
/** Applies a H matrix on a matrix (weighting functions).

    \retval k2      new matrix, k2 = h*k1
    \param  h       H matrix
    \param  k1      original matrix

    @exception logic_error The issparse flag is not 0 or 1.

    \author Patrick Eriksson 
    \date   2000-10-06 
*/
void h_apply (
              MATRIX&     k2, 
        const Hmatrix&    h,
        const MATRIX&     k1 )
{
  if ( h.issparse == 0 )
  {
    if ( h.full.dim(2) != k1.dim(1) )
      throw runtime_error("Sizes of h and k do not match.");
    k2 = h.full * k1;
  }

  else if ( h.issparse == 1 )
    throw runtime_error("H set to be SPARSE.");

  else
    throw logic_error("The isspare flag can only be 0 or 1.");
}
