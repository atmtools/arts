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
//   Math with H matrices
////////////////////////////////////////////////////////////////////////////

// The same code is found in both fuctions (vector and matrix). This to
// be more efficient. Accordingly, a change should be done for both
// functions.

//// h_apply (on vector) ///////////////////////////////////////////////////
//
/** Core function to apply a H matrix on a vector (spectra).

    \retval y2      new vector, y2 = h*y1
    \param  h       H matrix
    \param  y1      original vector

    @exception logic_error    The issparse flag is not 0 or 1.
    @exception runtime_error  Size of h and length of y do not match.

    \author Patrick Eriksson 
    \date   2000-10-22 
*/
void h_apply (
              VECTOR&     y2, 
        const Hmatrix&    h,
        const VECTOR&     y1 )
{
  if ( h.issparse == 0 )
  {
    if ( h.full.dim(2) != y1.dim() )
      throw  runtime_error("Size of h and length of y do not match.");
    y2 = h.full * y1;
  }

  else if ( h.issparse == 1 )
    throw runtime_error("H set to be SPARSE.");

  else
    throw logic_error("The isspare flag can only be 0 or 1.");
}



//// h_apply (on matrix) ///////////////////////////////////////////////////
//
/** Core function to apply a H matrix on a matrix (weighting functions).

    \retval k2      new matrix, k2 = h*k1
    \param  h       H matrix
    \param  k1      original matrix

    @exception logic_error The issparse flag is not 0 or 1.
    @exception runtime_error  Sizes of h and k do not match.

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



//// h_diff /////////////////////////////////////////////////////////////////
//
/** Calculates the difference between two H matrices.

    Both matrices must either be full or sparse.

    \retval hd      hd = h2 -h1
    \param  h2      a H matrix
    \param  k1      another H matrix

    @exception logic_error    The issparse flag is not 0 or 1.
    @exception runtime_error  The sizes of the H matrices do not match.
    @exception logic_error    Both matrices must either be sparse or full.
    \author Patrick Eriksson 
    \date   2000-10-22
*/
void h_diff (
              Hmatrix&    hd, 
        const Hmatrix&    h2,
        const Hmatrix&    h1 )
{
  if ( (h1.issparse==0) && (h2.issparse==0) )
  {
    if ( (h1.full.dim(1)!=h2.full.dim(1)) || (h1.full.dim(2)!=h2.full.dim(2)) )
      throw runtime_error("Sizes of h and k do not match.");
    hd.issparse = 0;
    hd.full = h2.full - h1.full;
  }

  else if ( (h1.issparse==1) && (h2.issparse==1) )
    throw runtime_error("H set to be SPARSE.");

  else
    throw logic_error("Both matrices must either be sparse or full.");
}



////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////


//// VectorApplyH //////////////////////////////////////////////////////////
//
/** Applies a H matrix on a vector (spectra).

    \retval y2        new vector, y2 = h*y1
    \param  name_y2   variable name of y2
    \param  h         H matrix
    \param  y1        original vector
    \param  name_y1   variable name of y1
    \param  name_h    variable name of h

    \author Patrick Eriksson 
    \date   2000-10-22 
*/
void VectorApplyH(// WS Generic Output:
              VECTOR&     y2, 
        const string&     name_y2,
        const Hmatrix&    h,
        const VECTOR&     y1,
        const string&     name_h,
        const string&     name_y1 )
{
  out2<<"  "<<name_y2<<"="<<name_h<<"*"<<name_y1<<"\n";

  if ( name_y1 == name_y2 )
  {
    VECTOR y;
    h_apply( y, h, y1 );
    y2 = y;
  }
  else
    h_apply( y2, h, y1 ); 
}



//// MatrixApplyH //////////////////////////////////////////////////////////
//
/** Applies a H matrix on a matrix (WFs).

    \retval k2        new matrix, k2 = h*k1
    \param  name_k2   variable name of k2
    \param  h         H matrix
    \param  k1        original matrix
    \param  name_k1   variable name of k1
    \param  name_h    variable name of h

    \author Patrick Eriksson 
    \date   2000-10-22 
*/
void MatrixApplyH(// WS Generic Output:
              MATRIX&     k2, 
        const string&     name_k2,
        const Hmatrix&    h,
        const MATRIX&     k1,
        const string&     name_h,
        const string&     name_k1 )
{
  out2<<"  "<<name_k2<<"="<<name_h<<"*"<<name_k1<<"\n";

  if ( name_k1 == name_k2 )
  {
    MATRIX k;
    h_apply( k, h, k1 );
    k2 = k;
  }
  else
    h_apply( k2, h, k1 ); 
}
