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
#include "auto_wsv.h"          
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
              Vector&     y2, 
        const Hmatrix&    h,
        const Vector&     y1 )
{
  if ( h.issparse == 0 )
  {
    if ( h.full.ncols() != y1.size() )
      throw  runtime_error("Size of h and length of y do not match.");
    if ( h.full.nrows() != y2.size() )
      {
	ostringstream os;
	os << "Size of h and length of output y do not match.\n"
	   << "H.nrows() = " << h.full.nrows() << "\n"
	   << "y.size()  = " << y2.size();
	throw  runtime_error(os.str());
      }
    mult( h.full, y1, y2 ); 	    //    y2 = h.full * y1;
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
              Matrix&     k2, 
        const Hmatrix&    h,
        const Matrix&     k1 )
{
  if ( h.issparse == 0 )
  {
    if ( h.full.ncols() != k1.nrows() )
      throw runtime_error("Sizes of h and k do not match.");
    if ( h.full.nrows() != k2.nrows() )
      throw  runtime_error("Sizes of h and output k do not match.");
    mult( h.full, k1, k2 );	    //    k2 = h.full * k1;
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

  // FIXME: There probably also should be a safety check here that hd
  // has the correct dimension.

  if ( (h1.issparse==0) && (h2.issparse==0) )
  {
    if ( (h1.full.nrows()!=h2.full.nrows()) || (h1.full.ncols()!=h2.full.ncols()) )
      throw runtime_error("Sizes of h and k do not match.");
    hd.issparse = 0;
    //    hd.full = h2.full - h1.full;
    copy( scaled(h1.full,-1), hd.full );
    add ( h2.full, hd.full );
  }

  else if ( (h1.issparse==1) && (h2.issparse==1) )
    throw runtime_error("H set to be SPARSE.");

  else
    throw logic_error("Both matrices must either be sparse or full.");
}



////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////


void VectorApplyH(
              Vector&     y2, 
        const String&     name_y2,
        const Hmatrix&    h,
        const Vector&     y1,
        const String&     name_h,
        const String&     name_y1 )
{
  out2<<"  "<<name_y2<<"="<<name_h<<"*"<<name_y1<<"\n";

  // Get output dimension.
  // FIXME: This has to be adapted for sparse matrices!
  size_t n = h.full.nrows();

  // For efficiency we distinguish between the case that input and
  // output argument are the same and the case that they are
  // different. In the latter case the output vector does not have to
  // be copied.
  if ( name_y1 == name_y2 )
  {
    Vector y(n);
    h_apply( y, h, y1 );
    y2 = y;
  }
  else
    {
      resize(  y2, n );
      h_apply( y2, h, y1 ); 
    }
}



void MatrixApplyH(
              Matrix&     k2, 
        const String&     name_k2,
        const Hmatrix&    h,
        const Matrix&     k1,
        const String&     name_h,
        const String&     name_k1 )
{
  out2<<"  "<<name_k2<<"="<<name_h<<"*"<<name_k1<<"\n";

  if ( name_k1 == name_k2 )
  {
    Matrix k(h.full.nrows(),k1.ncols());
    // FIXME: This has to be adapted for sparse matrices!
    h_apply( k, h, k1 );
    k2 = k;
  }
  else
    h_apply( k2, h, k1 ); 
}



