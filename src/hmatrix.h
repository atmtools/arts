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
   \file   hmatrix.h

   Stuff releated to H matrices.

   A matrix is a transfer matrix describing sensor characteristics and/or
   data reduction.

   This file contains the definition of:

     1. the structure holding H matrices
  
     2. functions to apply H matrices on data

     3. core functions to set up a H matrix  (To be done)

   \author Patrick Eriksson
   \date 2000-10-06
*/



#ifndef hmatrix_h
#define hmatrix_h


#include "matpackI.h"

////////////////////////////////////////////////////////////////////////////
//   The H matrix structure
////////////////////////////////////////////////////////////////////////////

/** Sensor and data reduction transfer matrix. 

    The H matrix structure has the fields:
    \verbatim
       Index         issparse
       Matrix         full
       SPARSE         sparse
    where 
       issparse  a flag to indicate if the matrix is sparse or "full"
       full      field used to store H when the matrix is full
       sparse    field used to store H when the matrix is sparse       
    \endverbatim

    One of the fields FULL and SPARSE should be empty. The field ISSPARSE
    gives the structure of H and witch of the data fields that is used.
    ISSPARSE=0 means that the matrix is full, and ISSPARSE=1 means that the
    matrix is sparse.

    \author Patrick Eriksson 
    \date 2000-10-06
*/
struct Hmatrix {
  Index         issparse;
  Matrix         full;
  // Here we should have a field for sparse matrix
};



////////////////////////////////////////////////////////////////////////////
//   Math with H matrices
////////////////////////////////////////////////////////////////////////////

void h_apply (
              Vector&     y2, 
        const Hmatrix&    h,
        const Vector&     y1 );

void h_apply (
              Matrix&     k2, 
        const Hmatrix&    h,
        const Matrix&     k1 );

void h_diff (
              Hmatrix&    hd, 
        const Hmatrix&    h2,
        const Hmatrix&    h1 );

#endif  // hmatrix_h
