/* Copyright (C) 2001-2007
   Stefan Buehler <sbuehler@ltu.se>
   Mattias Ekstroem <ekstrom@rss.chalmers.se>

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
  \file   matpackII.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
  
  \brief  Header file for sparse matrices.
  
  Notes:

  There are two different ways to index: 
  S.rw(3,4) = 1;                // Read and write
  cout << S.ro(3,4);            // Read only

  This distinction is necessary, because rw() creates elements if they
  don't already exist.

  The normal index operator "()" correspondes to ro, so "S(3,4)" is
  the same as S.ro(3,4).
*/

#ifndef matpackII_h
#define matpackII_h

#include "matpackI.h"

//! The Sparse class. 
/*! 
    The chosen storage format is the `compressed column' format. This
    is the same format used by Matlab. See Matlab User Guide for
    a description.
  
    \author Stefan Buehler <sbuehler@ltu.se>
    \date   Tue Jul 15 15:05:40 2003
*/

class Sparse {
public:
  // Constructors:
  Sparse();
  Sparse(Index r, Index c);
  Sparse(const Sparse& m);
  
  // Insert functions
  void insert_row(Index r, Vector v);

  // Make identity function:
  void make_I( Index r, Index c);

  // Resize function:
  void resize(Index r, Index c);

  // Destructor:
  ~Sparse();

  // Member functions:
  Index nrows() const;
  Index ncols() const;
  Index nnz()   const;

  const vector<Numeric> * const data () const {return mdata;};
  const vector<Index> * const rowind () const {return mrowind;};
  const vector<Index> * const colptr () const {return mcolptr;};

  // Index Operators:
  Numeric& rw(Index r, Index c);
  Numeric  ro(Index r, Index c) const;
  Numeric  operator() (Index r, Index c) const;

  // Assignment Operator:
  Sparse& operator=(const Sparse& m);

  // Friends:
  friend ostream& operator<<(ostream& os, const Sparse& v);
  friend void abs (Sparse& A, const Sparse& B );
  friend void mult (VectorView y, const Sparse& M, ConstVectorView x );
  friend void mult (MatrixView A, const Sparse& B, ConstMatrixView C );
  friend void mult (Sparse& A, const Sparse& B, const Sparse& C );
  friend void transpose (Sparse& A, const Sparse& B );

private:
  //! The actual data values.
  vector<Numeric> *mdata;
  //! Row indices.
  vector<Index> *mrowind;
  //! Pointers to first data element for each column.
  vector<Index> *mcolptr;
  //! Number of rows in the sparse matrix.
  Index mrr;
  //! Number of rows in the sparse matrix.
  Index mcr;
};


// Functions for general matrix operations
void abs(       Sparse& A,
          const Sparse& B );

void mult( VectorView y,
           const Sparse& M,
           ConstVectorView x );

void mult( MatrixView A,
           const Sparse& B,
           ConstMatrixView C );

void mult( Sparse& A,
           const Sparse& B,
           const Sparse& C );

void transpose( Sparse& A,
                const Sparse& B );

#endif
