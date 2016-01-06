/* Copyright (C) 2001-2016
   Stefan Buehler <sbuehler@ltu.se>
   Mattias Ekstroem <ekstrom@rss.chalmers.se>
   Simon Pfreundschuh <simonpf@chalmers.se>

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

#include "array.h"
#include "matpackI.h"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include <iostream>

//! The Sparse class.
/*!

  Wrapper class for Eigen sparse matrices. The storage format used
  is compressed row storage. Thus inserting of elements in a single
  row ordered by column index is performed in constant time, if the
  rows are themselves inserted in increasing order.

*/

class Sparse {
public:
    // Constructors:
    Sparse();
    Sparse(Index r, Index c);

    // Insert functions
    void insert_row(Index r, Vector v);

    // Make identity function:
    void identity();

    // Resize function:
    void resize(Index r, Index c);

    // Member functions:
    bool empty() const;
    Index nrows() const;
    Index ncols() const;
    Index nnz()   const;

    // Index Operators:
    Numeric& rw(Index r, Index c);
    Numeric  ro(Index r, Index c) const;
    Numeric  operator() (Index r, Index c) const;

    // Conversion to Dense Matrix:
    operator Matrix();

    // Matrix data access
    void list_elements( Vector &values,
                        ArrayOfIndex &row_indices,
                        ArrayOfIndex &column_indices ) const;

    // Friends:
    friend std::ostream& operator<<(std::ostream& os, const Sparse& v);
    friend void abs (Sparse& A, const Sparse& B );
    friend void mult (VectorView y, const Sparse& M, ConstVectorView x );
    friend void mult (MatrixView A, const Sparse& B, const ConstMatrixView& C );
    friend void mult (MatrixView A, const ConstMatrixView& B, const Sparse& C );
    friend void mult (Sparse& A, const Sparse& B, const Sparse& C );
    friend void add (Sparse& A, const Sparse& B, const Sparse& C );
    friend void sub (Sparse& A, const Sparse& B, const Sparse& C );
    friend void transpose (Sparse& A, const Sparse& B );

private:

    //! The actual matrix.
    Eigen::SparseMatrix<Numeric, Eigen::RowMajor> matrix;

};

// Functions for general matrix operations
void abs(       Sparse& A,
          const Sparse& B );

void mult( VectorView y,
           const Sparse& M,
           ConstVectorView x );

void mult( MatrixView A,
           const Sparse& B,
           const ConstMatrixView& C );

void mult( MatrixView A,
           const ConstMatrixView& B,
           const Sparse& C );

void mult( Sparse& A,
           const Sparse& B,
           const Sparse& C );

void add( Sparse& A,
          const Sparse& B,
          const Sparse& C );

void sub( Sparse& A,
          const Sparse& B,
          const Sparse& C );

void transpose( Sparse& A,
                const Sparse& B );

#endif
