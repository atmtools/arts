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

#include <iostream>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#pragma GCC diagnostic pop

#include "array.h"
#include "matpackI.h"

//! The Sparse class.
/*!

  Wrapper class for Eigen sparse matrices. The storage format used
  is compressed row storage. Thus inserting of elements in a single
  row ordered by column index is performed in constant time, if the
  rows are themselves inserted in increasing order.

*/

struct Sparse {
  // Constructors:
  Sparse();
  Sparse(Index r, Index c);

  void split(Index offset, Index nrows);

  // Insert functions
  void insert_row(Index r, Vector v);
  void insert_elements(Index nnz,
                       const ArrayOfIndex& rowind,
                       const ArrayOfIndex& colind,
                       ConstVectorView data);

  // Resize function:
  void resize(Index r, Index c);

  // Member functions:
  bool empty() const;
  Index nrows() const;
  Index ncols() const;
  Index nnz() const;

  /** Create a sparse matrix from a vector.
     *
     * @param v vector containing the diagonal elements.
     * @return Sparse matrix with the elements of the given vector
     *     on the diagonal.
     */
  static Sparse diagonal(ConstVectorView v);

  /** Diagonal elements as vector
     *
     * Extracts the diagonal elements from the sparse matrix.
     * matrix.
     *
     * @return A vector containing the diagonal elements.
     */
  Vector diagonal() const;

  // Index Operators:
  Numeric& rw(Index r, Index c);
  Numeric ro(Index r, Index c) const;
  Numeric operator()(Index r, Index c) const;

  // Arithmetic operators:
  Sparse& operator+=(const Sparse& x);
  Sparse& operator-=(const Sparse& x);

  // Scaling operators:
  Sparse& operator*=(Numeric x);
  Sparse& operator/=(Numeric x);

  // Conversion to Dense Matrix:
  explicit operator Matrix() const;

  // Matrix data access
  void list_elements(Vector& values,
                     ArrayOfIndex& row_indices,
                     ArrayOfIndex& column_indices) const;

  Numeric* get_element_pointer() { return matrix.valuePtr(); }
  int* get_column_index_pointer() { return matrix.innerIndexPtr(); }
  int* get_row_start_pointer() { return matrix.outerIndexPtr(); }

  // Friends:
  friend std::ostream& operator<<(std::ostream& os, const Sparse& v);
  friend void abs(Sparse& A, const Sparse& B);
  friend void mult(VectorView y, const Sparse& M, ConstVectorView x);
  friend void transpose_mult(VectorView y, const Sparse& M, ConstVectorView x);
  friend void mult(MatrixView A, const Sparse& B, const ConstMatrixView& C);
  friend void mult(MatrixView A, const ConstMatrixView& B, const Sparse& C);
  friend void mult(Sparse& A, const Sparse& B, const Sparse& C);
  friend void add(Sparse& A, const Sparse& B, const Sparse& C);
  friend void sub(Sparse& A, const Sparse& B, const Sparse& C);
  friend void transpose(Sparse& A, const Sparse& B);
  friend void id_mat(Sparse& A);
  
  //! The actual matrix.
  Eigen::SparseMatrix<Numeric, Eigen::RowMajor> matrix;
};

// Functions for general matrix operations
void abs(Sparse& A, const Sparse& B);

void mult(VectorView y, const Sparse& M, ConstVectorView x);

void transpose_mult(VectorView y, const Sparse& M, ConstVectorView x);

void mult(MatrixView A, const Sparse& B, const ConstMatrixView& C);

void mult(MatrixView A, const ConstMatrixView& B, const Sparse& C);

void mult(Sparse& A, const Sparse& B, const Sparse& C);

void add(Sparse& A, const Sparse& B, const Sparse& C);

void sub(Sparse& A, const Sparse& B, const Sparse& C);

void transpose(Sparse& A, const Sparse& B);

void id_mat(Sparse& A);

#endif
