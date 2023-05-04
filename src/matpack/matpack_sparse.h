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
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
#else
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#pragma GCC diagnostic pop

#include "array.h"
#include "matpack_data.h"

//! The Sparse class.
/*!

  Wrapper class for Eigen sparse matrices. The storage format used
  is compressed row storage. Thus inserting of elements in a single
  row ordered by column index is performed in constant time, if the
  rows are themselves inserted in increasing order.

*/

struct Sparse {
  static constexpr bool matpack_type{true};

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
  [[nodiscard]] bool empty() const;
  [[nodiscard]] Index nrows() const;
  [[nodiscard]] Index ncols() const;
  [[nodiscard]] Index nnz() const;

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
  [[nodiscard]] Vector diagonal() const;

  // Index Operators:
  [[nodiscard]] Numeric& rw(Index r, Index c);
  [[nodiscard]] Numeric ro(Index r, Index c) const;
  [[nodiscard]] Numeric operator()(Index r, Index c) const;

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

  Numeric* get_element_pointer();
  int* get_column_index_pointer();
  int* get_row_start_pointer();

  // Friends:
  friend std::ostream& operator<<(std::ostream& os, const Sparse& v);
  
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

/** An array of sparse matrices. */
using ArrayOfSparse = Array<Sparse>;

#endif
