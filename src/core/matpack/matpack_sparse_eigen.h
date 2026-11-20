#pragma once

/** Sparse as an Eigen sparse matrix
 *
 * Include this only where Eigen is needed for a sparse matrix.  Sparse stores the
 * compressed rows of an Eigen row-major sparse matrix, so it maps as one without a
 * copy, and an Eigen result is copied into a Sparse.
 */

#ifndef _MSC_VER
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
#endif

#include <Eigen/Core>
#include <Eigen/SparseCore>

#ifndef _MSC_VER
#pragma GCC diagnostic pop
#endif

#include "matpack_sparse.h"

namespace matpack {
//! The Eigen sparse matrix whose storage Sparse has
using eigen_sparse = Eigen::SparseMatrix<Numeric, Eigen::RowMajor, int>;

//! The matrix as an Eigen sparse matrix, without a copy, valid while the structure of the matrix is unchanged
inline Eigen::Map<const eigen_sparse> eigen_map(const Sparse& x) {
  return Eigen::Map<const eigen_sparse>(
      x.nrows(), x.ncols(), x.nnz(), x.row_starts().data(), x.column_indices().data(), x.values().data());
}

//! The matrix of an Eigen sparse matrix
Sparse from_eigen(eigen_sparse x);
}  // namespace matpack
