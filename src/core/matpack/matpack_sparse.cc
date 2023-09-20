/*!
  \file   matpackII.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \author Mattias Ekstroem
  \author Oliver Lemke
  \author Simon Pfreundschuh <simonpf@chalmers.se> (Eigen interface)
  \date   Wed Jan 06 10:17:15 2016

  \brief  Implementation of sparse matrices.

  Notes:

  There are two different ways to index:
  S.rw(3,4) = 1;                // Read and write
  cout << S.ro(3,4);            // Read only

  This distinction is necessary, because rw() creates elements if they
  don't already exist.

  The normal index operator "()" correspondes to ro, so "S(3,4)" is
  the same as S.ro(3,4).
*/

#include "matpack_eigen.h"
#include "matpack_sparse.h"

#include <algorithm>
#include <cmath>
#include <iostream>  // For debugging.
#include <iterator>
#include <set>
#include <vector>


// Simple member Functions
// ----------------

//! Returns true if variable size is zero.
bool Sparse::empty() const {
  return (matrix.rows() == 0 || matrix.cols() == 0);
}

//! Returns the number of rows.
Index Sparse::nrows() const { return matrix.rows(); }

//! Returns the number of columns.
Index Sparse::ncols() const { return matrix.cols(); }

//! Returns the number of nonzero elements.
Index Sparse::nnz() const { return matrix.nonZeros(); }

// Index Operators
// ---------------

//! Read and write index operator.
/*!
  This has to correctly handle two cases:

  1. The data element exists. In this case the operator acts similar
     to the const index operator in that the element is returned.

  2. The element does not exist. In this case it is created.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.
*/
Numeric& Sparse::rw(Index r, Index c) {
  // Check if indices are valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < nrows());
  ARTS_ASSERT(c < ncols());

  return matrix.coeffRef((int)r, (int)c);
}

//! Plain index operator.
/*!
  This is the same as the .ro index operator.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.

*/
Numeric Sparse::operator()(Index r, Index c) const {
  return matrix.coeff((int)r, (int)c);
}

//! Read only index operator.
/*!
  This has to correctly handle two cases:

  1. The data element exists. In this case the element is returned.

  2. The element does not exist. In this case the value 0 is returned.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices, or zero.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003 */
Numeric Sparse::ro(Index r, Index c) const {
  // Check if indices are valid:
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);
  ARTS_ASSERT(r < nrows());
  ARTS_ASSERT(c < ncols());

  return matrix.coeff((int)r, (int)c);
}

// Constructors
// ------------

//! Default constructor.
/*!
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
Sparse::Sparse() : matrix(0, 0) {
  // Nothing to do here
}

//! Constructor setting size.
/*!
  \param r Row dimension of new sparse matrix.
  \param c Column dimension of new sparse matrix.
*/
Sparse::Sparse(Index r, Index c) : matrix((int)r, (int)c) {
  // Nothing to do here.
}

Sparse Sparse::diagonal(ConstVectorView v) {
  Index n = v.nelem();
  ArrayOfIndex indices(n);
  for (Index i = 0; i < n; ++i) {
    indices[i] = i;
  }
  Sparse m(n, n);
  m.insert_elements(n, indices, indices, v);
  return m;
}

Vector Sparse::diagonal() const {
  Index m = std::min(nrows(), ncols());
  Vector diag(m);

  for (int i = 0; i < m; i++) {
    diag[i] = matrix.coeff(i, i);
  }
  return diag;
}

//! Convert to dense matrix.
/*!
  Converts a given sparse matrix to a dense matrix. Intended mainly
  for testing purposes.

  \return The dense representation of the given sparse matrix.
*/
Sparse::operator Matrix() const {
  Index m, n;
  m = nrows();
  n = ncols();
  Matrix A(m, n);
  A = 0.0;

  for (int i = 0; i < m; i++) {
    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, i);
    for (; it; ++it) {
      A(it.row(), it.col()) = it.value();
    }
  }
  return A;
}

//! Add sparse matrix.
/*!

  Element-wise, in-situ sum of two sparse matrices. Matrix must have the
  same dimensions.

  \param A The sparse matrix to add to the Sparse object.

  \return Reference to the resulting sum.
*/
Sparse& Sparse::operator+=(const Sparse& A) {
  matrix += A.matrix;
  return *this;
}

//! Subtract sparse matrix.
/*!

  Element-wise, in-situ difference of two sparse matrices. Matrix must have the
  same dimensions.

  \param A The sparse matrix to add to the Sparse object.

  \return Reference to the resulting sum.
*/

Sparse& Sparse::operator-=(const Sparse& A) {
  matrix -= A.matrix;
  return *this;
}

//! Scale matrix.
/*!

  \param x Scalar factor to scale the matrix with.

  \return Reference to the scaled sparse matrix object.
*/
Sparse& Sparse::operator*=(Numeric x) {
  matrix = matrix * x;
  return *this;
}

//! Scale matrix by reciprocal.
/*!

  \param x Reciprocal of the factor to scale the matrix with.

  \return Reference to the scaled sparse matrix object.
*/
Sparse& Sparse::operator/=(Numeric x) {
  matrix = matrix / x;
  return *this;
}

//! List elements in matrix.
/*!

  Writes the values in the matrix into the given Vector values and the row and
  column indices into the given ArrayOfIndex objects.

  \param values The values contained in the matrix.
  \param row_indices The row indices corresponding to the values.
  \param column_indices The column indices corresponding to the values.
*/
void Sparse::list_elements(Vector& values,
                           ArrayOfIndex& row_indices,
                           ArrayOfIndex& column_indices) const {
  Index m, n;
  m = nrows();
  n = ncols();
  Matrix A(m, n);
  A = 0.0;

  values.resize(nnz());
  row_indices.resize(nnz());
  column_indices.resize(nnz());

  Index j = 0;
  for (int i = 0; i < m; i++) {
    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, i);
    for (; it; ++it) {
      values[j] = it.value();
      row_indices[j] = it.row();
      column_indices[j] = it.col();
      j++;
    }
  }
}

Numeric* Sparse::get_element_pointer() { return matrix.valuePtr(); }
int* Sparse::get_column_index_pointer() { return matrix.innerIndexPtr(); }
int* Sparse::get_row_start_pointer() { return matrix.outerIndexPtr(); }

//! Reduce matrix to the row range [offset, offset + nrows]
void Sparse::split(Index offset, Index nrows_block) {
  Eigen::SparseMatrix<Numeric, Eigen::RowMajor> block_copy(
      matrix.block((int)offset, 0, (int)nrows_block, (int)ncols()));
  matrix.swap(block_copy);
}

//! Insert row function
/*!
  Inserts a Vector as row of elements at the given position.

  The row index must agree with the size of the matrix. This
  function can not be used to expand the matrix.
  Only non-zero values will be stored. If the destination row
  already exist it will be overwritten.

  \param r Where to insert the row
  \param v Vector to be inserted.
*/
void Sparse::insert_row(Index r, Vector v) {
  // Check if the row index and the Vector length are valid
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(r < nrows());
  ARTS_ASSERT(v.nelem() == ncols());

  for (int i = 0; i < ncols(); i++) {
    if (v[i] != 0) matrix.coeffRef((int)r, i) = v[i];
  }
}

//! Insert vector of elements with given row and column indices.
/*!
  Efficient inserting of a vector of elements into the sparse matrix.
  Overwrites elements currently in the matrix. The complexity is linear
  in the number of elements and should therefore be the preferred way
  of inserting elements into the sparse matrix.

  \param nelems The number of elements to insert.
  \param rowind A vector containing the row indices.
  \param colind A vector containing the column indices.
  \param data The vector containing the elements.
*/
void Sparse::insert_elements(Index nelems,
                             const ArrayOfIndex& rowind,
                             const ArrayOfIndex& colind,
                             ConstVectorView data) {
  typedef Eigen::Triplet<Numeric> T;
  std::vector<T> tripletList(nelems);

  for (Index i = 0; i < nelems; i++) {
    tripletList[i] = T((int)rowind[i], (int)colind[i], data[i]);
  }

  matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

//! Resize function.
/*!
  If the size is already correct this function does nothing.

  All data is lost after resizing! The new Sparse is not initialised
  so it will be empty.

  \param r New row dimension.
  \param c New column dimension.
*/
void Sparse::resize(Index r, Index c) {
  ARTS_ASSERT(0 <= r);
  ARTS_ASSERT(0 <= c);

  matrix.resize((int)r, (int)c);
}

//! Output operator for Sparse.
/*!
  \param os Output stream.
  \param v Sparse matrix to print.

  \return Output stream.
*/
std::ostream& operator<<(std::ostream& os, const Sparse& M) {
  os << M.matrix;
  return os;
}

// General matrix functions

//! Absolute value of sparse matrix elements
/*!
  Computes the absolute values of the elements in sparse matrix B.

  The output matrix A must have been initialized with the correct size.

  \param A Output: Absolute value matrix.
  \param B Original matrix.

  \author Mattias Ekstrom
  \date   2005-03-21
*/
void abs(Sparse& A, const Sparse& B) {
  // Check dimensions
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == B.ncols());

  A.matrix = B.matrix.cwiseAbs();
}

//! Sparse matrix - Vector multiplication.
/*!
  This calculates the product

  y = M*x, where M is sparse.

  Output comes first!

  Dimensions of y, M, and x must match. No memory reallocation takes
  place, only the data is copied.

  \param y Output: The multiplication result.
  \param M Matrix for multiplication (sparse).
  \param x Vector for multiplication.
*/
void mult(VectorView y, const Sparse& M, ConstVectorView x) {
  // Check dimensions:
  ARTS_ASSERT(y.nelem() == M.nrows());
  ARTS_ASSERT(M.ncols() == x.nelem());

    // Typedefs for Eigen interface
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1, Eigen::ColMajor>
      EigenColumnVector;
  typedef Eigen::Stride<1, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenColumnVector, 0, Stride> ColumnMap;

  Numeric* data;
  data = x.unsafe_data_handle();
  ColumnMap x_map(data, x.nelem(), Stride(1, x.stride(0)));
  data = y.unsafe_data_handle();
  ColumnMap y_map(data, y.nelem(), Stride(1, y.stride(0)));

  y_map = M.matrix * x_map;
}

//! Sparse matrix - Vector multiplication.
/*!
  This calculates the product

  y = M*x, where M is sparse.

  Output comes first!

  Dimensions of y, M, and x must match. No memory reallocation takes
  place, only the data is copied.

  \param y Output: The multiplication result.
  \param M Matrix for multiplication (sparse).
  \param x Vector for multiplication.
*/
void transpose_mult(VectorView y, const Sparse& M, ConstVectorView x) {
  // Check dimensions:
  ARTS_ASSERT(y.nelem() == M.ncols());
  ARTS_ASSERT(M.nrows() == x.nelem());

    // Typedefs for Eigen interface
  typedef Eigen::Matrix<Numeric, 1, Eigen::Dynamic, Eigen::RowMajor>
      EigenColumnVector;
  typedef Eigen::Stride<1, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenColumnVector, 0, Stride> ColumnMap;

  Numeric* data;
  data = x.unsafe_data_handle();
  ColumnMap x_map(data, x.nelem(), Stride(1, x.stride(0)));
  data = y.unsafe_data_handle();
  ColumnMap y_map(data, y.nelem(), Stride(1, y.stride(0)));

  y_map = x_map * M.matrix;
}

//! SparseMatrix - Matrix multiplication.
/*!
  Calculates the matrix product:

  A = B*C, where B is sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix (full).
  \param B First matrix to multiply (sparse).
  \param C Second matrix to multiply (full).

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
void mult(MatrixView A, const Sparse& B, const ConstMatrixView& C) {
  // Check dimensions:
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  
  // Typedefs for Eigen interface
  typedef Eigen::
      Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          EigenMatrix;
  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenMatrix, 0, Stride> MatrixMap;

  Index row_stride, column_stride;
  row_stride = C.stride(0);
  column_stride = C.stride(1);

  Numeric* data;
  data = C.unsafe_data_handle();
  Stride c_stride(row_stride, column_stride);
  MatrixMap C_map(data, C.nrows(), C.ncols(), c_stride);

  row_stride = A.stride(0);
  column_stride = A.stride(1);
  data = A.unsafe_data_handle();
  Stride a_stride(row_stride, column_stride);
  MatrixMap A_map(data, A.nrows(), A.ncols(), a_stride);

  A_map = B.matrix * C_map;
}

//! Matrix - SparseMatrix multiplication.
/*!
  Calculates the matrix product:

  A = B*C, where C is sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix (full).
  \param B First matrix to multiply (sparse).
  \param C Second matrix to multiply (full).

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
void mult(MatrixView A, const ConstMatrixView& B, const Sparse& C) {
  // Check dimensions:
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

    // Typedefs for Eigen interface.
  typedef Eigen::
      Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          EigenMatrix;
  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenMatrix, 0, Stride> MatrixMap;

  Index row_stride, column_stride;
  row_stride = B.stride(0);
  column_stride = B.stride(1);

  Numeric* data;
  data = B.unsafe_data_handle();
  Stride b_stride(row_stride, column_stride);
  MatrixMap B_map(data, B.nrows(), B.ncols(), b_stride);

  row_stride = A.stride(0);
  column_stride = A.stride(1);
  data = A.unsafe_data_handle();
  Stride a_stride(row_stride, column_stride);
  MatrixMap A_map(data, A.nrows(), A.ncols(), a_stride);

  A_map = B_map * C.matrix;
}

//! Transpose of sparse matrix
/*!
  Computes the transpose of the sparse matrix B.

  The output matrix A must have been initialized with the correct size.

  \param A Output: Transposed matrix.
  \param B Original matrix.

  \author Mattias Ekstroem
  \date   2003-08-05
*/
void transpose(Sparse& A, const Sparse& B) {
  // Check dimensions
  ARTS_ASSERT(A.nrows() == B.ncols());
  ARTS_ASSERT(A.ncols() == B.nrows());

  A.matrix = B.matrix.transpose();
}

//! Sparse - Sparse multiplication.
/*!
  Calculates A = B*C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix.
  \param B First product matrix.
  \param C Second product matrix.

  \author Mattias Ekstroem
  \date   2003-08-06
*/
void mult(Sparse& A, const Sparse& B, const Sparse& C) {
  // Check dimensions and make sure that A is empty
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  A.matrix = B.matrix * C.matrix;
}

//! Sparse - Sparse addition.
/*!
  Calculates A = B+C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of B, and C must match.  A will be resized.

  \param A Output: Result matrix.
  \param B First summand matrix.
  \param C Second summand matrix.

  \author Oliver Lemke
  \date   2009-09-03
*/
void add(Sparse& A, const Sparse& B, const Sparse& C) {
  // Check dimensions
  ARTS_ASSERT(B.ncols() == C.ncols());
  ARTS_ASSERT(B.nrows() == C.nrows());

  A.resize(B.nrows(), B.ncols());
  A.matrix = B.matrix + C.matrix;
}

//! Sparse identity matrix
/*!

  Set the given Sparse matrix object to the identity matrix. The matrix must
  be square.

  \param A The matrix to be set to the identity matrix.
*/
void id_mat(Sparse& A) {
  ARTS_ASSERT(A.ncols() == A.nrows());
  A.matrix.setIdentity();
}

//! Sparse - Sparse subtraction.
/*!
  Calculates A = B-C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of B, and C must match.  A will be resized.

  \param A Output: Result matrix.
  \param B First subtrahend matrix.
  \param C Second subtrahend matrix.

  \author Oliver Lemke
  \date   2009-09-03
*/
void sub(Sparse& A, const Sparse& B, const Sparse& C) {
  A.resize(B.nrows(), B.ncols());

  // Check dimensions
  ARTS_ASSERT(B.ncols() == C.ncols());
  ARTS_ASSERT(B.nrows() == C.nrows());

  A.matrix = B.matrix - C.matrix;
}

Range get_rowindex_for_mblock(const Sparse& sensor_response,
                              const Index& mblock_index) {
  const Index n1y = sensor_response.nrows();
  return Range(n1y * mblock_index, n1y);
}

std::ostream& operator<<(std::ostream& os, const ArrayOfSparse& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
