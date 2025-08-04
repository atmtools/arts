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

#include "matpack_sparse.h"

#include <double_imanip.h>
#include <xml_io_base.h>

#include <algorithm>
#include <cmath>
#include <istream>
#include <ostream>
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
  assert(0 <= r);
  assert(0 <= c);
  assert(r < nrows());
  assert(c < ncols());

  return matrix.coeffRef((int)r, (int)c);
}

//! Plain index operator.
/*!
  This is the same as the .ro index operator.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.

*/
Numeric Sparse::operator[](Index r, Index c) const {
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
  assert(0 <= r);
  assert(0 <= c);
  assert(r < nrows());
  assert(c < ncols());

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

Sparse::Sparse(const Sparse& other)                = default;
Sparse::Sparse(Sparse&& other) noexcept            = default;
Sparse& Sparse::operator=(const Sparse& other)     = default;
Sparse& Sparse::operator=(Sparse&& other) noexcept = default;

//! Constructor setting size.
/*!
  \param r Row dimension of new sparse matrix.
  \param c Column dimension of new sparse matrix.
*/
Sparse::Sparse(Index r, Index c) : matrix((int)r, (int)c) {
  // Nothing to do here.
}

Sparse Sparse::diagonal(StridedConstVectorView v) {
  Index n = v.size();
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
      A[it.row(), it.col()] = it.value();
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
      values[j]         = it.value();
      row_indices[j]    = it.row();
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
  assert(0 <= r);
  assert(r < nrows());
  assert(v.size() == static_cast<Size>(ncols()));

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
                             StridedConstVectorView data) {
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
  assert(0 <= r);
  assert(0 <= c);

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
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == B.ncols());

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
void mult(StridedVectorView y, const Sparse& M, StridedConstVectorView x) {
  // Check dimensions:
  assert(y.size() == static_cast<Size>(M.nrows()));
  assert(static_cast<Size>(M.ncols()) == x.size());

  // Typedefs for Eigen interface
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1, Eigen::ColMajor>
      EigenColumnVector;
  typedef Eigen::Stride<1, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenColumnVector, 0, Stride> ColumnMap;

  Numeric* data;
  data = const_cast<Numeric*>(x.data_handle());
  ColumnMap x_map(data, x.size(), Stride(1, x.stride(0)));
  data = const_cast<Numeric*>(y.data_handle());
  ColumnMap y_map(data, y.size(), Stride(1, y.stride(0)));

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
void transpose_mult(StridedVectorView y,
                    const Sparse& M,
                    StridedConstVectorView x) {
  // Check dimensions:
  assert(y.size() == static_cast<Size>(M.ncols()));
  assert(static_cast<Size>(M.nrows()) == x.size());

  // Typedefs for Eigen interface
  typedef Eigen::Matrix<Numeric, 1, Eigen::Dynamic, Eigen::RowMajor>
      EigenColumnVector;
  typedef Eigen::Stride<1, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenColumnVector, 0, Stride> ColumnMap;

  Numeric* data;
  data = const_cast<Numeric*>(x.data_handle());
  ColumnMap x_map(data, x.size(), Stride(1, x.stride(0)));
  data = const_cast<Numeric*>(y.data_handle());
  ColumnMap y_map(data, y.size(), Stride(1, y.stride(0)));

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
void mult(StridedMatrixView A,
          const Sparse& B,
          const StridedConstMatrixView& C) {
  // Check dimensions:
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

  // Typedefs for Eigen interface
  typedef Eigen::
      Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          EigenMatrix;
  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenMatrix, 0, Stride> MatrixMap;

  Index row_stride, column_stride;
  row_stride    = C.stride(0);
  column_stride = C.stride(1);

  Numeric* data;
  data = const_cast<Numeric*>(C.data_handle());
  Stride c_stride(row_stride, column_stride);
  MatrixMap C_map(data, C.nrows(), C.ncols(), c_stride);

  row_stride    = A.stride(0);
  column_stride = A.stride(1);
  data          = A.data_handle();
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
void mult(StridedMatrixView A,
          const StridedConstMatrixView& B,
          const Sparse& C) {
  // Check dimensions:
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

  // Typedefs for Eigen interface.
  typedef Eigen::
      Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
          EigenMatrix;
  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Stride;
  typedef Eigen::Map<EigenMatrix, 0, Stride> MatrixMap;

  Index row_stride, column_stride;
  row_stride    = B.stride(0);
  column_stride = B.stride(1);

  Numeric* data;
  data = const_cast<Numeric*>(B.data_handle());
  Stride b_stride(row_stride, column_stride);
  MatrixMap B_map(data, B.nrows(), B.ncols(), b_stride);

  row_stride    = A.stride(0);
  column_stride = A.stride(1);
  data          = A.data_handle();
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
  assert(A.nrows() == B.ncols());
  assert(A.ncols() == B.nrows());

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
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

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
  assert(B.ncols() == C.ncols());
  assert(B.nrows() == C.nrows());

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
  assert(A.ncols() == A.nrows());
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
  assert(B.ncols() == C.ncols());
  assert(B.nrows() == C.nrows());

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

void xml_io_stream<Sparse>::write(std::ostream& os_xml,
                                  const Sparse& sparse,
                                  bofstream* pbofs,
                                  std::string_view name) {
  XMLTag sparse_tag;
  XMLTag row_tag;
  XMLTag col_tag;
  XMLTag data_tag;
  XMLTag close_tag;

  sparse_tag.name = type_name;
  if (name.length()) sparse_tag.add_attribute("name", name);
  sparse_tag.add_attribute("nrows", sparse.nrows());
  sparse_tag.add_attribute("ncols", sparse.ncols());
  //sparse_tag.add_attribute ("nnz", sparse.nnz());
  row_tag.name = "RowIndex";
  row_tag.add_attribute("nelem", sparse.nnz());
  col_tag.name = "ColIndex";
  col_tag.add_attribute("nelem", sparse.nnz());
  data_tag.name = "SparseData";
  data_tag.add_attribute("nelem", sparse.nnz());

  sparse_tag.write_to_stream(os_xml);
  std::println(os_xml);

  row_tag.write_to_stream(os_xml);
  std::println(os_xml);

  ArrayOfIndex rowind(sparse.nnz()), colind(sparse.nnz());
  Vector data(sparse.nnz());
  sparse.list_elements(data, rowind, colind);

  // Write row indices.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      //FIXME: It should be the longer lines
      *pbofs << rowind[i];
    else
      std::println(os_xml, "{}", rowind[i]);
  }

  close_tag.name = "/RowIndex";
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);

  col_tag.write_to_stream(os_xml);
  std::println(os_xml);

  // Write column indices.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      //FIXME: It should be the longer lines
      *pbofs << colind[i];
    else
      std::println(os_xml, "{}", colind[i]);
  }

  close_tag.name = "/ColIndex";
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);

  data_tag.write_to_stream(os_xml);
  std::println(os_xml);
  xml_set_stream_precision(os_xml);

  // Write data.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      *pbofs << data[i];
    else
      std::print(os_xml, "{} ", data[i]);
  }
  std::println(os_xml);
  close_tag.name = "/SparseData";
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);

  close_tag.name = type_name;
  close_tag.write_to_end_stream(os_xml);

  std::println(os_xml);
}

void xml_io_stream<Sparse>::read(std::istream& is_xml,
                                 Sparse& sparse,
                                 bifstream* pbifs) {
  XMLTag tag;
  Index nrows, ncols, nnz;

  tag.read_from_stream(is_xml);
  tag.check_name("Sparse");

  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  sparse.resize(nrows, ncols);

  tag.read_from_stream(is_xml);
  tag.check_name("RowIndex");
  tag.get_attribute_value("nelem", nnz);

  ArrayOfIndex rowind(nnz), colind(nnz);
  Vector data(nnz);

  for (Index i = 0; i < nnz; i++) {
    if (pbifs) {
      *pbifs >> rowind[i];
      if (pbifs->fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Row index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    } else {
      is_xml >> rowind[i];
      if (is_xml.fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Row index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/RowIndex");

  tag.read_from_stream(is_xml);
  tag.check_name("ColIndex");

  for (Index i = 0; i < nnz; i++) {
    if (pbifs) {
      *pbifs >> colind[i];
      if (pbifs->fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Column index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    } else {
      is_xml >> colind[i];
      if (is_xml.fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Column index: " << i;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/ColIndex");

  tag.read_from_stream(is_xml);
  tag.check_name("SparseData");

  if (pbifs) {
    pbifs->readDoubleArray(data.data_handle(), nnz);
  } else {
    for (Index i = 0; i < nnz; i++) {
      is_xml >> double_imanip() >> data[i];
      if (is_xml.fail()) {
        std::ostringstream os;
        os << " near "
           << "\n  Data element: " << i;
        xml_data_parse_error(tag, os.str());
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/SparseData");

  tag.read_from_stream(is_xml);
  tag.check_name("/Sparse");

  sparse.insert_elements(nnz, rowind, colind, data);
}
