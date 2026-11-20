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
#include <numeric>
#include <ostream>
#include <vector>

#include "matpack_sparse_eigen.h"

// Simple member Functions
// ----------------

//! Returns true if variable size is zero.
bool Sparse::empty() const { return (nr == 0 || nc == 0); }

//! Returns the number of rows.
Index Sparse::nrows() const { return nr; }

//! Returns the number of columns.
Index Sparse::ncols() const { return nc; }

//! Returns the number of nonzero elements.
Index Sparse::nnz() const { return static_cast<Index>(vals.size()); }

std::span<const int>     Sparse::row_starts() const { return row_start; }
std::span<const int>     Sparse::column_indices() const { return col_index; }
std::span<const Numeric> Sparse::values() const { return vals; }
std::span<Numeric>       Sparse::values() { return vals; }

namespace {
//! The position of column c in row r, and whether the element is stored there
std::pair<Index, bool> find_element(std::span<const int> row_start, std::span<const int> col_index, Index r, Index c) {
  const auto first = col_index.begin() + row_start[r], last = col_index.begin() + row_start[r + 1];
  const auto it = std::lower_bound(first, last, static_cast<int>(c));
  return {it - col_index.begin(), it != last and *it == c};
}
}  // namespace

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

  const auto [k, found] = find_element(row_start, col_index, r, c);
  if (not found) {
    col_index.insert(col_index.begin() + k, static_cast<int>(c));
    vals.insert(vals.begin() + k, 0.0);
    for (Index i = r + 1; i <= nr; i++) row_start[i]++;
  }
  return vals[k];
}

//! Plain index operator.
/*!
  This is the same as the .ro index operator.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.

*/
Numeric Sparse::operator[](Index r, Index c) const { return ro(r, c); }

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

  const auto [k, found] = find_element(row_start, col_index, r, c);
  return found ? vals[k] : 0.0;
}

// Constructors
// ------------

//! Default constructor.
/*!
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
Sparse::Sparse() = default;

Sparse::Sparse(const Sparse& other)                = default;
Sparse::Sparse(Sparse&& other) noexcept            = default;
Sparse& Sparse::operator=(const Sparse& other)     = default;
Sparse& Sparse::operator=(Sparse&& other) noexcept = default;

//! Constructor setting size.
/*!
  \param r Row dimension of new sparse matrix.
  \param c Column dimension of new sparse matrix.
*/
Sparse::Sparse(Index r, Index c) : nr(r), nc(c), row_start(r + 1, 0) {}

Sparse::Sparse(
    Index r, Index c, std::vector<int> row_starts, std::vector<int> column_indices, std::vector<Numeric> values)
    : nr(r), nc(c), row_start(std::move(row_starts)), col_index(std::move(column_indices)), vals(std::move(values)) {
  assert(static_cast<Index>(row_start.size()) == nr + 1);
  assert(row_start.front() == 0);
  assert(static_cast<Size>(row_start.back()) == col_index.size());
  assert(col_index.size() == vals.size());
}

Sparse Sparse::diagonal(StridedConstVectorView v) {
  Index        n = v.size();
  ArrayOfIndex indices(n);
  for (Index i = 0; i < n; ++i) { indices[i] = i; }
  Sparse m(n, n);
  m.insert_elements(n, indices, indices, v);
  return m;
}

Vector Sparse::diagonal() const {
  Index  m = std::min(nrows(), ncols());
  Vector diag(m);

  for (Index i = 0; i < m; i++) { diag[i] = ro(i, i); }
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

  for (const auto [row, col, value] : *this | by_elem) A[row, col] = value;
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
  *this = matpack::from_eigen(matpack::eigen_map(*this) + matpack::eigen_map(A));
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
  *this = matpack::from_eigen(matpack::eigen_map(*this) - matpack::eigen_map(A));
  return *this;
}

//! Scale matrix.
/*!

  \param x Scalar factor to scale the matrix with.

  \return Reference to the scaled sparse matrix object.
*/
Sparse& Sparse::operator*=(Numeric x) {
  for (Numeric& v : vals) v *= x;
  return *this;
}

//! Scale matrix by reciprocal.
/*!

  \param x Reciprocal of the factor to scale the matrix with.

  \return Reference to the scaled sparse matrix object.
*/
Sparse& Sparse::operator/=(Numeric x) {
  for (Numeric& v : vals) v /= x;
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
void Sparse::list_elements(Vector& element_values, ArrayOfIndex& row_indices, ArrayOfIndex& column_indices) const {
  element_values.resize(nnz());
  row_indices.resize(nnz());
  column_indices.resize(nnz());

  Index j = 0;
  for (const auto [row, col, value] : *this | by_elem) {
    element_values[j] = value;
    row_indices[j]    = row;
    column_indices[j] = col;
    j++;
  }
}

Numeric* Sparse::get_element_pointer() { return vals.data(); }
int*     Sparse::get_column_index_pointer() { return col_index.data(); }
int*     Sparse::get_row_start_pointer() { return row_start.data(); }

//! Reduce matrix to the row range [offset, offset + nrows]
void Sparse::split(Index offset, Index nrows_block) {
  assert(0 <= offset and offset + nrows_block <= nr);

  const int first = row_start[offset], last = row_start[offset + nrows_block];

  std::vector<int> starts(nrows_block + 1);
  for (Index i = 0; i <= nrows_block; i++) starts[i] = row_start[offset + i] - first;

  row_start = std::move(starts);
  col_index = std::vector<int>(col_index.begin() + first, col_index.begin() + last);
  vals      = std::vector<Numeric>(vals.begin() + first, vals.begin() + last);
  nr        = nrows_block;
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

  // The stored elements of the row, overwritten or extended by the non-zero elements of v
  std::vector<int>     cols;
  std::vector<Numeric> row_vals;
  Index                k = row_start[r];
  for (Index i = 0; i < nc; i++) {
    const bool stored = k < row_start[r + 1] and col_index[k] == i;
    if (v[i] != 0 or stored) {
      cols.push_back(static_cast<int>(i));
      row_vals.push_back(v[i] != 0 ? v[i] : vals[k]);
    }
    if (stored) k++;
  }

  const int first = row_start[r], last = row_start[r + 1];
  col_index.erase(col_index.begin() + first, col_index.begin() + last);
  col_index.insert(col_index.begin() + first, cols.begin(), cols.end());
  vals.erase(vals.begin() + first, vals.begin() + last);
  vals.insert(vals.begin() + first, row_vals.begin(), row_vals.end());

  const int change = static_cast<int>(cols.size()) - (last - first);
  for (Index i = r + 1; i <= nr; i++) row_start[i] += change;
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
void Sparse::insert_elements(Index                  nelems,
                             const ArrayOfIndex&    rowind,
                             const ArrayOfIndex&    colind,
                             StridedConstVectorView data) {
  // As Eigen's setFromTriplets: replaces the elements, and sums duplicates in their order
  std::vector<Index> order(nelems);
  std::iota(order.begin(), order.end(), Index{0});
  std::ranges::stable_sort(
      order, [&](Index a, Index b) { return std::pair{rowind[a], colind[a]} < std::pair{rowind[b], colind[b]}; });

  row_start.assign(nr + 1, 0);
  col_index.clear();
  vals.clear();
  for (Index n = 0; n < nelems; n++) {
    const Index i = order[n];
    assert(0 <= rowind[i] and rowind[i] < nr and 0 <= colind[i] and colind[i] < nc);
    if (n > 0 and rowind[order[n - 1]] == rowind[i] and colind[order[n - 1]] == colind[i]) {
      vals.back() += data[i];
    } else {
      col_index.push_back(static_cast<int>(colind[i]));
      vals.push_back(data[i]);
      row_start[rowind[i] + 1]++;
    }
  }
  std::partial_sum(row_start.begin(), row_start.end(), row_start.begin());
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

  nr = r;
  nc = c;
  row_start.assign(r + 1, 0);
  col_index.clear();
  vals.clear();
}

//! Output operator for Sparse.
/*!
  \param os Output stream.
  \param v Sparse matrix to print.

  \return Output stream.
*/
std::ostream& operator<<(std::ostream& os, const Sparse& M) {
  arts::print(os, "{}", M);
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

  A = B;
  for (Numeric& v : A.values()) v = std::abs(v);
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
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1, Eigen::ColMajor> EigenColumnVector;
  typedef Eigen::Stride<1, Eigen::Dynamic>                           Stride;
  typedef Eigen::Map<EigenColumnVector, 0, Stride>                   ColumnMap;

  Numeric* data;
  data = const_cast<Numeric*>(x.data_handle());
  ColumnMap x_map(data, x.size(), Stride(1, x.stride(0)));
  data = const_cast<Numeric*>(y.data_handle());
  ColumnMap y_map(data, y.size(), Stride(1, y.stride(0)));

  y_map = matpack::eigen_map(M) * x_map;
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
void transpose_mult(StridedVectorView y, const Sparse& M, StridedConstVectorView x) {
  // Check dimensions:
  assert(y.size() == static_cast<Size>(M.ncols()));
  assert(static_cast<Size>(M.nrows()) == x.size());

  // Typedefs for Eigen interface
  typedef Eigen::Matrix<Numeric, 1, Eigen::Dynamic, Eigen::RowMajor> EigenColumnVector;
  typedef Eigen::Stride<1, Eigen::Dynamic>                           Stride;
  typedef Eigen::Map<EigenColumnVector, 0, Stride>                   ColumnMap;

  Numeric* data;
  data = const_cast<Numeric*>(x.data_handle());
  ColumnMap x_map(data, x.size(), Stride(1, x.stride(0)));
  data = const_cast<Numeric*>(y.data_handle());
  ColumnMap y_map(data, y.size(), Stride(1, y.stride(0)));

  y_map = x_map * matpack::eigen_map(M);
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
void mult(StridedMatrixView A, const Sparse& B, const StridedConstMatrixView& C) {
  // Check dimensions:
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

  // Typedefs for Eigen interface
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrix;
  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>                           Stride;
  typedef Eigen::Map<EigenMatrix, 0, Stride>                                      MatrixMap;

  Index row_stride, column_stride;
  row_stride    = C.stride(0);
  column_stride = C.stride(1);

  Numeric* data;
  data = const_cast<Numeric*>(C.data_handle());
  Stride    c_stride(row_stride, column_stride);
  MatrixMap C_map(data, C.nrows(), C.ncols(), c_stride);

  row_stride    = A.stride(0);
  column_stride = A.stride(1);
  data          = A.data_handle();
  Stride    a_stride(row_stride, column_stride);
  MatrixMap A_map(data, A.nrows(), A.ncols(), a_stride);

  A_map = matpack::eigen_map(B) * C_map;
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
void mult(StridedMatrixView A, const StridedConstMatrixView& B, const Sparse& C) {
  // Check dimensions:
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

  // Typedefs for Eigen interface.
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrix;
  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>                           Stride;
  typedef Eigen::Map<EigenMatrix, 0, Stride>                                      MatrixMap;

  Index row_stride, column_stride;
  row_stride    = B.stride(0);
  column_stride = B.stride(1);

  Numeric* data;
  data = const_cast<Numeric*>(B.data_handle());
  Stride    b_stride(row_stride, column_stride);
  MatrixMap B_map(data, B.nrows(), B.ncols(), b_stride);

  row_stride    = A.stride(0);
  column_stride = A.stride(1);
  data          = A.data_handle();
  Stride    a_stride(row_stride, column_stride);
  MatrixMap A_map(data, A.nrows(), A.ncols(), a_stride);

  A_map = B_map * matpack::eigen_map(C);
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

  A = matpack::from_eigen(matpack::eigen_map(B).transpose());
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

  A = matpack::from_eigen(matpack::eigen_map(B) * matpack::eigen_map(C));
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

  A = matpack::from_eigen(matpack::eigen_map(B) + matpack::eigen_map(C));
}

//! Sparse identity matrix
/*!

  Set the given Sparse matrix object to the identity matrix. The matrix must
  be square.

  \param A The matrix to be set to the identity matrix.
*/
void id_mat(Sparse& A) {
  assert(A.ncols() == A.nrows());
  const Index      n = A.nrows();
  std::vector<int> starts(n + 1), cols(n);
  std::iota(starts.begin(), starts.end(), 0);
  std::iota(cols.begin(), cols.end(), 0);
  A = Sparse(n, n, std::move(starts), std::move(cols), std::vector<Numeric>(n, 1.0));
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
  // Check dimensions
  assert(B.ncols() == C.ncols());
  assert(B.nrows() == C.nrows());

  A = matpack::from_eigen(matpack::eigen_map(B) - matpack::eigen_map(C));
}

Sparse matpack::from_eigen(eigen_sparse x) {
  x.makeCompressed();
  const Index nnz = x.nonZeros();
  return Sparse(x.rows(),
                x.cols(),
                std::vector<int>(x.outerIndexPtr(), x.outerIndexPtr() + x.rows() + 1),
                std::vector<int>(x.innerIndexPtr(), x.innerIndexPtr() + nnz),
                std::vector<Numeric>(x.valuePtr(), x.valuePtr() + nnz));
}

Range get_rowindex_for_mblock(const Sparse& sensor_response, const Index& mblock_index) {
  const Index n1y = sensor_response.nrows();
  return Range(n1y * mblock_index, n1y);
}

std::ostream& operator<<(std::ostream& os, const ArrayOfSparse& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

void xml_io_stream<Sparse>::write(std::ostream& os_xml, const Sparse& sparse, bofstream* pbofs, std::string_view name) {
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
  arts::println(os_xml);

  row_tag.write_to_stream(os_xml);
  arts::println(os_xml);

  ArrayOfIndex rowind(sparse.nnz()), colind(sparse.nnz());
  Vector       data(sparse.nnz());
  sparse.list_elements(data, rowind, colind);

  // Write row indices.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      //FIXME: It should be the longer lines
      *pbofs << rowind[i];
    else
      arts::println(os_xml, "{}", rowind[i]);
  }

  close_tag.name = "/RowIndex";
  close_tag.write_to_stream(os_xml);
  arts::println(os_xml);

  col_tag.write_to_stream(os_xml);
  arts::println(os_xml);

  // Write column indices.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      //FIXME: It should be the longer lines
      *pbofs << colind[i];
    else
      arts::println(os_xml, "{}", colind[i]);
  }

  close_tag.name = "/ColIndex";
  close_tag.write_to_stream(os_xml);
  arts::println(os_xml);

  data_tag.write_to_stream(os_xml);
  arts::println(os_xml);
  xml_set_stream_precision(os_xml);

  // Write data.

  for (Index i = 0; i < sparse.nnz(); i++) {
    if (pbofs)
      *pbofs << data[i];
    else
      arts::print(os_xml, "{} ", data[i]);
  }
  arts::println(os_xml);
  close_tag.name = "/SparseData";
  close_tag.write_to_stream(os_xml);
  arts::println(os_xml);

  close_tag.name = type_name;
  close_tag.write_to_end_stream(os_xml);

  arts::println(os_xml);
}

void xml_io_stream<Sparse>::read(std::istream& is_xml, Sparse& sparse, bifstream* pbifs) {
  XMLTag tag;
  Index  nrows, ncols, nnz;

  tag.read_from_stream(is_xml);
  tag.check_name("Sparse");

  tag.get_attribute_value("nrows", nrows);
  tag.get_attribute_value("ncols", ncols);
  sparse.resize(nrows, ncols);

  tag.read_from_stream(is_xml);
  tag.check_name("RowIndex");
  tag.get_attribute_value("nelem", nnz);

  ArrayOfIndex rowind(nnz), colind(nnz);
  Vector       data(nnz);

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
