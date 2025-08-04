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

#include <iosfwd>
#include <utility>

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

#include <xml.h>
#include <xml_io_stream.h>

#include "array.h"
#include "matpack_mdspan_data_t.h"

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
  Sparse(const Sparse& other);
  Sparse(Sparse&& other) noexcept;
  Sparse& operator=(const Sparse& other);
  Sparse& operator=(Sparse&& other) noexcept;

  void split(Index offset, Index nrows);

  // Insert functions
  void insert_row(Index r, Vector v);
  void insert_elements(Index nnz,
                       const ArrayOfIndex& rowind,
                       const ArrayOfIndex& colind,
                       StridedConstVectorView data);

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
  static Sparse diagonal(StridedConstVectorView v);

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
  [[nodiscard]] Numeric operator[](Index r, Index c) const;

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

void mult(StridedVectorView y, const Sparse& M, StridedConstVectorView x);

void transpose_mult(StridedVectorView y,
                    const Sparse& M,
                    StridedConstVectorView x);

void mult(StridedMatrixView A,
          const Sparse& B,
          const StridedConstMatrixView& C);

void mult(StridedMatrixView A,
          const StridedConstMatrixView& B,
          const Sparse& C);

void mult(Sparse& A, const Sparse& B, const Sparse& C);

void add(Sparse& A, const Sparse& B, const Sparse& C);

void sub(Sparse& A, const Sparse& B, const Sparse& C);

void transpose(Sparse& A, const Sparse& B);

void id_mat(Sparse& A);

/** Returns the "range" of *y* corresponding to a measurement block

    @param[in]   sensor_response    As the WSV.
    @param[in]   mblock_index       Index of the measurement block.

    @return  The range.

    @author Patrick Eriksson 
    @date   2009-10-16
 */
Range get_rowindex_for_mblock(const Sparse& sensor_response,
                              const Index& imblock);

/** An array of sparse matrices. */
using ArrayOfSparse = Array<Sparse>;

std::ostream& operator<<(std::ostream& os, const ArrayOfSparse& a);

template <>
struct std::formatter<Sparse> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Sparse& v, FmtContext& ctx) const {
    using Iter = Eigen::SparseMatrix<Numeric, Eigen::RowMajor>::InnerIterator;

    std::string_view first = tags.sep();
    std::string_view sep   = ""sv;

    tags.add_if_bracket(ctx, '[');

    for (int k = 0; k < v.matrix.outerSize(); ++k) {
      for (Iter it(v.matrix, k); it; ++it) {
        tags.format(ctx, std::exchange(sep, first));
        tags.add_if_bracket(ctx, '[');
        tags.format(ctx, it.row(), sep, it.col(), sep, it.value());
        tags.add_if_bracket(ctx, ']');
      }
    }

    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct xml_io_stream<Sparse> {
  static constexpr std::string_view type_name = "Sparse"sv;

  static void write(std::ostream& os,
                    const Sparse& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, Sparse& x, bifstream* pbifs = nullptr);
};

#endif
