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
#include <optional>
#include <ranges>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

#include <xml.h>
#include <xml_io_stream.h>

#include "array.h"
#include "matpack_mdspan_data_t.h"

//! The Sparse class.
/*!

  A sparse matrix in compressed row storage.  The elements of row r are at
  [row_starts()[r], row_starts()[r + 1]) of column_indices() and values(), with
  ascending and unique columns in each row.  This is the storage of an Eigen
  row-major sparse matrix, and matpack_sparse_eigen.h maps it as one where Eigen
  is needed, so that Eigen is not needed where the matrix is merely used.

*/

struct Sparse {
  static constexpr bool matpack_type{true};

  // Constructors:
  Sparse();
  Sparse(Index r, Index c);

  /** A matrix from its compressed row storage
   *
   * @param r The number of rows
   * @param c The number of columns
   * @param row_starts The start of each row in the other two, and their size last, so r + 1 elements
   * @param column_indices The ascending and unique columns of the elements of each row
   * @param values The values of the elements
   */
  Sparse(Index r, Index c, std::vector<int> row_starts, std::vector<int> column_indices, std::vector<Numeric> values);

  Sparse(const Sparse& other);
  Sparse(Sparse&& other) noexcept;
  Sparse& operator=(const Sparse& other);
  Sparse& operator=(Sparse&& other) noexcept;

  void split(Index offset, Index nrows);

  // Insert functions
  void insert_row(Index r, Vector v);
  void insert_elements(Index nnz, const ArrayOfIndex& rowind, const ArrayOfIndex& colind, StridedConstVectorView data);

  // Resize function:
  void resize(Index r, Index c);

  // Member functions:
  [[nodiscard]] bool  empty() const;
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
  [[nodiscard]] Numeric  ro(Index r, Index c) const;
  [[nodiscard]] Numeric  operator[](Index r, Index c) const;

  // Arithmetic operators:
  Sparse& operator+=(const Sparse& x);
  Sparse& operator-=(const Sparse& x);

  // Scaling operators:
  Sparse& operator*=(Numeric x);
  Sparse& operator/=(Numeric x);

  // Conversion to Dense Matrix:
  explicit operator Matrix() const;

  // Matrix data access
  void list_elements(Vector& element_values, ArrayOfIndex& row_indices, ArrayOfIndex& column_indices) const;

  Numeric* get_element_pointer();
  int*     get_column_index_pointer();
  int*     get_row_start_pointer();

  //! The start of each row in column_indices() and values(), and their size last
  [[nodiscard]] std::span<const int> row_starts() const;

  //! The ascending and unique columns of the elements of each row
  [[nodiscard]] std::span<const int> column_indices() const;

  //! The values of the elements
  [[nodiscard]] std::span<const Numeric> values() const;
  [[nodiscard]] std::span<Numeric>       values();

  // Friends:
  friend std::ostream& operator<<(std::ostream& os, const Sparse& v);

 private:
  Index                nr{0};
  Index                nc{0};
  std::vector<int>     row_start{0};
  std::vector<int>     col_index;
  std::vector<Numeric> vals;
};

// Functions for general matrix operations
void abs(Sparse& A, const Sparse& B);

void mult(StridedVectorView y, const Sparse& M, StridedConstVectorView x);

void transpose_mult(StridedVectorView y, const Sparse& M, StridedConstVectorView x);

void mult(StridedMatrixView A, const Sparse& B, const StridedConstMatrixView& C);

void mult(StridedMatrixView A, const StridedConstMatrixView& B, const Sparse& C);

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
Range get_rowindex_for_mblock(const Sparse& sensor_response, const Index& imblock);

/** An array of sparse matrices. */
using ArrayOfSparse = Array<Sparse>;

std::ostream& operator<<(std::ostream& os, const ArrayOfSparse& a);

// Stored entries only, including explicit zeros. Structural changes invalidate
// iterators; changing an existing value through the range is allowed.
namespace matpack {
template <bool Const> class sparse_element_range : public std::ranges::view_interface<sparse_element_range<Const>> {
  using source_type = std::conditional_t<Const, const Sparse, Sparse>;
  source_type* source;

 public:
  explicit sparse_element_range(source_type& value) : source(&value) {}
  // The storage is fetched once, so that iterating does not call the out-of-line accessors
  class iterator {
    using value_pointer = std::conditional_t<Const, const Numeric*, Numeric*>;
    using reference     = std::conditional_t<Const, const Numeric&, Numeric&>;

    const int*    starts{};
    const int*    cols{};
    value_pointer vals{};
    Index         nrows{};
    Index         row{};
    Index         pos{};

    void skip_finished_rows() {
      while (row < nrows and pos >= starts[row + 1]) ++row;
    }

   public:
    using difference_type  = std::ptrdiff_t;
    using value_type       = std::tuple<Index, Index, Numeric>;
    using iterator_concept = std::input_iterator_tag;
    iterator()             = default;
    explicit iterator(source_type& value)
        : starts(value.row_starts().data()),
          cols(value.column_indices().data()),
          vals(value.values().data()),
          nrows(value.nrows()) {
      skip_finished_rows();
    }
    auto operator*() const { return std::tuple<Index, Index, reference>{row, cols[pos], vals[pos]}; }
    iterator& operator++() {
      ++pos;
      skip_finished_rows();
      return *this;
    }
    void operator++(int) { ++*this; }
    bool operator==(std::default_sentinel_t) const { return row >= nrows; }
  };
  iterator                begin() const { return iterator{*source}; }
  std::default_sentinel_t end() const { return {}; }
};
}  // namespace matpack

inline auto operator|(Sparse& value, matpack::elemwise_r) { return matpack::sparse_element_range<false>{value}; }
inline auto operator|(const Sparse& value, matpack::elemwise_r) { return matpack::sparse_element_range<true>{value}; }
void        operator|(Sparse&&, matpack::elemwise_r)       = delete;
void        operator|(const Sparse&&, matpack::elemwise_r) = delete;

template <> struct std::formatter<Sparse> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const Sparse& v, FmtContext& ctx) const {
    std::string_view first = tags.sep();
    std::string_view sep   = ""sv;

    tags.add_if_bracket(ctx, "["sv);

    for (const auto [row, col, value] : v | by_elem) {
      tags.format(ctx, std::exchange(sep, first));
      tags.add_if_bracket(ctx, "["sv);
      tags.format(ctx, row, sep, col, sep, value);
      tags.add_if_bracket(ctx, "]"sv);
    }

    tags.add_if_bracket(ctx, "]"sv);
    return ctx.out();
  }
};


template <> struct xml_io_stream<Sparse> {
  static constexpr std::string_view type_name = "Sparse"sv;

  static void write(std::ostream& os, const Sparse& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, Sparse& x, bifstream* pbifs = nullptr);
};

#endif
