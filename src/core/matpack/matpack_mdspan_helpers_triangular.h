#pragma once

#include <array>
#include <cassert>
#include <concepts>
#include <type_traits>

#include "matpack_mdspan_common_types.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_view_t.h"

namespace matpack {
//! The number of elements on and below the diagonal of an N x N matrix
[[nodiscard]] constexpr Size triangular_size(Size N) { return N * (N + 1) / 2; }

//! The position of element (i, j), with 0 <= j <= i, in a lower triangle packed row by row
[[nodiscard]] constexpr Index triangular_index(Index i, Index j) {
  assert(0 <= j and j <= i);
  return i * (i + 1) / 2 + j;
}

/** An mdspan layout of the lower triangle of a square matrix, packed row by row
 *
 * Element (i, j) is at triangular_index(i, j).  This is the mapping of
 * std::linalg::layout_blas_packed<lower_triangle_t, row_major_t>, except that
 * only the elements on and below the diagonal, j <= i, may be accessed.
 */
struct layout_triangular {
  template <class Extents> class mapping {
   public:
    using extents_type = Extents;
    using index_type   = typename extents_type::index_type;
    using size_type    = typename extents_type::size_type;
    using rank_type    = typename extents_type::rank_type;
    using layout_type  = layout_triangular;

    static_assert(extents_type::rank() == 2, "A triangular matrix has two dimensions");

    constexpr mapping() noexcept = default;

    constexpr mapping(const extents_type& e) noexcept : exts(e) { assert(e.extent(0) == e.extent(1)); }

    template <class OtherExtents> requires(std::is_constructible_v<extents_type, OtherExtents>)
    constexpr explicit(not std::is_convertible_v<OtherExtents, extents_type>)
        mapping(const mapping<OtherExtents>& other) noexcept
        : mapping(extents_type(other.extents())) {}

    [[nodiscard]] constexpr const extents_type& extents() const noexcept { return exts; }

    [[nodiscard]] constexpr index_type required_span_size() const noexcept {
      return static_cast<index_type>(triangular_size(static_cast<Size>(exts.extent(0))));
    }

    template <std::integral I, std::integral J> [[nodiscard]] constexpr index_type operator()(I i, J j) const noexcept {
      assert(static_cast<index_type>(i) < exts.extent(0));
      return static_cast<index_type>(triangular_index(static_cast<Index>(i), static_cast<Index>(j)));
    }

    [[nodiscard]] static constexpr bool is_always_unique() noexcept { return true; }
    [[nodiscard]] static constexpr bool is_always_exhaustive() noexcept { return true; }
    [[nodiscard]] static constexpr bool is_always_strided() noexcept { return false; }
    [[nodiscard]] static constexpr bool is_unique() noexcept { return true; }
    [[nodiscard]] static constexpr bool is_exhaustive() noexcept { return true; }
    [[nodiscard]] static constexpr bool is_strided() noexcept { return false; }

    template <class OtherExtents>
    [[nodiscard]] friend constexpr bool operator==(const mapping& a, const mapping<OtherExtents>& b) noexcept {
      return a.extents() == b.extents();
    }

   private:
    extents_type exts{};
  };
};

//! A view of a lower triangular matrix
template <typename T> using triangular_view_t = mdspan<T, dextents<Index, 2>, layout_triangular>;

//! A lower triangular matrix of runtime size
template <typename T> class triangular_data_t {
  Index        n{0};
  data_t<T, 1> x{};

 public:
  triangular_data_t() = default;

  //! An N x N lower triangular matrix with all elements set to init
  explicit triangular_data_t(Index N, const T& init = T{})
      : n(N), x(static_cast<Index>(triangular_size(static_cast<Size>(N))), init) {}

  //! The number of rows, and columns
  [[nodiscard]] Index nrows() const { return n; }

  [[nodiscard]] T& operator[](Index i, Index j) {
    assert(i < n);
    return x[triangular_index(i, j)];
  }

  [[nodiscard]] const T& operator[](Index i, Index j) const {
    assert(i < n);
    return x[triangular_index(i, j)];
  }

  [[nodiscard]] triangular_view_t<T>       view() { return triangular_view_t<T>{x.data_handle(), n, n}; }
  [[nodiscard]] triangular_view_t<const T> view() const { return triangular_view_t<const T>{x.data_handle(), n, n}; }

  operator triangular_view_t<T>() { return view(); }
  operator triangular_view_t<const T>() const { return view(); }

  //! The elements as a vector, in the order of triangular_index()
  [[nodiscard]] view_t<T, 1>       flat() { return x; }
  [[nodiscard]] view_t<const T, 1> flat() const { return x; }
};

//! A lower triangular N x N matrix of compile-time size
template <typename T, Size N> struct triangular_cdata_t {
  static constexpr Size ndata = triangular_size(N);

  std::array<T, ndata> data{};

  //! The number of rows, and columns
  [[nodiscard]] static constexpr Index nrows() { return static_cast<Index>(N); }

  [[nodiscard]] constexpr T& operator[](Index i, Index j) {
    assert(i < nrows());
    return data[static_cast<Size>(triangular_index(i, j))];
  }

  [[nodiscard]] constexpr const T& operator[](Index i, Index j) const {
    assert(i < nrows());
    return data[static_cast<Size>(triangular_index(i, j))];
  }

  [[nodiscard]] constexpr triangular_view_t<T> view() { return triangular_view_t<T>{data.data(), nrows(), nrows()}; }

  [[nodiscard]] constexpr triangular_view_t<const T> view() const {
    return triangular_view_t<const T>{data.data(), nrows(), nrows()};
  }

  constexpr operator triangular_view_t<T>() { return view(); }
  constexpr operator triangular_view_t<const T>() const { return view(); }

  //! The elements as a vector, in the order of triangular_index()
  [[nodiscard]] constexpr view_t<T, 1> flat() {
    return view_t<T, 1>{mdview_t<T, 1>{data.data(), std::array{static_cast<Index>(ndata)}}};
  }

  [[nodiscard]] constexpr view_t<const T, 1> flat() const {
    return view_t<const T, 1>{mdview_t<T, 1>{const_cast<T*>(data.data()), std::array{static_cast<Index>(ndata)}}};
  }
};
}  // namespace matpack

using TriangularMatrix          = matpack::triangular_data_t<Numeric>;
using TriangularMatrixView      = matpack::triangular_view_t<Numeric>;
using ConstTriangularMatrixView = matpack::triangular_view_t<const Numeric>;
