#pragma once

#include <complex>
#include <concepts>
#include <type_traits>

namespace matpack {
template <typename T>
concept has_nelem = requires(T a) {
  { a.nelem() } -> std::integral;
};

template <typename T>
concept has_ncols = requires(T a) {
  { a.ncols() } -> std::integral;
};

template <typename T>
concept has_nrows = requires(T a) {
  { a.nrows() } -> std::integral;
};

template <typename T>
concept has_npages = requires(T a) {
  { a.npages() } -> std::integral;
};

template <typename T>
concept has_nbooks = requires(T a) {
  { a.nbooks() } -> std::integral;
};

template <typename T>
concept has_nshelves = requires(T a) {
  { a.nshelves() } -> std::integral;
};

template <typename T>
concept has_nvitrines = requires(T a) {
  { a.nvitrines() } -> std::integral;
};

template <typename T>
concept has_nlibraries = requires(T a) {
  { a.nlibraries() } -> std::integral;
};

//! For external class interoperability
namespace external_class {
//! Eigen uses cols() for column index
template <typename T>
concept has_cols = requires(T a) {
  { a.cols() } -> std::integral;
};

//! Eigen uses rows() for column index
template <typename T>
concept has_rows = requires(T a) {
  { a.rows() } -> std::integral;
};
}  // namespace external_class

//! Checks if the type has any accepted types of columns
template <typename T> concept column_keeper = has_ncols<T> or external_class::has_cols<T>;

//! Checks if the type has any accepted types of rows
template <typename T> concept row_keeper = has_nrows<T> or external_class::has_rows<T>;

//! Get a column size from x
constexpr auto column_size(column_keeper auto&& x) {
  using internal_type = decltype(x);
  if constexpr (external_class::has_cols<internal_type>) return x.cols();
  else return x.ncols();
}

//! Get a row size from x
constexpr auto row_size(row_keeper auto&& x) {
  using internal_type = decltype(x);
  if constexpr (external_class::has_rows<internal_type>) return x.rows();
  else return x.nrows();
}

//! A concept overload to remove non std::complex<> from list
template <typename T>
struct is_complex : std::false_type {};

//! A concept whose ::value member is true if this is a complex type
template <std::floating_point T>
struct is_complex<std::complex<T>> : std::true_type {};

//! A concept to state if the type is a floating point or a floating point complex
template <typename T>
concept complex_or_real = 
    std::floating_point<std::remove_cvref_t<T>> or
    is_complex<std::remove_cvref_t<T>>::value;

//! A concept for an Arts vector-like type with access operations
template <typename T>
concept vector_like = requires(T a) {
  { a.size() } -> std::integral;
  { a[0] } -> complex_or_real;
};

//! A concept for any of the Arts vector types
template <typename T>
concept vector = vector_like<T> and requires(T a) {
  { a.nelem() } -> std::integral;
  { a.delem() } -> std::integral;
  { a.selem() } -> std::integral;
  { *a.get_c_array() } -> complex_or_real;
};

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept matrix_like = requires(T a) {
  { column_size(a) } -> std::integral;
  { row_size(a) } -> std::integral;
  { a.size() } -> std::integral;
  { a(0, 0) } -> complex_or_real;
};

//! A concept for any of the Arts matrix types
template <typename T>
concept matrix = matrix_like<T> and has_ncols<T> and has_nrows<T> and requires(T a) {
  { a.drows() } -> std::integral;
  { a.dcols() } -> std::integral;
  { a.selem() } -> std::integral;
  { *a.get_c_array() } -> complex_or_real;
};
}  // namespace matpack
