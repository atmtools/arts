#pragma once

#include "debug.h"

#include <array>
#include <complex>
#include <concepts>
#include <type_traits>

namespace matpack {
template <typename T>
concept matpack_type = std::remove_cvref_t<T>::matpack_type;

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

// For instance vector and array uses size()
template <typename T>
concept has_size = requires(T a) {
  { a.size() } -> std::integral;
} and not matpack_type<T>;
}  // namespace external_class

//! Checks if the type has any accepted types of columns
template <typename T> concept column_keeper = has_ncols<T> or has_nelem<T> or external_class::has_size<T> or external_class::has_cols<T>;

//! Get a column size from x
constexpr auto column_size(column_keeper auto&& x) {
  using internal_type = decltype(x);
  if constexpr (external_class::has_cols<internal_type>) return std::forward<internal_type>(x).cols();
  else if constexpr (external_class::has_size<internal_type>) return std::forward<internal_type>(x).size();
  else if constexpr (has_nelem<internal_type>) return std::forward<internal_type>(x).nelem();
  else return std::forward<internal_type>(x).ncols();
}

//! Checks if the type has any accepted types of rows as well as previous sizes
template <typename T> concept row_keeper = column_keeper<T> and (has_nrows<T> or external_class::has_rows<T>);

//! Get a row size from x
constexpr auto row_size(row_keeper auto&& x) {
  using internal_type = decltype(x);
  if constexpr (external_class::has_rows<internal_type>) return std::forward<internal_type>(x).rows();
  else return std::forward<internal_type>(x).nrows();
}

//! Checks if the type has any accepted types of pages as well as previous sizes
template <typename T> concept page_keeper = row_keeper<T> and (has_npages<T>);

//! Get a page size from x
constexpr auto page_size(page_keeper auto&& x) {
  using internal_type = decltype(x);
  return std::forward<internal_type>(x).npages();
}

//! Checks if the type has any accepted types of books as well as previous sizes
template <typename T> concept book_keeper = page_keeper<T> and (has_nbooks<T>);

//! Get a book size from x
constexpr auto book_size(book_keeper auto&& x) {
  using internal_type = decltype(x);
  return std::forward<internal_type>(x).nbooks();
}

//! Checks if the type has any accepted types of shelves as well as previous sizes
template <typename T> concept shelf_keeper = book_keeper<T> and (has_nshelves<T>);

//! Get a shelf size from x
constexpr auto shelf_size(shelf_keeper auto&& x) {
  using internal_type = decltype(x);
  return std::forward<internal_type>(x).nshelves();
}

//! Checks if the type has any accepted types of vitrines as well as previous sizes
template <typename T> concept vitrine_keeper = shelf_keeper<T> and (has_nvitrines<T>);

//! Get a vitrine size from x
constexpr auto vitrine_size(vitrine_keeper auto&& x) {
  using internal_type = decltype(x);
  return std::forward<internal_type>(x).nvitrines();
}

//! Checks if the type has any accepted types of libraries as well as previous sizes
template <typename T> concept library_keeper = vitrine_keeper<T> and (has_nlibraries<T>);

//! Get a library size from x
constexpr auto library_size(library_keeper auto&& x) {
  using internal_type = decltype(x);
  return std::forward<internal_type>(x).nlibraries();
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
  { column_size(a) } -> std::integral;
  { a[0] } -> complex_or_real;
};

//! A concept for any of the Arts vector types
template <typename T>
concept vector = matpack_type<T> and vector_like<T> and requires(T a) {
  { a.nelem() } -> std::integral;
  { a.delem() } -> std::integral;
  { a.selem() } -> std::integral;
  { *a.get_c_array() } -> complex_or_real;
};

//! A concept precluding Arts vector objects but allowing things similar to vectors
template <typename T>
concept vector_like_not_vector = vector_like<T> and not matpack_type<T>;

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept matrix_like = requires(T a) {
  { row_size(a) } -> std::integral;
  { a(0, 0) } -> complex_or_real;
};

//! A concept for any of the Arts matrix types
template <typename T>
concept matrix = matpack_type<T> and matrix_like<T> and has_ncols<T> and has_nrows<T> and requires(T a) {
  { a.drows() } -> std::integral;
  { a.dcols() } -> std::integral;
  { a.selem() } -> std::integral;
  { *a.get_c_array() } -> complex_or_real;
};

//! A concept precluding Arts matrix objects but allowing things similar to matrices
template <typename T>
concept matrix_like_not_matrix = matrix_like<T> and not matpack_type<T>;

//! Matrix or vector
template <typename T>
concept matrix_or_vector = matrix<T> or vector<T>;

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept tensor3_like = requires(T a) {
  { page_size(a) } -> std::integral;
  { a(0, 0, 0) } -> complex_or_real;
};

//! A concept precluding Arts types but allowing the tensor-like object
template <typename T>
concept tensor3_like_not_tensor3 = tensor3_like<T> and not matpack_type<T>;

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept tensor4_like = requires(T a) {
  { book_size(a) } -> std::integral;
  { a(0, 0, 0, 0) } -> complex_or_real;
};

//! A concept precluding Arts types but allowing the tensor-like object
template <typename T>
concept tensor4_like_not_tensor4 = tensor4_like<T> and not matpack_type<T>;

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept tensor5_like = requires(T a) {
  { shelf_size(a) } -> std::integral;
  { a(0, 0, 0, 0, 0) } -> complex_or_real;
};

//! A concept precluding Arts types but allowing the tensor-like object
template <typename T>
concept tensor5_like_not_tensor5 = tensor5_like<T> and not matpack_type<T>;

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept tensor6_like = requires(T a) {
  { vitrine_size(a) } -> std::integral;
  { a(0, 0, 0, 0, 0, 0) } -> complex_or_real;
};

//! A concept precluding Arts types but allowing the tensor-like object
template <typename T>
concept tensor6_like_not_tensor6 = tensor6_like<T> and not matpack_type<T>;

//! A concept for an Arts matrix-like type with access operations
template <typename T>
concept tensor7_like = requires(T a) {
  { library_size(a) } -> std::integral;
  { a(0, 0, 0, 0, 0, 0, 0) } -> complex_or_real;
};

//! A concept precluding Arts types but allowing the tensor-like object
template <typename T>
concept tensor7_like_not_tensor7 = tensor7_like<T> and not matpack_type<T>;

template <typename T>
concept matpack_like =
    tensor7_like<T> or tensor6_like<T> or tensor5_like<T> or tensor4_like<T> or
    tensor3_like<T> or matrix_like<T> or vector_like<T>;

template <matpack_like T> constexpr std::size_t dim() {
  return library_keeper<T> + vitrine_keeper<T> + shelf_keeper<T> +
         book_keeper<T> + page_keeper<T> + row_keeper<T> + column_keeper<T>;
}

//! Creates a shape array (note that this might create compile time errors if N is larger than x allows for)
template <class IndexType, std::size_t N>
constexpr std::array<IndexType, N> shape(matpack_like auto &&x) {
  using T = decltype(x);
  static_assert(N > 0 and N < 8, "Out of range");

  static_assert((N == 1 and dim<T>() == 2) or (N == dim<T>()),
                "Dimensions must agree, or a matrix-like that is runtime "
                "defined as vector be used");

  if constexpr (N == 7)
    return {static_cast<IndexType>(library_size(std::forward<T>(x))),
            static_cast<IndexType>(vitrine_size(std::forward<T>(x))),
            static_cast<IndexType>(shelf_size(std::forward<T>(x))),
            static_cast<IndexType>(book_size(std::forward<T>(x))),
            static_cast<IndexType>(page_size(std::forward<T>(x))),
            static_cast<IndexType>(row_size(std::forward<T>(x))),
            static_cast<IndexType>(column_size(std::forward<T>(x)))};
  else if constexpr (N == 6)
    return {static_cast<IndexType>(vitrine_size(std::forward<T>(x))),
            static_cast<IndexType>(shelf_size(std::forward<T>(x))),
            static_cast<IndexType>(book_size(std::forward<T>(x))),
            static_cast<IndexType>(page_size(std::forward<T>(x))),
            static_cast<IndexType>(row_size(std::forward<T>(x))),
            static_cast<IndexType>(column_size(std::forward<T>(x)))};
  else if constexpr (N == 5)
    return {static_cast<IndexType>(shelf_size(std::forward<T>(x))),
            static_cast<IndexType>(book_size(std::forward<T>(x))),
            static_cast<IndexType>(page_size(std::forward<T>(x))),
            static_cast<IndexType>(row_size(std::forward<T>(x))),
            static_cast<IndexType>(column_size(std::forward<T>(x)))};
  else if constexpr (N == 4)
    return {static_cast<IndexType>(book_size(std::forward<T>(x))),
            static_cast<IndexType>(page_size(std::forward<T>(x))),
            static_cast<IndexType>(row_size(std::forward<T>(x))),
            static_cast<IndexType>(column_size(std::forward<T>(x)))};
  else if constexpr (N == 3)
    return {static_cast<IndexType>(page_size(std::forward<T>(x))),
            static_cast<IndexType>(row_size(std::forward<T>(x))),
            static_cast<IndexType>(column_size(std::forward<T>(x)))};
  else if constexpr (N == 2)
    return {static_cast<IndexType>(row_size(std::forward<T>(x))),
            static_cast<IndexType>(column_size(std::forward<T>(x)))};
  else if constexpr (N == 1) {
    if constexpr (dim<T>() == 2) {
      ARTS_USER_ERROR_IF(
          static_cast<IndexType>(std::min(column_size(std::forward<T>(x)),
                                          row_size(std::forward<T>(x)))) not_eq
              static_cast<IndexType>(1),
          "Cannot perform vector conversion for matrix: ", x)
      return {static_cast<IndexType>(std::max(column_size(std::forward<T>(x)),
                                              row_size(std::forward<T>(x))))};
    } else
      return {static_cast<IndexType>(column_size(std::forward<T>(x)))};
  }
}
}  // namespace matpack
