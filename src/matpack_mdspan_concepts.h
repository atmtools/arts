#pragma once

#include "matpack.h"
#include "matpack_concepts.h"

#include <mdspan/include/experimental/mdspan>

#include <concepts>
#include <type_traits>

#pragma GCC diagnostic push
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-comma-subscript"
#else
#pragma GCC diagnostic ignored "-Wcomma-subscript"
#endif
namespace matpack::md {
//! The type has basic arithmetic properties
template <typename T>
concept arithmetic = std::is_arithmetic_v<std::remove_cvref_t<T>>;

//! The type is complex
template <typename T>
concept complex_type = requires(T a) {
  { a.real() } -> arithmetic;
  { a.imag() } -> arithmetic;
};

//! The type is arithmetic or complex (for interface reasons, this cannot be const, volatile, reference, or any combination)
template <typename T>
concept math_type = (arithmetic<T> or complex_type<T>) and 
    std::same_as<T, std::remove_cvref_t<T>>;

//! The core data owning type
template <math_type T, Index N> class simple_data;

//! A continous view of the data held by simple_data
template <math_type T, Index N, bool constant> class simple_view;

//! A strided view of the data help by simple_data
template <math_type T, Index N, bool constant> class strided_view;

//! A type that denotes all values should be accessed
using Joker = std::experimental::full_extent_t;

//! A value to denote that all values should be accessed
inline constexpr Joker joker = std::experimental::full_extent;

//! A type that denotes some values should be accessed
struct strided_access;

//! The type is an integer
template <typename T>
concept integral = std::integral<std::remove_cvref_t<T>>;

//! A type that is useful as an access operation
template <typename T>
concept access_operator =
    integral<T> or std::is_same_v<std::remove_cvref_t<T>, Joker> or
    std::is_same_v<std::remove_cvref_t<T>, strided_access>;

//! A view that can be mutated
template <typename U, typename T, Index N>
concept mutable_view =
    std::same_as<std::remove_cvref_t<U>, simple_view<T, N, false>> or
    std::same_as<std::remove_cvref_t<U>, strided_view<T, N, false>>;

//! A view or a non-const of the original data
template <typename U, typename T, Index N>
concept any_mutable =
    mutable_view<U, T, N> or
    std::same_as<std::remove_volatile_t<std::remove_reference_t<U>>,
                 simple_data<T, N>>;

//! A view that cannot be mutated
template <typename U, typename T, Index N>
concept const_view =
    std::same_as<std::remove_cvref_t<U>, simple_view<T, N, true>> or
    std::same_as<std::remove_cvref_t<U>, strided_view<T, N, true>>;

//! Any view of the data of a specific type and size
template <typename U, typename T, Index N>
concept any_view = mutable_view<U, T, N> or const_view<U, T, N>;

//! Any matpack data type of a specific type of a specific size
template <typename U, typename T, Index N>
concept any_md = any_view<U, T, N> or
                 std::same_as<std::remove_cvref_t<U>, simple_data<T, N>>;

//! Any view of the data with a derived type of a specific size
template <typename U, Index N>
concept any_same_view = any_view<U, typename U::value_type, N>;

//! Any matpack data type with a derived type of a specific size
template <typename U, Index N>
concept any_same_md = any_md<U, typename U::value_type, N>;

//! Any strided data with a derived type of a specific size
template <typename U, Index N>
concept any_strided = std::same_as<std::remove_cvref_t<U>, strided_view<typename U::value_type, N, true>> or 
                      std::same_as<std::remove_cvref_t<U>, strided_view<typename U::value_type, N, false>>;

//! Any continous view with a derived type of a specific size
template <typename U, Index N>
concept any_exhaustive = any_md<U, typename U::value_type, N> and not any_strided<U, N>;

//! Helper to get the subtype of a complex type (e.g., Numeric from std::complex<Numeric>)
template <complex_type T>
using complex_subtype = std::remove_cvref_t<decltype(T{}.real())>;

//! A vector or view thereof
template <typename U> concept matpack_vector = any_same_md<U, 1>;

//! A matrix or view thereof
template <typename U> concept matpack_matrix = any_same_md<U, 2>;

//! A matrix or vector, or view of either
template <typename U> concept matpack_matrix_or_vector = matpack_matrix<U> or matpack_vector<U>;

//! Checks if the type has some extent
template <typename T>
concept has_extent = requires(T a) {
  { a.extent(0) } -> integral;
};

//! Checks if the type has some extent
template <typename T>
concept has_dimension = requires(T a) {
  { a.dimension(0) } -> integral;
};

//! Checks if the type has some nelem
template <typename T>
concept has_nelem = requires(T a) {
  { a.nelem() } -> integral;
};

//! Checks if the type has some ncols
template <typename T>
concept has_ncols = requires(T a) {
  { a.ncols() } -> integral;
};

//! Checks if the type has some nrows
template <typename T>
concept has_nrows = requires(T a) {
  { a.nrows() } -> integral;
};

//! Checks if the type has some npages
template <typename T>
concept has_npages = requires(T a) {
  { a.npages() } -> integral;
};

//! Checks if the type has some nbooks
template <typename T>
concept has_nbooks = requires(T a) {
  { a.nbooks() } -> integral;
};

//! Checks if the type has some nshelves
template <typename T>
concept has_nshelves = requires(T a) {
  { a.nshelves() } -> integral;
};

//! Checks if the type has some nvitrines
template <typename T>
concept has_nvitrines = requires(T a) {
  { a.nvitrines() } -> integral;
};

//! Checks if the type has some nlibraries
template <typename T>
concept has_nlibraries = requires(T a) {
  { a.nlibraries() } -> integral;
};

//! Eigen uses cols() for column index [exclude ncols and nelem]
template <typename T>
concept has_cols = requires(T a) {
  { a.cols() } ->  integral;
} and not has_ncols<T> and not has_nelem<T>;

//! Eigen uses rows() for column index [exclude nrows]
template <typename T>
concept has_rows = requires(T a) {
  { a.rows() } -> integral;
} and not has_nrows<T>;

// For instance vector and array uses size() [exclude ncols and nelem]
template <typename T>
concept has_size = requires(T a) {
  { a.size() } -> integral;
} and not has_ncols<T> and not has_nelem<T> and not has_cols<T>;

//! The multidimensional array has a rank
template <typename T>
concept has_rank = requires(T) {
  { T::rank() } -> integral;
};

//! The multidimensional array has a NumIndices
template <typename T>
concept has_NumIndices = requires(T) {
  { T::NumIndices } -> integral;
};

//! The multidimensional array has IsVectorAtCompileTime to indicate (if true)
//! if it is a vector or (if false) if it is a matrix
template <typename T>
concept has_IsVectorAtCompileTime = requires(T) {
  { T::IsVectorAtCompileTime } -> std::convertible_to<Index>;
};

//! Checks if the type has any accepted types of columns
template <typename T> concept column_keeper = has_ncols<T> or has_nelem<T> or has_size<T> or has_cols<T>;

//! Get a column size from x
template<column_keeper U>
constexpr auto column_size(const U& x) {
  if constexpr (has_cols<U>) {
      // Special Eigen::Matrix workaround for vectors
    if constexpr (has_IsVectorAtCompileTime<U>) {
      if constexpr (U::IsVectorAtCompileTime) {
        return std::max(x.rows(), x.cols());
      }
    }
    return x.cols();
  }
  else if constexpr (has_size<U>) return x.size();
  else if constexpr (has_nelem<U>) return x.nelem();
  else return x.ncols();
}

//! Checks if the type has any accepted types of rows as well as previous sizes
template <typename T> concept row_keeper = column_keeper<T> and (has_nrows<T> or has_rows<T>);

//! Get a row size from x
template<row_keeper U>
constexpr auto row_size(const U& x) {
  if constexpr (has_rows<U>) return x.rows();
  else return x.nrows();
}

//! Checks if the type has any accepted types of pages as well as previous sizes
template <typename T> concept page_keeper = row_keeper<T> and (has_npages<T>);

//! Get a page size from x
template<page_keeper U>
constexpr auto page_size(const U& x) {
  return x.npages();
}

//! Checks if the type has any accepted types of books as well as previous sizes
template <typename T> concept book_keeper = page_keeper<T> and (has_nbooks<T>);

//! Get a book size from x
template<book_keeper U>
constexpr auto book_size(const U&& x) {
  return x.nbooks();
}

//! Checks if the type has any accepted types of shelves as well as previous sizes
template <typename T> concept shelf_keeper = book_keeper<T> and (has_nshelves<T>);

//! Get a shelf size from x
template<shelf_keeper U>
constexpr auto shelf_size(const U& x) {
  return x.nshelves();
}

//! Checks if the type has any accepted types of vitrines as well as previous sizes
template <typename T> concept vitrine_keeper = shelf_keeper<T> and (has_nvitrines<T>);

//! Get a vitrine size from x
template<shelf_keeper U>
constexpr auto vitrine_size(const U& x) {
  return x.nvitrines();
}

//! Checks if the type has any accepted types of libraries as well as previous sizes
template <typename T> concept library_keeper = vitrine_keeper<T> and (has_nlibraries<T>);

//! Get a library size from x
template<library_keeper U>
constexpr auto library_size(const U& x) {
  return x.nlibraries();
}

//! A rankable multidimensional array
template <typename T> concept rankable = has_rank<T> or has_NumIndices<T> or has_IsVectorAtCompileTime<T>;

//! Get the rank of the multidimensional array at compile time
template <rankable T>
consteval auto rank() {
  if constexpr (has_NumIndices<T>) return T::NumIndices;
  else if constexpr (has_IsVectorAtCompileTime<T>) return 2 - T::IsVectorAtCompileTime;
  else if constexpr (has_size<T>) return 1;
  else return T::rank();
}

template <rankable T, auto dim = rank<T>()> Index dimsize(const T& v, integral auto ind) {
  if constexpr (has_extent<T>) {
    return v.extent(ind);
  } else if constexpr (has_dimension<T>) {
    return v.dimension(ind);
  } else {
    switch (dim - ind - 1) {
      case 0: if constexpr (dim > 0) return column_size(v);
      case 1: if constexpr (dim > 1) return row_size(v);
      case 2: if constexpr (dim > 2) return page_size(v);
      case 3: if constexpr (dim > 3) return book_size(v);
      case 4: if constexpr (dim > 4) return shelf_size(v);
      case 5: if constexpr (dim > 5) return vitrine_size(v);
      case 6: if constexpr (dim > 6) return library_size(v);
      default: {}
    }
    return 0;
  }
}

//! Ensure that T is of the type std::array<integral, N>
template <typename T, Index N>
concept integral_array = 
  requires(T a) { {a[0] } -> integral; } and
  std::same_as<std::array<std::remove_cvref_t<decltype(T{}[0])>, N>,
               std::remove_cvref_t<T>>;

//! Checks if the type has a shape
template <typename T, Index N>
concept has_shape = rankable<T> and requires(T a) {
  { a.shape() } -> integral_array<rank<T>()>;
};

//! Get the shape
template <rankable T, auto dim = rank<T>()>
constexpr std::array<Index, dim> mdshape(const T& v) {
  if constexpr (has_shape<T, dim>) {
    return v.shape();
  } else {
    std::array<Index, dim> out;
    for (Index i=0; i<dim; i++) out[i] = static_cast<Index>(dimsize(v, i));
    return out;
  }
}

//! Checks that a type is a pointer to a math_type
template <typename T>
concept math_type_ptr = math_type<std::remove_cvref_t<std::remove_pointer_t<T>>> and std::is_pointer_v<T>;

//! Checks that a type has a data_handle
template <typename T>
concept has_data_handle = requires (T a) {
  { a.data_handle() } -> math_type_ptr;
};

//! Checks that a type has a data
template <typename T>
concept has_data = requires (T a) {
  { a.data() } -> math_type_ptr;
};

//! Has some raw data pointer
template <typename T>
concept has_data_pointer = has_data_handle<T> or has_data<T>;

//! Get the data pointer
template <has_data_pointer T>
constexpr auto mddata(const T& v) {
  if constexpr (has_data<T>) return v.data();
  else return v.data_handle();
}

//! The type's stride can be had with stride()
template <typename T>
concept has_stride = requires(T a) {
  { a.stride(0) } -> integral;
};

//! The type's stride can be had with stride()
template <typename T>
concept has_innerStride = requires(T a) {
  { a.innerStride() } -> integral;
  { a.outerStride() } -> integral;
};

//! Get the stride sizes of some index
template <rankable T, auto dim = rank<T>()>
constexpr Index stridesize(const T& v, integral auto ind) {
  if constexpr (has_stride<T>) return v.stride(ind);
  else if constexpr (has_size<T>) return 1;
  else if constexpr (has_innerStride<T>) {
    static_assert(dim < 3 and dim > 0);
    
    if constexpr (dim == 1) {
      return v.innerStride();
    } else {
      const bool t = true; // TRANSPOSE AND false
      if (ind == 1) {
        return t ? v.innerStride() : v.outerStride();
      } else {
        return t ? v.outerStride() : v.innerStride();
      }
    }
  }
  else return 0;
}

template <rankable T, auto dim = rank<T>()>
constexpr std::array<Index, dim> mdstride(const T& v) {
  std::array<Index, dim> out;
  for (decltype(dim) i=0; i<dim; i++) out[i] = static_cast<Index>(stridesize(v, i));
  return out;
}

template <typename T>
concept matpack_convertible = rankable<T> and requires (T a) {
  { mddata(a) } -> math_type_ptr;
  { mdshape(a) } -> integral_array<rank<T>()>;
} and not any_same_md<T, rank<T>()>;
}  // namespace matpack::md
#pragma GCC diagnostic pop
