#pragma once

#include "config.h"
#include "debug.h"

#include <experimental/mdspan>

#include <array>
#include <complex>
#include <concepts>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace stdx = std::experimental;

//! The base number type of matpack
using Numeric = NUMERIC;

//! The base integer type of matpack
using Index = INDEX;

using Complex = std::complex<Numeric>;

//! A type that denotes all values should be accessed
using Joker = stdx::full_extent_t;

//! A value to denote that all values should be accessed
inline constexpr Joker joker = stdx::full_extent;

//! Contains all the matpack core concepts and types
namespace matpack {

//! A type that denotes some values should be accessed
struct matpack_strided_access;

//! The basic view type
template <typename T, Index N, bool constant, bool strided>
class matpack_view;

//! The basic data type
template <typename T, Index N>
class matpack_data;

//! The constexpr view type
template <typename T, bool constant, Index... alldim> struct matpack_constant_view;

//! The constexpr data type
template <typename T, Index... alldim> struct matpack_constant_data;

//! Helper bool
template <typename T>
inline constexpr bool is_matpack_strided_access =
    std::is_same_v<std::remove_cvref_t<T>, matpack_strided_access>;

//! The type is an integer
template <typename T>
concept integral = std::integral<std::remove_cvref_t<T>>;

//! A type that is useful as an access operation
template <typename T>
concept access_operator =
    integral<T> or std::is_same_v<std::remove_cvref_t<T>, Joker> or
    is_matpack_strided_access<T>;

//! Holds true for mdspan(s) that are always continuous in memory, i.e., exhaustive_mdspan and not strided_mdspan
template <typename T>
concept is_always_exhaustive_v = std::remove_cvref_t<T>::is_always_exhaustive();

//! The type has basic arithmetic properties
template <typename T>
concept arithmetic = std::is_arithmetic_v<std::remove_cvref_t<T>>;

//! Checks that the type is a pure arithmetic complex type
template <typename T>
concept complex_type =
    requires(T a) {
      { a.real() } -> arithmetic;
      { a.imag() } -> arithmetic;
    } and sizeof(std::remove_cvref_t<decltype(T{}.real())>) * 2 == sizeof(T) and
    sizeof(std::remove_cvref_t<decltype(T{}.imag())>) * 2 == sizeof(T);

//! Helper to get the subtype of a complex type (e.g., Numeric from std::complex<Numeric>)
template <complex_type T>
using complex_subtype = std::remove_cvref_t<decltype(T{}.real())>;

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
  { std::remove_cvref_t<T>::rank() } -> integral;
};

//! The multidimensional array has a NumIndices
template <typename T>
concept has_NumIndices = requires(T) {
  { std::remove_cvref_t<T>::NumIndices } -> integral;
};

//! The multidimensional array has IsVectorAtCompileTime to indicate (if true)
//! if it is a vector or (if false) if it is a matrix
template <typename T>
concept has_IsVectorAtCompileTime = requires(T) {
  { std::remove_cvref_t<T>::IsVectorAtCompileTime } -> std::convertible_to<Index>;
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

//! Returns the 1D size of a given shape
template <Index N>
constexpr Index mdsize(const std::array<Index, N>& shape) {
  return std::reduce(shape.begin(), shape.end(), Index{1}, std::multiplies<>());
}

//! A rankable multidimensional array
template <typename T> concept rankable = has_rank<T> or has_NumIndices<T> or has_IsVectorAtCompileTime<T> or has_size<T> or has_nelem<T>;

//! Get the rank of the multidimensional array at compile time
template <rankable T>
constexpr auto rank() {
  if constexpr (has_NumIndices<T>) return std::remove_cvref_t<T>::NumIndices;
  else if constexpr (has_IsVectorAtCompileTime<T>) return 2 - std::remove_cvref_t<T>::IsVectorAtCompileTime;
  else if constexpr (has_rank<T>) return std::remove_cvref_t<T>::rank();
  else if constexpr (has_size<T>) return 1;
  else if constexpr (has_nelem<T>) return 1;
  else return -1;
}

//! Gets the dimension size of some extent
template <rankable T, auto dim = rank<T>()> constexpr Index dimsize(const T& v, integral auto ind) {
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

//! The type can call dimsize
template <typename T>
concept has_dimsize = requires(T a) {
  { dimsize(a, 0) } -> std::same_as<Index>;
};

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

//! Test that the object can have a shape
template <typename T>
concept has_mdshape = rankable<T> and requires(T a) {
  { mdshape(a) } -> std::same_as<std::array<Index, rank<T>()>>;
};

//! Test if the object can be accessed by an Index
template <typename T>
concept has_index_access = requires(T a) {
  a[Index{}];
};

//! Thest if the object can be iterated over
template <typename T>
concept is_iterable = requires(T a) {
  std::begin(a);
  std::begin(a) + Index{};
  std::end(a);
  std::size(a);
};

//! Get a positional value
template <rankable T, auto dim = rank<T>()>
constexpr auto mdvalue(const T& v, const std::array<Index, dim>& pos) {
  if constexpr (dim == 1) {
    if constexpr (has_index_access<T>) return v[pos[0]];
    else if constexpr (is_iterable<T>) return *(std::begin(v)+pos[0]);
    else {/**/}
  } else return std::apply([&v](auto... inds){return v(inds...);}, pos);
}

//! The type is arithmetic or complex (for interface reasons, this cannot be
//! const, volatile, reference, or any combination)
template <typename T>
concept math_type = arithmetic<T> or complex_type<T>;

//! Test that you can convert the underlying value type of an object to a
//! specific type
template <typename U, typename T>
concept mdvalue_type_compatible = requires(U a) {
                                    {
                                      mdvalue(a, mdshape(a))
                                      } -> std::convertible_to<T>;
                                  };

//! The underlying value type of the type T
template <typename T>
using matpack_value_type = typename std::remove_cvref_t<T>::value_type;

//! Helper struct to test that all the types that define it has the same value type
template <typename first, typename second, typename... rest>
class same_value_type {
/** Test the types that they are the same
 * 
 * @return true if all of first, second, rest..., have the same matpack_value_type
 * @return false otherwise
 */
static constexpr bool same() {
  using L = matpack_value_type<first>;
  using R = matpack_value_type<second>;
  if constexpr (sizeof...(rest) == 0) return std::same_as<L, R>;
  else return std::same_as<L, R> and same_value_type<second, rest...>::value;
}

public:
  static constexpr bool value = same();
};

//! Helper bool for ::value on same_value_type
template <typename first, typename second, typename... rest>
inline constexpr bool same_value_type_v = same_value_type<first, second, rest...>::value;


//! Test that the type U is exactly matpack_view<T, N, false/true, false/true> (w/o qualifiers) for any constant, strided
template <typename U, typename T, Index N>
concept ranked_matpack_view = std::same_as<std::remove_cvref_t<U>, matpack_view<T, N, false, false>> or
                              std::same_as<std::remove_cvref_t<U>, matpack_view<T, N, false, true>> or
                              std::same_as<std::remove_cvref_t<U>, matpack_view<T, N, true, false>> or
                              std::same_as<std::remove_cvref_t<U>, matpack_view<T, N, true, true>>;

//! Test that the type U is exactly matpack_view<U::value_type, N, constant, strided> (w/o qualifiers) for any constant, strided
template <typename U, Index N>
concept strict_rank_matpack_view = ranked_matpack_view<U, matpack_value_type<U>, N>;

//! Test that the type U is any matpack_view<...> (w/o qualifiers)
template <typename U>
concept any_matpack_view = rankable<U> and ranked_matpack_view<U, matpack_value_type<U>, rank<U>()>;


//! Test that the type U is exactly matpack_data<T, N> (w/o qualifiers)
template <typename U, typename T, Index N>
concept ranked_matpack_data = std::same_as<std::remove_cvref_t<U>, matpack_data<T, N>>;

//! Test that the type U is exactly matpack_data<U::value_type, N> (w/o qualifiers)
template <typename U, Index N>
concept strict_rank_matpack_data = ranked_matpack_data<U, matpack_value_type<U>, N>;

//! Test that the type U is any matpack_data<...> (w/o qualifiers)
template <typename U>
concept any_matpack_data = rankable<U> and ranked_matpack_data<U, matpack_value_type<U>, rank<U>()>;


/** Test that the type U is a matpack_constant_view<T, true/false, ...> of a given rank and type (w/o qualifiers) using magic function */
template <typename U, typename T, Index N>
concept ranked_matpack_constant_view =
  std::remove_cvref_t<U>::template magic_test_of_matpack_constant_view<T, N>();

/** Test that the type U is a matpack_constant_view<T, true/false, ...> of a given rank (w/o qualifiers) */
template <typename U, Index N>
concept strict_rank_matpack_constant_view = ranked_matpack_constant_view<U, matpack_value_type<U>, N>;

/** Test that the type U is a matpack_constant_view<T, true/false, ...> (w/o qualifiers) */
template <typename U>
concept any_matpack_constant_view = ranked_matpack_constant_view<U, matpack_value_type<U>, rank<U>()>;


/** Test that the type U is a matpack_constant_data<T, ...> of a given rank and type (w/o qualifiers) using magic function */
template <typename U, typename T, Index N>
concept ranked_matpack_constant_data =
  std::remove_cvref_t<U>::template magic_test_of_matpack_constant_data<T, N>();

/** Test that the type U is a matpack_constant_data<T, ...> of a given rank (w/o qualifiers) */
template <typename U, Index N>
concept strict_rank_matpack_constant_data = ranked_matpack_constant_data<U, matpack_value_type<U>, N>;

/** Test that the type U is a matpack_constant_data<T, ...> (w/o qualifiers) */
template <typename U>
concept any_matpack_constant_data = ranked_matpack_constant_data<U, matpack_value_type<U>, rank<U>()>;


//! Test that the type U is exactly matpack_data<T, N> or matpack_view<T, N, ...> or their constant friends (w/o qualifiers)
template <typename U, typename T, Index N>
concept ranked_matpack_type = ranked_matpack_data<U, T, N> or
                              ranked_matpack_view<U, T, N> or 
                              ranked_matpack_constant_data<U, T, N> or
                              ranked_matpack_constant_view<U, T, N>;

//! Test that the type U is strictly ranked matpack_data or matpack_view (w/o qualifiers)
template<typename U, Index N>
concept strict_rank_matpack_type = ranked_matpack_type<U, matpack_value_type<U>, N>;

//! Test that the type is any matpack_view or matpack_data
template <typename U>
concept any_matpack_type = ranked_matpack_type<U, matpack_value_type<U>, rank<U>()>;


//! Test if the type can be converted to a ranked matpack_data
template <typename U, typename T, Index N>
concept matpack_convertible = not any_matpack_type<U> and has_mdshape<U> and mdvalue_type_compatible<U, T> and rank<U>() == N;


//! Helper to rank a mdspan
template <Index N> using ranking = stdx::dextents<Index, N>;

//! Helper to define a strided layout mdspan
using strided = stdx::layout_stride;

//! An exhaustive layout mdspan in matpack style
template <typename T, Index N>
using exhaustive_mdspan = stdx::mdspan<T, ranking<N>>;

//! A strided mdspan in matpack style
template <typename T, Index N>
using strided_mdspan = stdx::mdspan<T, ranking<N>, strided>;
}  // namespace matpack
