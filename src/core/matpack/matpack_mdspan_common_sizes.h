#pragma once

#include <configtypes.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdlib>
#include <tuple>
#include <type_traits>

namespace matpack {
//! The type is an integer
template <typename T>
concept integral = std::integral<std::remove_cvref_t<T>>;

template <typename T>
concept has_extent = requires(T a) {
  { a.extent(0) } -> integral;
};

//! Checks if the type has some extent
template <typename T>
concept has_dimension = requires(T a) {
  { a.dimension(0) } -> integral;
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
  { a.cols() } -> integral;
} and not has_ncols<T>;

//! Eigen uses rows() for column index [exclude nrows]
template <typename T>
concept has_rows = requires(T a) {
  { a.rows() } -> integral;
} and not has_nrows<T>;

// For instance vector and array uses size() [exclude ncols and nelem]
template <typename T>
concept has_size = requires(T a) {
  { a.size() } -> integral;
} and not has_ncols<T> and not has_cols<T>;

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
  {
    std::remove_cvref_t<T>::IsVectorAtCompileTime
  } -> std::convertible_to<Index>;
};

//! Checks if the type has any accepted types of columns
template <typename T>
concept column_keeper =
    has_ncols<T> or has_size<T> or has_cols<T>;

//! Get a column size from x
template <column_keeper U>
constexpr auto column_size(const U& x) {
  if constexpr (has_cols<U>) {
    // Special Eigen::Matrix workaround for vectors
    if constexpr (has_IsVectorAtCompileTime<U>) {
      if constexpr (U::IsVectorAtCompileTime) {
        return std::max(x.rows(), x.cols());
      }
    }
    return x.cols();
  } else if constexpr (has_size<U>)
    return x.size();
  else
    return x.ncols();
}

//! Checks if the type has any accepted types of rows as well as previous sizes
template <typename T>
concept row_keeper = column_keeper<T> and (has_nrows<T> or has_rows<T>);

//! Get a row size from x
template <row_keeper U>
constexpr auto row_size(const U& x) {
  if constexpr (has_rows<U>)
    return x.rows();
  else
    return x.nrows();
}

//! Checks if the type has any accepted types of pages as well as previous sizes
template <typename T>
concept page_keeper = row_keeper<T> and (has_npages<T>);

//! Get a page size from x
template <page_keeper U>
constexpr auto page_size(const U& x) {
  return x.npages();
}

//! Checks if the type has any accepted types of books as well as previous sizes
template <typename T>
concept book_keeper = page_keeper<T> and (has_nbooks<T>);

//! Get a book size from x
template <book_keeper U>
constexpr auto book_size(const U&& x) {
  return x.nbooks();
}

//! Checks if the type has any accepted types of shelves as well as previous sizes
template <typename T>
concept shelf_keeper = book_keeper<T> and (has_nshelves<T>);

//! Get a shelf size from x
template <shelf_keeper U>
constexpr auto shelf_size(const U& x) {
  return x.nshelves();
}

//! Checks if the type has any accepted types of vitrines as well as previous sizes
template <typename T>
concept vitrine_keeper = shelf_keeper<T> and (has_nvitrines<T>);

//! Get a vitrine size from x
template <shelf_keeper U>
constexpr auto vitrine_size(const U& x) {
  return x.nvitrines();
}

//! Checks if the type has any accepted types of libraries as well as previous sizes
template <typename T>
concept library_keeper = vitrine_keeper<T> and (has_nlibraries<T>);

//! Get a library size from x
template <library_keeper U>
constexpr auto library_size(const U& x) {
  return x.nlibraries();
}

//! A rankable multidimensional array
template <typename T>
concept rankable = has_rank<T> or has_NumIndices<T> or
                   has_IsVectorAtCompileTime<T> or has_size<T>;

//! Get the rank of the multidimensional array at compile time
template <rankable T>
constexpr Size rank() {
  if constexpr (has_NumIndices<T>)
    return std::remove_cvref_t<T>::NumIndices;
  else if constexpr (has_IsVectorAtCompileTime<T>)
    return 2 - std::remove_cvref_t<T>::IsVectorAtCompileTime;
  else if constexpr (has_rank<T>)
    return std::remove_cvref_t<T>::rank();
  else if constexpr (has_size<T>)
    return 1;
  else
    return -1;
}

//! Gets the dimension size of some extent
template <rankable T, Size dim = rank<T>()>
constexpr Size dimsize(const T& v, integral auto ind) {
  if constexpr (has_extent<T>) {
    return v.extent(ind);
  } else if constexpr (has_dimension<T>) {
    return v.dimension(ind);
  } else {
    switch (dim - ind - 1) {
      case 0:
        if constexpr (dim > 0) return column_size(v);
      case 1:
        if constexpr (dim > 1) return row_size(v);
      case 2:
        if constexpr (dim > 2) return page_size(v);
      case 3:
        if constexpr (dim > 3) return book_size(v);
      case 4:
        if constexpr (dim > 4) return shelf_size(v);
      case 5:
        if constexpr (dim > 5) return vitrine_size(v);
      case 6:
        if constexpr (dim > 6) return library_size(v);
      default: {
      }
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
template <typename T, Size N>
concept integral_array =
    requires(T a) {
      { a[0] } -> integral;
    } and std::same_as<std::array<std::remove_cvref_t<decltype(T{}[0])>, N>,
                       std::remove_cvref_t<T>>;

//! Checks if the type has a shape
template <typename T, Size N>
concept has_shape = rankable<T> and requires(T a) {
  { a.shape() } -> integral_array<rank<T>()>;
};

//! Get the shape
template <rankable T, Size dim = rank<T>()>
constexpr std::array<Index, dim> mdshape(const T& v) {
  if constexpr (has_shape<T, dim>) {
    return v.shape();
  } else {
    std::array<Index, dim> out;
    for (Size i = 0; i < dim; i++) out[i] = static_cast<Index>(dimsize(v, i));
    return out;
  }
}

template <typename T>
concept has_mdshape = rankable<T> and requires(T a) {
  { mdshape(a) } -> std::same_as<std::array<Index, rank<T>()>>;
};

//! Test if the object can be accessed by an Index
template <typename T>
concept has_index_access = requires(T a) { a[Index{}]; };

//! Thest if the object can be iterated over
template <typename T>
concept is_iterable = requires(T a) {
  std::begin(a);
  std::begin(a) + Index{};
  std::end(a);
  std::size(a);
};

//! Get a positional value
template <rankable T, Size dim = rank<T>()>
constexpr auto mdvalue(const T& v, const std::array<Index, dim>& pos) {
  if constexpr (dim == 1) {
    if constexpr (has_index_access<T>)
      return v[pos[0]];
    else
      return *(std::begin(v) + pos[0]);
  } else {
    return std::apply(
        [&v](auto... inds) {
          if constexpr (requires { v(inds...); }) {
            return v(inds...);
          } else {
            return v[inds...];
          }
        },
        pos);
  }
}

//! Test that you can convert the underlying value type of an object to a
//! specific type
template <typename U, typename T>
concept mdvalue_type_compatible = requires(U a) {
  { mdvalue(a, mdshape(a)) } -> std::convertible_to<T>;
};

template <Size N>
constexpr Size mdsize(const std::array<Index, N>& shape) {
  Size size = 1;
  for (Size i = 0; i < N; i++) size *= shape[i];
  return size;
}

template <Size N>
constexpr std::array<Index, N> mdpos(const std::array<Index, N>& shape,
                                     Index i) {
  std::array<Index, N> pos{};

  Index n = mdsize<N>(shape);
  for (Size j = 0; j < N; j++) {
    n              /= shape[j];
    const auto res  = std::div(i, n);
    pos[j]          = res.quot;
    i               = res.rem;
  }

  return pos;
}

template <rankable T, Size dim = rank<T>()>
constexpr auto mdvalue(const T& v,
                       const std::array<Index, dim>& shape,
                       const Index i) {
  return mdvalue(v, mdpos(shape, i));
}
}  // namespace matpack
