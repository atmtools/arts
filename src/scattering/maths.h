/** \file scattering/maths.h
 *
 * Type aliases and helper function for math operations on
 * scattering data.
 *
 * @author Simon Pfreundschuh, 2020-2022
 */
#pragma once

#include <Eigen/CXX11/Tensor>
#include <Eigen/Core>
#include <Eigen/MatrixFunctions>
#include <iostream>
#include <memory>
#include <type_traits>
#include <complex>

#include <arts.h>
#include <matpackI.h>
#include <matpackVII.h>
#include <double_imanip.h>

namespace scattering {
namespace math {

using Index = Eigen::Index;

template <int rank>
using IndexArray = std::array<Eigen::DenseIndex, rank>;

//
// Vectors
//

/** Variable-length vector.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
template <typename Scalar>
using Vector = Eigen::Matrix<Scalar, 1, -1, Eigen::RowMajor>;
/** Fixed-length vector.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
template <typename Scalar, long int N>
using VectorFixedSize = Eigen::Matrix<Scalar, 1, N, Eigen::RowMajor>;
/** Variable-length vector that doesn't own its data.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
template <typename Scalar>
using VectorPtr = std::shared_ptr<Vector<Scalar>>;
template <typename Scalar>
using ConstVectorPtr = std::shared_ptr<const Vector<Scalar>>;
template <typename Scalar>
using VectorMap = Eigen::Map<Vector<Scalar>>;
template <typename Scalar>
using VectorMapDynamic =
    Eigen::Map<Vector<Scalar>, 0, Eigen::Stride<1, Eigen::Dynamic>>;
template <typename Scalar>
using ConstVectorMap = Eigen::Map<const Vector<Scalar>>;
template <typename Scalar>
using ConstVectorMapDynamic =
    Eigen::Map<const Vector<Scalar>, 0, Eigen::Stride<1, Eigen::Dynamic>>;
template <typename Scalar>
using VectorRef = Eigen::Ref<Vector<Scalar>>;
template <typename Scalar>
using ConstVectorRef = Eigen::Ref<const Vector<Scalar>>;

//
// Matrices
//

/** A variable-size matrix containing coefficients of type Scalar.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
template <typename Scalar>
using Matrix = Eigen::Matrix<Scalar, -1, -1, Eigen::RowMajor>;
template <typename Scalar, long int N>
/** Matrix with fixed number of rows.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
using MatrixFixedRows =
    Eigen::Matrix<Scalar, -1, N, (N > 1) ? Eigen::RowMajor : Eigen::ColMajor>;
template <typename Scalar, long int M, long int N>
/** Matrix with fixed number of rows.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
using MatrixFixedSize =
    Eigen::Matrix<Scalar, M, N, (N > 1) ? Eigen::RowMajor : Eigen::ColMajor>;
/** A matrix that doesn't own its data.
 * @tparam Scalar The type used to represent coefficients of the matrix.
 */
template <typename Scalar>
using MatrixMap = Eigen::Map<Matrix<Scalar>>;
template <typename Scalar>
using MatrixMapDynamic =
    Eigen::Map<Matrix<Scalar>,
               Eigen::RowMajor,
               Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
template <typename Scalar>
using ConstMatrixMap = Eigen::Map<const Matrix<Scalar>>;
template <typename Scalar>
using ConstMatrixMapDynamic =
    Eigen::Map<const Matrix<Scalar>,
               Eigen::RowMajor,
               Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
template <typename Scalar>
using MatrixRef = Eigen::Ref<Matrix<Scalar>>;
template <typename Scalar>
using ConstMatrixRef = Eigen::Ref<const Matrix<Scalar>>;

//
// Tensors
//

template <typename Scalar, int rank>
using Tensor = Eigen::Tensor<Scalar, rank, Eigen::RowMajor>;
template <typename Scalar, int rank>
using TensorPtr = std::shared_ptr<Eigen::Tensor<Scalar, rank, Eigen::RowMajor>>;
template <typename Scalar, int rank>
using TensorMap =
    Eigen::TensorMap<Eigen::Tensor<Scalar, rank, Eigen::RowMajor>>;
template <typename Scalar, int rank>
using ConstTensorMap =
    Eigen::TensorMap<const Eigen::Tensor<Scalar, rank, Eigen::RowMajor>>;

//
// General map type.
//

template <typename Tensor>
struct Map {
  using type = Eigen::TensorMap<Tensor>;
};

template <typename Scalar, int rows, int cols, int options>
struct Map<Eigen::Matrix<Scalar, rows, cols, options>> {
  using type = Eigen::Map<Eigen::Matrix<Scalar, rows, cols, options>>;
};

template <typename Scalar, int rows, int cols, int options>
struct Map<const Eigen::Matrix<Scalar, rows, cols, options>> {
  using type = Eigen::Map<const Eigen::Matrix<Scalar, rows, cols, options>>;
};

template <typename Derived>
struct Map<Eigen::TensorMap<Derived>> {
  using type = typename Map<Derived>::type;
};

template <typename Derived>
struct Map<Eigen::Map<Derived>> {
  using type = typename Map<Derived>::type;
};

template <typename Derived>
struct Map<const Eigen::TensorMap<Derived>> {
  using type = typename Map<Derived>::type;
};

template <typename Derived>
struct Map<const Eigen::Map<Derived>> {
  using type = typename Map<Derived>::type;
};

////////////////////////////////////////////////////////////////////////////////
// Tensor indexing
////////////////////////////////////////////////////////////////////////////////

//
// Helper trait to reduce rank.
//

template <typename T, int n_indices>
struct IndexResult {
  using CoeffType = decltype(((T *)1)->operator()({}));
  using Scalar = typename std::remove_reference<
      typename std::remove_cv<CoeffType>::type>::type;

  static constexpr int rank = T::NumIndices;
  static constexpr bool is_const =
      std::is_const<CoeffType>::value || std::is_const<T>::value;

  using NonConstReturnType = TensorMap<Scalar, rank - n_indices>;
  using ConstReturnType = ConstTensorMap<Scalar, rank - n_indices>;
  using ReturnType = typename std::
      conditional<is_const, ConstReturnType, NonConstReturnType>::type;

  using type =
      typename std::conditional<(n_indices < rank), ReturnType, Scalar>::type;
};

//
// Template helper to create map to sub tensor.
//

template <typename T, int rank_in, int n_indices>
struct TensorIndexer {
  using Tensor = typename std::remove_reference<T>::type;
  using TensorIndex = typename Tensor::Index;
  using ReturnType = typename IndexResult<Tensor, n_indices>::type;
  using Scalar = typename Tensor::Scalar;
  static constexpr int rank_out = rank_in - n_indices;

  __attribute__((always_inline)) static inline ReturnType get(
      T &t, std::array<TensorIndex, n_indices> index_array) {
    auto dimensions_in = t.dimensions();
    std::array<TensorIndex, rank_in> index{};
    for (int i = 0; i < n_indices; ++i) {
      index[i] = index_array[i];
    }
    Eigen::DSizes<TensorIndex, rank_out> dimensions_out{};
    for (int i = 0; i < rank_out; ++i) {
      dimensions_out[i] = dimensions_in[n_indices + i];
    }
    auto offset = dimensions_in.IndexOfRowMajor(index);
    return ReturnType(t.data() + offset, dimensions_out);
  }
};

template <typename T, int rank_in>
struct TensorIndexer<T, rank_in, rank_in> {
  using Tensor = typename std::remove_reference<T>::type;
  using TensorIndex = typename Tensor::Index;
  using Scalar = typename Tensor::Scalar;

  __attribute__((always_inline)) static inline Scalar get(
      T &t, std::array<TensorIndex, rank_in> index_array) {
    return t(index_array);
  };
};

// pxx :: export
// pxx :: instance(["4", "scattering::eigen::Tensor<double, 4>"])
// pxx :: instance(["3", "scattering::eigen::Tensor<double, 4>"])
// pxx :: instance(["2", "scattering::eigen::Tensor<double, 4>"])
// pxx :: instance(["1", "scattering::eigen::Tensor<double, 4>"])
// pxx :: instance(["3", "scattering::eigen::Tensor<double, 3>"])
// pxx :: instance(["2", "scattering::eigen::Tensor<double, 3>"])
// pxx :: instance(["1", "scattering::eigen::Tensor<double, 3>"])
// pxx :: instance(["2", "scattering::eigen::Tensor<double, 2>"])
// pxx :: instance(["1", "scattering::eigen::Tensor<double, 2>"])
// pxx :: instance(["1", "scattering::eigen::Tensor<double, 1>"])
template <size_t N, typename T>
__attribute__((always_inline)) inline auto tensor_index(
    T &t, std::array<typename T::Index, N> indices) ->
    typename IndexResult<T, N>::type {
  return TensorIndexer<T, T::NumIndices, N>::get(t, indices);
}

//
// Tensor map
//

namespace detail {

template <int N, int i = 0>
struct MapOverDimensionsImpl {
  template <typename TensorTypeOut,
            typename TensorTypeIn,
            typename F,
            typename... Indices>
  inline static void run(TensorTypeIn &&out,
                         TensorTypeOut &&in,
                         F f,
                         Indices... indices) {
    for (math::Index j = 0; j < in.dimension(i); ++j) {
      MapOverDimensionsImpl<N, i + 1>::run(std::forward<TensorTypeIn>(out),
                                           std::forward<TensorTypeOut>(in),
                                           f,
                                           indices...,
                                           j);
    }
  }
};

template <int N>
struct MapOverDimensionsImpl<N, N> {
  template <typename TensorTypeOut,
            typename TensorTypeIn,
            typename F,
            typename... Indices>
  inline static void run(TensorTypeOut &&out,
                         TensorTypeIn &&in,
                         F f,
                         Indices... indices) {
    for (math::Index j = 0; j < in.dimension(N - 1); ++j) {
      f(tensor_index<sizeof...(Indices)>(out, {indices...}),
        tensor_index<sizeof...(Indices)>(in, {indices...}));
    }
  }
};

}  // namespace detail

template <int N, typename TensorTypeOut, typename TensorTypeIn, typename F>
void map_over_dimensions(TensorTypeOut &&out, TensorTypeIn &&in, F f) {
  detail::MapOverDimensionsImpl<N>::run(
      std::forward<TensorTypeOut>(out), std::forward<TensorTypeIn>(in), f);
}

/** Create matrix map onto tensor.
 *
 * @param T The scalar type of the tensor.
 * @param t The rank-2 tensor to create a matrix map of.
 * @return A MatrixMap onto the memory of the tensor.
 */
template <typename T>
auto to_matrix_map(T &t) -> MatrixMap<typename T::Scalar> {
  using Scalar = typename T::Scalar;
  static_assert(T::NumIndices == 2,
                "Tensor must be of rank 2 to be convertible to matrix.");
  return MatrixMap<Scalar>(t.data(), t.dimension(0), t.dimension(1));
}

/** Create const matrix map onto tensor.
 *
 * @param T The scalar type of the tensor.
 * @param t The rank-2 tensor to create a matrix map of.
 * @return A ConstMatrixMap onto the memory of the tensor.
 */
template <typename T>
auto to_matrix_map(const T &t) -> ConstMatrixMap<typename T::Scalar> {
  using Scalar = typename T::Scalar;
  static_assert(T::NumIndices == 2,
                "Tensor must be of rank 2 to be convertible to matrix.");
  return ConstMatrixMap<Scalar>(t.data(), t.dimension(0), t.dimension(1));
}

/** Create vector map onto tensor.
 *
 * @param T The scalar type of the tensor.
 * @param t The rank-1 tensor to create a vector map of.
 * @return A vector onto the memory of the tensor.
 */
template <typename T>
auto to_vector_map(T &t) -> VectorMap<typename T::Scalar> {
  using Scalar = typename T::Scalar;
  static_assert(T::NumIndices == 1,
                "Tensor must be of rank 1 to be convertible to vector.");
  return VectorMap<Scalar>(t.data(), t.dimension(0));
}

/** Create vector map onto tensor.
 *
 * @param T The scalar type of the tensor.
 * @param t The rank-1 tensor to create a vector map of.
 * @return A vector onto the memory of the tensor.
 */
template <typename T>
auto to_vector_map(const T &t) -> ConstVectorMap<typename T::Scalar> {
  using Scalar = typename T::Scalar;
  static_assert(T::NumIndices == 1,
                "Tensor must be of rank 1 to be convertible to vector.");
  return ConstVectorMap<Scalar>(t.data(), t.dimension(0));
}

//
// Access to sub-matrix of tensor.
//

namespace detail {
template <typename TensorType>
auto calculate_strides(const TensorType &t) {
  constexpr int rank = TensorType::NumIndices;
  auto dimensions = t.dimensions();
  decltype(dimensions) strides{};
  Index stride = 1;
  for (int i = 0; i < rank; ++i) {
    strides[rank - i - 1] = stride;
    stride *= dimensions[rank - i - 1];
  }
  return strides;
}

}  // namespace detail

// pxx :: export
// pxx :: instance(["1", "3", "scattering::eigen::Tensor<double, 5>", "std::array<int, 3>"])
template <int m,
          int n,
          typename TensorType,
          typename IndexArray = std::array<typename TensorType::Index,
                                           TensorType::NumIndices - 2>>
auto inline get_submatrix(TensorType &t, IndexArray matrix_index) ->
    typename std::conditional<
        !std::is_const<decltype(*(std::declval<TensorType>().data()))>::value &&
            !std::is_const<TensorType>::value,
        MatrixMapDynamic<typename TensorType::Scalar>,
        ConstMatrixMapDynamic<typename TensorType::Scalar>>::type {
  using CoeffType = decltype(*(std::declval<TensorType>().data()));
  using ResultType = typename std::conditional<
      !std::is_const<CoeffType>::value && !std::is_const<TensorType>::value,
      MatrixMapDynamic<typename TensorType::Scalar>,
      ConstMatrixMapDynamic<typename TensorType::Scalar>>::type;

  using TensorIndex = typename TensorType::Index;
  constexpr int rank = TensorType::NumIndices;

  // Extend matrix dimensions to tensor dimension
  // to calculate offset.
  auto dimensions_in = t.dimensions();
  std::array<TensorIndex, rank> index{0};
  int dimension_index = 0;
  for (int i = 0; i < rank; ++i) {
    if ((i != m) && (i != n)) {
      index[i] = matrix_index[dimension_index];
      dimension_index++;
    }
  }
  auto offset = dimensions_in.IndexOfRowMajor(index);
  auto strides = detail::calculate_strides(t);

  Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> matrix_strides{strides[m],
                                                               strides[n]};
  auto map = ResultType(
      t.data() + offset, dimensions_in[m], dimensions_in[n], matrix_strides);
  return map;
}

// pxx :: export
// pxx :: instance(["3", "scattering::eigen::Tensor<double, 5>", "std::array<int, 4>"])
template <int m,
          typename TensorType,
          typename IndexArray = std::array<typename TensorType::Index,
                                           TensorType::NumIndices - 1>>
auto inline get_subvector(TensorType &t, IndexArray vector_index) ->
    typename std::conditional<
        !std::is_const<decltype(*(std::declval<TensorType>().data()))>::value &&
            !std::is_const<TensorType>::value,
        VectorMapDynamic<typename TensorType::Scalar>,
        ConstVectorMapDynamic<typename TensorType::Scalar>>::type {
  using CoeffType = decltype(*(std::declval<TensorType>().data()));
  using ResultType = typename std::conditional<
      !std::is_const<CoeffType>::value && !std::is_const<TensorType>::value,
      VectorMapDynamic<typename TensorType::Scalar>,
      ConstVectorMapDynamic<typename TensorType::Scalar>>::type;

  using TensorIndex = typename TensorType::Index;
  constexpr int rank = TensorType::NumIndices;

  // Extend vector dimensions to tensor dimension
  // to calculate offset.
  auto dimensions_in = t.dimensions();
  std::array<TensorIndex, rank> index{0};
  int dimension_index = 0;
  for (int i = 0; i < rank; ++i) {
    if (i != m) {
      index[i] = vector_index[dimension_index];
      dimension_index++;
    }
  }
  auto offset = dimensions_in.IndexOfRowMajor(index);
  auto strides = detail::calculate_strides(t);

  Eigen::Stride<1, Eigen::Dynamic> vector_strides{1, strides[m]};
  auto map = ResultType(t.data() + offset, 1, dimensions_in[m], vector_strides);
  return map;
}

////////////////////////////////////////////////////////////////////////////////
// Dimension counter
////////////////////////////////////////////////////////////////////////////////

/** Tensor index counter to loop over elements.
 *
 * The dimensions counter is a helper class that allows looping over the
 * indices in a tensor in a single loop. It loops over all elements in the
 * tensor in row-major order, meaning that the last index increased with
 * every increment of the counter. The current tensor index can be
 * accessed through the coordinates member.
 */
template <int rank>
struct DimensionCounter {
  /** Create counter
   * @param dims Array containing the dimension of the tensor to loop over.
   */
  DimensionCounter(std::array<Eigen::DenseIndex, rank> dims) {
    dimensions = dims;
  }



  /// Increment the counter.
  DimensionCounter &operator++() {
    for (int i = rank - 1; i >= 0; i--) {
      coordinates[i]++;
      if (coordinates[i] == dimensions[i]) {
        coordinates[i] = 0;
        if (i == 0) {
          exhausted_ = true;
        }
      } else {
        break;
      }
    }
    return *this;
  }

  /// Have all elements been looped over?
  operator bool() { return !exhausted_; }

  bool exhausted_ = false;
  /// The index of the current tensor element.
  std::array<Eigen::DenseIndex, rank> coordinates{0};
  std::array<Eigen::DenseIndex, rank> dimensions{0};
};

template <int rank>
std::ostream &operator<<(std::ostream &out, const DimensionCounter<rank> &c) {
  out << "Dimension counter: ";
  for (int i = 0; i < rank - 1; ++i) {
    out << c.coordinates[i] << ", ";
  }
  out << c.coordinates[rank - 1] << std::endl;
  return out;
}

//////////////////////////////////////////////////////////////////////
// Adaptive copy
//////////////////////////////////////////////////////////////////////

/// Create a zero tensor.
template <typename Scalar, typename... Types>
Tensor<Scalar, sizeof...(Types)> zeros(Types... dimensions) {
  constexpr int rank = sizeof...(Types);
  std::array<Index, rank> dimension_array({dimensions...});
  return Tensor<Scalar, sizeof...(Types)>(dimension_array).setZero();
}

/// Add dummy dimension into tensor.
template <int... dims, typename Scalar, int rank>
Tensor<Scalar, rank + sizeof...(dims)> unsqueeze(
    const Tensor<Scalar, rank> &tensor) {
  constexpr int new_rank = rank + sizeof...(dims);
  auto dimensions = tensor.dimensions();
  std::array<Index, new_rank> new_dimensions;
  std::array<Index, sizeof...(dims)> trivial_dimensions{dims...};
  std::sort(trivial_dimensions.begin(), trivial_dimensions.end());
  auto dimension_iterator = dimensions.begin();
  auto trivial_dimension_iterator = trivial_dimensions.begin();

  for (int i = 0; i < new_rank; ++i) {
    if (i == *trivial_dimension_iterator) {
      new_dimensions[i] = 1;
      trivial_dimension_iterator++;
    } else {
      new_dimensions[i] = *dimension_iterator;
      dimension_iterator++;
    }
  }
  return tensor.reshape(new_dimensions);
}

/// Cycle dimension that exceed tensor size.
template <typename Scalar, int rank>
Tensor<Scalar, rank> cycle_dimensions(const Tensor<Scalar, rank> &t) {
  std::array<Index, rank> dimensions = {};
  for (int i = 0; i < rank; ++i) {
    dimensions[i] = (i + 1) % rank;
  }
  return t.shuffle(dimensions);
}

/** Size-adaptive copy generator.
 *
 * An Eigen tensor generator that copies a tensor into a tensor of arbitrary
 * size and fill excess elements with 0.0.
 */
template <typename TensorType>
struct CopyGenerator {
  static constexpr int rank = TensorType::NumIndices;
  using Scalar = typename TensorType::Scalar;
  CopyGenerator(const TensorType &from_) : from(from_){};

  Scalar operator()(const std::array<Index, rank> &coordinates) const {
    for (size_t i = 0; i < rank; ++i) {
      if (coordinates[i] >= from.dimension(i)) {
        return 0.0;
      }
    }
    return from(coordinates);
  }

  const TensorType &from;
};

/** Copy content of a tensor into another one.
 *
 * Copies a tensor into another, possibly differently-sized, tensor.
 * Missing elements are filled with zeros.
 * @dest The tensor into which to copy the elements into.
 * @source The tensor from which to copy the elements.
 */
template <typename TensorType1, typename TensorType2>
void copy(TensorType1 &dest, const TensorType2 &source) {
  dest = dest.generate(CopyGenerator<TensorType2>(source));
}

//////////////////////////////////////////////////////////////////////
// Misc. helper functions
//////////////////////////////////////////////////////////////////////

template <typename VectorType1, typename VectorType2>
bool equal(const VectorType1 &left, const VectorType2 &right) {
  if (left.size() != right.size()) {
    return false;
  }
  return left == right;
}

template <typename Scalar>
auto colatitudes(const Vector<Scalar> &input) {
  return input.unaryExpr([](Scalar x) { return -1.0 * cos(x); });
}

template <Index n, typename TensorType>
std::array<Index, n> get_dimensions(const TensorType &t) {
  std::array<Index, n> result{};
  for (Index i = 0; i < n; ++i) {
    if (i < TensorType::NumIndices) {
      result[i] = t.dimension(i);
    } else {
      result[i] = 0;
    }
  }
  return result;
}

template <typename Scalar>
bool equal(Scalar a, Scalar b, Scalar epsilon = 1e-6) {
  return std::abs(a - b) <=
         ((std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a)) * epsilon);
}


template <typename Scalar>
bool small(Scalar a, Scalar epsilon = 1e-6) {
    return std::abs(a) < epsilon;
}

template <typename Scalar>
Scalar save_acos(Scalar a, Scalar epsilon = 1e-6) {
  if (equal(a, 1.0, epsilon)) {
    return 0.0;
  }
  if (equal(a, -1.0, epsilon)) {
    return M_PI;
  }
  return acos(a);
}

/*** Indirect sorting of a vector.
 *
 * Indirectly sorts elements in a vector and returns a vector of
 * sorted indices.
 *
 * @param v: The vector to sort.
 * @return A vector containing the indices that sort the given
 *    vector into ascending order.
 */
template <typename Scalar>
    math::Vector<Eigen::Index> indirect_sort(const math::Vector<Scalar>& v) {
    math::Vector<Eigen::Index> indices;
    indices.setLinSpaced(v.size(), 0, v.size() - 1);

    auto comp = [&v](size_t i, size_t j) { return v[i] < v[j]; };
    std::sort(indices.begin(), indices.end(), comp);
    return indices;
}

/*** Write Eigen::Vector to output stream.
 *
 * @param output The output stream to write the vector to.
 * @param vector The Eigen::Vector to serialize.
 * @return reference to the output stream.
 */
template<typename Scalar>
    std::ostream& serialize(std::ostream &output, const Vector<Scalar> &input)  {
    Index size = input.size();
    output.write(reinterpret_cast<char const*>(&size), sizeof(Index));
    output.write(reinterpret_cast<char const*>(input.data()), input.size() * sizeof(Scalar));
    return output;
}

/*** Read Eigen::Vector from output stream.
 *
 * @param input The input stream from which to read the vector.
 * @param output The vector object in which to store the results.
 * @return reference to the input stream.
 */
template<typename Scalar>
    std::istream& deserialize(std::istream &input, Vector<Scalar> &output)  {
    Index size;
    input.read(reinterpret_cast<char *>(&size), sizeof(Index));
    output.resize(1, size);
    input.read(reinterpret_cast<char *>(output.data()), output.size() * sizeof(Scalar));
    return input;
}

/*** Write Eigen::Matrix to output stream.
 *
 * @param output The output stream to write the matrix to.
 * @param input The Eigen::Matrix to serialize.
 * @return reference to the output stream.
 */
template<typename Scalar>
    std::ostream& serialize(std::ostream &output, const Matrix<Scalar> &input)  {
    Index rows{input.rows()}, cols{input.cols()};
    output.write(reinterpret_cast<char const*>(&rows), sizeof(Index));
    output.write(reinterpret_cast<char const*>(&cols), sizeof(Index));
    output.write(reinterpret_cast<char const*>(input.data()), input.size() * sizeof(Scalar));
    return output;
}

/*** Read Eigen::Matrix from input stream.
 *
 * @param input The input stream from which to read the matrix.
 * @param output The matrix object in which to store the results.
 * @return reference to the input stream.
 */
template<typename Scalar>
    std::istream& deserialize(std::istream &input, Matrix<Scalar> &output)  {
    Index rows, cols;
    input.read(reinterpret_cast<char *>(&rows), sizeof(Index));
    input.read(reinterpret_cast<char *>(&cols), sizeof(Index));
    output.resize(rows, cols);
    input.read(reinterpret_cast<char *>(output.data()), output.size() * sizeof(Scalar));
    return input;
}

template<typename Scalar, int rank>
    std::ostream& serialize(std::ostream &output, const Tensor<Scalar, rank> &input)  {
    std::array<Index, rank> dims = input.dimensions();
    output.write(reinterpret_cast<char const*>(&dims), rank * sizeof(Index));
    output.write(reinterpret_cast<char const*>(input.data()), input.size() * sizeof(Scalar));
    return output;
}

/*** Read Eigen::Tensor from input stream.
 *
 * @param input The input stream from which to read the tensor.
 * @param output The tensor object in which to store the results.
 * @return reference to the input stream.
 */
template<typename Scalar, int rank>
std::istream& deserialize(std::istream &input, Tensor<Scalar, rank> &output)  {
    std::array<Index, rank> dims;
    input.read(reinterpret_cast<char *>(&dims), rank * sizeof(Index));
    output.resize(dims);
    input.read(reinterpret_cast<char *>(output.data()), output.size() * sizeof(Scalar));
    return input;
}

}  // namespace math

////////////////////////////////////////////////////////////////////////////////
// Conversion between ARTS and Eigen
////////////////////////////////////////////////////////////////////////////////

using EigenVector = Eigen::Matrix<Numeric, 1, -1, Eigen::RowMajor>;
using EigenVectorPtr = std::shared_ptr<EigenVector>;
using EigenVectorMap = Eigen::Map<EigenVector>;
using EigenConstVectorMap = Eigen::Map<const EigenVector>;

inline EigenConstVectorMap to_eigen(const Vector &vector) {
  return EigenConstVectorMap{vector.get_c_array(), vector.nelem()};
}

inline Vector to_arts(const EigenVector &vector) {
  auto n = vector.size();
  Vector result(n);
  std::copy(vector.data(), vector.data() + n, result.get_c_array());
  return result;
}

using EigenMatrix = Eigen::Matrix<Numeric, -1, -1, Eigen::RowMajor>;
using EigenMatrixPtr = std::shared_ptr<EigenMatrix>;

inline Matrix to_arts(const EigenMatrix &matrix) {
    auto m = matrix.rows();
    auto n = matrix.cols();
    Matrix result(m, n);
    std::copy(matrix.data(), matrix.data() + m * n, result.get_c_array());
    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Tensors
////////////////////////////////////////////////////////////////////////////////

template <int rank>
using EigenTensor = Eigen::Tensor<Numeric, rank, Eigen::RowMajor>;
template <int rank>
using EigenComplexTensor = Eigen::Tensor<std::complex<Numeric>, rank, Eigen::RowMajor>;
template <int rank>
using EigenConstTensorMap = Eigen::TensorMap<const EigenTensor<rank>>;

inline EigenConstTensorMap<7> to_eigen(const Tensor7 &tensor) {
  std::array<Eigen::Index, 7> dimensions = {tensor.nlibraries(),
                                            tensor.nvitrines(),
                                            tensor.nshelves(),
                                            tensor.nbooks(),
                                            tensor.npages(),
                                            tensor.nrows(),
                                            tensor.ncols()};
  return EigenConstTensorMap<7>{tensor.get_c_array(), dimensions};
}

inline EigenConstTensorMap<6> to_eigen(const Tensor6 &tensor) {
    std::array<Eigen::Index, 6> dimensions = {tensor.nvitrines(),
                                              tensor.nshelves(),
                                              tensor.nbooks(),
                                              tensor.npages(),
                                              tensor.nrows(),
                                              tensor.ncols()};
    return EigenConstTensorMap<6>{tensor.get_c_array(), dimensions};
}

inline Tensor6 to_arts(const EigenTensor<6>& input) {
    auto dimensions = input.dimensions();
    Tensor6 result{dimensions[0], dimensions[1], dimensions[2],
            dimensions[3], dimensions[4], dimensions[5],
            };
    std::copy(input.data(), input.data() + input.size(), result.get_c_array());
    return result;
}

inline Tensor7 to_arts(const EigenTensor<7>& input) {
    auto dimensions = input.dimensions();
    Tensor7 result{dimensions[0], dimensions[1], dimensions[2],
            dimensions[3], dimensions[4], dimensions[5], dimensions[6]
            };
    std::copy(input.data(), input.data() + input.size(), result.get_c_array());
    return result;
}


inline EigenConstTensorMap<7> to_arts(const Tensor7 &tensor) {
    std::array<Eigen::Index, 7> dimensions = {tensor.nlibraries(),
                                              tensor.nvitrines(),
                                              tensor.nshelves(),
                                              tensor.nbooks(),
                                              tensor.npages(),
                                              tensor.nrows(),
                                              tensor.ncols()};
    return EigenConstTensorMap<7>{tensor.get_c_array(), dimensions};
}

inline EigenConstTensorMap<5> const to_eigen(const Tensor5 &tensor) {
  std::array<Eigen::Index, 5> dimensions = {tensor.nshelves(),
                                            tensor.nbooks(),
                                            tensor.npages(),
                                            tensor.nrows(),
                                            tensor.ncols()};
  return EigenConstTensorMap<5>{tensor.get_c_array(), dimensions};
}
}  // namespace scattering
