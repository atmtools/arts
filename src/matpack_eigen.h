#pragma once

#include <type_traits>

#include "matpackI.h"
#include "matpack_complex.h"
#include "matpack_concepts.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
#if !defined(__clang__)
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#include <Eigen/Dense>

#pragma GCC diagnostic pop

namespace matpack {

// Eigen library interactions numeric:
using StrideType = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
using MatrixType = Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixViewMap = Eigen::Map<MatrixType, 0, StrideType>;
using ConstMatrixViewMap = Eigen::Map<const MatrixType, 0, StrideType>;

// Eigen library interactions complex:
using ComplexMatrixType = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ComplexMatrixViewMap = Eigen::Map<ComplexMatrixType, 0, StrideType>;
using ConstComplexMatrixViewMap = Eigen::Map<const ComplexMatrixType, 0, StrideType>;

auto MapToEigenCol(matpack::vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;
  constexpr bool complex_type = matpack::is_complex<internal_type>::value;

  auto* data = x.get_c_array() + x.selem();
  auto stride = StrideType(1, x.delem());
  const Index row_size = 1;
  const Index col_size = x.nelem();

  if constexpr (complex_type and constant_type) {
    return ConstComplexMatrixViewMap(data, row_size, col_size, stride);
  } else if constexpr (not complex_type and constant_type) {
    return ConstMatrixViewMap(data, row_size, col_size, stride);
  } else if constexpr (complex_type and not constant_type) {
    return ComplexMatrixViewMap(data, row_size, col_size, stride);
  } else {
    return MatrixViewMap(data, row_size, col_size, stride);
  }
}

auto MapToEigenRow(matpack::vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;
  constexpr bool complex_type = matpack::is_complex<internal_type>::value;

  auto* data = x.get_c_array() + x.selem();
  auto stride = StrideType(x.delem(), 1);
  const Index row_size = x.nelem();
  const Index col_size = 1;

  if constexpr (complex_type and constant_type) {
    return ConstComplexMatrixViewMap(data, row_size, col_size, stride);
  } else if constexpr (not complex_type and constant_type) {
    return ConstMatrixViewMap(data, row_size, col_size, stride);
  } else if constexpr (complex_type and not constant_type) {
    return ComplexMatrixViewMap(data, row_size, col_size, stride);
  } else {
    return MatrixViewMap(data, row_size, col_size, stride);
  }
}

auto MapToEigen(matpack::matrix auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;
  constexpr bool complex_type = matpack::is_complex<internal_type>::value;

  auto* data = x.get_c_array() + x.selem();
  auto stride = StrideType(x.drows(), x.dcols());
  const auto row_size = x.nrows();
  const auto col_size = x.ncols();

  if constexpr (complex_type and constant_type) {
    return ConstComplexMatrixViewMap(data, row_size, col_size, stride);
  } else if constexpr (not complex_type and constant_type) {
    return ConstMatrixViewMap(data, row_size, col_size, stride);
  } else if constexpr (complex_type and not constant_type) {
    return ComplexMatrixViewMap(data, row_size, col_size, stride);
  } else {
    return MatrixViewMap(data, row_size, col_size, stride);
  }
}

auto MapToEigen(matpack::vector auto&& x) {return MapToEigenRow(x);}

template <typename T>
concept standard_vector = requires(T a) {
  a.data();
  { a.size() } -> std::integral;
  { a[0] } -> complex_or_real;
};

auto ConstantEigenColumnVector(standard_vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.data())>>;
  constexpr bool complex_type = matpack::is_complex<internal_type>::value;

  auto* data = x.data();
  auto stride = StrideType(1, 1);
  const Index row_size = 1;
  const Index col_size = x.size();

  if constexpr (complex_type) {
    return ConstComplexMatrixViewMap(data, row_size, col_size, stride);
  } else  {
    return ConstMatrixViewMap(data, row_size, col_size, stride);
  } 
}

auto ConstantEigenRowVector(standard_vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.data())>>;
  constexpr bool complex_type = matpack::is_complex<internal_type>::value;

  auto* data = x.data();
  auto stride = StrideType(1, 1);
  const Index col_size = 1;
  const Index row_size = x.size();

  if constexpr (complex_type) {
    return ConstComplexMatrixViewMap(data, row_size, col_size, stride);
  } else  {
    return ConstMatrixViewMap(data, row_size, col_size, stride);
  } 
}
}  // namespace matpack
