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
template<typename internal_type, bool constant_type=true>
struct eigen {
  using stride_type = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
  using matrix_type = Eigen::Matrix<internal_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using constant_matrix_type = std::add_const_t<matrix_type>;
  using mutable_matrix_type = std::remove_const_t<matrix_type>;
  using map = Eigen::Map<std::conditional_t<constant_type, constant_matrix_type, mutable_matrix_type>, Eigen::Unaligned, stride_type>;
};

template<typename internal_type, bool constant_type=true>
using eigen_map = typename eigen<internal_type, constant_type>::map;

template<typename internal_type, bool constant_type=true>
using eigen_stride = typename eigen<internal_type, constant_type>::stride_type;

auto MapToEigenCol(matpack::vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;

  using stride_type = eigen_stride<internal_type, constant_type>;
  using matrix_map = eigen_map<internal_type, constant_type>;
  
  return matrix_map(x.get_c_array() + x.selem(), 1, x.nelem(), stride_type(1, x.delem()));
}

auto MapToEigenRow(matpack::vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;

  using stride_type = eigen_stride<internal_type, constant_type>;
  using matrix_map = eigen_map<internal_type, constant_type>;
  
  return matrix_map(x.get_c_array() + x.selem(), x.nelem(), 1, stride_type(x.delem(), 1));
}

auto MapToEigen(matpack::matrix auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;

  using stride_type = eigen_stride<internal_type, constant_type>;
  using matrix_map = eigen_map<internal_type, constant_type>;
  
  return matrix_map(x.get_c_array() + x.selem(), x.nrows(), x.ncols(), stride_type(x.drows(), x.dcols()));
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

  using stride_type = eigen_stride<internal_type>;
  using matrix_map = eigen_map<internal_type>;

  return matrix_map(x.data(), 1, x.size(), stride_type(1, 1));
}

auto ConstantEigenRowVector(standard_vector auto&& x) {
  using internal_type = std::remove_cvref_t<std::remove_pointer_t<decltype(x.data())>>;

  using stride_type = eigen_stride<internal_type>;
  using matrix_map = eigen_map<internal_type>;

  return matrix_map(x.data(), x.size(), 1, stride_type(1, 1)); 
}
}  // namespace matpack
