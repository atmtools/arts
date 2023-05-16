#pragma once

#include "matpack_concepts.h"
#include "matpack_constexpr.h"

#include "enums.h"

#include <concepts>
#include <type_traits>

namespace matpack::lazy {
enum class op : char { plus, minus, times, divided_by };

template <typename T, Index n>
struct constexpr_vec_data {
  using value_type = T;
  static constexpr Index size() {return n;};
  using vec_t = matpack_constant_data<value_type, size()>;

  const vec_t& x;  // SPEEDTEST: const& WAS faster on GCC than by-value.  Clang, reverse but not so strong

  constexpr T operator[](Index i) const {return x[i];}
  constexpr operator vec_t() && {return x;}
};

//! Lvalue initialization does not owns the data
template <strict_rank_matpack_constant_data<1> Data>
constexpr_vec_data(const Data &)
    -> constexpr_vec_data<matpack_value_type<Data>,
                          std::remove_cvref_t<Data>::size()>;

//! A type compatible with the constexpr_vec_data type
template <typename T>
concept constexpr_vec_data_like =
    std::same_as<std::remove_cvref_t<T>,
                 constexpr_vec_data<typename std::remove_cvref_t<T>::value_type,
                                    std::remove_cvref_t<T>::size()>> or
    std::same_as<
        typename std::remove_cvref_t<T>::vec_t,
        matpack_constant_data<typename std::remove_cvref_t<T>::value_type,
                              std::remove_cvref_t<T>::size()>>;

template <constexpr_vec_data_like LHS, constexpr_vec_data_like RHS, op oper>
struct vop {
  using value_type = std::common_type_t<typename LHS::value_type, typename RHS::value_type>;
  static constexpr Index size() {return RHS::size();};
  static_assert(LHS::size() == size());
  using vec_t = matpack_constant_data<value_type, size()>;

  LHS lhs;
  RHS rhs;

  constexpr vop(const LHS& a, const RHS& b) : lhs(a), rhs(b) {}
  
  constexpr value_type operator[](Index i) const {
    if constexpr (oper == op::plus) return static_cast<value_type>(lhs[i]) + static_cast<value_type>(rhs[i]);
    if constexpr (oper == op::minus) return static_cast<value_type>(lhs[i]) - static_cast<value_type>(rhs[i]);
    if constexpr (oper == op::times) return static_cast<value_type>(lhs[i]) * static_cast<value_type>(rhs[i]);
    if constexpr (oper == op::divided_by) return static_cast<value_type>(lhs[i]) / static_cast<value_type>(rhs[i]);
  }

  constexpr operator vec_t() && {
    vec_t out;
    for (Index i=0; i<size(); i++) out[i] = this->operator[](i);
    return out;
  }

  constexpr vop<vop, constexpr_vec_data<value_type, size()>, op::plus> operator+(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<vop, constexpr_vec_data<value_type, size()>, op::minus> operator-(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<vop, constexpr_vec_data<value_type, size()>, op::times> operator*(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<vop, constexpr_vec_data<value_type, size()>, op::divided_by> operator/(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}

  friend constexpr vop<constexpr_vec_data<value_type, size()>, vop, op::plus> operator+(const vec_t& l, const vop& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, vop, op::minus> operator-(const vec_t& l, const vop& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, vop, op::times> operator*(const vec_t& l, const vop& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, vop, op::divided_by> operator/(const vec_t& l, const vop& r) {return {constexpr_vec_data{l}, r};}

  static_assert(constexpr_vec_data_like<vop>);
};

template <constexpr_vec_data_like LHS, constexpr_vec_data_like RHS> constexpr vop<LHS, RHS, op::plus> operator+(const LHS& l, const RHS& r) {return {l, r};}
template <constexpr_vec_data_like LHS, constexpr_vec_data_like RHS> constexpr vop<LHS, RHS, op::minus> operator-(const LHS& l, const RHS& r) {return {l, r};}
template <constexpr_vec_data_like LHS, constexpr_vec_data_like RHS> constexpr vop<LHS, RHS, op::times> operator*(const LHS& l, const RHS& r) {return {l, r};}
template <constexpr_vec_data_like LHS, constexpr_vec_data_like RHS> constexpr vop<LHS, RHS, op::divided_by> operator/(const LHS& l, const RHS& r) {return {l, r};}

template <constexpr_vec_data_like Vec>
struct vscl {
  using value_type = typename Vec::value_type;
  static constexpr Index size() {return Vec::size();};
  using vec_t = matpack_constant_data<value_type, size()>;

  value_type scl;
  Vec vec;

  constexpr vscl(const Vec& b) : scl(1.0), vec(b) {}
  constexpr vscl(value_type a, const Vec& b) : scl(a), vec(b) {}
  
  constexpr value_type operator[](Index i) const {
    return scl * vec[i];
  }

  constexpr operator vec_t() && {
    vec_t out;
    for (Index i=0; i<size(); i++) out[i] = this->operator[](i);
    return out;
  }

  constexpr vop<vscl, constexpr_vec_data<value_type, size()>, op::plus> operator+(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<vscl, constexpr_vec_data<value_type, size()>, op::minus> operator-(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<vscl, constexpr_vec_data<value_type, size()>, op::times> operator*(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<vscl, constexpr_vec_data<value_type, size()>, op::divided_by> operator/(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}

  friend constexpr vop<constexpr_vec_data<value_type, size()>, vscl, op::plus> operator+(const vec_t& l, const vscl& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, vscl, op::minus> operator-(const vec_t& l, const vscl& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, vscl, op::times> operator*(const vec_t& l, const vscl& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, vscl, op::divided_by> operator/(const vec_t& l, const vscl& r) {return {constexpr_vec_data{l}, r};}

  static_assert(constexpr_vec_data_like<vscl>);
};

template <constexpr_vec_data_like Vec> constexpr vscl<Vec> operator*(std::convertible_to<typename Vec::value_type> auto a, const Vec& b) {return {static_cast<typename Vec::value_type>(a), b};}






template <typename T, Index n>
struct constexpr_smat_data {
  using value_type = T;
  static constexpr Index size() {return n;};
  using smat_t = matpack_constant_data<value_type, size(), size()>;

  const smat_t& x;  // SPEEDTEST: const& WAS faster on GCC than by-value.  Clang, reverse but not so strong

  constexpr T operator()(Index r, Index c) const {return x(r, c);}
  constexpr operator smat_t() && {return x;}
};

//! Lvalue initialization does not owns the data
template <strict_rank_matpack_constant_data<2> Data>
constexpr_smat_data(const Data &)
    -> constexpr_smat_data<matpack_value_type<Data>,
                          std::remove_cvref_t<Data>::shape()[0]>;

//! A type compatible with the constexpr_smat_data type
template <typename T>
concept constexpr_smat_data_like =
    std::same_as<std::remove_cvref_t<T>,
                 constexpr_smat_data<typename std::remove_cvref_t<T>::value_type,
                                    std::remove_cvref_t<T>::size()>> or
    std::same_as<
        typename std::remove_cvref_t<T>::smat_t,
        matpack_constant_data<typename std::remove_cvref_t<T>::value_type,
                              std::remove_cvref_t<T>::size(),
                              std::remove_cvref_t<T>::size()>>;

template <constexpr_smat_data_like LHS, constexpr_smat_data_like RHS, op oper>
struct smop {
  static_assert(oper not_eq op::divided_by);

  using value_type = std::common_type_t<typename LHS::value_type, typename RHS::value_type>;
  static constexpr Index size() {return RHS::size();};
  static_assert(LHS::size() == size());
  using smat_t = matpack_constant_data<value_type, size(), size()>;

  LHS lhs;
  RHS rhs;

  constexpr smop(const LHS& a, const RHS& b) : lhs(a), rhs(b) {}
  
  constexpr value_type operator()(Index r, Index c) const {
    if constexpr (oper == op::plus) return static_cast<value_type>(lhs(r, c)) + static_cast<value_type>(rhs(r, c));
    if constexpr (oper == op::minus) return static_cast<value_type>(lhs(r, c)) - static_cast<value_type>(rhs(r, c));
    if constexpr (oper == op::times) {
      value_type out{0};
      for (Index i=0; i<size(); i++) out += static_cast<value_type>(lhs(r, i)) * static_cast<value_type>(rhs(i, c));
      return out;
    }
  }

  constexpr operator smat_t() && {
    smat_t out;
    if constexpr (oper == op::times)  // mini-optimization
      for (Index i = 0; i < size(); i++)
        for (Index k = 0; k < size(); k++)
          for (Index j = 0; j < size(); j++)
            out(i, j) += static_cast<value_type>(lhs(i, k)) *
                         static_cast<value_type>(rhs(k, j));
    else
      for (Index i = 0; i < size(); i++)
        for (Index j = 0; j < size(); j++)
          out(i, j) = this->operator()(i, j);
    return out;
  }

  constexpr smop<smop, constexpr_smat_data<value_type, size()>, op::plus> operator+(const smat_t& r) const {return {*this, constexpr_smat_data{r}};}
  constexpr smop<smop, constexpr_smat_data<value_type, size()>, op::minus> operator-(const smat_t& r) const {return {*this, constexpr_smat_data{r}};}
  constexpr smop<smop, constexpr_smat_data<value_type, size()>, op::times> operator*(const smat_t& r) const {return {*this, constexpr_smat_data{r}};}

  friend constexpr smop<constexpr_smat_data<value_type, size()>, smop, op::plus> operator+(const smat_t& l, const smop& r) {return {constexpr_smat_data{l}, r};}
  friend constexpr smop<constexpr_smat_data<value_type, size()>, smop, op::minus> operator-(const smat_t& l, const smop& r) {return {constexpr_smat_data{l}, r};}
  friend constexpr smop<constexpr_smat_data<value_type, size()>, smop, op::times> operator*(const smat_t& l, const smop& r) {return {constexpr_smat_data{l}, r};}

  static_assert(constexpr_smat_data_like<smop>);
};

template <constexpr_smat_data_like LHS, constexpr_smat_data_like RHS> constexpr smop<LHS, RHS, op::plus> operator+(const LHS& l, const RHS& r) {return {l, r};}
template <constexpr_smat_data_like LHS, constexpr_smat_data_like RHS> constexpr smop<LHS, RHS, op::minus> operator-(const LHS& l, const RHS& r) {return {l, r};}
template <constexpr_smat_data_like LHS, constexpr_smat_data_like RHS> constexpr smop<LHS, RHS, op::times> operator*(const LHS& l, const RHS& r) {return {l, r};}


template <constexpr_smat_data_like Mat>
struct smscl {
  using value_type = typename Mat::value_type;
  static constexpr Index size() {return Mat::size();};
  using smat_t = matpack_constant_data<value_type, size(), size()>;

  value_type scl;
  Mat mat;

  constexpr smscl(const Mat& b) : scl(1.0), mat(b) {}
  constexpr smscl(value_type a, const Mat& b) : scl(a), mat(b) {}
  
  constexpr value_type operator()(Index i, Index j) const {
    return scl * mat(i, j);
  }

  constexpr operator smat_t() && {
    smat_t out;
    for (Index i=0; i<size(); i++) for (Index j=0; j<size(); j++) out(i, j) = this->operator()(i, j);
    return out;
  }

  constexpr smop<smscl, constexpr_smat_data<value_type, size()>, op::plus> operator+(const smat_t& r) const {return {*this, constexpr_smat_data{r}};}
  constexpr smop<smscl, constexpr_smat_data<value_type, size()>, op::minus> operator-(const smat_t& r) const {return {*this, constexpr_smat_data{r}};}
  constexpr smop<smscl, constexpr_smat_data<value_type, size()>, op::times> operator*(const smat_t& r) const {return {*this, constexpr_smat_data{r}};}

  friend constexpr smop<constexpr_smat_data<value_type, size()>, smscl, op::plus> operator+(const smat_t& l, const smscl& r) {return {constexpr_smat_data{l}, r};}
  friend constexpr smop<constexpr_smat_data<value_type, size()>, smscl, op::minus> operator-(const smat_t& l, const smscl& r) {return {constexpr_smat_data{l}, r};}
  friend constexpr smop<constexpr_smat_data<value_type, size()>, smscl, op::times> operator*(const smat_t& l, const smscl& r) {return {constexpr_smat_data{l}, r};}

  static_assert(constexpr_smat_data_like<smscl>);
};

template <constexpr_smat_data_like Mat> constexpr vscl<Mat> operator*(std::convertible_to<typename Mat::value_type> auto a, const Mat& b) {return {static_cast<typename Mat::value_type>(a), b};}



template <constexpr_smat_data_like Mat, constexpr_vec_data_like Vec>
struct smvmul {
  using value_type = std::common_type_t<typename Mat::value_type, typename Vec::value_type>;
  static constexpr Index size() {return Vec::size();};
  static_assert(Mat::size() == size());

using vec_t = matpack_constant_data<value_type, size()>;

  Mat A;
  Vec x;

  constexpr smvmul(const Mat& a, const Vec& b) : A(a), x(b) {}

  constexpr value_type operator[](Index i) const {
    value_type out{};
    for (Index j = 0; j < size(); j++) out += A(i, j) * x[j];
    return out;
  }

constexpr operator vec_t() const {
  vec_t out;
  for (Index i=0; i<size(); i++) out[i] = this->operator[](i);
  return out;
}

  constexpr vop<smvmul, constexpr_vec_data<value_type, size()>, op::plus> operator+(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<smvmul, constexpr_vec_data<value_type, size()>, op::minus> operator-(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<smvmul, constexpr_vec_data<value_type, size()>, op::times> operator*(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}
  constexpr vop<smvmul, constexpr_vec_data<value_type, size()>, op::divided_by> operator/(const vec_t& r) const {return {*this, constexpr_vec_data{r}};}

  friend constexpr vop<constexpr_vec_data<value_type, size()>, smvmul, op::plus> operator+(const vec_t& l, const smvmul& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, smvmul, op::minus> operator-(const vec_t& l, const smvmul& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, smvmul, op::times> operator*(const vec_t& l, const smvmul& r) {return {constexpr_vec_data{l}, r};}
  friend constexpr vop<constexpr_vec_data<value_type, size()>, smvmul, op::divided_by> operator/(const vec_t& l, const smvmul& r) {return {constexpr_vec_data{l}, r};}

  static_assert(constexpr_vec_data_like<smvmul>);
};
} // namespace matpack::lazy