#pragma once

#include "matpack_concepts.h"

#include <complex>
#include <type_traits>

template <typename T> struct ComplexLayout {
  T real, imag;
};

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> &real_val(T &c) noexcept {
  return reinterpret_cast<ComplexLayout<matpack::complex_subtype<T>>(&)>(c)
      .real;
}

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> &imag_val(T &c) noexcept {
  return reinterpret_cast<ComplexLayout<matpack::complex_subtype<T>>(&)>(c)
      .imag;
}

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> real_val(const T &c) noexcept {
  return reinterpret_cast<const ComplexLayout<matpack::complex_subtype<T>>(&)>(
             c)
      .real;
}

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> imag_val(const T &c) noexcept {
  return reinterpret_cast<const ComplexLayout<matpack::complex_subtype<T>>(&)>(
             c)
      .imag;
}

template <matpack::complex_type T>
constexpr T operator+(matpack::complex_subtype<T> x, T c) noexcept {
  return x + c;
}

template <matpack::complex_type T>
constexpr T operator-(matpack::complex_subtype<T> x, T c) noexcept {
  return x - c;
}

template <matpack::complex_type T>
constexpr T operator*(matpack::complex_subtype<T> x, T c) noexcept {
  return x * c;
}

template <matpack::complex_type T>
constexpr T operator/(matpack::complex_subtype<T> x, T c) noexcept {
  return x / c;
}
template <matpack::complex_type T>
constexpr T operator+(T x, matpack::complex_subtype<T> c) noexcept {
  return x + c;
}

template <matpack::complex_type T>
constexpr T operator-(T x, matpack::complex_subtype<T> c) noexcept {
  return x - c;
}

template <matpack::complex_type T>
constexpr T operator*(T x, matpack::complex_subtype<T> c) noexcept {
  return x * c;
}
template <matpack::complex_type T>
constexpr T operator/(T x, matpack::complex_subtype<T> c) noexcept {
  return x / c;
}

template <matpack::complex_type T> constexpr T operator+(T x, T c) {
  return x + c;
}

template <matpack::complex_type T> constexpr T operator-(T x, T c) {
  return x - c;
}

template <matpack::complex_type T> constexpr T operator*(T x, T c) {
  return x * c;
}

template <matpack::complex_type T> constexpr T operator/(T x, T c) {
  return x / c;
}

template <matpack::complex_type T> constexpr T operator-(T c) { return -c; }

template <matpack::complex_type T> constexpr T operator+(T c) { return c; }

template <matpack::complex_type T, matpack::complex_type U>
constexpr auto operator+(T x, U c) noexcept
  requires(not std::same_as<T, U>)
{
  using K = std::common_type_t<matpack::complex_subtype<T>,
                               matpack::complex_subtype<U>>;
  return std::complex<K>{static_cast<K>(x.real()), static_cast<K>(x.imag())} +
         std::complex<K>{static_cast<K>(c.real()), static_cast<K>(c.imag())};
}

template <matpack::complex_type T, matpack::complex_type U>
constexpr auto operator-(T x, U c) noexcept
  requires(not std::same_as<T, U>)
{
  using K = std::common_type_t<matpack::complex_subtype<T>,
                               matpack::complex_subtype<U>>;
  return std::complex<K>{static_cast<K>(x.real()), static_cast<K>(x.imag())} -
         std::complex<K>{static_cast<K>(c.real()), static_cast<K>(c.imag())};
}

template <matpack::complex_type T, matpack::complex_type U>
constexpr auto operator*(T x, U c) noexcept
  requires(not std::same_as<T, U>)
{
  using K = std::common_type_t<matpack::complex_subtype<T>,
                               matpack::complex_subtype<U>>;
  return std::complex<K>{static_cast<K>(x.real()), static_cast<K>(x.imag())} *
         std::complex<K>{static_cast<K>(c.real()), static_cast<K>(c.imag())};
}

template <matpack::complex_type T, matpack::complex_type U>
constexpr auto operator/(T x, U c) noexcept
  requires(not std::same_as<T, U>)
{
  using K = std::common_type_t<matpack::complex_subtype<T>,
                               matpack::complex_subtype<U>>;
  return std::complex<K>{static_cast<K>(x.real()), static_cast<K>(x.imag())} /
         std::complex<K>{static_cast<K>(c.real()), static_cast<K>(c.imag())};
}

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> abs2(T c) {
  return c.real() * c.real() + c.imag() * c.imag();
}

template <matpack::complex_type T> constexpr T conj(T c) {
  return {c.real(), -c.imag()};
}

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> real(T c) {
  return c.real();
}

template <matpack::complex_type T>
constexpr matpack::complex_subtype<T> imag(T c) {
  return c.imag();
}
