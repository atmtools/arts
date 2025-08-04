#pragma once

#include <configtypes.h>
#include <format_tags.h>

#include <complex>

using Complex = std::complex<Numeric>;

namespace matpack {
template <typename T>
concept arithmetic = std::is_arithmetic_v<std::remove_cvref_t<T>>;

template <typename T>
concept complex_type =
    requires(T a) {
      { a.real() } -> arithmetic;
      { a.imag() } -> arithmetic;
    } and sizeof(std::remove_cvref_t<decltype(T{}.real())>) * 2 == sizeof(T) and
    sizeof(std::remove_cvref_t<decltype(T{}.imag())>) * 2 == sizeof(T);

template <complex_type T>
struct complex_subtype {
  using type = std::remove_cvref_t<typename T::value_type>;
};

template <complex_type T>
using complex_subtype_t = typename complex_subtype<T>::type;
}  // namespace matpack

template <typename T>
struct std::formatter<std::complex<T>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const std::complex<T>& v, FmtContext& ctx) const {
    if (tags.io) {
      tags.format(ctx, v.real(), ' ', v.imag());
    } else if (tags.comma) {
      tags.format(ctx, '(', v.real(), ' ', v.imag(), "j)"sv);
    } else {
      if (v.imag() < 0) {
        tags.format(ctx, v.real(), v.imag(), 'j');
      } else {
        tags.format(ctx, v.real(), '+', v.imag(), 'j');
      }
    }

    return ctx.out();
  }
};

template <typename T>
struct ComplexLayout {
  T real, imag;
};

template <typename T>
constexpr T& real_val(std::complex<T>& c) {
  return reinterpret_cast<ComplexLayout<T>(&)>(c).real;
}

template <typename T>
constexpr const T& real_val(const std::complex<T>& c) {
  return reinterpret_cast<const ComplexLayout<T>(&)>(c).real;
}

template <typename T>
constexpr T& imag_val(std::complex<T>& c) {
  return reinterpret_cast<ComplexLayout<T>(&)>(c).imag;
}

template <typename T>
constexpr const T& imag_val(const std::complex<T>& c) {
  return reinterpret_cast<const ComplexLayout<T>(&)>(c).imag;
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator+(std::complex<T> x, U c) {
  return x + static_cast<T>(c);
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator-(std::complex<T> x, U c) {
  return x - static_cast<T>(c);
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator*(std::complex<T> x, U c) {
  return x * static_cast<T>(c);
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator/(std::complex<T> x, U c) {
  return x / static_cast<T>(c);
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator+(U c, std::complex<T> x) {
  return static_cast<T>(c) + x;
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator-(U c, std::complex<T> x) {
  return static_cast<T>(c) - x;
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator*(U c, std::complex<T> x) {
  return static_cast<T>(c) * x;
}

template <typename T, std::convertible_to<T> U>
  requires(not std::same_as<T, U>)
constexpr std::complex<T> operator/(U c, std::complex<T> x) {
  return static_cast<T>(c) / x;
}

/** Compute |c|^2
 *  
 * @param c A complex type, a + ib
 * @return a^2 + b^2
 */
template <matpack::complex_type T>
constexpr matpack::complex_subtype_t<T> abs2(T c) {
  return c.real() * c.real() + c.imag() * c.imag();
}

/** Returns the conjugate
 * 
 * @param c A complex type, a + ib
 * @return Its conjugate, a - ib
 */
template <matpack::complex_type T>
constexpr T conj(T c) {
  return {c.real(), -c.imag()};
}
