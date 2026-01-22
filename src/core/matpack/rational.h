/**
 * @file   rational.h
 * @author Richard Larsson
 * @date   2012-10-31
 * 
 * @brief  Contains the rational class definition
 **/

#ifndef rational_h
#define rational_h

#include <array.h>
#include <bifstream.h>
#include <bofstream.h>
#include <configtypes.h>
#include <format_tags.h>
#include <mystring.h>
#include <nonstd.h>
#include <xml.h>

#include <cmath>
#include <compare>
#include <concepts>
#include <format>
#include <iosfwd>
#include <numeric>
#include <type_traits>
#include <utility>

#include "xml_io_stream.h"

struct Rational;  // Forward declaration

template <typename T>
concept RationalFriend =
    not std::is_same_v<std::remove_cvref_t<T>, Rational> and
    std::constructible_from<Rational, T>;

/** Implements rational numbers to work with other ARTS types */
struct Rational {
  Index numer{};
  Index denom{1};

  static constexpr Rational reduce_by_gcd(Rational a) noexcept {
    const Index div = std::gcd(a.numer, a.denom);
    if (div) return {a.numer / div, a.denom / div};
    return a;
  }

  /** Initialization call for integer values
   * 
   * @param[in] nom Numerator
   */
  explicit constexpr Rational(const Index n) noexcept : numer(n) {}

  /** Initialization call
   * 
   * @param[in] nom Numerator
   * @param[in] denom Denominator
   */
  constexpr Rational(const Index n, const Index d) noexcept
      : numer(d ? n / std::gcd(n, d) : 0), denom(d ? d / std::gcd(n, d) : 0) {}

  /** Initialization call
   * 
   * Sets the rational from the string. Formats accepted are
   * 
   * Numeric:  1.234567890
   * Fraction: 12345/67890
   * Index:    1234567890
   * 
   * Note that overflow is possible and we do not care to capture it
   * 
   * @param[in] s String of the value
   */
  explicit Rational(const std::string_view s);

  constexpr Rational()                               = default;
  constexpr Rational(Rational&&) noexcept            = default;
  constexpr Rational(const Rational&)                = default;
  constexpr Rational& operator=(Rational&&) noexcept = default;
  constexpr Rational& operator=(const Rational&)     = default;

  template <RationalFriend T>
  constexpr Rational& operator=(T&& v) noexcept {
    return *this = Rational(std::forward<decltype(v)>(v));
  }

  [[nodiscard]] constexpr Index toIndex(Index n = 1) const noexcept {
    return (n * numer) / denom;
  }

  [[nodiscard]] constexpr Numeric toNumeric() const noexcept {
    return Numeric(numer) / Numeric(denom);
  }

  [[nodiscard]] constexpr int toInt(Index n = 1) const noexcept {
    return static_cast<int>(toIndex(n));
  }

  constexpr Rational& operator+=(const Rational& a) noexcept {
    numer  = numer * a.denom + a.numer * denom;
    denom *= a.denom;
    *this  = reduce_by_gcd(*this);
    return *this;
  }

  constexpr Rational& operator-=(const Rational& a) noexcept {
    numer  = numer * a.denom - a.numer * denom;
    denom *= a.denom;
    *this  = reduce_by_gcd(*this);
    return *this;
  }

  constexpr Rational& operator/=(const Rational& a) noexcept {
    numer *= a.denom;
    denom *= a.numer;
    *this  = reduce_by_gcd(*this);
    return *this;
  }

  constexpr Rational& operator*=(const Rational& a) noexcept {
    numer *= a.numer;
    denom *= a.denom;
    *this  = reduce_by_gcd(*this);
    return *this;
  }

  constexpr Rational operator+(const Rational& a) const noexcept {
    Rational r{*this};
    r += a;
    return r;
  }

  constexpr Rational operator-(const Rational& a) const noexcept {
    Rational r{*this};
    r -= a;
    return r;
  }

  constexpr Rational operator/(const Rational& a) const noexcept {
    Rational r{*this};
    r /= a;
    return r;
  }

  constexpr Rational operator*(const Rational& a) const noexcept {
    Rational r{*this};
    r *= a;
    return r;
  }

  template <RationalFriend T>
  constexpr Rational& operator+=(T&& a) noexcept {
    return *this += Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational& operator-=(T&& a) noexcept {
    return *this -= Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational& operator/=(T&& a) noexcept {
    return *this /= Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational& operator*=(T&& a) noexcept {
    return *this *= Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational operator+(T&& a) const noexcept {
    return *this + Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational operator-(T&& a) const noexcept {
    return *this - Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational operator/(T&& a) const noexcept {
    return *this / Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  constexpr Rational operator*(T&& a) const noexcept {
    return *this * Rational(std::forward<T>(a));
  }

  template <RationalFriend T>
  friend constexpr Rational operator+(T&& a, const Rational& b) noexcept {
    return Rational(std::forward<T>(a)) + b;
  }

  template <RationalFriend T>
  friend constexpr Rational operator-(T&& a, const Rational& b) noexcept {
    return Rational(std::forward<T>(a)) - b;
  }

  template <RationalFriend T>
  friend constexpr Rational operator/(T&& a, const Rational& b) noexcept {
    return Rational(std::forward<T>(a)) / b;
  }

  template <RationalFriend T>
  friend constexpr Rational operator*(T&& a, const Rational& b) noexcept {
    return Rational(std::forward<T>(a)) * b;
  }

  constexpr Rational operator++(int) const noexcept {
    return {numer + denom, denom};
  }

  constexpr Rational operator--(int) const noexcept {
    return {numer - denom, denom};
  }

  constexpr Rational& operator++() noexcept {
    numer += denom;
    return *this;
  }

  constexpr Rational& operator--() noexcept {
    numer -= denom;
    return *this;
  }

  constexpr Rational operator-() const noexcept { return {-numer, denom}; }
  constexpr Rational operator+() const noexcept { return *this; }

  friend constexpr Rational operator%(const Rational a,
                                      const Rational b) noexcept {
    return (a.denom == b.denom)
               ? Rational(a.numer % b.numer, a.denom)
               : Rational((a.numer * b.denom) % (a.denom * b.numer),
                          a.denom * b.denom);
  }

  template <RationalFriend T>
  friend constexpr Rational operator%(const Rational a, T&& b) noexcept {
    return a % Rational(std::forward<T>(b));
  }

  template <RationalFriend T>
  friend constexpr Rational operator%(T&& a, const Rational b) noexcept {
    return Rational(std::forward<T>(a)) % b;
  }

  constexpr auto operator<=>(const Rational& b) const noexcept = default;

  constexpr bool operator==(const Rational& b) const noexcept {
    return numer == b.numer and denom == b.denom;
  }

  template <RationalFriend T>
  constexpr auto operator<=>(T&& b) const noexcept {
    return *this <=> Rational(std::forward<T>(b));
  }

  template <RationalFriend T>
  constexpr bool operator==(T&& b) const noexcept {
    return *this == Rational(std::forward<T>(b));
  }

  template <RationalFriend T>
  friend constexpr auto operator<=>(T&& b, const Rational& a) noexcept {
    return a <=> Rational{std::forward<T>(b)};
  }

  template <RationalFriend T>
  friend constexpr bool operator==(T&& b, const Rational& a) noexcept {
    return Rational{std::forward<T>(b)} == a;
  }

  explicit constexpr operator bool() const noexcept { return numer != 0; }
  explicit constexpr operator Numeric() const noexcept { return toNumeric(); }
  explicit constexpr operator Index() const noexcept { return toIndex(); }
  explicit constexpr operator int() const noexcept { return toInt(); }

  bifstream& read(bifstream& bif);
  bofstream& write(bofstream& bof) const;

  friend std::ostream& operator<<(std::ostream& os, const Rational& a);
  friend std::istream& operator>>(std::istream& is, Rational& a);
};  // Rational

/** Absolute
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational Absolute value of the Rational
 */
constexpr Rational abs(const Rational a) noexcept { return a < 0 ? -a : a; }

/** Maximum
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational Largest of a and b
 */
constexpr const Rational& maxr(const Rational& a, const Rational& b) noexcept {
  return a < b ? b : a;
}

/** Minimum
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational Smallest of a and b
 */
constexpr const Rational& minr(const Rational& a, const Rational& b) noexcept {
  return a < b ? a : b;
}

using ArrayOfRational = Array<Rational>;

/** Returns true if even integer
 * 
 * @param[in]  r Any rational
 * @return  true if r is even, otherwise false
 */
constexpr bool iseven(const Rational r) noexcept {
  return (r.numer % (2 * r.denom) == 0);
}

/** Multiplication with numeric */
constexpr Numeric operator*(Rational y, Numeric x) noexcept {
  return static_cast<Numeric>(y) * x;
}
constexpr Numeric operator*(Numeric x, Rational y) noexcept {
  return x * static_cast<Numeric>(y);
}

/** Division with numeric */
constexpr Numeric operator/(Rational y, Numeric x) noexcept {
  return static_cast<Numeric>(y) / x;
}
constexpr Numeric operator/(Numeric x, Rational y) noexcept {
  return x / static_cast<Numeric>(y);
}

/** Addition with numeric */
constexpr Numeric operator+(Rational y, Numeric x) noexcept {
  return static_cast<Numeric>(y) + x;
}
constexpr Numeric operator+(Numeric x, Rational y) noexcept {
  return x + static_cast<Numeric>(y);
}

/** Subtraction with numeric */
constexpr Numeric operator-(Rational y, Numeric x) noexcept {
  return static_cast<Numeric>(y) - x;
}
constexpr Numeric operator-(Numeric x, Rational y) noexcept {
  return x - static_cast<Numeric>(y);
}

inline Numeric sqrtr(const Rational r) {
  return std::sqrt(static_cast<Numeric>(r));
}

inline Numeric powr(const Rational base, Numeric exp) {
  return nonstd::pow(static_cast<Numeric>(base), exp);
}

inline Numeric powr(Numeric base, const Rational exp) {
  return nonstd::pow(base, static_cast<Numeric>(exp));
}

inline Numeric powr(const Rational base, const Rational exp) {
  return powr(base, static_cast<Numeric>(exp));
}

template <>
struct std::formatter<Rational> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Rational& v, FmtContext& ctx) const {
    if (v.denom == 1) return tags.format(ctx, v.numer);
    return tags.format(ctx, v.numer, '/', v.denom);
  }
};

template <>
struct xml_io_stream<Rational> {
  static constexpr std::string_view type_name = "Rational"sv;

  static void write(std::ostream& os,
                    const Rational& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, Rational& x, bifstream* pbifs = nullptr);
};

#endif  // rational_h
