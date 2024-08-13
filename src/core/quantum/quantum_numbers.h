#pragma once

#include <enums.h>

#include <algorithm>
#include <compare>
#include <cstddef>
#include <istream>
#include <limits>
#include <ostream>
#include <string_view>
#include <utility>
#include <vector>

#include "array.h"
#include "debug.h"
#include "isotopologues.h"
#include "nonstd.h"
#include "rational.h"

constexpr Index quantum_number_error_value = -999'999'999;

namespace Quantum::Number {
//! Three tags, S: str, I: index, H: half-index
enum class QuantumNumberValueType : char { S, I, H };
std::ostream& operator<<(std::ostream&, QuantumNumberValueType);

//! Holds string values but can only hold sizeof(Index) long values
struct StringValue {
  static constexpr std::size_t N = 16;
  std::array<char, N> x{'\0'};  // First is \0 to be "empty"

  //! Returns the value in such a way that no \0 remains in the view
  [[nodiscard]] constexpr std::string_view val() const noexcept {
    // Look for the first \0 or return everything
    for (std::size_t i = 0; i < N; i++)
      if (x[i] == '\0') return {x.data(), i};
    return {x.data(), N};
  }

  //! Default initializer to a zero-first solution
  constexpr StringValue() = default;

  //! Set to expected value from a view of something at most N-char long
  explicit constexpr StringValue(std::string_view s) {
    const std::size_t n = s.size();

    /** This means that we will need to either redesign quantum numbers
      or make an exception for the type, e.g., is \tilde{X} perhaps OK as ~X?
      We have no opinion at this moment */
    ARTS_USER_ERROR_IF(n > N,
                       "The value \"",
                       s,
                       "\" is too long.  Can be only ",
                       N,
                       " chars but is ",
                       n)

    // Fill with correct values or zero characters
    std::size_t i = 0;
    for (; i < n; i++) x[i] = s[i];
    if (i < N) x[i] = '\0';
  }

  constexpr std::strong_ordering operator<=>(const StringValue& sv) const {
    for (std::size_t i = 0; i < N; i++) {
      if (x[i] < sv.x[i]) return std::strong_ordering::less;
      if (sv.x[i] < x[i]) return std::strong_ordering::greater;
    }
    return std::strong_ordering::equal;
  }
};

//! Holds integer values
struct IntegerValue {
  Index x{std::numeric_limits<Index>::lowest()};

  //! Returns the value as a rational
  [[nodiscard]] constexpr Rational val() const noexcept { return x; }

  constexpr std::strong_ordering operator<=>(const IntegerValue& i) const {
    return x <=> i.x;
  }
};

//! Holds half integer values, but only its denominator
struct HalfIntegerValue {
  Index x{std::numeric_limits<Index>::lowest()};

  //! Returns the value as a Rational, keeping that this is a half-integer
  [[nodiscard]] constexpr Rational val() const noexcept {
    return Rational(x, 2);
  }

  constexpr std::strong_ordering operator<=>(const HalfIntegerValue& h) const {
    return x <=> h.x;
  }
};

/** Common value type of a given quantum number type 
 * 
 * Guards against integer/half-integer comparisons failing for bad reasons
 * and allows IO to not need some information 
 *
 * Keep new entries to a 1-to-1 key-to-key match, and sort them the same as QuantumNumberType for ease of reading
 *
 * @param type A quantum number type
 * @return constexpr QuantumNumberValueType The common type
 */
constexpr QuantumNumberValueType common_value_type(
    QuantumNumberType type) noexcept {
  switch (type) {
    case QuantumNumberType::alpha:
      return QuantumNumberValueType::S;
    case QuantumNumberType::config:
      return QuantumNumberValueType::S;
    case QuantumNumberType::ElecStateLabel:
      return QuantumNumberValueType::S;
    case QuantumNumberType::F:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F1:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F10:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F11:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F12:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F2:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F3:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F4:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F5:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F6:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F7:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F8:
      return QuantumNumberValueType::H;
    case QuantumNumberType::F9:
      return QuantumNumberValueType::H;
    case QuantumNumberType::I:
      return QuantumNumberValueType::I;
    case QuantumNumberType::J:
      return QuantumNumberValueType::H;
    case QuantumNumberType::K:
      return QuantumNumberValueType::I;
    case QuantumNumberType::Ka:
      return QuantumNumberValueType::I;
    case QuantumNumberType::Kc:
      return QuantumNumberValueType::I;
    case QuantumNumberType::L:
      return QuantumNumberValueType::I;
    case QuantumNumberType::Lambda:
      return QuantumNumberValueType::I;
    case QuantumNumberType::N:
      return QuantumNumberValueType::I;
    case QuantumNumberType::Omega:
      return QuantumNumberValueType::H;
    case QuantumNumberType::S:
      return QuantumNumberValueType::H;
    case QuantumNumberType::Sigma:
      return QuantumNumberValueType::H;
    case QuantumNumberType::SpinComponentLabel:
      return QuantumNumberValueType::I;
    case QuantumNumberType::asSym:
      return QuantumNumberValueType::S;
    case QuantumNumberType::elecInv:
      return QuantumNumberValueType::S;
    case QuantumNumberType::elecRefl:
      return QuantumNumberValueType::S;
    case QuantumNumberType::elecSym:
      return QuantumNumberValueType::S;
    case QuantumNumberType::kronigParity:
      return QuantumNumberValueType::S;
    case QuantumNumberType::l:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l1:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l10:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l11:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l12:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l2:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l3:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l4:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l5:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l6:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l7:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l8:
      return QuantumNumberValueType::I;
    case QuantumNumberType::l9:
      return QuantumNumberValueType::I;
    case QuantumNumberType::n:
      return QuantumNumberValueType::S;
    case QuantumNumberType::parity:
      return QuantumNumberValueType::S;
    case QuantumNumberType::r:
      return QuantumNumberValueType::I;
    case QuantumNumberType::rotSym:
      return QuantumNumberValueType::S;
    case QuantumNumberType::rovibSym:
      return QuantumNumberValueType::S;
    case QuantumNumberType::sym:
      return QuantumNumberValueType::S;
    case QuantumNumberType::tau:
      return QuantumNumberValueType::S;
    case QuantumNumberType::term:
      return QuantumNumberValueType::S;
    case QuantumNumberType::v:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v1:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v10:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v11:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v12:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v2:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v3:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v4:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v5:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v6:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v7:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v8:
      return QuantumNumberValueType::I;
    case QuantumNumberType::v9:
      return QuantumNumberValueType::I;
    case QuantumNumberType::vibInv:
      return QuantumNumberValueType::S;
    case QuantumNumberType::vibRefl:
      return QuantumNumberValueType::S;
    case QuantumNumberType::vibSym:
      return QuantumNumberValueType::S;
  }
  std::unreachable();
}

/** Return a common type between a and b
 * 
 * If they are the same, returns the value
 * If either is H and the other is I, returns H
 * Otherwise returns FINAL as an error flag
 *
 * @param a A value type
 * @param b Another value type
 * @return constexpr QuantumNumberValueType with FINAL as error flag
 */
constexpr QuantumNumberValueType common_value_type(
    QuantumNumberValueType a, QuantumNumberValueType b) noexcept {
  // Same is same, H is I, a has both:
  if (a == b or
      (a == QuantumNumberValueType::H and b == QuantumNumberValueType::I))
    return a;

  // H is I:
  if (b == QuantumNumberValueType::H and a == QuantumNumberValueType::I)
    return QuantumNumberValueType::H;

  return static_cast<QuantumNumberValueType>(-1);
}

//! A union of the three type of values we need to consider
union ValueHolder {
  StringValue s;
  IntegerValue i;
  HalfIntegerValue h;

  constexpr ValueHolder(QuantumNumberValueType t) noexcept
      : s(StringValue{"NODEF"}) {
    switch (t) {
      case QuantumNumberValueType::H:
        h = HalfIntegerValue{0};
        break;
      case QuantumNumberValueType::I:
        i = IntegerValue{0};
        break;
      default: {
      }
    }
  }

  constexpr ValueHolder(QuantumNumberType t) noexcept
      : ValueHolder(common_value_type(t)) {}
  constexpr ValueHolder(const ValueHolder&)                = default;
  constexpr ValueHolder(ValueHolder&&) noexcept            = default;
  constexpr ValueHolder& operator=(const ValueHolder&)     = default;
  constexpr ValueHolder& operator=(ValueHolder&&) noexcept = default;
};

/** A complete description of a value, its type and value
 * 
 * Intended to be returned from IO operations
 */
struct ValueDescription {
  QuantumNumberValueType type;
  ValueHolder val;

  constexpr ValueDescription(QuantumNumberValueType t) noexcept
      : type(t), val(t) {}
  constexpr ValueDescription(const ValueDescription&)                = default;
  constexpr ValueDescription(ValueDescription&&) noexcept            = default;
  constexpr ValueDescription& operator=(const ValueDescription&)     = default;
  constexpr ValueDescription& operator=(ValueDescription&&) noexcept = default;

  //! Debug output only
  friend std::ostream& operator<<(std::ostream& os, ValueDescription x);
};

constexpr std::strong_ordering cmp(std::strong_ordering&& a,
                                   std::strong_ordering&& b) {
  return a == std::strong_ordering::equal ? b : a;
}

/** The values of two levels
 * 
 * Its ValueDescription constructor ensures that we have two valid 
 * types and converts integer input in one to half-integer in case
 * there's a 'mismatch'
 */
struct TwoLevelValueHolder {
  ValueHolder upp;
  ValueHolder low;

  //! Constructor to ensure two ValueDescription have the same type, for IO purposes
  constexpr TwoLevelValueHolder(ValueDescription u,
                                ValueDescription l,
                                QuantumNumberType t)
      : upp(u.val), low(l.val) {
    auto ct = common_value_type(t);
    if (u.type not_eq ct) {
      ARTS_USER_ERROR_IF(u.type not_eq QuantumNumberValueType::I,
                         "Cannot convert from ",
                         u.type,
                         " to ",
                         ct)
      upp.h.x = 2 * u.val.i.x;
    }

    if (l.type not_eq ct) {
      ARTS_USER_ERROR_IF(l.type not_eq QuantumNumberValueType::I,
                         "Cannot convert from ",
                         l.type,
                         " to ",
                         ct)
      low.h.x = 2 * l.val.i.x;
    }
  }

  constexpr TwoLevelValueHolder(QuantumNumberType t) noexcept
      : upp(t), low(t) {}
  constexpr TwoLevelValueHolder(const TwoLevelValueHolder&)     = default;
  constexpr TwoLevelValueHolder(TwoLevelValueHolder&&) noexcept = default;
  constexpr TwoLevelValueHolder& operator=(const TwoLevelValueHolder&) =
      default;
  constexpr TwoLevelValueHolder& operator=(TwoLevelValueHolder&&) noexcept =
      default;

  [[nodiscard]] constexpr std::strong_ordering order(
      const TwoLevelValueHolder& tv, QuantumNumberValueType t) const {
    switch (t) {
      case QuantumNumberValueType::S:
        return cmp(upp.s <=> tv.upp.s, low.s <=> tv.low.s);
      case QuantumNumberValueType::I:
        return cmp(upp.i <=> tv.upp.i, low.i <=> tv.low.i);
      case QuantumNumberValueType::H:
        return cmp(upp.h <=> tv.upp.h, low.h <=> tv.low.h);
    }
    return std::strong_ordering::equal;
  }
};

/** Takes a rational and determine which type of quantum number it is,
 * returning this information or throwing a runtime error if there's
 * an error
 *
 * @param r_ A rational
 * @return constexpr ValueDescription 
 */
[[nodiscard]] constexpr ValueDescription value_holder(Rational r_) {
  const Rational r = reduce_by_gcd(r_);

  ValueDescription x{QuantumNumberValueType::H};

  // We must now have a half-integer or not
  if (r.denom == 2) {
    x.val.h.x = r.numer;
    return x;
  }

  if (r.denom == 1) {
    x.type    = QuantumNumberValueType::I;
    x.val.i.x = r.numer;
    return x;
  }

  x.type    = QuantumNumberValueType::I;
  x.val.i.x = quantum_number_error_value;
  return x;
}

[[nodiscard]] constexpr ValueDescription value_holder(std::string_view s) {
  ValueDescription x{QuantumNumberValueType::S};
  x.val.s = StringValue(s);
  return x;
}

/** Returns a rational if possible or RATIONAL_UNDEFINED otherwise
 * 
 * Spaces are not allowed and results in RATIONAL_UNDEFINED being returns
 *
 * Also only accepts rationals that match how quantum numbers are defined,
 * so half or full integers, in case a decimal string is given.  Note also
 * that only a single decimal value can be given
 *
 * If the rational is not a full integer or a half-integer, returns RATIONAL_UNDEFINED
 *
 * @param s Some view of a string
 * @return constexpr Rational 
 */
[[nodiscard]] constexpr Rational cast_qnrat(std::string_view s) noexcept {
  // Counts for divides, decimals, and existence
  int div = 0, dot = 0, any = 0, minus = false;
  const std::size_t n = s.size();

  // Counts relevant items
  for (std::size_t i = 0; i < n; i++) {
    auto x = s[i];
    if (x == '-') {
      minus = true;
      if (i) return RATIONAL_UNDEFINED;
    } else if (x == '+') {
      if (i) return RATIONAL_UNDEFINED;
    } else if (x == '/')
      div++;  // Count divs because we can have at most one
    else if (x == '.')
      dot++;  // Count dots for the same reason
    else if (not nonstd::isdigit(x))
      return RATIONAL_UNDEFINED;  // Error!

    // There is a value!
    any++;
  }

  // Can only have one of div or dot and need some data
  if ((div + dot) > 1 or any == 0) return RATIONAL_UNDEFINED;

  // We have a rational!  Lets see which we have got

  // We have a simple rational
  if (div) {
    Index num = 0, den = 0;
    std::size_t i = 0;

    // Numerator
    for (; s[i] not_eq '/'; ++i) {
      if (s[i] == '-' or s[i] == '+') continue;
      num *= 10;
      num += s[i] - '0';
    }

    // Denominator
    i++;
    for (; i < s.size(); ++i) {
      den *= 10;
      den += s[i] - '0';
    }

    // Guard for QN style rationals
    return Rational(minus ? -num : num, den);
  }

  // We have a decimal number
  if (dot) {
    Index f = 0, d = 0;
    std::size_t i = 0;

    // The integer part
    for (; s[i] not_eq '.'; ++i) {
      if (s[i] == '-' or s[i] == '+') continue;
      f *= 10;
      f += s[i] - '0';
    }

    // The decimal part
    i++;
    for (; i < s.size(); ++i) {
      d *= 10;
      d += s[i] - '0';
    }

    if (d == 0) return minus ? -f : f;
    if (d == 5) return Rational((minus ? -1 : 1) * (2 * f + 1), 2);
    return RATIONAL_UNDEFINED;
  }

  std::size_t num = 0;
  for (auto x : s) {
    if (x == '-' or x == '+') continue;
    num *= 10;
    num += x - '0';
  }
  return minus ? -num : num;
}

/** Returns a value description for the quantum number
 * 
 * Note that several branches can throw as the input is assumed to be from
 * the user
 *
 * @param s Some view of a string
 * @return constexpr ValueDescription 
 */
[[nodiscard]] constexpr ValueDescription value_holder(std::string_view s,
                                                      QuantumNumberType t) {
  switch (common_value_type(t)) {
    case QuantumNumberValueType::I:
    case QuantumNumberValueType::H:
      return value_holder(cast_qnrat(s));
    case QuantumNumberValueType::S:
      return value_holder(s);
  }
  std::unreachable();
}

//! Struct that converts to bool automatically but allows checking both energy levels matching status
struct LevelMatch {
  bool upp{true};
  bool low{true};

  //! Convert automatically to bool so exact matches are easy, non-explicit by design
  constexpr operator bool() const noexcept { return upp and low; }
};

/** Count all space-separated items in s
 * 
 * Example: "X 1 3123 1/3 1,,,,,2 " returns 5
 *
 * @param s Any set of characters
 * @return constexpr Index The number of space-separated items in s
 */
constexpr Index count_items(std::string_view s) noexcept {
  // Checks if we are in-between items, we start true as we are inbetween items
  bool last_space = true;

  Index count = 0;
  for (auto& x : s) {
    const bool this_space = nonstd::isspace(x);

    // If we had a space and now no longer do, we are in an item
    if (last_space and not this_space) count++;

    // The current state must be remembere
    last_space = this_space;
  }
  return count;
}

/** Strips spaces at the end of x before returning it
 * 
 * @param x any string view
 * @return constexpr std::string_view stripped
 */
constexpr std::string_view rstrip(std::string_view x) {
  while (not x.empty() and nonstd::isspace(x.back())) x.remove_suffix(1);
  return x;
}

/** Strips spaces at the beginning x before returning it
 * 
 * @param x any string view
 * @return constexpr std::string_view stripped
 */
constexpr std::string_view lstrip(std::string_view x) {
  while (not x.empty() and nonstd::isspace(x.front())) x.remove_prefix(1);
  return x;
}

/** Strips spaces at the beginning and end of x before returning it
 * 
 * @param x any string view
 * @return constexpr std::string_view stripped
 */
constexpr std::string_view strip(std::string_view x) {
  return rstrip(lstrip(x));
}

/** Get a view of a number of space-separated items from the list
 * 
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=0, n=3 returns "X    1   3123"
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=2, n=3 returns "3123 1/3 1,,,,,2"
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=4, n=3 returns "1,,,,,2"
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=5, n=3 returns ""
 *
 * @tparam n Number of items
 * @param s Any set of characters
 * @param i The first item from the original list in the list of n items
 * @return constexpr std::string_view 
 */
template <std::size_t n = 1>
constexpr std::string_view items(std::string_view s, std::size_t i) noexcept {
  static_assert(n > 0, "Must want some items");
  bool last_space = true;

  std::size_t beg = 0, count = 0, end = s.size();
  if (end == 0) return s;

  for (std::size_t ind = 0; ind < end; ind++) {
    const bool this_space = nonstd::isspace(s[ind]);

    // Return when we find the end of the final item
    if (this_space and count == i + n) return {&s[beg], ind - beg};

    // If we had a space and now no longer do, we are in an item
    if (last_space and not this_space) {
      count++;

      // If that is our first item, we are good!
      if (count - 1 == i) beg = ind;
    }

    // Count up the beginning until we have found the first item
    if (count - 1 < i) beg = ind;

    // The current state must be remembere
    last_space = this_space;
  }

  // Remove spaces at the end to be sure
  while (nonstd::isspace(s[end - 1]) and end > beg) end--;
  return {&s[beg], end - beg};
}

//! A complete quantum number value with type information
struct Value {
  QuantumNumberType type;
  TwoLevelValueHolder qn;

  constexpr std::strong_ordering operator<=>(const Value& v) const {
    if (type < v.type) return std::strong_ordering::less;
    if (v.type < type) return std::strong_ordering::greater;
    return qn.order(v.qn, common_value_type(type));
  }

  constexpr Value(QuantumNumberType t = QuantumNumberType::term)
      : type(t), qn(type) {}
  Value(const Value&)                = default;
  Value(Value&&) noexcept            = default;
  Value& operator=(const Value&)     = default;
  Value& operator=(Value&&) noexcept = default;

  constexpr Value(QuantumNumberType t, Rational upp_, Rational low_)
      : Value(t) {
    Rational upp = reduce_by_gcd(upp_), low = reduce_by_gcd(low_);

    if (common_value_type(type) == QuantumNumberValueType::H) {
      ARTS_ASSERT(upp.denom <= 2 and low.denom <= 2)
      if (upp.denom not_eq 2) upp *= 2;
      if (low.denom not_eq 2) low *= 2;
      qn.upp.h.x = upp.numer;
      qn.low.h.x = low.numer;
    } else if (common_value_type(type) == QuantumNumberValueType::I) {
      ARTS_ASSERT(upp.denom == 1 and low.denom == 1)
      qn.upp.i.x = upp.numer;
      qn.low.i.x = low.numer;
    } else {
      ARTS_USER_ERROR(
          t, " is a string-type, so cannot be constructed from rationals")
    }
  }

  //! Default constructor from some string of values
  constexpr Value(std::string_view s)
      : Value(to<QuantumNumberType>(items(s, 0))) {
    ARTS_USER_ERROR_IF(count_items(s) not_eq 3,
                       "Must have ' TYPE UPPNUM LOWNUM ' but got: '",
                       s,
                       '\'')

    // Get values and ensure they are good types
    auto upv = value_holder(items(s, 1), type);
    auto lov = value_holder(items(s, 2), type);

    // Deal with errors while setting the level values
    qn = TwoLevelValueHolder(upv, lov, type);
  }

  //! Returns the upper quantum number rational if it exists or an undefined
  [[nodiscard]] constexpr Rational upp() const noexcept {
    switch (common_value_type(type)) {
      case QuantumNumberValueType::I:
        return qn.upp.i.val();
      case QuantumNumberValueType::H:
        return qn.upp.h.val();
      default: {
      }
    }
    return RATIONAL_UNDEFINED;
  }

  //! Returns the lower quantum number rational if it exists or an undefined
  [[nodiscard]] constexpr Rational low() const noexcept {
    switch (common_value_type(type)) {
      case QuantumNumberValueType::I:
        return qn.low.i.val();
      case QuantumNumberValueType::H:
        return qn.low.h.val();
      default: {
      }
    }
    return RATIONAL_UNDEFINED;
  }

  //! Returns the upper quantum number string copy
  [[nodiscard]] String str_upp() const noexcept;

  //! Returns the lower quantum number string copy
  [[nodiscard]] String str_low() const noexcept;

  //! Legacy way to swap the values between two Values
  void swap_values(Value& x);

  //! Set level value
  constexpr void set(std::string_view s, bool upp) {
    const ValueDescription v = value_holder(s, type);
    const TwoLevelValueHolder nqn(v, v, type);
    if (upp) {
      qn.upp = nqn.upp;
    } else {
      qn.low = nqn.low;
    }
  }

  /** Returns a description of whether both levels match
   * 
   * Note that the LevelMatch type should automatically transform to bool
   * so there's no need for extra work if this is your target question
   *
   * @param other Another value
   * @return constexpr LevelMatch
   */
  [[nodiscard]] constexpr LevelMatch level_match(Value other) const noexcept {
    if (type == other.type) {
      switch (common_value_type(type)) {
        case QuantumNumberValueType::I:
        case QuantumNumberValueType::H:
          return {upp() == other.upp(), low() == other.low()};
        case QuantumNumberValueType::S:
          return {qn.upp.s.x == other.qn.upp.s.x,
                  qn.low.s.x == other.qn.low.s.x};
      }
    }
    return {false, false};
  }

  //! Standard output
  friend std::ostream& operator<<(std::ostream& os, Value x);

  //! Standard input
  friend std::istream& operator>>(std::istream& is, Value& x);

  bofstream& write(bofstream& bof) const;

  bifstream& read(bifstream& bif);

  [[nodiscard]] constexpr bool good() const { return level_match(*this); }
};

//! Status of comparing two lists that are supposedly of some type
enum class CheckValue : char { Full, AinB, BinA, Miss };

//! Level-by-level version of CheckValue
struct CheckMatch {
  CheckValue upp{CheckValue::Full};
  CheckValue low{CheckValue::Full};

  //! Convert automatically to bool so exact matches are easy, non-explicit by design
  constexpr operator bool() const noexcept {
    return upp == CheckValue::Full and low == CheckValue::Full;
  }
};

//! Updates old by what a new check says it should be
constexpr CheckValue update(CheckValue val, CheckValue res) noexcept {
  if (val == CheckValue::Miss) return val;
  if (val == CheckValue::Full) return res;
  if (res == CheckValue::Miss) return res;
  if (res == CheckValue::Full) return val;
  if (val == res) return val;
  return CheckValue::Miss;
}

//! Updates old by what a new check says it should be
constexpr CheckMatch update(CheckMatch val, CheckValue res) noexcept {
  return {update(val.upp, res), update(val.low, res)};
}

//! Updates old by what a new check says it should be
constexpr CheckMatch update(CheckMatch val, CheckMatch res) noexcept {
  return {update(val.upp, res.upp), update(val.low, res.low)};
}

/** Checks if an array of types is sorted
 * 
 * @tparam N Number of types
 * @param types Array of types
 * @return true if it is sorted
 * @return false if it is not sorted
 */
template <size_t N>
constexpr bool is_sorted(
    const std::array<QuantumNumberType, N>& types) noexcept {
  for (size_t i = 1; i < N; i++)
    if (not(types[i - 1] < types[i])) return false;
  return true;
}

//! A list of many quantum numbers.  Should always remain sorted
struct ValueList {
  Array<Value> values;

  //! Internal sort function.  Should be called whenever new items are created
  void sort_by_type();

  //! Internal check function.  Remember to sort by type before calling this
  [[nodiscard]] bool has_unique_increasing_types() const;

  //! From text
  explicit ValueList(std::string_view s, bool legacy = false);

  //! From legacy text
  ValueList(std::string_view upp, std::string_view low);

  std::strong_ordering operator<=>(const ValueList& v) const;

  //! From values (resorted)
  explicit ValueList(Array<Value> values_) : values(std::move(values_)) {
    finalize();
  }

  //! Empty
  ValueList() : values(0) {}

  //! For iterators
  Array<Value>::iterator begin() { return values.begin(); }
  Array<Value>::iterator end() { return values.end(); }
  Array<Value>::iterator cbegin() { return values.begin(); }
  Array<Value>::iterator cend() { return values.end(); }
  [[nodiscard]] Array<Value>::const_iterator begin() const {
    return values.begin();
  }
  [[nodiscard]] Array<Value>::const_iterator end() const {
    return values.end();
  }
  [[nodiscard]] Array<Value>::const_iterator cbegin() const {
    return values.cbegin();
  }
  [[nodiscard]] Array<Value>::const_iterator cend() const {
    return values.cend();
  }

  //! Should always be called before this object is handed to another user
  void finalize();

  //! Return number of quantum numbers
  [[nodiscard]] Index size() const ARTS_NOEXCEPT { return values.size(); }

  //! Finds whether two ValueList describe completely different sets of quantum numbers (e.g., local vs global)
  [[nodiscard]] bool perpendicular(const ValueList& that) const ARTS_NOEXCEPT;

  //! Returns whether all the Types are part of the list, the types must be sorted
  template <typename... Types>
  [[nodiscard]] bool has(Types... ts) const ARTS_NOEXCEPT {
    static_assert(sizeof...(Types) > 0);

    ARTS_ASSERT(is_sorted(std::array{QuantumNumberType(ts)...}))

    auto ptr = cbegin();
    auto end = cend();
    for (QuantumNumberType t : {QuantumNumberType(ts)...}) {
      ptr = std::find_if(ptr, end, [t](auto& x) { return x.type == t; });
      if (ptr == end) return false;
    }
    return true;
  }

  //! Returns the value of the Type (assumes it exist)
  const Value& operator[](QuantumNumberType t) const ARTS_NOEXCEPT;

  //! Legacy manipulation operator access
  Value& operator[](Index i) { return values.at(i); }

  //! Add for manipulation
  Value& add(QuantumNumberType t);

  //! Add for manipulation
  Value& add(Value v);

  //! Sets the value if it exists or adds it otherwise
  void set(Value v);

  //! Set a value in value list
  void set(Index i, std::string_view upp, std::string_view low);

  //! Returns upper and lower matching status
  [[nodiscard]] CheckMatch check_match(const ValueList& other) const
      ARTS_NOEXCEPT;

  //! ouptut stream if all values
  friend std::ostream& operator<<(std::ostream& os, const ValueList& vl);

  //! input stream must have pre-set size
  friend std::istream& operator>>(std::istream& is, ValueList& vl);

  //! Add a type without sorting (WARNING, many things might break if you don't sort in the end)
  void add_type_wo_sort(QuantumNumberType);

  [[nodiscard]] bool good() const;

  void reserve(Size);

  template <typename... Ts>
  auto& emplace_back(Ts&&... xs) {
    return values.emplace_back(std::forward<Ts>(xs)...);
  }
};

ValueList from_hitran(std::string_view upp, std::string_view low);

//! A logical struct for local quantum numbers
struct LocalState {
  ValueList val{};

  std::strong_ordering operator<=>(const LocalState& l) const {
    return val <=> l.val;
  }
  auto operator==(const LocalState& l) const {
    return std::strong_ordering::equal == (*this <=> l);
  }
  auto operator!=(const LocalState& l) const {
    return std::strong_ordering::equal != (*this <=> l);
  }

  LocalState() = default;

  template <typename... Values>
  LocalState(Values... vals) : val(Array<Value>{Value(vals)...}) {}

  void set_unsorted_qns(const Array<QuantumNumberType>& vals);

  [[nodiscard]] String keys() const;

  [[nodiscard]] String values() const;

  [[nodiscard]] bool same_types_as(const LocalState& that) const;

  //! ouptut stream if all values
  friend std::ostream& operator<<(std::ostream& os, const LocalState& vl);

  //! input stream must have pre-set size
  friend std::istream& operator>>(std::istream& is, LocalState& vl);

  //! Test if there are bad quantum numbers (undefined ones)
  [[nodiscard]] bool good() const;
};

struct LevelTest {
  bool upp{true}, low{true};
};

//! A logical struct for global quantum numbers with species identifiers
struct GlobalState {
  static constexpr Index version = 1;  // Second version of quantum identifiers

  Index isotopologue_index{"Ar-8"_isot_index};
  ValueList val{};

  GlobalState() = default;

  explicit GlobalState(Index i, ValueList v = {})
      : isotopologue_index(i), val(std::move(v)) {}

  explicit GlobalState(const SpeciesIsotope& ir)
      : isotopologue_index(Species::find_species_index(ir)) {}

  explicit GlobalState(std::string_view s, Index v = version);

  [[nodiscard]] SpeciesIsotope Isotopologue() const noexcept;
  [[nodiscard]] SpeciesEnum Species() const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const GlobalState& gs);

  friend std::istream& operator>>(std::istream& is, GlobalState& gs);

  [[nodiscard]] GlobalState LowerLevel() const;
  [[nodiscard]] GlobalState UpperLevel() const;

  std::strong_ordering operator<=>(const GlobalState& g) const {
    if (isotopologue_index < g.isotopologue_index)
      return std::strong_ordering::less;
    if (g.isotopologue_index < isotopologue_index)
      return std::strong_ordering::greater;
    return val <=> g.val;
  }
  auto operator==(const GlobalState& g) const {
    return std::strong_ordering::equal == (*this <=> g);
  }
  auto operator!=(const GlobalState& g) const {
    return std::strong_ordering::equal != (*this <=> g);
  }

  //! Checks whether all of the LHS is part of RHS
  [[nodiscard]] bool part_of(const GlobalState& other) const;

  //! Checks whether all of the RHS is part of LHS (uses reverse part_of)
  [[nodiscard]] bool may_be(const GlobalState& other) const;

  //! Checks whether all of the LHS is part of any of the RHS
  [[nodiscard]] LevelTest part_of(const GlobalState& g,
                                  const LocalState& l) const;

  //! Test if there are bad quantum numbers (undefined ones) or if the isotopologue is not a normal target
  [[nodiscard]] bool good() const;
};

//! StateMatchType operates so that a check less than a level should be 'better', bar None
enum class StateMatchType : char { Full, Level, Isotopologue, Species, None };

//! StateMatch, where you need to check the upp and low values when manipulating energy levels
struct StateMatch {
  StateMatchType type{StateMatchType::None};
  bool upp{false}, low{false};

  constexpr StateMatch() = default;

  StateMatch(const GlobalState& target,
             const LocalState& local,
             const GlobalState& global);

  StateMatch(const GlobalState& target, const GlobalState& key);

  //! It is of the desired type if it is less than the value, bar None
  constexpr bool operator==(StateMatchType x) const noexcept {
    return x == type;
  }

  constexpr bool operator!=(StateMatchType x) const noexcept {
    return not((*this) == x);
  }
};

//! VAMDC classes of quantum number cases
enum class VAMDC : char {
  asymcs,  // Schema for specifying the quantum numbers of closed-shell asymmetric top molecules
  asymos,  // Schema for specifying the quantum numbers of open-shell asymmetric top molecules
  dcs,  // Schema for specifying the quantum numbers of closed-shell, diatomic molecules
  hunda,  // Schema for specifying the quantum numbers for Hund's case (a) diatomic molecules
  hundb,  // Schema for specifying the quantum numbers for Hund's case (b) diatomic molecules
  lpcs,  // Schema for specifying the quantum numbers of closed-shell linear polyatomic molecules
  lpos,  // Schema for specifying the quantum numbers of open-shell linear polyatomic molecules
  ltcs,  // Schema for specifying the quantum numbers of closed-shell linear triatomic molecules
  ltos,  // Schema for specifying the quantum numbers of open-shell linear triatomic molecules
  nltcs,  // Schema for specifying the quantum numbers of closed-shell non-linear triatomic molecules
  nltos,  // Schema for specifying the quantum numbers of open-shell non-linear triatomic molecules
  sphcs,  // Schema for specifying the quantum numbers of closed-shell spherical top molecules
  sphos,  // Schema for specifying the quantum numbers of closed-shell spherical top molecules
  stcs  // Schema for specifying the quantum numbers of closed-shell, symmetric top molecules
};

/** Checks if a ValueList can belong to a given VAMDC type by
 * ensuring it cannot have some quantum numbers
 * 
 * @param l A list of values
 * @param type A type of VAMDC molecular model
 * @return true If it can belong to the VAMDC type
 * @return false If it cannot belong to the VAMDC type
 */
bool vamdcCheck(const ValueList& l, VAMDC type) ARTS_NOEXCEPT;

//! A default state of global quantum numbers
[[maybe_unused]] inline constexpr std::array global_types{
    QuantumNumberType::alpha,
    QuantumNumberType::config,
    QuantumNumberType::ElecStateLabel,
    QuantumNumberType::L,
    QuantumNumberType::Lambda,
    QuantumNumberType::Omega,
    QuantumNumberType::S,
    QuantumNumberType::Sigma,
    QuantumNumberType::SpinComponentLabel,
    QuantumNumberType::asSym,
    QuantumNumberType::elecInv,
    QuantumNumberType::elecRefl,
    QuantumNumberType::elecSym,
    QuantumNumberType::kronigParity,
    QuantumNumberType::l,
    QuantumNumberType::l1,
    QuantumNumberType::l10,
    QuantumNumberType::l11,
    QuantumNumberType::l12,
    QuantumNumberType::l2,
    QuantumNumberType::l3,
    QuantumNumberType::l4,
    QuantumNumberType::l5,
    QuantumNumberType::l6,
    QuantumNumberType::l7,
    QuantumNumberType::l8,
    QuantumNumberType::l9,
    QuantumNumberType::n,
    QuantumNumberType::parity,
    QuantumNumberType::r,
    QuantumNumberType::rotSym,
    QuantumNumberType::rovibSym,
    QuantumNumberType::sym,
    QuantumNumberType::tau,
    QuantumNumberType::term,
    QuantumNumberType::v,
    QuantumNumberType::v1,
    QuantumNumberType::v10,
    QuantumNumberType::v11,
    QuantumNumberType::v12,
    QuantumNumberType::v2,
    QuantumNumberType::v3,
    QuantumNumberType::v4,
    QuantumNumberType::v5,
    QuantumNumberType::v6,
    QuantumNumberType::v7,
    QuantumNumberType::v8,
    QuantumNumberType::v9,
    QuantumNumberType::vibInv,
    QuantumNumberType::vibRefl,
    QuantumNumberType::vibSym};

//! A default state of local quantum numbers
[[maybe_unused]] inline constexpr std::array local_types{QuantumNumberType::F,
                                                         QuantumNumberType::F1,
                                                         QuantumNumberType::F10,
                                                         QuantumNumberType::F11,
                                                         QuantumNumberType::F12,
                                                         QuantumNumberType::F2,
                                                         QuantumNumberType::F3,
                                                         QuantumNumberType::F4,
                                                         QuantumNumberType::F5,
                                                         QuantumNumberType::F6,
                                                         QuantumNumberType::F7,
                                                         QuantumNumberType::F8,
                                                         QuantumNumberType::F9,
                                                         QuantumNumberType::I,
                                                         QuantumNumberType::J,
                                                         QuantumNumberType::K,
                                                         QuantumNumberType::Ka,
                                                         QuantumNumberType::Kc,
                                                         QuantumNumberType::N};

std::ostream& operator<<(std::ostream& os, const Array<GlobalState>& a);
}  // namespace Quantum::Number

using QuantumNumberValue       = Quantum::Number::Value;
using QuantumNumberValueList   = Quantum::Number::ValueList;
using QuantumNumberLocalState  = Quantum::Number::LocalState;
using QuantumIdentifier        = Quantum::Number::GlobalState;
using ArrayOfQuantumIdentifier = Array<QuantumIdentifier>;

std::ostream& operator<<(std::ostream& os, const Array<QuantumNumberType>& a);

namespace std {
//! Allow Quantum::Number::TwoLevelValueHolder to be used in hashes
template <>
struct hash<Quantum::Number::TwoLevelValueHolder> {
  std::size_t operator()(const Quantum::Number::TwoLevelValueHolder& g) const {
    return std::hash<std::string_view>{}(g.upp.s.val()) ^
           (std::hash<std::string_view>{}(g.low.s.val()) << 1);
  }
};

//! Allow QuantumNumberValueList to be used in hashes
template <>
struct hash<QuantumNumberValueList> {
  std::size_t operator()(const QuantumNumberValueList& g) const {
    std::size_t out = 0;
    std::size_t i   = 1;
    for (auto& x : g) {
      out ^= (std::hash<QuantumNumberType>{}(x.type) ^
              (std::hash<Quantum::Number::TwoLevelValueHolder>{}(x.qn) << 1))
             << i;
      i++;
    }
    return out;
  }
};

//! Allow QuantumNumberLocalState to be used in hashes
struct LocalStateHash {
  std::size_t operator()(const QuantumNumberLocalState& g) const {
    return std::hash<QuantumNumberValueList>{}(g.val);
  }
};

//! Allow QuantumIdentifier to be used in hashes
template <>
struct hash<QuantumIdentifier> {
  std::size_t operator()(const QuantumIdentifier& g) const {
    return static_cast<std::size_t>(g.isotopologue_index) ^
           (std::hash<QuantumNumberValueList>{}(g.val) << 1);
  }
};
}  // namespace std

template <>
struct std::formatter<QuantumNumberValue> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumNumberValue& v,
                              FmtContext& ctx) const {
    const auto quote = tags.quote();
    const auto sep   = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.type,
                sep,
                quote,
                v.str_upp(),
                quote,
                sep,
                quote,
                v.str_low(),
                quote);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<QuantumNumberValueList> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumNumberValueList& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.values);
  }
};

template <>
struct std::formatter<QuantumNumberLocalState> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumNumberLocalState& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.val);
  }
};

template <>
struct std::formatter<QuantumIdentifier> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumIdentifier& v,
                              FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.Isotopologue().FullName(), tags.sep(), v.val);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<Quantum::Number::QuantumNumberValueType> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Quantum::Number::QuantumNumberValueType& v,
                              FmtContext& ctx) const {
    const auto quote = tags.quote();

    switch (v) {
      case Quantum::Number::QuantumNumberValueType::I:
        return format_to(ctx.out(), "{}I{}", quote, quote);
      case Quantum::Number::QuantumNumberValueType::S:
        return format_to(ctx.out(), "{}S{}", quote, quote);
      case Quantum::Number::QuantumNumberValueType::H:
        return format_to(ctx.out(), "{}H{}", quote, quote);
    }

    return format_to(ctx.out(), "{}bad-value{}", quote, quote);
  }
};
