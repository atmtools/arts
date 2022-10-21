#ifndef quantun_numbers_h
#define quantun_numbers_h

#include <algorithm>
#include <compare>
#include <cstddef>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "debug.h"
#include "enums.h"
#include "isotopologues.h"
#include "matpack.h"
#include "nonstd.h"
#include "rational.h"

constexpr Index quantum_number_error_value = -999'999'999;

namespace Quantum::Number {

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
    for (std::size_t i=0; i<N; i++) {
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

  constexpr std::strong_ordering operator<=>(const IntegerValue& i) const {return x<=> i.x;}
};

//! Holds half integer values, but only its denominator
struct HalfIntegerValue {
  Index x{std::numeric_limits<Index>::lowest()};

  //! Returns the value as a Rational, keeping that this is a half-integer
  [[nodiscard]] constexpr Rational val() const noexcept {
    return Rational(x, 2);
  }

  constexpr std::strong_ordering operator<=>(const HalfIntegerValue& h) const {return x<=> h.x;}
};

//! Three tags, S: str, I: index, H: half-index
ENUMCLASS(ValueType, char, S, I, H)

//! Different types of quantum numbers according to VAMDC with some unbounded counted up to 12 times
ENUMCLASS(Type,
          char,
          alpha,  // FIXME: Not in VAMDC
          config,  // FIXME: Not in VAMDC
          ElecStateLabel,
          F,
          F1,
          F10,
          F11,
          F12,
          F2,
          F3,
          F4,
          F5,
          F6,
          F7,
          F8,
          F9,
          I,
          J,
          K,
          Ka,
          Kc,
          L,  // FIXME: Not in VAMDC
          Lambda,
          N,
          Omega,
          S,
          Sigma,
          SpinComponentLabel,
          asSym,
          elecInv,
          elecRefl,
          elecSym,
          kronigParity,
          l,
          l1,
          l10,
          l11,
          l12,
          l2,
          l3,
          l4,
          l5,
          l6,
          l7,
          l8,
          l9,
          n,  // FIXME: Not in VAMDC
          parity,
          r,
          rotSym,
          rovibSym,
          sym,
          tau,  // FIXME: Not in VAMDC
          term,  // FIXME: Not in VAMDC
          v,
          v1,
          v10,
          v11,
          v12,
          v2,
          v3,
          v4,
          v5,
          v6,
          v7,
          v8,
          v9,
          vibInv,
          vibRefl,
          vibSym)

/** Common value type of a given quantum number type 
 * 
 * Guards against integer/half-integer comparisons failing for bad reasons
 * and allows IO to not need some information 
 *
 * Keep new entries to a 1-to-1 key-to-key match, and sort them the same as Type for ease of reading
 *
 * @param type A quantum number type
 * @return constexpr ValueType The common type
 */
constexpr ValueType common_value_type(Type type) noexcept {
  switch (type) {
    case Type::alpha:
      return ValueType::S;
    case Type::config:
      return ValueType::S;
    case Type::ElecStateLabel:
      return ValueType::S;
    case Type::F:
      return ValueType::H;
    case Type::F1:
      return ValueType::H;
    case Type::F10:
      return ValueType::H;
    case Type::F11:
      return ValueType::H;
    case Type::F12:
      return ValueType::H;
    case Type::F2:
      return ValueType::H;
    case Type::F3:
      return ValueType::H;
    case Type::F4:
      return ValueType::H;
    case Type::F5:
      return ValueType::H;
    case Type::F6:
      return ValueType::H;
    case Type::F7:
      return ValueType::H;
    case Type::F8:
      return ValueType::H;
    case Type::F9:
      return ValueType::H;
    case Type::I:
      return ValueType::I;
    case Type::J:
      return ValueType::H;
    case Type::K:
      return ValueType::I;
    case Type::Ka:
      return ValueType::I;
    case Type::Kc:
      return ValueType::I;
    case Type::L:
      return ValueType::I;
    case Type::Lambda:
      return ValueType::I;
    case Type::N:
      return ValueType::I;
    case Type::Omega:
      return ValueType::H;
    case Type::S:
      return ValueType::H;
    case Type::Sigma:
      return ValueType::H;
    case Type::SpinComponentLabel:
      return ValueType::I;
    case Type::asSym:
      return ValueType::S;
    case Type::elecInv:
      return ValueType::S;
    case Type::elecRefl:
      return ValueType::S;
    case Type::elecSym:
      return ValueType::S;
    case Type::kronigParity:
      return ValueType::S;
    case Type::l:
      return ValueType::I;
    case Type::l1:
      return ValueType::I;
    case Type::l10:
      return ValueType::I;
    case Type::l11:
      return ValueType::I;
    case Type::l12:
      return ValueType::I;
    case Type::l2:
      return ValueType::I;
    case Type::l3:
      return ValueType::I;
    case Type::l4:
      return ValueType::I;
    case Type::l5:
      return ValueType::I;
    case Type::l6:
      return ValueType::I;
    case Type::l7:
      return ValueType::I;
    case Type::l8:
      return ValueType::I;
    case Type::l9:
      return ValueType::I;
    case Type::n:
      return ValueType::S;
    case Type::parity:
      return ValueType::S;
    case Type::r:
      return ValueType::I;
    case Type::rotSym:
      return ValueType::S;
    case Type::rovibSym:
      return ValueType::S;
    case Type::sym:
      return ValueType::S;
    case Type::tau:
      return ValueType::S;
    case Type::term:
      return ValueType::S;
    case Type::v:
      return ValueType::I;
    case Type::v1:
      return ValueType::I;
    case Type::v10:
      return ValueType::I;
    case Type::v11:
      return ValueType::I;
    case Type::v12:
      return ValueType::I;
    case Type::v2:
      return ValueType::I;
    case Type::v3:
      return ValueType::I;
    case Type::v4:
      return ValueType::I;
    case Type::v5:
      return ValueType::I;
    case Type::v6:
      return ValueType::I;
    case Type::v7:
      return ValueType::I;
    case Type::v8:
      return ValueType::I;
    case Type::v9:
      return ValueType::I;
    case Type::vibInv:
      return ValueType::S;
    case Type::vibRefl:
      return ValueType::S;
    case Type::vibSym:
      return ValueType::S;
    case Type::FINAL: {
    }
  }
  return ValueType::FINAL;
}

/** Return a common type between a and b
 * 
 * If they are the same, returns the value
 * If either is H and the other is I, returns H
 * Otherwise returns FINAL as an error flag
 *
 * @param a A value type
 * @param b Another value type
 * @return constexpr ValueType with FINAL as error flag
 */
constexpr ValueType common_value_type(ValueType a, ValueType b) noexcept {
  // Same is same, H is I, a has both:
  if (a == b or (a == ValueType::H and b == ValueType::I)) return a;

  // H is I:
  if (b == ValueType::H and a == ValueType::I) return ValueType::H;

  // Something has gone very wrong:
  return ValueType::FINAL;
}

//! A union of the three type of values we need to consider
union ValueHolder {
  StringValue s;
  IntegerValue i;
  HalfIntegerValue h;

  constexpr ValueHolder(ValueType t) noexcept : s(StringValue{"NODEF"}) {
    switch (t) {
      case ValueType::H:
        h = HalfIntegerValue{0};
        break;
      case ValueType::I:
        i = IntegerValue{0};
        break;
      default: {
      }
    }
  }

  constexpr ValueHolder(Type t) noexcept : ValueHolder(common_value_type(t)) {}
  constexpr ValueHolder(const ValueHolder&) = default;
  constexpr ValueHolder(ValueHolder&&) noexcept = default;
  constexpr ValueHolder& operator=(const ValueHolder&) = default;
  constexpr ValueHolder& operator=(ValueHolder&&) noexcept = default;
};

/** A complete description of a value, its type and value
 * 
 * Intended to be returned from IO operations
 */
struct ValueDescription {
  ValueType type;
  ValueHolder val;

  constexpr ValueDescription(ValueType t) noexcept : type(t), val(t) {}
  constexpr ValueDescription(const ValueDescription&) = default;
  constexpr ValueDescription(ValueDescription&&) noexcept = default;
  constexpr ValueDescription& operator=(const ValueDescription&) = default;
  constexpr ValueDescription& operator=(ValueDescription&&) noexcept = default;

  //! Debug output only
  friend std::ostream& operator<<(std::ostream& os, ValueDescription x);
};

constexpr std::strong_ordering cmp(std::strong_ordering&& a, std::strong_ordering&& b) {
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
  constexpr TwoLevelValueHolder(ValueDescription u, ValueDescription l, Type t)
      : upp(u.val), low(l.val) {
    auto ct = common_value_type(t);
    if (u.type not_eq ct) {
      ARTS_USER_ERROR_IF(u.type not_eq ValueType::I,
                         "Cannot convert from ",
                         u.type,
                         " to ",
                         ct)
      upp.h.x = 2 * u.val.i.x;
    }

    if (l.type not_eq ct) {
      ARTS_USER_ERROR_IF(l.type not_eq ValueType::I,
                         "Cannot convert from ",
                         l.type,
                         " to ",
                         ct)
      low.h.x = 2 * l.val.i.x;
    }
  }

  constexpr TwoLevelValueHolder(Type t) noexcept : upp(t), low(t) {}
  constexpr TwoLevelValueHolder(const TwoLevelValueHolder&) = default;
  constexpr TwoLevelValueHolder(TwoLevelValueHolder&&) noexcept = default;
  constexpr TwoLevelValueHolder& operator=(const TwoLevelValueHolder&) = default;
  constexpr TwoLevelValueHolder& operator=(TwoLevelValueHolder&&) noexcept = default;

  [[nodiscard]] constexpr std::strong_ordering order(const TwoLevelValueHolder& tv, ValueType t) const {
    switch (t) {
    case ValueType::S:
      return cmp(upp.s <=> tv.upp.s, low.s <=> tv.low.s);
    case ValueType::I:
      return cmp(upp.i <=> tv.upp.i, low.i <=> tv.low.i);
    case ValueType::H:
      return cmp(upp.h <=> tv.upp.h, low.h <=> tv.low.h);
    case ValueType::FINAL: break;
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

  // We must now have a half-integer or not
  if (r.denom == 2) {
    ValueDescription x{ValueType::H};
    x.val.h.x = r.numer;
    return x;
  }
  
  if (r.denom == 1) {
    ValueDescription x{ValueType::I};
    x.val.i.x = r.numer;
    return x;
  }

  ValueDescription x{ValueType::I};
  x.val.i.x = quantum_number_error_value;
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
  std::size_t const n = s.size();

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
                                                      Type t) {
  switch (common_value_type(t)) {
    case ValueType::I: 
    case ValueType::H: {
      const Rational r = cast_qnrat(s);
      return value_holder(r);
    }
    case ValueType::S: {
      ValueDescription x{ValueType::S};
      x.val.s = StringValue(s);
      return x;
    }
    case ValueType::FINAL: {
    }
  }
  
  return ValueType::FINAL;
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
    bool const this_space = nonstd::isspace(x);

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
    bool const this_space = nonstd::isspace(s[ind]);

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
  Type type;
  TwoLevelValueHolder qn;

  constexpr std::strong_ordering operator<=>(const Value& v) const {
    if (type < v.type) return std::strong_ordering::less;
    if (v.type < type) return std::strong_ordering::greater;
    return qn.order(v.qn, common_value_type(type));
  }

  constexpr Value(Type t=Type::FINAL) : type(t), qn(type) {}
  Value(const Value&) = default;
  Value(Value&&) noexcept = default;
  Value& operator=(const Value&) = default;
  Value& operator=(Value&&) noexcept = default;

  constexpr Value(Type t, Rational upp_, Rational low_) : Value(t) {
    Rational upp = reduce_by_gcd(upp_), low = reduce_by_gcd(low_);

    if (common_value_type(type) == ValueType::H) {
      ARTS_ASSERT(upp.denom <= 2 and low.denom <= 2)
      if (upp.denom not_eq 2) upp *= 2;
      if (low.denom not_eq 2) low *= 2;
      qn.upp.h.x = upp.numer;
      qn.low.h.x = low.numer;
    } else if (common_value_type(type) == ValueType::I) {
      ARTS_ASSERT(upp.denom == 1 and low.denom == 1)
      qn.upp.i.x = upp.numer;
      qn.low.i.x = low.numer;
    } else {
      ARTS_USER_ERROR(
          t, " is a string-type, so cannot be constructed from rationals")
    }
  }

  //! Default constructor from some string of values
  constexpr Value(std::string_view s) : Value(toTypeOrThrow(items(s, 0))) {
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
      case ValueType::I:
        return qn.upp.i.val();
      case ValueType::H:
        return qn.upp.h.val();
      default: {
      }
    }
    return RATIONAL_UNDEFINED;
  }

  //! Returns the lower quantum number rational if it exists or an undefined
  [[nodiscard]] constexpr Rational low() const noexcept {
    switch (common_value_type(type)) {
      case ValueType::I:
        return qn.low.i.val();
      case ValueType::H:
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
    ValueDescription const v = value_holder(s, type);
    TwoLevelValueHolder const nqn(v, v, type);
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
        case ValueType::I:
        case ValueType::H:
          return {upp() == other.upp(), low() == other.low()};
        case ValueType::S:
          return {qn.upp.s.x == other.qn.upp.s.x,
                  qn.low.s.x == other.qn.low.s.x};
        case ValueType::FINAL: {
        }
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

  [[nodiscard]] constexpr bool good() const {return level_match(*this);}
};

//! Status of comparing two lists that are supposedly of some type
ENUMCLASS(CheckValue, char, Full, AinB, BinA, Miss)

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
constexpr bool is_sorted(const std::array<Type, N>& types) noexcept {
  for (size_t i = 1; i < N; i++)
    if (not(types[i - 1] < types[i])) return false;
  return true;
}

//! A list of many quantum numbers.  Should always remain sorted
class ValueList {
  Array<Value> values;

  //! Internal sort function.  Should be called whenever new items are created
  void sort_by_type();

  //! Internal check function.  Remember to sort by type before calling this
  [[nodiscard]] bool has_unique_increasing_types() const;

 public:
  //! From text
  explicit ValueList(std::string_view s, bool legacy = false);

  //! From legacy text
  ValueList(std::string_view upp, std::string_view low);

  std::strong_ordering operator<=>(const ValueList& v) const {
    const std::size_t n = std::min(values.size(), v.values.size());
    for (std::size_t i=0; i<n; i++) {
      if (values[i] < v.values[i]) return std::strong_ordering::less;
      if (v.values[i] < values[i]) return std::strong_ordering::greater;
    }
    return values.size() <=> v.values.size();
  }

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
  [[nodiscard]] Index nelem() const { return values.nelem(); }

  //! Finds whether two ValueList describe completely different sets of quantum numbers (e.g., local vs global)
  [[nodiscard]] bool perpendicular(const ValueList& that) const ARTS_NOEXCEPT;

  //! Returns whether all the Types are part of the list, the types must be sorted
  template <typename... Types>
  [[nodiscard]] bool has(Types... ts) const ARTS_NOEXCEPT {
    static_assert(sizeof...(Types) > 0);

    ARTS_ASSERT(is_sorted(std::array{Type(ts)...}))

    auto ptr = cbegin();
    auto end = cend();
    for (Type t : {Type(ts)...}) {
      ptr = std::find_if(ptr, end, [t](auto& x) { return x.type == t; });
      if (ptr == end) return false;
    }
    return true;
  }

  //! Returns the value of the Type (assumes it exist)
  const Value& operator[](Type t) const ARTS_NOEXCEPT;

  //! Legacy manipulation operator access
  Value& operator[](Index i) { return values.at(i); }

  //! Add for manipulation
  Value& add(Type t);

  //! Add for manipulation
  Value& add(Value v);

  //! Sets the value if it exists or adds it otherwise
  void set(Value v);

  //! Set a value in value list
  void set(Index i, std::string_view upp, std::string_view low);

  //! Returns upper and lower matching status
  [[nodiscard]] CheckMatch check_match(const ValueList& other) const noexcept;

  //! ouptut stream if all values
  friend std::ostream& operator<<(std::ostream& os, const ValueList& vl);

  //! input stream must have pre-set size
  friend std::istream& operator>>(std::istream& is, ValueList& vl);

  //! Add a type without sorting (WARNING, many things might break if you don't sort in the end)
  void add_type_wo_sort(Type);

  [[nodiscard]] bool good() const;
};

ValueList from_hitran(std::string_view upp, std::string_view low);

//! A logical struct for local quantum numbers
struct LocalState {
  ValueList val{};

  std::strong_ordering operator<=>(const LocalState& l) const {
    return val <=> l.val;
  }
  auto operator==(const LocalState& l) const {return std::strong_ordering::equal == (*this <=> l);}
  auto operator!=(const LocalState& l) const {return std::strong_ordering::equal != (*this <=> l);}

  LocalState() = default;

  template <typename... Values>
  LocalState(Values... vals) : val(Array<Value>{Value(vals)...}) {}

  void set_unsorted_qns(const Array<Type>& vals);

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

struct LevelTest {bool upp{true}, low{true};};

//! A logical struct for global quantum numbers with species identifiers
struct GlobalState {
  static constexpr Index version = 1;  // Second version of quantum identifiers

  Index isotopologue_index{-1};
  ValueList val{};

  GlobalState() = default;

  explicit GlobalState(Index i, ValueList v = {})
      : isotopologue_index(i), val(std::move(v)) {}

  explicit GlobalState(const Species::IsotopeRecord& ir)
      : isotopologue_index(Species::find_species_index(ir)) {}

  explicit GlobalState(std::string_view s, Index v = version);

  [[nodiscard]] Species::IsotopeRecord Isotopologue() const noexcept;
  [[nodiscard]] Species::Species Species() const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const GlobalState& gs);

  friend std::istream& operator>>(std::istream& is, GlobalState& gs);

  [[nodiscard]] GlobalState LowerLevel() const;
  [[nodiscard]] GlobalState UpperLevel() const;

  std::strong_ordering operator<=>(const GlobalState& g) const {
    if (isotopologue_index < g.isotopologue_index) return std::strong_ordering::less;
    if (g.isotopologue_index < isotopologue_index) return std::strong_ordering::greater;
    return val <=> g.val;
  }
  auto operator==(const GlobalState& g) const {return std::strong_ordering::equal == (*this <=> g);}
  auto operator!=(const GlobalState& g) const {return std::strong_ordering::equal != (*this <=> g);}

  //! Checks wheter all of the LHS is part of RHS
  [[nodiscard]] bool part_of(const GlobalState& other) const;

  //! Checks wheter all of the LHS is part of any of the RHS
  [[nodiscard]] LevelTest part_of(const GlobalState& g, const LocalState& l) const;

  //! Test if there are bad quantum numbers (undefined ones) or if the isotopologue is not a normal target
  [[nodiscard]] bool good() const;
};

//! StateMatchType operates so that a check less than a level should be 'better', bar None
ENUMCLASS(StateMatchType, char, Full, Level, Isotopologue, Species, None)

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

  constexpr bool operator!=(StateMatchType x) const noexcept { return not((*this) == x); }
};

//! VAMDC classes of quantum number cases
ENUMCLASS(
    VAMDC,
    char,
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
)

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
[[maybe_unused]] constexpr std::array global_types{Type::alpha,
                                                   Type::config,
                                                   Type::ElecStateLabel,
                                                   Type::L,
                                                   Type::Lambda,
                                                   Type::Omega,
                                                   Type::S,
                                                   Type::Sigma,
                                                   Type::SpinComponentLabel,
                                                   Type::asSym,
                                                   Type::elecInv,
                                                   Type::elecRefl,
                                                   Type::elecSym,
                                                   Type::kronigParity,
                                                   Type::l,
                                                   Type::l1,
                                                   Type::l10,
                                                   Type::l11,
                                                   Type::l12,
                                                   Type::l2,
                                                   Type::l3,
                                                   Type::l4,
                                                   Type::l5,
                                                   Type::l6,
                                                   Type::l7,
                                                   Type::l8,
                                                   Type::l9,
                                                   Type::n,
                                                   Type::parity,
                                                   Type::r,
                                                   Type::rotSym,
                                                   Type::rovibSym,
                                                   Type::sym,
                                                   Type::tau,
                                                   Type::term,
                                                   Type::v,
                                                   Type::v1,
                                                   Type::v10,
                                                   Type::v11,
                                                   Type::v12,
                                                   Type::v2,
                                                   Type::v3,
                                                   Type::v4,
                                                   Type::v5,
                                                   Type::v6,
                                                   Type::v7,
                                                   Type::v8,
                                                   Type::v9,
                                                   Type::vibInv,
                                                   Type::vibRefl,
                                                   Type::vibSym};

//! A default state of local quantum numbers
[[maybe_unused]] constexpr std::array local_types{Type::F,
                                                  Type::F1,
                                                  Type::F10,
                                                  Type::F11,
                                                  Type::F12,
                                                  Type::F2,
                                                  Type::F3,
                                                  Type::F4,
                                                  Type::F5,
                                                  Type::F6,
                                                  Type::F7,
                                                  Type::F8,
                                                  Type::F9,
                                                  Type::I,
                                                  Type::J,
                                                  Type::K,
                                                  Type::Ka,
                                                  Type::Kc,
                                                  Type::N};
}  // namespace Quantum::Number

using QuantumNumberType = Quantum::Number::Type;
using QuantumNumberValue = Quantum::Number::Value;
using QuantumNumberValueList = Quantum::Number::ValueList;
using QuantumNumberLocalState = Quantum::Number::LocalState;
using QuantumIdentifier = Quantum::Number::GlobalState;
using ArrayOfQuantumIdentifier = Array<QuantumIdentifier>;

#endif  // quantun_numbers_h
