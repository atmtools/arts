#ifndef quantun_numbers_h
#define quantun_numbers_h

#include <algorithm>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "debug.h"
#include "enums.h"
#include "matpack.h"
#include "nonstd.h"
#include "rational.h"

namespace QuantumNumber {

//! Holds string values but can only hold sizeof(Index) long values
struct StringValue {
  static constexpr std::size_t N = sizeof(Index);
  std::array<char, N> x;

  //! Returns the value in such a way that no \0 remains in the view
  [[nodiscard]] constexpr std::string_view val() const noexcept {
    for (std::size_t i = 0; i < N; i++)
      if (x[i] == '\0') return {x.data(), i};
    return {x.data(), N};
  }

  constexpr StringValue() : x() {}

  explicit constexpr StringValue(std::string_view s) : x() {
    const std::size_t n = s.size();
    ARTS_USER_ERROR_IF(n > N,
                       "The value \"",
                       s,
                       "\" is too long.  Can be only ",
                       N,
                       " chars but is ",
                       n)

    // Fill with correct values or zero characters
    for (std::size_t i = 0; i < n; i++) x[i] = s[i];
    for (std::size_t i = n; i < N; i++) x[i] = '\0';
  }
};

//! Holds integer values
struct IntegerValue {
  Index x;
  [[nodiscard]] constexpr Rational val() const noexcept { return x; }
};

//! Holds half integer values, but only its denominator
struct HalfIntegerValue {
  Index x;
  [[nodiscard]] constexpr Rational val() const noexcept {
    return Rational(x, 2);
  }
};

//! Value type tags, S: str, I: index, H: half-index
ENUMCLASS(ValueType, char, S, I, H)

//! Different types of quantum numbers according to VAMDC with some unbounded counted up to 12 times
ENUMCLASS(Type,
          char,
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
          parity,
          r,
          rotSym,
          rovibSym,
          sym,
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

//! Common value type of a given quantum number type
constexpr ValueType common_value_type(Type type) noexcept {
  switch (type) {
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
  if (a == b or (a == ValueType::H and b == ValueType::I)) return a;
  if (b == ValueType::H and a == ValueType::I) return b;
  return ValueType::FINAL;
}

//! A union of the three type of values we need to consider
union ValueHolder {
  StringValue s;
  IntegerValue i;
  HalfIntegerValue h;
  constexpr ValueHolder() noexcept : i() {}
};

/** A complete description of a value, its type and value
 * 
 * Intended to be returned from IO operations
 */
struct ValueDescription {
  ValueType type;
  ValueHolder val;

  //! Debug output only
  friend std::ostream& operator<<(std::ostream& os, ValueDescription x) {
    os << x.type << ' ';
    switch (x.type) {
      case ValueType::S:
        os << x.val.s.val();
        break;
      case ValueType::I:
        os << x.val.i.x;
        break;
      case ValueType::H:
        os << Rational(x.val.h.x, 2);
        break;
      case ValueType::FINAL: {
      }
    }
    return os;
  }
};

/** The values of two levels
 * 
 * Its ValueDescription constructor ensures that we have two valid 
 * types and converts integer input in one to half-integer in case
 * there's a 'mismatch'
 */
struct TwoLevelValueHolder {
  ValueHolder upp, low;
  constexpr TwoLevelValueHolder(ValueDescription u, ValueDescription l)
      : upp(u.val), low(l.val) {
    if (u.type not_eq l.type) {
      ARTS_USER_ERROR_IF(
          u.type == ValueType::S,
          "The upper quantum number is a string but not the lower.  Value: ",
          upp.s.val())
      ARTS_USER_ERROR_IF(
          l.type == ValueType::S,
          "The lower quantum number is a string but not the upper.  Value: ",
          low.s.val())

      // either u or l is half
      if (u.type == ValueType::H)
        low.h.x = 2 * l.val.i.x;
      else
        upp.h.x = 2 * u.val.i.x;
    }
  }

  TwoLevelValueHolder() = default;
};

/** Takes a rational and determine which type of quantum number it is,
 * returning this information or throwing a runtime error if there's
 * an error
 *
 * @param r_ A rational
 * @return constexpr ValueDescription 
 */
[[nodiscard]] constexpr ValueDescription value_holder(Rational r_) {
  ValueDescription x{};

  const Rational r = reduce_by_gcd(r_);

  if (r.Denom() == 2) {
    x.type = ValueType::H;
    x.val.h.x = r.Nom();
  } else if (r.Denom() == 1) {
    x.type = ValueType::I;
    x.val.i.x = r.Nom();
  } else {
    ARTS_USER_ERROR(
        "Cannot convert ", r, " to half-integer or to full integer value")
  }

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
 * @param s Some view of a string
 * @return constexpr Rational 
 */
[[nodiscard]] constexpr Rational cast_as_qn_rational(
    std::string_view s) noexcept {
  int div = 0, dot = 0, any = 0;
  for (auto x : s) {
    if (x == '/')
      div++;  // Count divs because we can have at most one
    else if (x == '.')
      dot++;  // Count dots for the same reason
    else if (not nonstd::isdigit(x))
      return RATIONAL_UNDEFINED;
    any++;
  }

  // Can only have one of div or dot but is otherwise a Rational
  if ((div + dot) > 1 or any == 0) return RATIONAL_UNDEFINED;

  if (div) {
    Index n = 0, d = 0;
    std::size_t i = 0;
    for (; s[i] not_eq '/'; ++i) {
      n *= 10;
      n += s[i] - '0';
    }
    for (; i < s.size(); ++i) {
      d *= 10;
      d += s[i] - '0';
    }
    return Rational(n, d);
  }

  if (dot) {
    Index f = 0;
    std::size_t i = 0;
    std::size_t n = s.size();

    for (; s[i] not_eq '.'; ++i) {
      f *= 10;
      f += s[i] - '0';
    }

    if (i + 1 > n) return RATIONAL_UNDEFINED;
    if (i + 1 == n or s.back() == '0') return f;
    if (s.back() == '5') return Rational(2 * f + 1, 2);
    return RATIONAL_UNDEFINED;
  }

  Index n = 0;
  for (auto x : s) {
    n *= 10;
    n += x - '0';
  }
  return n;
}

/** Returns a value description for the quantum number
 * 
 * Note that several branches can throw as the input is assumed to be from
 * the user
 *
 * @param s Some view of a string
 * @return constexpr ValueDescription 
 */
[[nodiscard]] constexpr ValueDescription value_holder(std::string_view s) {
  if (const Rational r = cast_as_qn_rational(s); r.isDefined())
    return value_holder(r);

  ValueDescription x{};
  x.type = ValueType::S;
  x.val.s = StringValue(s);
  return x;
}

//! Struct that converts to bool automatically but allows checking both energy levels matching status
struct LevelMatch {
  bool upp{true}, low{true};
  constexpr LevelMatch operator&&(LevelMatch that) const noexcept {
    return {upp and that.upp, low and that.low};
  }
  constexpr operator bool() const noexcept { return upp and low; }
};

//! Count all space-separated items in s
constexpr Index count_items(std::string_view s) noexcept {
  bool currently_space = true;
  Index count = 0;
  for (auto& x : s) {
    bool space = nonstd::isspace(x);
    if (currently_space and not space) count++;
    currently_space = space;
  }
  return count;
}

//! Access items by index (note, first item is 1, not 0)
template<std::size_t n=1>
constexpr std::string_view items(std::string_view s,
                                 std::size_t i) noexcept {
  static_assert(n > 0, "Must want some items");
  bool currently_space = true, found = false;
  std::size_t beg = 0, count = 0, end = s.size();

  for (std::size_t ind = 0; ind < end; ind++) {
    bool space = nonstd::isspace(s[ind]);
    if (space and count == i + n - 1) return {&s[beg], ind - beg};
    if (currently_space and not space) count++;
    if (count == i and currently_space) {
      beg = ind;
      found = true;
    }
    currently_space = space;
  }
  if (found) return {&s[beg], end - beg};
  return {&s[0], 0};
}

//! A complete quantum number value with type information
struct Value {
  Type type;
  ValueType value_type;
  TwoLevelValueHolder qn;

  Value() = default;

  constexpr Value(std::string_view s)
      : type(toTypeOrThrow(items(s, 1))), value_type(), qn() {
    ARTS_USER_ERROR_IF(count_items(s) not_eq 3,
                       "Must have ' TYPE UPPNUM LOWNUM ' but got: '",
                       s,
                       "'")

    auto upv = value_holder(items(s, 2));
    auto lov = value_holder(items(s, 3));
    value_type = common_value_type(upv.type, lov.type);  // Ignore error flag
    qn = TwoLevelValueHolder(upv, lov);
  }

  //! Returns the upper quantum number rational if it exists or an undefined
  [[nodiscard]] constexpr Rational upp() const noexcept {
    switch (value_type) {
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
    switch (value_type) {
      case ValueType::I:
        return qn.low.i.val();
      case ValueType::H:
        return qn.low.h.val();
      default: {
      }
    }
    return RATIONAL_UNDEFINED;
  }

  //! Returns the upper quantum number string view if it exist or throws
  [[nodiscard]] String str_upp() const noexcept {
    if (ValueType::S == value_type) return qn.upp.s.val();
    return var_string(upp());
  }

  //! Returns the lower quantum number string view if it exist or throws
  [[nodiscard]] String str_low() const noexcept {
    if (ValueType::S == value_type) return qn.low.s.val();
    return var_string(low());
  }

  //! Checks if this value is good.  Assumed true after construction
  [[nodiscard]] constexpr bool good() const {
    const auto common = common_value_type(type);

    // Check that the common value type for this quantum number type is
    // also common with the stored value here
    return common == common_value_type(common, value_type);
  }

  /** Returns a description of whether both levels match
   * 
   * Note that the LevelMatch type should automatically transform to bool
   * so there's no need for extra work if this is your target question
   *
   * @param other Another value
   * @return constexpr LevelMatch
   */
  constexpr LevelMatch operator==(Value other) const noexcept {
    if (type == other.type and
        good_enum(common_value_type(value_type, other.value_type))) {
      switch (value_type) {
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
  friend std::ostream& operator<<(std::ostream& os, Value x) {
    return os << x.type << ' ' << x.str_upp() << ' ' << x.str_low();
  }

  //! Standard input
  friend std::istream& operator>>(std::istream& is, Value& x) {
    std::string upp, low;
    is >> x.type >> upp >> low;

    // Let the x.qn constructor deal with type errors and value_holder deal with read errors
    auto upv = value_holder(upp);
    auto lov = value_holder(low);
    x.value_type = common_value_type(upv.type, lov.type);  // Ignore error flag
    x.qn = TwoLevelValueHolder(upv, lov);

    ARTS_USER_ERROR_IF(not x.good(),
                       "Cannot understand: ",
                       x.type,
                       ' ',
                       upp,
                       ' ',
                       low,
                       '\n',
                       "Expected type is: ",
                       common_value_type(x.type),
                       " but got:\n"
                       "\tupper type: ",
                       upv,
                       "\n\tlower type: ",
                       lov)
    return is;
  }
};

//! Status of comparing two lists that are supposedly of some type
ENUMCLASS(CheckValue, char, Full, AinB, BinA, Miss)

//! Level-by-level version of CheckValue
struct CheckMatch {
  CheckValue upp;
  CheckValue low;
  constexpr operator bool() const noexcept { return upp == CheckValue::Full and low== CheckValue::Full; }
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

//! A list of many quantum numbers.  Should always remain sorted
class ValueList {
  Array<Value> values;

  //! Internal sort function.  Should be called whenever new items are created
  void sort_by_type() {
    std::sort(values.begin(), values.end(), [](auto& a, auto& b) {
      return a.type < b.type;
    });
  }

  //! Internal check function.  Remember to sort by type before calling this
  [[nodiscard]] bool has_unique_increasing_types() const {
    return std::adjacent_find(
               values.begin(), values.end(), [](auto& a, auto& b) {
                 return a.type >= b.type;
               }) == values.end();
  }

 public:
  ValueList(std::string_view s) : values(0) {
    Index n = count_items(s);
    ARTS_USER_ERROR_IF(
        n % 3, "Must have multiple of three items, got ", n, " in:\n", s)
    for (Index i = 1; i < n; i += 3) values.emplace_back(Value(items<3>(s, i)));
    sort_by_type();
    ARTS_USER_ERROR_IF(not has_unique_increasing_types(),
                       "Bad Value List (contains copies or cannot sort):\n",
                       *this)
  }

  //! Should always be called before this object is handed to another user
  void finalize() {
    sort_by_type();
    ARTS_USER_ERROR_IF(not has_unique_increasing_types(),
                       "The quantum number list: [",
                       *this,
                       "] contains copies of types")
  }

  //! Finds whether two ValueList describe completely different sets of quantum numbers (e.g., local vs global)
  [[nodiscard]] bool perpendicular(const ValueList& that) const {
    ARTS_ASSERT(has_unique_increasing_types())
    ARTS_ASSERT(that.has_unique_increasing_types())

    auto this_val = values.begin();
    auto that_val = that.values.begin();
    while (that_val not_eq that.values.end() and this_val not_eq values.end()) {
      if (that_val->type < this_val->type)
        that_val++;
      else if (this_val->type < that_val->type)
        this_val++;
      else
        return false;
    }
    return true;
  }

  //! Returns whether the Type is part of the list
  [[nodiscard]] bool has(Type t) const noexcept {
    return std::find_if(values.begin(), values.end(), [t](auto& x) {
             return x.type == t;
           }) not_eq values.end();
  }

  //! Returns the value of the Type (assumes it exist)
  const Value& operator[](Type t) const ARTS_NOEXCEPT {
    auto val = std::find_if(
        values.begin(), values.end(), [t](auto& x) { return x.type == t; });
    ARTS_ASSERT(val not_eq values.end())
    return *val;
  }

  friend std::ostream& operator<<(std::ostream& os, const ValueList& vl) {
    for (Index i = 0; i < vl.values.nelem(); i++) {
      if (i) os << ' ';
      os << vl.values[i];
    }
    return os;
  }

  //! input stream must have pre-set size
  friend std::istream& operator>>(std::istream& is, ValueList& vl) {
    for (auto& x : vl.values) is >> x;
    vl.sort_by_type();
    ARTS_USER_ERROR_IF(not vl.has_unique_increasing_types(),
                       "Bad Value List (contains copies or cannot sort):\n",
                       vl)
    return is;
  }
};

/** Check if two lists of quantum numbers agree with eachother
 *
 * @param[in] A a list of quantum numbers
 * @return Full: perfect match;                                A is B
 * @return AinB: all A == B but there are missing Types in A;  A in B
 * @return BinA: all B == A but there are missing Types in B;  B in A
 * @return Miss: no match;                                     A != B
 */
inline CheckMatch checkValueList(const ValueList& a, const ValueList& b) {
  CheckMatch status = {CheckValue::Full, CheckValue::Full};

  for (Type t : enumtyps::TypeTypes) {
    const bool ahas = a.has(t), bhas = b.has(t); 

    if (ahas and bhas) {
      const LevelMatch levels = a[t] == b[t];
      if (not levels.upp) status.upp = CheckValue::Miss;
      if (not levels.low) status.low = CheckValue::Miss;
    } else if (ahas and not bhas) {
      status = {update(status.upp, CheckValue::BinA), update(status.low, CheckValue::BinA)};
    } else if (not ahas and bhas) {
      status = {update(status.upp, CheckValue::AinB), update(status.low, CheckValue::AinB)};
    }
  }
  return status;
}
}  // namespace QuantumNumber

using QuantumNumberValue = QuantumNumber::Value;
using ArrayOfQuantumNumberValue = Array<QuantumNumber::Value>;

#endif  // quantun_numbers_h
