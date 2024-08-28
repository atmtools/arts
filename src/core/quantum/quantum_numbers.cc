#include "quantum_numbers.h"

#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "debug.h"
#include "isotopologues.h"
#include "species.h"

namespace Quantum::Number {
std::ostream& operator<<(std::ostream& os, QuantumNumberValueType x) {
  switch (x) {
    case QuantumNumberValueType::S:
      return os << "S";
    case QuantumNumberValueType::I:
      return os << "I";
    case QuantumNumberValueType::H:
      return os << "H";
  }
  return os;
}

std::string_view StringValue::val() const noexcept {
  // Look for the first \0 or return everything
  for (std::size_t i = 0; i < N; i++)
    if (x[i] == '\0') return {x.data(), i};
  return {x.data(), N};
}

StringValue::StringValue(std::string_view s) {
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

std::strong_ordering StringValue::operator<=>(const StringValue& sv) const {
  for (std::size_t i = 0; i < N; i++) {
    if (x[i] < sv.x[i]) return std::strong_ordering::less;
    if (sv.x[i] < x[i]) return std::strong_ordering::greater;
  }
  return std::strong_ordering::equal;
}

std::ostream& operator<<(std::ostream& os, ValueDescription x) {
  os << x.type << ' ';
  switch (x.type) {
    case QuantumNumberValueType::S:
      return os << x.val.s.val();
    case QuantumNumberValueType::I:
      return os << x.val.i.x;
    case QuantumNumberValueType::H:
      return os << Rational(x.val.h.x, 2);
  }
  return os;
}

QuantumNumberValueType common_value_type(QuantumNumberType type) noexcept {
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

QuantumNumberValueType common_value_type(QuantumNumberValueType a,
                                         QuantumNumberValueType b) noexcept {
  // Same is same, H is I, a has both:
  if (a == b or
      (a == QuantumNumberValueType::H and b == QuantumNumberValueType::I))
    return a;

  // H is I:
  if (b == QuantumNumberValueType::H and a == QuantumNumberValueType::I)
    return QuantumNumberValueType::H;

  return static_cast<QuantumNumberValueType>(-1);
}

ValueHolder::ValueHolder(QuantumNumberValueType t) noexcept
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

ValueHolder::ValueHolder(QuantumNumberType t) noexcept
    : ValueHolder(common_value_type(t)) {}

TwoLevelValueHolder::TwoLevelValueHolder(ValueDescription u,
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

std::strong_ordering TwoLevelValueHolder::order(
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

ValueDescription value_holder(Rational r_) {
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
ValueDescription value_holder(std::string_view s) {
  ValueDescription x{QuantumNumberValueType::S};
  x.val.s = StringValue(s);
  return x;
}

Rational cast_qnrat(std::string_view s) noexcept {
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

ValueDescription value_holder(std::string_view s, QuantumNumberType t) {
  switch (common_value_type(t)) {
    case QuantumNumberValueType::I:
    case QuantumNumberValueType::H:
      return value_holder(cast_qnrat(s));
    case QuantumNumberValueType::S:
      return value_holder(s);
  }
  std::unreachable();
}

Index count_items(std::string_view s) noexcept {
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

std::string_view rstrip(std::string_view x) {
  while (not x.empty() and nonstd::isspace(x.back())) x.remove_suffix(1);
  return x;
}

std::string_view lstrip(std::string_view x) {
  while (not x.empty() and nonstd::isspace(x.front())) x.remove_prefix(1);
  return x;
}

std::string_view strip(std::string_view x) { return rstrip(lstrip(x)); }

std::string_view items(std::string_view s, std::size_t i, std::size_t n) noexcept {  
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

std::strong_ordering Value::operator<=>(const Value& v) const {
  if (type < v.type) return std::strong_ordering::less;
  if (v.type < type) return std::strong_ordering::greater;
  return qn.order(v.qn, common_value_type(type));
}

Value::Value(QuantumNumberType t, Rational upp_, Rational low_) : Value(t) {
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

Value::Value(std::string_view s) : Value(to<QuantumNumberType>(items(s, 0))) {
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

Rational Value::upp() const noexcept {
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

Rational Value::low() const noexcept {
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

void Value::set(std::string_view s, bool upp) {
  const ValueDescription v = value_holder(s, type);
  const TwoLevelValueHolder nqn(v, v, type);
  if (upp) {
    qn.upp = nqn.upp;
  } else {
    qn.low = nqn.low;
  }
}

LevelMatch Value::level_match(Value other) const noexcept {
  if (type == other.type) {
    switch (common_value_type(type)) {
      case QuantumNumberValueType::I:
      case QuantumNumberValueType::H:
        return {upp() == other.upp(), low() == other.low()};
      case QuantumNumberValueType::S:
        return {qn.upp.s.x == other.qn.upp.s.x, qn.low.s.x == other.qn.low.s.x};
    }
  }
  return {false, false};
}

bool Value::good() const { return level_match(*this); }

CheckValue update(CheckValue val, CheckValue res) noexcept {
  if (val == CheckValue::Miss) return val;
  if (val == CheckValue::Full) return res;
  if (res == CheckValue::Miss) return res;
  if (res == CheckValue::Full) return val;
  if (val == res) return val;
  return CheckValue::Miss;
}

CheckMatch update(CheckMatch val, CheckValue res) noexcept {
  return {update(val.upp, res), update(val.low, res)};
}

CheckMatch update(CheckMatch val, CheckMatch res) noexcept {
  return {update(val.upp, res.upp), update(val.low, res.low)};
}

ValueList::ValueList(Array<Value> values_) : values(std::move(values_)) {
  finalize();
}

std::strong_ordering LocalState::operator<=>(const LocalState& l) const {
  return val <=> l.val;
}

bool LocalState::operator==(const LocalState& l) const {
  return std::strong_ordering::equal == (*this <=> l);
}

bool LocalState::operator!=(const LocalState& l) const {
  return std::strong_ordering::equal != (*this <=> l);
}

std::strong_ordering GlobalState::operator<=>(const GlobalState& g) const {
  if (isotopologue_index < g.isotopologue_index)
    return std::strong_ordering::less;
  if (g.isotopologue_index < isotopologue_index)
    return std::strong_ordering::greater;
  return val <=> g.val;
}

bool GlobalState::operator==(const GlobalState& g) const {
  return std::strong_ordering::equal == (*this <=> g);
}

bool GlobalState::operator!=(const GlobalState& g) const {
  return std::strong_ordering::equal != (*this <=> g);
}

bool StateMatch::operator==(StateMatchType x) const noexcept {
  return x == type;
}

bool StateMatch::operator!=(StateMatchType x) const noexcept {
  return not((*this) == x);
}

String Value::str_upp() const noexcept {
  if (QuantumNumberValueType::S == common_value_type(type))
    return String{qn.upp.s.val()};
  return var_string(upp());
}
String Value::str_low() const noexcept {
  if (QuantumNumberValueType::S == common_value_type(type))
    return String{qn.low.s.val()};
  return var_string(low());
}

std::ostream& operator<<(std::ostream& os, Value x) {
  return os << x.type << ' ' << x.str_upp() << ' ' << x.str_low();
}

std::istream& operator>>(std::istream& is, Value& x) {
  std::string upp, low;
  is >> x.type >> upp >> low;

  x.qn = TwoLevelValueHolder(
      value_holder(upp, x.type), value_holder(low, x.type), x.type);

  return is;
}

void ValueList::sort_by_type() {
  std::sort(begin(), end(), [](auto& a, auto& b) { return a.type < b.type; });
}

bool ValueList::has_unique_increasing_types() const {
  return std::adjacent_find(begin(), end(), [](auto& a, auto& b) {
           return a.type >= b.type;
         }) == values.end();
}

//! Fix legacy catalog, where some values are rationals even though they shouldn't be
std::pair<std::string_view, String> fix_legacy(std::string_view key,
                                               std::string_view val) {
  if (key == "ElectronState") {
    if (val == "X") return {"ElecStateLabel", "X"};
    if (val == "a") return {"ElecStateLabel", "a"};
    if (val == "b") return {"ElecStateLabel", "b"};
    if (val == "c") return {"ElecStateLabel", "c"};
    if (val == "A") return {"ElecStateLabel", "A"};
    if (val == "'") return {"ElecStateLabel", "'"};
    if (val == "B") return {"ElecStateLabel", "B"};
    if (val == "88") return {"ElecStateLabel", "X"};
    if (val == "97") return {"ElecStateLabel", "a"};
    if (val == "98") return {"ElecStateLabel", "b"};
    if (val == "99") return {"ElecStateLabel", "c"};
    if (val == "65") return {"ElecStateLabel", "A"};
    if (val == "39") return {"ElecStateLabel", "'"};
    if (val == "66") return {"ElecStateLabel", "B"};
    goto error;
  }

  if (key == "parity") {
    if (val == "-") return {key, "-"};
    if (val == "+") return {key, "+"};
    if (val == "-1") return {key, "-"};
    if (val == "1") return {key, "+"};
    goto error;
  }

  if (key == "Hund") return {"config", String{val}};

  if (key == "kronigParity") {
    if (val == "e") return {key, "e"};
    if (val == "f") return {key, "f"};
    if (val == "101") return {key, "e"};
    if (val == "102") return {key, "f"};
    goto error;
  }

  return {key, String{val}};

error:
  ARTS_USER_ERROR("Cannot read combination ", key, ' ', val)
}

ValueList::ValueList(std::string_view s, bool legacy) : values(0) {
  const Index n = count_items(s);

  if (not legacy) {
    ARTS_USER_ERROR_IF(
        n % 3, "Must have multiple of three items, got ", n, " in:\n", s)
    for (Index i = 0; i < n; i += 3) values.emplace_back(Value(items(s, i, 3)));
  } else {
    ARTS_USER_ERROR_IF(n < 2, "Must have two items:\n", s)
    Index i       = 0;
    auto key_type = items(s, i);
    if (key_type == "ALL") return;
    if (key_type == "NONE") return;

    if (key_type == "TR") {
      // Transition type, will look like SPEC TR UP QN VAL QN VAL QN VAL LO QN VAL QN VAL QN VAL
      ARTS_USER_ERROR_IF(n < 4, "Must have at least four items")

      i++;
      ARTS_USER_ERROR_IF(
          items(s, i) not_eq "UP", "Bad legacy quantum numbers in:\n", s)
      i++;

      bool upp = true;
      for (; i < n; i += 2) {
        auto key = items(s, i);
        if (key == "LO") {
          i++;
          key = items(s, i);
          upp = false;
        }

        auto [k, val] = fix_legacy(key, items(s, i + 1));

        auto t = to<QuantumNumberType>(k);
        if (has(t)) {
          std::find_if(begin(), end(), [t](auto& x) {
            return x.type == t;
          })->set(val, upp);
        } else {
          Value value = add(t);
          value.set(val, upp);
          set(value);
        }
      }
      return;
    }

    if (key_type == "EN") {
      // Transition type, will look like SPEC EN QN VAL QN VAL QN VAL
      i++;

      for (; i < n; i += 2) {
        auto [key, val] = fix_legacy(items(s, i), items(s, i + 1));
        auto t          = to<QuantumNumberType>(key);

        if (has(t)) {
          auto valptr = std::find_if(
              begin(), end(), [t](auto& x) { return x.type == t; });
          valptr->set(val, true);
          valptr->set(val, false);
        } else {
          Value& value = add(t);
          value.set(val, true);
          value.set(val, false);
        }
      }
      return;
    }
  }
  finalize();
}

ValueList::ValueList(std::string_view upp, std::string_view low) {
  const Index nu = count_items(upp);
  const Index nl = count_items(low);

  ARTS_USER_ERROR_IF(
      nu % 2,
      "Uneven count of items for legacy upper quantum number list: ",
      upp)
  ARTS_USER_ERROR_IF(
      nl % 2,
      "Uneven count of items for legacy lower quantum number list: ",
      low)

  for (Index i = 0; i < nu; i += 2) {
    auto [k, val] = fix_legacy(items(upp, i), items(upp, i + 1));
    auto key      = to<QuantumNumberType>(k);

    auto ptr =
        std::find_if(begin(), end(), [key](auto& a) { return a.type == key; });
    if (ptr == end()) {
      add(key);
      ptr = std::find_if(
          begin(), end(), [key](auto& a) { return a.type == key; });
    }

    ptr->set(val, true);
  }

  for (Index i = 0; i < nl; i += 2) {
    auto [k, val] = fix_legacy(items(low, i), items(low, i + 1));
    auto key      = to<QuantumNumberType>(k);

    auto ptr =
        std::find_if(begin(), end(), [key](auto& a) { return a.type == key; });
    if (ptr == end()) {
      add(key);
      ptr = std::find_if(
          begin(), end(), [key](auto& a) { return a.type == key; });
    }

    ptr->set(val, false);
  }
}

/** Returns some input "ASDASDS=asdAS" as ["ASDASDS", "asdAS"] for Hitran online data
 * 
 * Note that there is a special exception for F#XYZ values
 *
 * @param x A string
 * @return constexpr std::pair<std::string_view, std::string_view> 
 */
std::pair<std::string_view, std::string_view> split_hitran_qn(
    std::string_view x) {
  auto eq = x.find('=');
  std::pair<std::string_view, std::string_view> out{strip(x.substr(0, eq)),
                                                    strip(x.substr(eq + 1))};

  if (x.size() > 1 and 'F' == out.first.front() and '#' == out.first[1])
    out.first = out.first.substr(0, 1);

  return out;
}

ValueList from_hitran(std::string_view upp, std::string_view low) {
  ValueList out;

  upp = strip(upp);
  while (not upp.empty()) {
    auto sep    = upp.find(';');
    auto [t, v] = split_hitran_qn(upp.substr(0, sep));
    auto type   = to<QuantumNumberType>(t);
    ARTS_USER_ERROR_IF(
        out.has(type),
        "Type ",
        t,
        " already exist, this is a problem, there should be only one per level!")
    out.add(type).set(v, true);

    if (sep == upp.npos) break;
    upp = upp.substr(sep + 1);
  }

  low = strip(low);
  while (not low.empty()) {
    auto sep    = low.find(';');
    auto [t, v] = split_hitran_qn(low.substr(0, sep));
    auto type   = to<QuantumNumberType>(t);
    if (out.has(type)) {
      Value val = out[type];
      val.set(v, false);
      out.set(val);
    } else {
      out.add(type).set(v, false);
    }

    if (sep == low.npos) break;
    low = low.substr(sep + 1);
  }

  return out;
}

void ValueList::finalize() {
  sort_by_type();
  ARTS_USER_ERROR_IF(not has_unique_increasing_types(),
                     "The quantum number list: [",
                     *this,
                     "] contains copies of types")
}

bool ValueList::perpendicular(const ValueList& that) const ARTS_NOEXCEPT {
  ARTS_ASSERT(has_unique_increasing_types())
  ARTS_ASSERT(that.has_unique_increasing_types())

  auto this_val = cbegin();
  auto that_val = that.cbegin();
  while (that_val not_eq that.cend() and this_val not_eq cend()) {
    if (that_val->type < this_val->type)
      that_val++;
    else if (this_val->type < that_val->type)
      this_val++;
    else
      return false;
  }
  return true;
}

CheckMatch ValueList::check_match(const ValueList& other) const ARTS_NOEXCEPT {
  CheckMatch status = {CheckValue::Full, CheckValue::Full};

  for (const QuantumNumberType t : enumtyps::QuantumNumberTypeTypes) {
    const bool ahas = has(t), bhas = other.has(t);

    if (ahas and bhas) {
      const LevelMatch levels = operator[](t).level_match(other[t]);
      if (not levels.upp) status.upp = CheckValue::Miss;
      if (not levels.low) status.low = CheckValue::Miss;
    } else if (ahas and not bhas) {
      status = update(status, CheckValue::BinA);
    } else if (not ahas and bhas) {
      status = update(status, CheckValue::AinB);
    }
  }
  return status;
}

const Value& ValueList::operator[](QuantumNumberType t) const ARTS_NOEXCEPT {
  auto val =
      std::find_if(cbegin(), cend(), [t](auto& x) { return x.type == t; });
  ARTS_ASSERT(val not_eq cend())
  return *val;
}

Value& ValueList::add(QuantumNumberType t) {
  Value v{t};
  values.push_back(v);
  finalize();

  // We have the value, it is unique, so no need to check this pointer (we still need to find it, since we sort things)
  return *std::find_if(begin(), end(), [t](auto& x) { return x.type == t; });
}

Value& ValueList::add(Value v) {
  values.push_back(v);
  finalize();

  // We have the value, it is unique, so no need to check this pointer (we still need to find it, since we sort things)
  return *std::find_if(
      begin(), end(), [t = v.type](auto& x) { return x.type == t; });
}

void ValueList::set(Value v) {
  if (has(v.type))
    *std::find_if(
        begin(), end(), [t = v.type](auto& x) { return x.type == t; }) = v;
  else
    add(v);
}

void ValueList::set(Index i, std::string_view upp, std::string_view low) {
  values[i] = Value(var_string(values[i].type, ' ', upp, ' ', low));
}

std::ostream& operator<<(std::ostream& os, const ValueList& vl) {
  for (Size i = 0; i < vl.values.size(); i++) {
    if (i) os << ' ';
    os << vl.values[i];
  }
  return os;
}

std::istream& operator>>(std::istream& is, ValueList& vl) {
  for (auto& x : vl.values) is >> x;
  vl.finalize();
  return is;
}

String LocalState::keys() const {
  std::ostringstream os;

  bool has = false;
  for (auto& qn : val) {
    if (has) os << ' ';
    has = true;
    os << qn.type;
  }
  return os.str();
}

String LocalState::values() const {
  std::ostringstream os;

  bool first = true;
  for (auto& x : val) {
    if (first)
      first = false;
    else
      os << ' ';
    os << x.str_low();
  }

  for (auto& x : val) {
    os << ' ' << x.str_upp();
  }

  return os.str();
}

SpeciesIsotope GlobalState::Isotopologue() const noexcept {
  return isotopologue_index < 0 ? SpeciesIsotope()
                                : Species::Isotopologues[isotopologue_index];
}

SpeciesEnum GlobalState::Species() const noexcept {
  return Isotopologue().spec;
}

std::strong_ordering ValueList::operator<=>(const ValueList& v) const {
  const std::size_t n = std::min(values.size(), v.values.size());

  if (auto f = values.size() <=> v.values.size();
      f == std::strong_ordering::equal) {
    return f;
  }

  for (std::size_t i = 0; i < n; i++) {
    if (auto q = values[i] <=> v.values[i]; q != std::strong_ordering::equal) {
      return q;
    }
  }

  return std::strong_ordering::equal;
}

std::ostream& operator<<(std::ostream& os, const GlobalState& gs) {
  os << gs.Isotopologue().FullName();
  for (auto& x : gs.val) {
    os << ' ' << x;
  }
  return os;
}

std::istream& operator>>(std::istream& is, GlobalState& gs) {
  String spec;
  is >> spec >> gs.val;
  gs.isotopologue_index = Species::find_species_index(spec);
  ARTS_USER_ERROR_IF(gs.isotopologue_index < 0, "Cannot read species: ", spec)
  return is;
}

bool vamdcCheck(const ValueList& l, VAMDC type) ARTS_NOEXCEPT {
  switch (type) {
    case VAMDC::asymcs:
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::asSym)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::kronigParity)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::asymos:
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::asSym)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::kronigParity)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::dcs:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F2)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v1)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v2)) return false;
      if (l.has(QuantumNumberType::v3)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::hunda:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F2)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v1)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v2)) return false;
      if (l.has(QuantumNumberType::v3)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::hundb:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F2)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v1)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v2)) return false;
      if (l.has(QuantumNumberType::v3)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::lpcs:
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      return true;
    case VAMDC::lpos:
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::ltcs:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::ltos:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::nltcs:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::nltos:
      if (l.has(QuantumNumberType::F10)) return false;
      if (l.has(QuantumNumberType::F11)) return false;
      if (l.has(QuantumNumberType::F12)) return false;
      if (l.has(QuantumNumberType::F3)) return false;
      if (l.has(QuantumNumberType::F4)) return false;
      if (l.has(QuantumNumberType::F5)) return false;
      if (l.has(QuantumNumberType::F6)) return false;
      if (l.has(QuantumNumberType::F7)) return false;
      if (l.has(QuantumNumberType::F8)) return false;
      if (l.has(QuantumNumberType::F9)) return false;
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::l1)) return false;
      if (l.has(QuantumNumberType::l10)) return false;
      if (l.has(QuantumNumberType::l11)) return false;
      if (l.has(QuantumNumberType::l12)) return false;
      if (l.has(QuantumNumberType::l2)) return false;
      if (l.has(QuantumNumberType::l3)) return false;
      if (l.has(QuantumNumberType::l4)) return false;
      if (l.has(QuantumNumberType::l5)) return false;
      if (l.has(QuantumNumberType::l6)) return false;
      if (l.has(QuantumNumberType::l7)) return false;
      if (l.has(QuantumNumberType::l8)) return false;
      if (l.has(QuantumNumberType::l9)) return false;
      if (l.has(QuantumNumberType::rotSym)) return false;
      if (l.has(QuantumNumberType::rovibSym)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::v10)) return false;
      if (l.has(QuantumNumberType::v11)) return false;
      if (l.has(QuantumNumberType::v12)) return false;
      if (l.has(QuantumNumberType::v4)) return false;
      if (l.has(QuantumNumberType::v5)) return false;
      if (l.has(QuantumNumberType::v6)) return false;
      if (l.has(QuantumNumberType::v7)) return false;
      if (l.has(QuantumNumberType::v8)) return false;
      if (l.has(QuantumNumberType::v9)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      if (l.has(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::sphcs:
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::asSym)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::kronigParity)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::sphos:
      if (l.has(QuantumNumberType::K)) return false;
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::asSym)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::kronigParity)) return false;
      if (l.has(QuantumNumberType::l)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::vibInv)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::stcs:
      if (l.has(QuantumNumberType::Ka)) return false;
      if (l.has(QuantumNumberType::Kc)) return false;
      if (l.has(QuantumNumberType::Lambda)) return false;
      if (l.has(QuantumNumberType::N)) return false;
      if (l.has(QuantumNumberType::Omega)) return false;
      if (l.has(QuantumNumberType::S)) return false;
      if (l.has(QuantumNumberType::Sigma)) return false;
      if (l.has(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.has(QuantumNumberType::asSym)) return false;
      if (l.has(QuantumNumberType::elecInv)) return false;
      if (l.has(QuantumNumberType::elecRefl)) return false;
      if (l.has(QuantumNumberType::elecSym)) return false;
      if (l.has(QuantumNumberType::kronigParity)) return false;
      if (l.has(QuantumNumberType::sym)) return false;
      if (l.has(QuantumNumberType::v)) return false;
      if (l.has(QuantumNumberType::vibRefl)) return false;
      return true;
  }
  return false;
}

StateMatch::StateMatch(const GlobalState& target,
                       const LocalState& local,
                       const GlobalState& global) {
  if (target.isotopologue_index == global.isotopologue_index)
    type = StateMatchType::Isotopologue;
  else if (target.Species() == global.Species())
    type = StateMatchType::Species;

  if (type == StateMatchType::Isotopologue) {
    auto g = target.val.check_match(global.val);
    auto l = target.val.check_match(local.val);

    bool ug = g.upp == CheckValue::Full or g.upp == CheckValue::BinA;
    bool ul = l.upp == CheckValue::Full or l.upp == CheckValue::BinA;
    bool lg = g.low == CheckValue::Full or g.low == CheckValue::BinA;
    bool ll = l.low == CheckValue::Full or l.low == CheckValue::BinA;

    upp = ug and ul;
    low = lg and ll;

    if (upp and low)
      type = StateMatchType::Full;
    else if (upp or low)
      type = StateMatchType::Level;
  }
}

StateMatch::StateMatch(const GlobalState& target, const GlobalState& key) {
  if (target.isotopologue_index == key.isotopologue_index)
    type = StateMatchType::Isotopologue;
  else if (target.Species() == key.Species())
    type = StateMatchType::Species;

  if (type == StateMatchType::Isotopologue) {
    auto m = target.val.check_match(key.val);

    upp = m.upp == CheckValue::Full or m.upp == CheckValue::BinA;
    low = m.low == CheckValue::Full or m.low == CheckValue::BinA;
    if (upp and low)
      type = StateMatchType::Full;
    else if (upp or low)
      type = StateMatchType::Level;
  }
}

bool GlobalState::part_of(const GlobalState& other) const {
  if (other.isotopologue_index not_eq isotopologue_index) return false;

  auto test = other.val.check_match(val);
  return (test.upp == CheckValue::Full or test.upp == CheckValue::AinB) and
         (test.low == CheckValue::Full or test.low == CheckValue::AinB);
}

bool GlobalState::may_be(const GlobalState& other) const {
  return other.part_of(*this);
}

std::ostream& operator<<(std::ostream& os, const LocalState& vl) {
  return os << vl.values();
}

std::istream& operator>>(std::istream& is, LocalState& vl) {
  String val;
  for (auto& v : vl.val) {
    is >> val;
    v.set(val, false);
  }
  for (auto& v : vl.val) {
    is >> val;
    v.set(val, true);
  }
  return is;
}

bool LocalState::same_types_as(const LocalState& that) const {
  return std::equal(val.begin(),
                    val.end(),
                    that.val.begin(),
                    that.val.end(),
                    [](auto& a, auto& b) { return a.type == b.type; });
}

GlobalState GlobalState::LowerLevel() const {
  GlobalState out = *this;
  for (auto& value : out.val) value.qn.upp = value.qn.low;
  return out;
}

GlobalState GlobalState::UpperLevel() const {
  GlobalState out = *this;
  for (auto& value : out.val) value.qn.low = value.qn.upp;
  return out;
}

void Value::swap_values(Value& x) {
  // Make copies
  Value _this = *this;
  Value _x    = x;

  // Copy values from the other part
  _this.qn = x.qn;
  _x.qn    = qn;

  // Assign by reinterpreting the data using standard operations
  *this = Value(var_string(_this));
  x     = Value(var_string(_x));
}

GlobalState::GlobalState(std::string_view s, Index v) {
  auto n             = count_items(s);
  auto specname      = items(s, 0);
  isotopologue_index = Species::find_species_index(specname);

  if (isotopologue_index < 0)
    throw std::runtime_error("Invalid isotopologue: " + std::string(specname) +
                             " from " + std::string(s));
  if (Species::Isotopologues[isotopologue_index].joker() or
      Species::is_predefined_model(
          Species::Isotopologues[isotopologue_index])) {
    throw std::runtime_error("Expects valid standard isotopologue, got: " +
                             std::string(specname) + " from " + std::string(s));
  }

  if (version == v) {
    if (n > 1) val = ValueList(s.substr(specname.length() + 1));
  } else if (v == 0 or v == 1) {
    val = ValueList(s.substr(specname.length() + 1), true);
  } else {
    ARTS_USER_ERROR("Unknown version: ", v)
  }
}

void ValueList::add_type_wo_sort(QuantumNumberType t) {
  values.emplace_back().type = t;
}

void LocalState::set_unsorted_qns(const Array<QuantumNumberType>& vals) {
  val = ValueList("");
  for (auto qn : vals) val.add_type_wo_sort(qn);
}

[[nodiscard]] LevelTest GlobalState::part_of(const GlobalState& g,
                                             const LocalState& l) const {
  bool upp = true, low = true;
  if (isotopologue_index not_eq g.isotopologue_index) return {false, false};

  for (auto& qn : val) {
    bool any = false;

    if (g.val.has(qn.type)) {
      auto& v = g.val[qn.type];

      upp = upp and v.str_upp() == qn.str_upp();
      low = low and v.str_low() == qn.str_low();
      any = true;
    }

    if (l.val.has(qn.type)) {
      auto& v = l.val[qn.type];

      upp = upp and v.str_upp() == qn.str_upp();
      low = low and v.str_low() == qn.str_low();
      ARTS_USER_ERROR_IF(any,
                         "Repeated type: ",
                         qn.type,
                         '\n',
                         "From original key: ",
                         *this,
                         '\n',
                         "With global key: ",
                         g,
                         '\n',
                         "With local key: ",
                         l.val)
      any = true;
    }

    if (not any) return {false, false};
  }
  return {upp, low};
}

bofstream& Quantum::Number::Value::write(bofstream& bof) const {
  switch (common_value_type(type)) {
    case QuantumNumberValueType::I:
      bof << qn.upp.i.x << qn.low.i.x;
      break;
    case QuantumNumberValueType::H:
      bof << qn.upp.h.x << qn.low.h.x;
      break;
    case QuantumNumberValueType::S:
      bof.writeString(qn.upp.s.x.data(), StringValue::N);
      bof.writeString(qn.low.s.x.data(), StringValue::N);
      break;
  }
  return bof;
}

bifstream& Quantum::Number::Value::read(bifstream& bif) {
  switch (common_value_type(type)) {
    case QuantumNumberValueType::I:
      bif >> qn.upp.i.x >> qn.low.i.x;
      break;
    case QuantumNumberValueType::H:
      bif >> qn.upp.h.x >> qn.low.h.x;
      break;
    case QuantumNumberValueType::S:
      bif.readString(qn.upp.s.x.data(), StringValue::N);
      bif.readString(qn.low.s.x.data(), StringValue::N);
      break;
  }
  return bif;
}

bool Quantum::Number::ValueList::good() const {
  {
    for (auto& qn : values)
      if (not qn.good()) {
        return false;
      }
  }
  return true;
}

bool Quantum::Number::LocalState::good() const { return val.good(); }

bool Quantum::Number::GlobalState::good() const {
  return Species::is_normal_isotopologue(Isotopologue()) and val.good();
}

std::ostream& operator<<(std::ostream& os, const Array<GlobalState>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

void ValueList::reserve(Size n) { values.reserve(n); }
}  // namespace Quantum::Number

std::ostream& operator<<(std::ostream& os, const Array<QuantumNumberType>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
