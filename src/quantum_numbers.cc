#include "quantum_numbers.h"

#include <algorithm>
#include <string_view>
#include <utility>

#include "debug.h"
#include "species.h"

namespace Quantum::Number {
std::ostream& operator<<(std::ostream& os, ValueDescription x) {
  os << x.type << ' ';
  switch (x.type) {
    case ValueType::S:
      return os << x.val.s.val();
    case ValueType::I:
      return os << x.val.i.x;
    case ValueType::H:
      return os << Rational(x.val.h.x, 2);
    case ValueType::FINAL: {
    }
  }
  return os;
}

String Value::str_upp() const noexcept {
  if (ValueType::S == common_value_type(type)) return qn.upp.s.val();
  return var_string(upp());
}
String Value::str_low() const noexcept {
  if (ValueType::S == common_value_type(type)) return qn.low.s.val();
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

  if (key == "Hund") return {"config", val};

  if (key == "kronigParity") {
    if (val == "e") return {key, "e"};
    if (val == "f") return {key, "f"};
    if (val == "101") return {key, "e"};
    if (val == "102") return {key, "f"};
    goto error;
  }

  return {key, val};

error:
  ARTS_USER_ERROR("Cannot read combination ", key, ' ', val)
}

ValueList::ValueList(std::string_view s, bool legacy) : values(0) {
  const Index n = count_items(s);

  if (not legacy) {
    ARTS_USER_ERROR_IF(
        n % 3, "Must have multiple of three items, got ", n, " in:\n", s)
    for (Index i = 0; i < n; i += 3) values.emplace_back(Value(items<3>(s, i)));
  } else {
    ARTS_USER_ERROR_IF(n < 2, "Must have two items:\n", s)
    Index i = 0;
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

        auto t = toTypeOrThrow(k);
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
        auto t = toTypeOrThrow(key);

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
    auto key = toTypeOrThrow(k);

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
    auto key = toTypeOrThrow(k);

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
constexpr std::pair<std::string_view, std::string_view> split_hitran_qn(
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
    auto sep = upp.find(';');
    auto [t, v] = split_hitran_qn(upp.substr(0, sep));
    auto type = toTypeOrThrow(t);
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
    auto sep = low.find(';');
    auto [t, v] = split_hitran_qn(low.substr(0, sep));
    auto type = toTypeOrThrow(t);
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

  for (Type const t : enumtyps::TypeTypes) {
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

const Value& ValueList::operator[](Type t) const ARTS_NOEXCEPT {
  auto val =
      std::find_if(cbegin(), cend(), [t](auto& x) { return x.type == t; });
  ARTS_ASSERT(val not_eq cend())
  return *val;
}

Value& ValueList::add(Type t) {
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
  for (Index i = 0; i < vl.values.nelem(); i++) {
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

Species::IsotopeRecord GlobalState::Isotopologue() const noexcept {
  return isotopologue_index < 0 ? Species::IsotopeRecord() : Species::Isotopologues[isotopologue_index];
}

Species::Species GlobalState::Species() const noexcept {
  return Isotopologue().spec;
}

std::ostream& operator<<(std::ostream& os, const GlobalState& gs) {
  os << gs.Isotopologue().FullName();
  for (auto& x: gs.val) {
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
      if (l.has(Type::K)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::asSym)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::kronigParity)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::vibRefl)) return false;
      return true;
    case VAMDC::asymos:
      if (l.has(Type::K)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::asSym)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::kronigParity)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::vibRefl)) return false;
      return true;
    case VAMDC::dcs:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F2)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v1)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v2)) return false;
      if (l.has(Type::v3)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::hunda:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F2)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v1)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v2)) return false;
      if (l.has(Type::v3)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::hundb:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F2)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v1)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v2)) return false;
      if (l.has(Type::v3)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::lpcs:
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      return true;
    case VAMDC::lpos:
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::ltcs:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::ltos:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::nltcs:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::nltos:
      if (l.has(Type::F10)) return false;
      if (l.has(Type::F11)) return false;
      if (l.has(Type::F12)) return false;
      if (l.has(Type::F3)) return false;
      if (l.has(Type::F4)) return false;
      if (l.has(Type::F5)) return false;
      if (l.has(Type::F6)) return false;
      if (l.has(Type::F7)) return false;
      if (l.has(Type::F8)) return false;
      if (l.has(Type::F9)) return false;
      if (l.has(Type::K)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::l1)) return false;
      if (l.has(Type::l10)) return false;
      if (l.has(Type::l11)) return false;
      if (l.has(Type::l12)) return false;
      if (l.has(Type::l2)) return false;
      if (l.has(Type::l3)) return false;
      if (l.has(Type::l4)) return false;
      if (l.has(Type::l5)) return false;
      if (l.has(Type::l6)) return false;
      if (l.has(Type::l7)) return false;
      if (l.has(Type::l8)) return false;
      if (l.has(Type::l9)) return false;
      if (l.has(Type::rotSym)) return false;
      if (l.has(Type::rovibSym)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::v10)) return false;
      if (l.has(Type::v11)) return false;
      if (l.has(Type::v12)) return false;
      if (l.has(Type::v4)) return false;
      if (l.has(Type::v5)) return false;
      if (l.has(Type::v6)) return false;
      if (l.has(Type::v7)) return false;
      if (l.has(Type::v8)) return false;
      if (l.has(Type::v9)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      if (l.has(Type::vibSym)) return false;
      return true;
    case VAMDC::sphcs:
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::asSym)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::kronigParity)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      return true;
    case VAMDC::sphos:
      if (l.has(Type::K)) return false;
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::asSym)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::kronigParity)) return false;
      if (l.has(Type::l)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::vibInv)) return false;
      if (l.has(Type::vibRefl)) return false;
      return true;
    case VAMDC::stcs:
      if (l.has(Type::Ka)) return false;
      if (l.has(Type::Kc)) return false;
      if (l.has(Type::Lambda)) return false;
      if (l.has(Type::N)) return false;
      if (l.has(Type::Omega)) return false;
      if (l.has(Type::S)) return false;
      if (l.has(Type::Sigma)) return false;
      if (l.has(Type::SpinComponentLabel)) return false;
      if (l.has(Type::asSym)) return false;
      if (l.has(Type::elecInv)) return false;
      if (l.has(Type::elecRefl)) return false;
      if (l.has(Type::elecSym)) return false;
      if (l.has(Type::kronigParity)) return false;
      if (l.has(Type::sym)) return false;
      if (l.has(Type::v)) return false;
      if (l.has(Type::vibRefl)) return false;
      return true;
    case VAMDC::FINAL: {
    }
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
  Value _x = x;

  // Copy values from the other part
  _this.qn = x.qn;
  _x.qn = qn;

  // Assign by reinterpreting the data using standard operations
  *this = Value(var_string(_this));
  x = Value(var_string(_x));
}

GlobalState::GlobalState(std::string_view s, Index v) {
  auto n = count_items(s);
  auto specname = items(s, 0);
  isotopologue_index = Species::find_species_index(specname);
  ARTS_USER_ERROR_IF(isotopologue_index < 0, "Bad species in: ", s)

  if (version == v) {
    if (n > 1) val = ValueList(s.substr(specname.length() + 1));
  } else if (v == 0 or v == 1) {
    val = ValueList(s.substr(specname.length() + 1), true);
  } else {
    ARTS_USER_ERROR("Unknown version: ", v)
  }
}

void ValueList::add_type_wo_sort(Type t) { values.emplace_back().type = t; }

void LocalState::set_unsorted_qns(const Array<Type>& vals) {
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
    case ValueType::I:
      bof << qn.upp.i.x << qn.low.i.x;
      break;
    case ValueType::H:
      bof << qn.upp.h.x << qn.low.h.x;
      break;
    case ValueType::S:
      bof.writeString(qn.upp.s.x.data(), StringValue::N);
      bof.writeString(qn.low.s.x.data(), StringValue::N);
      break;
    case ValueType::FINAL: {
    }
  }
  return bof;
}

bifstream& Quantum::Number::Value::read(bifstream& bif) {
  switch (common_value_type(type)) {
    case ValueType::I:
      bif >> qn.upp.i.x >> qn.low.i.x;
      break;
    case ValueType::H:
      bif >> qn.upp.h.x >> qn.low.h.x;
      break;
    case ValueType::S:
      bif.readString(qn.upp.s.x.data(), StringValue::N);
      bif.readString(qn.low.s.x.data(), StringValue::N);
      break;
    case ValueType::FINAL: {
    }
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
}  // namespace Quantum::Number
