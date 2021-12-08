#include "quantum_numbers.h"
#include "debug.h"
#include <utility>

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
  if (ValueType::S == value_type) return qn.upp.s.val();
  return var_string(upp());
}
String Value::str_low() const noexcept {
  if (ValueType::S == value_type) return qn.low.s.val();
  return var_string(low());
}

std::ostream& operator<<(std::ostream& os, Value x) {
  return os << x.type << ' ' << x.str_upp() << ' ' << x.str_low();
}

std::istream& operator>>(std::istream& is, Value& x) {
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

void ValueList::sort_by_type() {
  std::sort(values.begin(), values.end(), [](auto& a, auto& b) {
    return a.type < b.type;
  });
}

bool ValueList::has_unique_increasing_types() const {
  return std::adjacent_find(values.begin(), values.end(), [](auto& a, auto& b) {
           return a.type >= b.type;
         }) == values.end();
}

ValueList::ValueList(std::string_view s) : values(0) {
  Index n = count_items(s);
  ARTS_USER_ERROR_IF(
      n % 3, "Must have multiple of three items, got ", n, " in:\n", s)
  for (Index i = 0; i < n; i += 3) values.emplace_back(Value(items<3>(s, i)));
  sort_by_type();
  ARTS_USER_ERROR_IF(not has_unique_increasing_types(),
                     "Bad Value List (contains copies or cannot sort):\n",
                     *this)
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

CheckMatch ValueList::operator==(const ValueList& other) const noexcept {
  CheckMatch status = {CheckValue::Full, CheckValue::Full};

  for (Type t : enumtyps::TypeTypes) {
    const bool ahas = has(t), bhas = other.has(t);

    if (ahas and bhas) {
      const LevelMatch levels = operator[](t) == other[t];
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
  auto val = std::find_if(
      values.begin(), values.end(), [t](auto& x) { return x.type == t; });
  ARTS_ASSERT(val not_eq values.end())
  return *val;
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
  vl.sort_by_type();
  ARTS_USER_ERROR_IF(not vl.has_unique_increasing_types(),
                     "Bad Value List (contains copies or cannot sort):\n",
                     vl)
  return is;
}

String LocalState::keys() const {
  std::ostringstream os;

  bool has = false;
  for (Type qn : enumtyps::TypeTypes) {
    if (val.has(qn)) {
      if (has) os << ' ';
      has = true;
      os << qn;
    }
  }
  return os.str();
}

String LocalState::values() const {
  std::ostringstream os;

  bool has = false;
  for (Type qn : enumtyps::TypeTypes) {
    if (val.has(qn)) {
      if (has) os << ' ';
      has = true;
      auto& v = val[qn];
      os << v.str_upp() << ' ' << v.str_low();
    }
  }
  return os.str();
}

StateMatch checkLocalGlobal(const GlobalState& target,
                            const LocalState& local,
                            const GlobalState& global,
                            bool level) ARTS_NOEXCEPT {
  if (target.Isotopologue() == global.Isotopologue()) {}
  else if (target.Isotopologue().spec == global.Isotopologue().spec) return StateMatch::Species;
  else return StateMatch::None;

  ARTS_ASSERT(
      local.val.perpendicular(global.val),
      "The local and global quantum number lists must be independent\n\n",
      "Global: ",
      global.val,
      "\nLocal: ",
      local.val)

  auto global_match = target.val == global.val;
  auto local_match = target.val == local.val;

  bool lowl = local.val.nelem() and (local_match.low == CheckValue::Full or local_match.low == CheckValue::BinA);
  bool uppl = local.val.nelem() and (local_match.upp == CheckValue::Full or local_match.upp == CheckValue::BinA);
  bool lowg = (global_match.low == CheckValue::Full or global_match.low == CheckValue::BinA);
  bool uppg = (global_match.upp == CheckValue::Full or global_match.upp == CheckValue::BinA);

  if (lowl and uppl and lowg and uppg) return StateMatch::FullMatch;
  if (level and lowl and lowg) return StateMatch::FullLower;
  if (level and uppl and uppg) return StateMatch::FullUpper;
  if (level and lowg) return StateMatch::BandLower;
  if (level and uppg) return StateMatch::BandUpper;
  if (uppg and lowg) return StateMatch::BandMatch;
  return StateMatch::Isotopologue;
}

const Species::IsotopeRecord& GlobalState::Isotopologue() const noexcept {
  return Species::Isotopologues[isotopologue_index];
}

Species::Species GlobalState::Species() const noexcept {return Isotopologue().spec;}

std::ostream& operator<<(std::ostream& os, const GlobalState& gs) {
  return os << gs.Isotopologue().FullName() << ' ' << gs.val;
}

std::istream& operator<<(std::istream& is, GlobalState& gs) {
  String spec;
  is >> spec >> gs.val;
  gs.isotopologue_index = Species::find_species_index(spec);
  ARTS_USER_ERROR_IF(gs.isotopologue_index < 0, "Cannot read species: ", spec)
  return is;
}
}  // namespace Quantum::Number
