#include "quantum.h"

#include <debug.h>
#include <enumsQuantumNumberType.h>
#include <isotopologues.h>
#include <nonstd.h>
#include <rational.h>

#include <compare>
#include <exception>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <variant>

#include "quantum_term_symbol.h"

namespace Quantum {
using enum QuantumNumberType;

Value::Value(String x) : value(std::move(x)) {}

Value::Value(Index x) : value(Rational(x, 1)) {}

Value::Value(Rational x) : value(x) {}

Value::Value() : Value(Rational{0}) {}

Value::Value(QuantumNumberType type) {
  switch (type) {
    case alpha:              value = String{"*"}; return;
    case config:             value = String{"*"}; return;
    case ElecStateLabel:     value = String{"*"}; return;
    case F:                  value = Rational{0}; return;
    case F1:                 value = Rational{0}; return;
    case F10:                value = Rational{0}; return;
    case F11:                value = Rational{0}; return;
    case F12:                value = Rational{0}; return;
    case F2:                 value = Rational{0}; return;
    case F3:                 value = Rational{0}; return;
    case F4:                 value = Rational{0}; return;
    case F5:                 value = Rational{0}; return;
    case F6:                 value = Rational{0}; return;
    case F7:                 value = Rational{0}; return;
    case F8:                 value = Rational{0}; return;
    case F9:                 value = Rational{0}; return;
    case I:                  value = Rational{0}; return;
    case J:                  value = Rational{0}; return;
    case K:                  value = Rational{0}; return;
    case Ka:                 value = Rational{0}; return;
    case Kc:                 value = Rational{0}; return;
    case L:                  value = Rational{0}; return;
    case Lambda:             value = Rational{0}; return;
    case N:                  value = Rational{0}; return;
    case Omega:              value = Rational{0}; return;
    case S:                  value = Rational{0}; return;
    case Sigma:              value = Rational{0}; return;
    case SpinComponentLabel: value = Rational{0}; return;
    case asSym:              value = String{"*"}; return;
    case elecInv:            value = String{"*"}; return;
    case elecRefl:           value = String{"*"}; return;
    case elecSym:            value = String{"*"}; return;
    case kronigParity:       value = String{"*"}; return;
    case l:                  value = Rational{0}; return;
    case l1:                 value = Rational{0}; return;
    case l10:                value = Rational{0}; return;
    case l11:                value = Rational{0}; return;
    case l12:                value = Rational{0}; return;
    case l2:                 value = Rational{0}; return;
    case l3:                 value = Rational{0}; return;
    case l4:                 value = Rational{0}; return;
    case l5:                 value = Rational{0}; return;
    case l6:                 value = Rational{0}; return;
    case l7:                 value = Rational{0}; return;
    case l8:                 value = Rational{0}; return;
    case l9:                 value = Rational{0}; return;
    case n:                  value = String{"*"}; return;
    case parity:             value = String{"*"}; return;
    case r:                  value = Rational{0}; return;
    case rotSym:             value = String{"*"}; return;
    case rovibSym:           value = String{"*"}; return;
    case sym:                value = String{"*"}; return;
    case tau:                value = String{"*"}; return;
    case term:               value = String{"*"}; return;
    case v:                  value = Rational{0}; return;
    case v1:                 value = Rational{0}; return;
    case v10:                value = Rational{0}; return;
    case v11:                value = Rational{0}; return;
    case v12:                value = Rational{0}; return;
    case v2:                 value = Rational{0}; return;
    case v3:                 value = Rational{0}; return;
    case v4:                 value = Rational{0}; return;
    case v5:                 value = Rational{0}; return;
    case v6:                 value = Rational{0}; return;
    case v7:                 value = Rational{0}; return;
    case v8:                 value = Rational{0}; return;
    case v9:                 value = Rational{0}; return;
    case vibInv:             value = String{"*"}; return;
    case vibRefl:            value = String{"*"}; return;
    case vibSym:             value = String{"*"}; return;
  }
  std::unreachable();
}

std::istream& operator>>(std::istream& is, Value& v) {
  return std::visit([&is](auto& x) -> std::istream& { return is >> x; },
                    v.value);
}

Value::operator Rational() const { return get<Rational>(); }

Value::operator String() const { return get<String>(); }

std::istream& operator>>(std::istream& is, UpperLower& ul) {
  return is >> ul.upper >> ul.lower;
}

Level upper_level(const State& state) {
  Level lvl{};
  lvl.reserve(state.size());
  for (auto& [key, ul] : state) lvl.emplace(key, ul.upper);
  return lvl;
}

[[nodiscard]] Level lower_level(const State& state) {
  Level lvl{};
  lvl.reserve(state.size());
  for (auto& [key, ul] : state) lvl.emplace(key, ul.lower);
  return lvl;
}

std::istream& operator>>(std::istream& is, Level& level) {
  level.clear();
  Size n;
  is >> n;

  for (Size i = 0; i < n; i++) {
    QuantumNumberType qn;
    is >> qn;
    auto [it, _] = level.emplace(qn, Value{qn});
    is >> it->second;
  }

  return is;
}

std::istream& operator>>(std::istream& is, State& state) {
  state.clear();
  Size n;
  is >> n;

  for (Size i = 0; i < n; i++) {
    QuantumNumberType qn;
    is >> qn;
    auto [it, _] =
        state.emplace(qn, UpperLower{.upper = Value(qn), .lower = Value(qn)});
    is >> it->second;
  }

  return is;
}

namespace {
Size count(const std::string_view s) {
  if (s.empty()) return 0;

  bool is_space = nonstd::isspace(s.front());
  Size nelem    = not is_space;
  for (Size i = 1; i < s.size(); i++) {
    bool this_space = nonstd::isspace(s[i]);

    nelem += (is_space and not this_space);

    is_space = this_space;
  }

  return nelem;
}

std::string_view next(std::string_view& s) {
  Size len = 0, st = 0;
  while ((st + len) < s.size() and nonstd::isspace(s[st + len])) st++;
  while ((st + len) < s.size() and not nonstd::isspace(s[st + len])) len++;
  std::string_view out = s.substr(st, len);
  s.remove_prefix(st + len);
  return out;
}
}  // namespace

void Value::set(const std::string_view x) {
  std::visit([x]<typename T>(T& v) { v = T{x}; }, value);
}

State state_from(std::string_view s) {
  Size n = count(s);

  if (n % 3) throw std::runtime_error("Bad state, must be from 3N items");

  n /= 3;
  State state{n};

  for (Size i = 0; i < n; i++) {
    QuantumNumberType qn = to<QuantumNumberType>(next(s));
    auto [it, _] =
        state.emplace(qn, UpperLower{.upper = Value(qn), .lower = Value(qn)});
    it->second.upper.set(next(s));
    it->second.lower.set(next(s));
  }

  return state;
}

Identifier::Identifier(std::string_view s_) try {
  std::string_view s = s_;

  isot  = SpeciesIsotope(next(s));
  state = state_from(s);
} catch (std::exception& e) {
  throw std::runtime_error(std::format(
      "Cannot construct QuantumIdentifier from \"{}\", got error:\n{}",
      s_,
      e.what()));
}

Level level_from(std::string_view s) {
  Size n = count(s);

  if (n % 2) throw std::runtime_error("Bad level, must be from 2N items");

  n /= 2;
  Level state{n};

  for (Size i = 0; i < n; i++) {
    QuantumNumberType qn = to<QuantumNumberType>(next(s));
    auto [it, _]         = state.emplace(qn, Value(qn));
    it->second.set(next(s));
  }

  return state;
}

LevelIdentifier::LevelIdentifier(const std::string_view s_) try {
  std::string_view s = s_;
  isot               = SpeciesIsotope(next(s));
  state              = level_from(s);
} catch (std::exception& e) {
  throw std::runtime_error(std::format(
      "Cannot construct QuantumLevelIdentifier from \"{}\", got error:\n{}",
      s_,
      e.what()));
}

Identifier::Identifier(const SpeciesIsotope s) : isot(s), state{} {}

Identifier::Identifier(const SpeciesIsotope s, const State& q)
    : isot(s), state{q} {}

LevelIdentifier::LevelIdentifier(const SpeciesIsotope s) : isot(s), state{} {}

LevelIdentifier::LevelIdentifier(const SpeciesIsotope s, const Level& q)
    : isot(s), state{q} {}

Size Identifier::size() const { return state.size(); }

Size LevelIdentifier::size() const { return state.size(); }

LevelIdentifier Identifier::lower() const { return {isot, lower_level(state)}; }

LevelIdentifier Identifier::upper() const { return {isot, upper_level(state)}; }

bool LevelIdentifier::operator==(const LevelIdentifier& l) const {
  if (l.isot != isot or l.size() != size()) return false;

  return l.state == state;
}

bool LevelIdentifier::operator!=(const LevelIdentifier& l) const {
  return not this->operator==(l);
}

bool Identifier::operator==(const Identifier& l) const {
  if (l.isot != isot or l.size() != size()) return false;

  return l.state == state;
}

bool Identifier::operator!=(const Identifier& l) const {
  return not this->operator==(l);
}

std::strong_ordering Identifier::operator<=>(const Identifier& g) const {
  if (const auto test = isot <=> g.isot; test != std::strong_ordering::equal)
    return test;

  for (auto& qns : enumtyps::QuantumNumberTypeTypes) {
    auto thi = state.find(qns);
    auto tha = g.state.find(qns);

    if (thi == state.end() and tha != g.state.end())
      return std::strong_ordering::greater;
    if (thi != state.end() and tha == g.state.end())
      return std::strong_ordering::less;

    if (std::strong_ordering test = thi->second <=> tha->second;
        test != std::strong_ordering::equal)
      return test;
  }

  return std::strong_ordering::equal;
}

std::strong_ordering LevelIdentifier::operator<=>(
    const LevelIdentifier& g) const {
  if (const auto test = isot <=> g.isot; test != std::strong_ordering::equal)
    return test;

  for (auto& qns : enumtyps::QuantumNumberTypeTypes) {
    auto thi = state.find(qns);
    auto tha = g.state.find(qns);

    if (thi == state.end() and tha != g.state.end())
      return std::strong_ordering::greater;
    if (thi != state.end() and tha == g.state.end())
      return std::strong_ordering::less;

    if (std::strong_ordering test = thi->second <=> tha->second;
        test != std::strong_ordering::equal)
      return test;
  }

  return std::strong_ordering::equal;
}

namespace {
std::string_view rstrip(std::string_view x) {
  while (not x.empty() and nonstd::isspace(x.back())) x.remove_suffix(1);
  return x;
}

std::string_view lstrip(std::string_view x) {
  while (not x.empty() and nonstd::isspace(x.front())) x.remove_prefix(1);
  return x;
}

std::string_view strip(std::string_view x) { return rstrip(lstrip(x)); }
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
}  // namespace

State from_hitran(std::string_view upp, std::string_view low) {
  State out;

  upp = strip(upp);
  while (not upp.empty()) {
    auto sep    = upp.find(';');
    auto [t, v] = split_hitran_qn(upp.substr(0, sep));
    auto type   = to<QuantumNumberType>(t);

    auto [it, _] = out.emplace(
        type, UpperLower{.upper = Value{type}, .lower = Value{type}});

    //! FIMXE: when string view streams are available
    std::istringstream is{std::string{v}};
    is >> it->second.upper;

    if (sep == upp.npos) break;
    upp = upp.substr(sep + 1);
  }

  low = strip(low);
  while (not low.empty()) {
    auto sep    = low.find(';');
    auto [t, v] = split_hitran_qn(low.substr(0, sep));
    auto type   = to<QuantumNumberType>(t);

    //! FIMXE: when string view streams are available
    std::istringstream is{std::string{v}};

    if (out.contains(type)) {
      is >> out.at(type).lower;
    } else {
      auto [it, _] = out.emplace(
          type, UpperLower{.upper = Value{type}, .lower = Value{type}});
      is >> it->second.lower;
    }

    if (sep == low.npos) break;
    low = low.substr(sep + 1);
  }

  return out;
}

bool contains_any_of(const State& state,
                     std::initializer_list<QuantumNumberType> keys) {
  for (auto& key : keys)
    if (state.contains(key)) return true;
  return false;
}

bool contains_all_of(const State& state,
                     std::initializer_list<QuantumNumberType> keys) {
  for (auto& key : keys)
    if (not state.contains(key)) return false;
  return true;
}

bool vamdcCheck(const State& st, VAMDC type) {
  switch (type) {
    case VAMDC::asymcs:
      return not contains_any_of(st,
                                 {K,
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
                                  sym,
                                  v,
                                  vibRefl});
    case VAMDC::asymos:
      return not contains_any_of(st,
                                 {K,
                                  Lambda,
                                  Omega,
                                  Sigma,
                                  SpinComponentLabel,
                                  asSym,
                                  elecRefl,
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
                                  sym,
                                  v,
                                  vibRefl});
    case VAMDC::dcs:
      return not contains_any_of(
          st, {F10,     F11,      F12,      F2,    F3,
               F4,      F5,       F6,       F7,    F8,
               F9,      K,        Ka,       Kc,    Lambda,
               N,       Omega,    S,        Sigma, SpinComponentLabel,
               elecInv, elecRefl, elecSym,  l,     l1,
               l10,     l11,      l12,      l2,    l3,
               l4,      l5,       l6,       l7,    l8,
               l9,      rotSym,   rovibSym, sym,   v1,
               v10,     v11,      v12,      v2,    v3,
               v4,      v5,       v6,       v7,    v8,
               v9,      vibInv,   vibRefl,  vibSym});
    case VAMDC::hunda:
      return not contains_any_of(
          st, {F10,     F11, F12, F2,  F3,  F4,     F5,      F6,
               F7,      F8,  F9,  K,   Ka,  Kc,     N,       SpinComponentLabel,
               elecSym, l,   l1,  l10, l11, l12,    l2,      l3,
               l4,      l5,  l6,  l7,  l8,  l9,     rotSym,  rovibSym,
               sym,     v1,  v10, v11, v12, v2,     v3,      v4,
               v5,      v6,  v7,  v8,  v9,  vibInv, vibRefl, vibSym});
    case VAMDC::hundb:
      return not contains_any_of(
          st,
          {F10,    F11,      F12, F2, F3,    F4,     F5,      F6,    F7, F8,
           F9,     K,        Ka,  Kc, Omega, Sigma,  elecSym, l,     l1, l10,
           l11,    l12,      l2,  l3, l4,    l5,     l6,      l7,    l8, l9,
           rotSym, rovibSym, sym, v1, v10,   v11,    v12,     v2,    v3, v4,
           v5,     v6,       v7,  v8, v9,    vibInv, vibRefl, vibSym});
    case VAMDC::lpcs:
      return not contains_any_of(st,
                                 {K,
                                  Ka,
                                  Kc,
                                  Lambda,
                                  N,
                                  Omega,
                                  S,
                                  Sigma,
                                  SpinComponentLabel,
                                  elecInv,
                                  elecRefl,
                                  elecSym,
                                  rotSym,
                                  rovibSym,
                                  sym,
                                  v});
    case VAMDC::lpos:
      return not contains_any_of(st,
                                 {K,
                                  Ka,
                                  Kc,
                                  Omega,
                                  Sigma,
                                  SpinComponentLabel,
                                  elecSym,
                                  rotSym,
                                  rovibSym,
                                  sym,
                                  v,
                                  vibSym});
    case VAMDC::ltcs:
      return not contains_any_of(st,
                                 {F10,
                                  F11,
                                  F12,
                                  F3,
                                  F4,
                                  F5,
                                  F6,
                                  F7,
                                  F8,
                                  F9,
                                  K,
                                  Ka,
                                  Kc,
                                  Lambda,
                                  N,
                                  Omega,
                                  S,
                                  Sigma,
                                  SpinComponentLabel,
                                  elecInv,
                                  elecRefl,
                                  elecSym,
                                  l,
                                  l1,
                                  l10,
                                  l11,
                                  l12,
                                  l3,
                                  l4,
                                  l5,
                                  l6,
                                  l7,
                                  l8,
                                  l9,
                                  rotSym,
                                  rovibSym,
                                  sym,
                                  v,
                                  v10,
                                  v11,
                                  v12,
                                  v4,
                                  v5,
                                  v6,
                                  v7,
                                  v8,
                                  v9,
                                  vibInv,
                                  vibRefl,
                                  vibSym});
    case VAMDC::ltos:
      return not contains_any_of(st,
                                 {F10,     F11,    F12,      F3,
                                  F4,      F5,     F6,       F7,
                                  F8,      F9,     K,        Ka,
                                  Kc,      Omega,  Sigma,    SpinComponentLabel,
                                  elecSym, l,      l1,       l10,
                                  l11,     l12,    l3,       l4,
                                  l5,      l6,     l7,       l8,
                                  l9,      rotSym, rovibSym, sym,
                                  v,       v10,    v11,      v12,
                                  v4,      v5,     v6,       v7,
                                  v8,      v9,     vibInv,   vibRefl,
                                  vibSym});
    case VAMDC::nltcs:
      return not contains_any_of(st,
                                 {F10,
                                  F11,
                                  F12,
                                  F3,
                                  F4,
                                  F5,
                                  F6,
                                  F7,
                                  F8,
                                  F9,
                                  K,
                                  Lambda,
                                  N,
                                  Omega,
                                  S,
                                  Sigma,
                                  SpinComponentLabel,
                                  elecInv,
                                  elecRefl,
                                  elecSym,
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
                                  rotSym,
                                  rovibSym,
                                  sym,
                                  v,
                                  v10,
                                  v11,
                                  v12,
                                  v4,
                                  v5,
                                  v6,
                                  v7,
                                  v8,
                                  v9,
                                  vibInv,
                                  vibRefl,
                                  vibSym});
    case VAMDC::nltos:
      return not contains_any_of(
          st, {F10,     F11,      F12,   F3,     F4,
               F5,      F6,       F7,    F8,     F9,
               K,       Lambda,   Omega, Sigma,  SpinComponentLabel,
               elecInv, elecRefl, l,     l1,     l10,
               l11,     l12,      l2,    l3,     l4,
               l5,      l6,       l7,    l8,     l9,
               rotSym,  rovibSym, sym,   v,      v10,
               v11,     v12,      v4,    v5,     v6,
               v7,      v8,       v9,    vibInv, vibRefl,
               vibSym});
    case VAMDC::sphcs:
      return not contains_any_of(st,
                                 {K,
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
                                  v,
                                  vibInv,
                                  vibRefl});
    case VAMDC::sphos:
      return not contains_any_of(st,
                                 {K,
                                  Ka,
                                  Kc,
                                  Lambda,
                                  Omega,
                                  Sigma,
                                  SpinComponentLabel,
                                  asSym,
                                  elecRefl,
                                  kronigParity,
                                  l,
                                  v,
                                  vibInv,
                                  vibRefl});
    case VAMDC::stcs:
      return not contains_any_of(st,
                                 {Ka,
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
                                  sym,
                                  v,
                                  vibRefl});
  }
  std::unreachable();
}
}  // namespace Quantum

void xml_io_stream<QuantumLevel>::parse(std::span<QuantumLevel> qns,
                                        std::istream& is_xml) {
  std::string str;
  parse_xml_tag_content_as_string(is_xml, str);
  std::istringstream is{str};
  for (auto& x : qns) is >> x;
}

void xml_io_stream<QuantumLevel>::write(std::ostream& os_xml,
                                        const QuantumLevel& qns,
                                        bofstream*,
                                        std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os_xml);
  std::println(os_xml, "{:IO}", qns);
  tag.write_to_end_stream(os_xml);
}

void xml_io_stream<QuantumLevel>::read(std::istream& is_xml,
                                       QuantumLevel& x,
                                       bifstream*) try {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  parse({&x, 1}, is_xml);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error reading QuantumLevel:\n{}", e.what())
}

void xml_io_stream<QuantumState>::parse(std::span<QuantumState> qns,
                                        std::istream& is_xml) {
  std::string str;
  parse_xml_tag_content_as_string(is_xml, str);
  std::istringstream is{str};
  for (auto& x : qns) is >> x;
}

void xml_io_stream<QuantumState>::write(std::ostream& os_xml,
                                        const QuantumState& qns,
                                        bofstream*,
                                        std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os_xml);
  std::println(os_xml, "{:IO}", qns);
  tag.write_to_end_stream(os_xml);
}

void xml_io_stream<QuantumState>::read(std::istream& is_xml,
                                       QuantumState& x,
                                       bifstream*) try {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  parse({&x, 1}, is_xml);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error reading QuantumState:\n{}", e.what())
}

void xml_io_stream<QuantumIdentifier>::write(std::ostream& os_xml,
                                             const QuantumIdentifier& qns,
                                             bofstream*,
                                             std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os_xml);
  std::println(os_xml, "{:IO}", qns);
  tag.write_to_end_stream(os_xml);
}

void xml_io_stream<QuantumIdentifier>::read(std::istream& is_xml,
                                            QuantumIdentifier& x,
                                            bifstream*) try {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  std::string str;
  parse_xml_tag_content_as_string(is_xml, str);
  x = QuantumIdentifier(str);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error reading QuantumIdentifier:\n{}", e.what())
}

void xml_io_stream<QuantumLevelIdentifier>::write(
    std::ostream& os_xml,
    const QuantumLevelIdentifier& qns,
    bofstream*,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os_xml);
  std::println(os_xml, "{:IO}", qns);
  tag.write_to_end_stream(os_xml);
}

void xml_io_stream<QuantumLevelIdentifier>::read(std::istream& is_xml,
                                                 QuantumLevelIdentifier& x,
                                                 bifstream*) try {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  std::string str;
  parse_xml_tag_content_as_string(is_xml, str);
  x = QuantumLevelIdentifier(str);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error reading QuantumLevelIdentifier:\n{}", e.what())
}

std::string to_educational_string(const QuantumIdentifier& q) {
  return Quantum::Helpers::molecular_term_symbol(q);
}
