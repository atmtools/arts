#include "quantum.h"

#include <debug.h>
#include <enumsQuantumNumberType.h>
#include <isotopologues.h>
#include <nonstd.h>
#include <rational.h>

#include <boost/container_hash/hash.hpp>
#include <boost/functional/hash.hpp>
#include <compare>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <variant>

namespace Quantum {
NoSpaceString::NoSpaceString(String&& s) : str_(std::move(s)) {
  ARTS_USER_ERROR_IF(str_.empty(), "No empty strings")
  for (auto& c : str_) {
    if (nonstd::isspace(c)) c = '*';
  }
}

const String& NoSpaceString::str() const { return str_; }

std::istream& operator>>(std::istream& is, NoSpaceString& str) {
  String s;
  is >> s;
  str = NoSpaceString(std::move(s));
  return is;
}

Value::Value(String x) : value(std::move(x)) {}

Value::Value(Rational x) : value(x) {}

Value::Value() : Value(Rational{0, 0}) {}

Value::Value(QuantumNumberType type) {
  switch (type) {
    using enum QuantumNumberType;
    case alpha:              value = NoSpaceString{"*"}; return;
    case config:             value = NoSpaceString{"*"}; return;
    case ElecStateLabel:     value = NoSpaceString{"*"}; return;
    case F:                  value = Rational{}; return;
    case F1:                 value = Rational{}; return;
    case F10:                value = Rational{}; return;
    case F11:                value = Rational{}; return;
    case F12:                value = Rational{}; return;
    case F2:                 value = Rational{}; return;
    case F3:                 value = Rational{}; return;
    case F4:                 value = Rational{}; return;
    case F5:                 value = Rational{}; return;
    case F6:                 value = Rational{}; return;
    case F7:                 value = Rational{}; return;
    case F8:                 value = Rational{}; return;
    case F9:                 value = Rational{}; return;
    case I:                  value = Rational{}; return;
    case J:                  value = Rational{}; return;
    case K:                  value = Rational{}; return;
    case Ka:                 value = Rational{}; return;
    case Kc:                 value = Rational{}; return;
    case L:                  value = Rational{}; return;
    case Lambda:             value = Rational{}; return;
    case N:                  value = Rational{}; return;
    case Omega:              value = Rational{}; return;
    case S:                  value = Rational{}; return;
    case Sigma:              value = Rational{}; return;
    case SpinComponentLabel: value = Rational{}; return;
    case asSym:              value = NoSpaceString{"*"}; return;
    case elecInv:            value = NoSpaceString{"*"}; return;
    case elecRefl:           value = NoSpaceString{"*"}; return;
    case elecSym:            value = NoSpaceString{"*"}; return;
    case kronigParity:       value = NoSpaceString{"*"}; return;
    case l:                  value = Rational{}; return;
    case l1:                 value = Rational{}; return;
    case l10:                value = Rational{}; return;
    case l11:                value = Rational{}; return;
    case l12:                value = Rational{}; return;
    case l2:                 value = Rational{}; return;
    case l3:                 value = Rational{}; return;
    case l4:                 value = Rational{}; return;
    case l5:                 value = Rational{}; return;
    case l6:                 value = Rational{}; return;
    case l7:                 value = Rational{}; return;
    case l8:                 value = Rational{}; return;
    case l9:                 value = Rational{}; return;
    case n:                  value = NoSpaceString{"*"}; return;
    case parity:             value = NoSpaceString{"*"}; return;
    case r:                  value = Rational{}; return;
    case rotSym:             value = NoSpaceString{"*"}; return;
    case rovibSym:           value = NoSpaceString{"*"}; return;
    case sym:                value = NoSpaceString{"*"}; return;
    case tau:                value = NoSpaceString{"*"}; return;
    case term:               value = NoSpaceString{"*"}; return;
    case v:                  value = Rational{}; return;
    case v1:                 value = Rational{}; return;
    case v10:                value = Rational{}; return;
    case v11:                value = Rational{}; return;
    case v12:                value = Rational{}; return;
    case v2:                 value = Rational{}; return;
    case v3:                 value = Rational{}; return;
    case v4:                 value = Rational{}; return;
    case v5:                 value = Rational{}; return;
    case v6:                 value = Rational{}; return;
    case v7:                 value = Rational{}; return;
    case v8:                 value = Rational{}; return;
    case v9:                 value = Rational{}; return;
    case vibInv:             value = NoSpaceString{"*"}; return;
    case vibRefl:            value = NoSpaceString{"*"}; return;
    case vibSym:             value = NoSpaceString{"*"}; return;
  }
  std::unreachable();
}

std::istream& operator>>(std::istream& is, Value& v) {
  if (auto* ptr = std::get_if<Rational>(&v.value)) return is >> *ptr;
  if (auto* ptr = std::get_if<NoSpaceString>(&v.value)) return is >> *ptr;
  return is;
}

template <>
const String& Value::get() const {
  return std::get<NoSpaceString>(value).str();
}

template <>
[[nodiscard]] const Rational& Value::get() const {
  return std::get<Rational>(value);
}

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
Size count_items(std::string_view s) noexcept {
  // Checks if we are in-between items, we start true as we are inbetween items
  bool last_space = true;

  Size count = 0;
  for (auto& x : s) {
    const bool this_space = nonstd::isspace(x);

    // If we had a space and now no longer do, we are in an item
    if (last_space and not this_space) count++;

    // The current state must be remembere
    last_space = this_space;
  }
  return count;
}
}  // namespace

Identifier::Identifier(const std::string_view s) try {
  std::istringstream is(String{s});  //! Fixme when view-streams are a thing

  String str;

  is >> str;
  isot = SpeciesIsotope(str);

  Size n = count_items(s) - 1;
  if (n % 3) throw std::runtime_error("Bad count of items, must be 1 + 3N");
  n /= 3;

  state.reserve(n);

  QuantumNumberType qn;
  for (Size i = 0; i < n; i++) {
    is >> qn;
    auto [it, _] =
        state.emplace(qn, UpperLower{.upper = Value(qn), .lower = Value(qn)});
    is >> it->second;
  }
} catch (std::exception& e) {
  throw std::runtime_error(std::format(
      "Cannot construct QuantumIdentifier from \"{}\", got error:\n{}",
      s,
      e.what()));
}

LevelIdentifier::LevelIdentifier(const std::string_view s) {
  std::istringstream is(String{s});  //! Fixme when view-streams are a thing

  String str;

  is >> str;
  isot = SpeciesIsotope(str);

  Size n = count_items(s) - 1;
  if (n % 2) throw std::runtime_error("Bad count of items, must be 1 + 2N");
  n /= 2;

  state.reserve(n);

  QuantumNumberType qn;
  for (Size i = 0; i < n; i++) {
    is >> qn;
    auto [it, _] = state.emplace(qn, Value(qn));
    is >> it->second;
  }
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

bool vamdcCheck(const State& l, VAMDC type) {
  switch (type) {
    case VAMDC::asymcs:
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::asSym)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::kronigParity)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::asymos:
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::asSym)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::kronigParity)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::dcs:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F2)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v1)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v2)) return false;
      if (l.contains(QuantumNumberType::v3)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::hunda:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F2)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v1)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v2)) return false;
      if (l.contains(QuantumNumberType::v3)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::hundb:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F2)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v1)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v2)) return false;
      if (l.contains(QuantumNumberType::v3)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::lpcs:
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      return true;
    case VAMDC::lpos:
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::ltcs:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::ltos:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::nltcs:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::nltos:
      if (l.contains(QuantumNumberType::F10)) return false;
      if (l.contains(QuantumNumberType::F11)) return false;
      if (l.contains(QuantumNumberType::F12)) return false;
      if (l.contains(QuantumNumberType::F3)) return false;
      if (l.contains(QuantumNumberType::F4)) return false;
      if (l.contains(QuantumNumberType::F5)) return false;
      if (l.contains(QuantumNumberType::F6)) return false;
      if (l.contains(QuantumNumberType::F7)) return false;
      if (l.contains(QuantumNumberType::F8)) return false;
      if (l.contains(QuantumNumberType::F9)) return false;
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::l1)) return false;
      if (l.contains(QuantumNumberType::l10)) return false;
      if (l.contains(QuantumNumberType::l11)) return false;
      if (l.contains(QuantumNumberType::l12)) return false;
      if (l.contains(QuantumNumberType::l2)) return false;
      if (l.contains(QuantumNumberType::l3)) return false;
      if (l.contains(QuantumNumberType::l4)) return false;
      if (l.contains(QuantumNumberType::l5)) return false;
      if (l.contains(QuantumNumberType::l6)) return false;
      if (l.contains(QuantumNumberType::l7)) return false;
      if (l.contains(QuantumNumberType::l8)) return false;
      if (l.contains(QuantumNumberType::l9)) return false;
      if (l.contains(QuantumNumberType::rotSym)) return false;
      if (l.contains(QuantumNumberType::rovibSym)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::v10)) return false;
      if (l.contains(QuantumNumberType::v11)) return false;
      if (l.contains(QuantumNumberType::v12)) return false;
      if (l.contains(QuantumNumberType::v4)) return false;
      if (l.contains(QuantumNumberType::v5)) return false;
      if (l.contains(QuantumNumberType::v6)) return false;
      if (l.contains(QuantumNumberType::v7)) return false;
      if (l.contains(QuantumNumberType::v8)) return false;
      if (l.contains(QuantumNumberType::v9)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      if (l.contains(QuantumNumberType::vibSym)) return false;
      return true;
    case VAMDC::sphcs:
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::asSym)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::kronigParity)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::sphos:
      if (l.contains(QuantumNumberType::K)) return false;
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::asSym)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::kronigParity)) return false;
      if (l.contains(QuantumNumberType::l)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::vibInv)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      return true;
    case VAMDC::stcs:
      if (l.contains(QuantumNumberType::Ka)) return false;
      if (l.contains(QuantumNumberType::Kc)) return false;
      if (l.contains(QuantumNumberType::Lambda)) return false;
      if (l.contains(QuantumNumberType::N)) return false;
      if (l.contains(QuantumNumberType::Omega)) return false;
      if (l.contains(QuantumNumberType::S)) return false;
      if (l.contains(QuantumNumberType::Sigma)) return false;
      if (l.contains(QuantumNumberType::SpinComponentLabel)) return false;
      if (l.contains(QuantumNumberType::asSym)) return false;
      if (l.contains(QuantumNumberType::elecInv)) return false;
      if (l.contains(QuantumNumberType::elecRefl)) return false;
      if (l.contains(QuantumNumberType::elecSym)) return false;
      if (l.contains(QuantumNumberType::kronigParity)) return false;
      if (l.contains(QuantumNumberType::sym)) return false;
      if (l.contains(QuantumNumberType::v)) return false;
      if (l.contains(QuantumNumberType::vibRefl)) return false;
      return true;
  }
  return false;
}
}  // namespace Quantum

template <>
struct std::hash<Quantum::Value> {
  std::size_t operator()(const Quantum::Value& g) const {
    if (auto* ptr = std::get_if<Quantum::NoSpaceString>(&g.value))
      return std::hash<String>{}(ptr->str());
    if (auto* ptr = std::get_if<Rational>(&g.value)) {
      std::size_t seed{};

      boost::hash_combine(seed, ptr->numer);
      boost::hash_combine(seed, ptr->denom);

      return seed;
    }
    return 0;
  }
};

template <>
struct std::hash<Quantum::UpperLower> {
  std::size_t operator()(const Quantum::UpperLower& g) const {
    std::size_t seed{};

    boost::hash_combine(seed, std::hash<Quantum::Value>{}(g.upper));
    boost::hash_combine(seed, std::hash<Quantum::Value>{}(g.lower));

    return seed;
  }
};

template <>
struct std::hash<QuantumLevel> {
  std::size_t operator()(const QuantumLevel& g) const {
    std::size_t seed{};

    // Has to be ordered or it will affect the hash
    for (auto& qns : enumtyps::QuantumNumberTypeTypes) {
      if (auto iter = g.find(qns); iter != g.end()) {
        boost::hash_combine(seed, iter->first);
        boost::hash_combine(seed, std::hash<Quantum::Value>{}(iter->second));
      }
    }

    return seed;
  }
};

template <>
struct std::hash<QuantumState> {
  std::size_t operator()(const QuantumState& g) const {
    std::size_t seed{};

    // Has to be ordered or it will affect the hash
    for (auto& qns : enumtyps::QuantumNumberTypeTypes) {
      if (auto iter = g.find(qns); iter != g.end()) {
        boost::hash_combine(seed, iter->first);
        boost::hash_combine(seed,
                            std::hash<Quantum::UpperLower>{}(iter->second));
      }
    }

    return seed;
  }
};

std::size_t std::hash<QuantumLevelIdentifier>::operator()(
    const QuantumLevelIdentifier& g) const {
  std::size_t seed{};

  boost::hash_combine(seed, g.isot.spec);
  boost::hash_combine(seed, g.isot.isotname);
  boost::hash_combine(seed, std::hash<QuantumLevel>{}(g.state));

  return seed;
}

std::size_t std::hash<QuantumIdentifier>::operator()(
    const QuantumIdentifier& g) const {
  std::size_t seed{};

  boost::hash_combine(seed, g.isot.spec);
  boost::hash_combine(seed, g.isot.isotname);
  boost::hash_combine(seed, std::hash<QuantumState>{}(g.state));

  return seed;
}

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
  std::println(os_xml,
               R"(<{0} name="{1}">
{2:IO}
</{0}>)",
               type_name,
               name,
               qns);
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
  std::println(os_xml,
               R"(<{0} name="{1}">
{2:IO}
</{0}>)",
               type_name,
               name,
               qns);
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
  std::println(os_xml,
               R"(<{0} name="{1}">
{2}
</{0}>)",
               type_name,
               name,
               qns);
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
  std::println(os_xml,
               R"(<{0} name="{1}">
{2}
</{0}>)",
               type_name,
               name,
               qns);
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
