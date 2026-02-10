#pragma once

#include <format_tags.h>
#include <isotopologues.h>
#include <mystring.h>
#include <rational.h>
#include <xml.h>

#include <boost/container_hash/hash.hpp>
#include <concepts>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>

namespace Quantum {
template <typename T>
concept quantum_value_holder =
    std::same_as<T, String> or std::same_as<T, Rational>;

struct Value {
  std::variant<Rational, String> value;

  Value();

  Value(String);             // Default initialize by string
  Value(Index);              // Default initialize by index
  Value(Rational);           // Default initialize by rational
  Value(QuantumNumberType);  // Default initialize by type
  Value(const Value&)                = default;
  Value(Value&&) noexcept            = default;
  Value& operator=(const Value&)     = default;
  Value& operator=(Value&&) noexcept = default;

  template <quantum_value_holder T>
  T get() const {
    return std::get<T>(value);
  }

  operator Rational() const;
  operator String() const;

  void set(const std::string_view);

  std::strong_ordering operator<=>(const Value& l) const = default;

  friend std::istream& operator>>(std::istream&, Value&);
};

struct UpperLower {
  Value upper;
  Value lower;

  std::strong_ordering operator<=>(const UpperLower& g) const = default;

  friend std::istream& operator>>(std::istream&, UpperLower&);
};

using Level = std::unordered_map<QuantumNumberType, Value>;

using State = std::unordered_map<QuantumNumberType, UpperLower>;

Level level_from(const std::string_view);

State state_from(const std::string_view);

[[nodiscard]] Level upper_level(const State& state);
[[nodiscard]] Level lower_level(const State& state);

std::istream& operator>>(std::istream&, Level&);
std::istream& operator>>(std::istream&, State&);

struct LevelIdentifier {
  SpeciesIsotope isot;
  Level state;

  LevelIdentifier(SpeciesIsotope, const Level&);
  explicit LevelIdentifier(const std::string_view);
  explicit LevelIdentifier(const SpeciesIsotope);
  LevelIdentifier()                                      = default;
  LevelIdentifier(const LevelIdentifier&)                = default;
  LevelIdentifier(LevelIdentifier&&) noexcept            = default;
  LevelIdentifier& operator=(const LevelIdentifier&)     = default;
  LevelIdentifier& operator=(LevelIdentifier&&) noexcept = default;

  [[nodiscard]] Size size() const;

  bool operator==(const LevelIdentifier& l) const;
  bool operator!=(const LevelIdentifier& l) const;
  std::strong_ordering operator<=>(const LevelIdentifier& g) const;
};

struct Identifier {
  SpeciesIsotope isot;
  State state;

  Identifier(SpeciesIsotope, const State&);
  explicit Identifier(const std::string_view);
  explicit Identifier(const SpeciesIsotope);
  Identifier()                                 = default;
  Identifier(const Identifier&)                = default;
  Identifier(Identifier&&) noexcept            = default;
  Identifier& operator=(const Identifier&)     = default;
  Identifier& operator=(Identifier&&) noexcept = default;

  [[nodiscard]] LevelIdentifier lower() const;
  [[nodiscard]] LevelIdentifier upper() const;

  [[nodiscard]] Size size() const;

  bool operator==(const Identifier& l) const;
  bool operator!=(const Identifier& l) const;
  std::strong_ordering operator<=>(const Identifier& g) const;
};

State from_hitran(std::string_view upp, std::string_view low);

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
/** Selects the global state
 *
 * @param qns Quantum numbers to select, must just be iterable for template to work
 * @param qid State to select from
 * @return State of all qns in qid
 */
template <typename list_type>
[[nodiscard]] Identifier global_state(const list_type& qns,
                                      const Identifier& qid) {
  Identifier out(qid.isot);
  for (auto qn : qns) {
    if (qid.state.contains(qn)) out.state.emplace(qn, qid.state.at(qn));
  }
  return out;
}
/** Selects the global state
 *
 * @param qns Quantum numbers to select, must just be iterable for template to work
 * @param qid State to select from
 * @return State of all qns in qid
 */
template <typename list_type>
[[nodiscard]] State local_state(const list_type& qns, const Identifier& qid) {
  State out;
  for (auto qn : qns) {
    if (qid.state.contains(qn)) out.emplace(qn, qid.state.at(qn));
  }
  return out;
}

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

bool vamdcCheck(const State& l, VAMDC type);

bool contains_any_of(const State& state,
                     std::initializer_list<QuantumNumberType> keys);

bool contains_all_of(const State& state,
                     std::initializer_list<QuantumNumberType> keys);
}  // namespace Quantum

using QuantumIdentifier      = Quantum::Identifier;
using QuantumLevelIdentifier = Quantum::LevelIdentifier;
using QuantumState           = Quantum::State;
using QuantumLevel           = Quantum::Level;

using ArrayOfQuantumIdentifier      = Array<QuantumIdentifier>;
using ArrayOfQuantumLevelIdentifier = Array<QuantumLevelIdentifier>;

template <>
struct std::hash<Quantum::Value> {
  static std::size_t operator()(const Quantum::Value& g) {
    if (auto* ptr = std::get_if<String>(&g.value))
      return std::hash<String>{}(*ptr);
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
  static std::size_t operator()(const Quantum::UpperLower& g) {
    std::size_t seed{};

    boost::hash_combine(seed, std::hash<Quantum::Value>{}(g.upper));
    boost::hash_combine(seed, std::hash<Quantum::Value>{}(g.lower));

    return seed;
  }
};

template <>
struct std::hash<QuantumLevel> {
  static std::size_t operator()(const QuantumLevel& g) {
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
  static std::size_t operator()(const QuantumState& g) {
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

template <>
struct std::hash<QuantumLevelIdentifier> {
  static std::size_t operator()(const QuantumLevelIdentifier& g) {
    std::size_t seed{};

    boost::hash_combine(seed, g.isot.spec);
    boost::hash_combine(seed, g.isot.isotname);
    boost::hash_combine(seed, std::hash<QuantumLevel>{}(g.state));

    return seed;
  }
};

template <>
struct std::hash<QuantumIdentifier> {
  static std::size_t operator()(const QuantumIdentifier& g) {
    std::size_t seed{};

    boost::hash_combine(seed, g.isot.spec);
    boost::hash_combine(seed, g.isot.isotname);
    boost::hash_combine(seed, std::hash<QuantumState>{}(g.state));

    return seed;
  }
};

template <>
struct std::formatter<Quantum::Value> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Quantum::Value& q, FmtContext& ctx) const {
    return tags.format(ctx, q.value);
  }
};

template <>
struct std::formatter<Quantum::UpperLower> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Quantum::UpperLower& q,
                              FmtContext& ctx) const {
    return tags.format(ctx, q.upper, ' ', q.lower);
  }
};

template <>
struct std::formatter<QuantumLevel> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumLevel& q, FmtContext& ctx) const {
    std::string_view first = ""sv;

    if (tags.io) tags.format(ctx, q.size(), std::exchange(first, " "sv));

    for (auto& v : q) {
      tags.format(ctx, std::exchange(first, " "sv), v.first, " "sv, v.second);
    }

    return ctx.out();
  }
};

template <>
struct std::formatter<QuantumState> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumState& q, FmtContext& ctx) const {
    std::string_view first = ""sv;

    if (tags.io) tags.format(ctx, q.size(), std::exchange(first, " "sv));

    for (auto& v : q) {
      tags.format(ctx, std::exchange(first, " "sv), v.first, " "sv, v.second);
    }

    return ctx.out();
  }
};

std::string to_educational_string(const QuantumIdentifier&);

template <>
struct std::formatter<QuantumIdentifier> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumIdentifier& q,
                              FmtContext& ctx) const {
    if (tags.help) {
      tags.format(ctx, "Species: "sv, q.isot, " Quantum Numbers:"sv);
      for (auto& [k, v] : q.state) {
        tags.format(ctx, ' ', k, ": "sv, v, ',');
      }
      return ctx.out();
    }

    if (tags.depth > 0) {
      return tags.format(ctx, to_educational_string(q));
    }

    tags.format(ctx, q.isot);

    for (auto v : enumtyps::QuantumNumberTypeTypes) {
      if (q.state.contains(v)) {
        tags.format(ctx, " "sv, v, " "sv, q.state.at(v));
      }
    }

    return ctx.out();
  }
};

template <>
struct std::formatter<QuantumLevelIdentifier> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumLevelIdentifier& q,
                              FmtContext& ctx) const {
    tags.format(ctx, q.isot);

    for (auto v : enumtyps::QuantumNumberTypeTypes) {
      if (q.state.contains(v)) {
        tags.format(ctx, " "sv, v, " "sv, q.state.at(v));
      }
    }

    return ctx.out();
  }
};

template <>
struct xml_io_stream<QuantumLevel> {
  constexpr static std::string_view type_name = "QuantumLevel"sv;

  static void parse(std::span<QuantumLevel>, std::istream& is);

  static void write(std::ostream& os,
                    const QuantumLevel& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   QuantumLevel& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<QuantumState> {
  constexpr static std::string_view type_name = "QuantumState"sv;

  static void parse(std::span<QuantumState>, std::istream& is);

  static void write(std::ostream& os,
                    const QuantumState& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   QuantumState& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<QuantumIdentifier> {
  constexpr static std::string_view type_name = "QuantumIdentifier"sv;

  static void write(std::ostream& os,
                    const QuantumIdentifier& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   QuantumIdentifier& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<QuantumLevelIdentifier> {
  constexpr static std::string_view type_name = "QuantumLevelIdentifier"sv;

  static void write(std::ostream& os,
                    const QuantumLevelIdentifier& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   QuantumLevelIdentifier& x,
                   bifstream* pbifs = nullptr);
};
