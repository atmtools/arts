#pragma once

#include <format_tags.h>
#include <isotopologues.h>
#include <mystring.h>
#include <rational.h>
#include <xml.h>

#include <concepts>
#include <unordered_map>
#include <utility>
#include <variant>

namespace Quantum {
class NoSpaceString {
  String str_;

 public:
  NoSpaceString(String&&);  // Replaces ' ' with '*'
  NoSpaceString(const NoSpaceString&)                = default;
  NoSpaceString(NoSpaceString&&) noexcept            = default;
  NoSpaceString& operator=(const NoSpaceString&)     = default;
  NoSpaceString& operator=(NoSpaceString&&) noexcept = default;

  //! Get the value, you cannot modify it
  [[nodiscard]] const String& str() const;

  std::strong_ordering operator<=>(const NoSpaceString& l) const = default;

  friend std::istream& operator>>(std::istream&, NoSpaceString&);
};

struct Value {
  std::variant<Rational, NoSpaceString> value;

  Value();

  explicit Value(String);    // Default initialize by string
  explicit Value(Rational);  // Default initialize by rational
  explicit Value(const std::integral auto& x) : value(Rational(x, 1)) {}
  explicit Value(QuantumNumberType);  // Default initialize by type
  Value(const Value&)                = default;
  Value(Value&&) noexcept            = default;
  Value& operator=(const Value&)     = default;
  Value& operator=(Value&&) noexcept = default;

  template <typename T>
  const T& get() const = delete;

  template <>
  [[nodiscard]] const String& get() const;

  template <>
  [[nodiscard]] const Rational& get() const;

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

bool contains_any_of(const State& state, std::initializer_list<QuantumNumberType> keys);

bool contains_all_of(const State& state, std::initializer_list<QuantumNumberType> keys);
}  // namespace Quantum

using QuantumIdentifier      = Quantum::Identifier;
using QuantumLevelIdentifier = Quantum::LevelIdentifier;
using QuantumState           = Quantum::State;
using QuantumLevel           = Quantum::Level;

using ArrayOfQuantumIdentifier      = Array<QuantumIdentifier>;
using ArrayOfQuantumLevelIdentifier = Array<QuantumLevelIdentifier>;

template <>
struct std::hash<QuantumLevelIdentifier> {
  std::size_t operator()(const QuantumLevelIdentifier& g) const;
};

template <>
struct std::hash<QuantumIdentifier> {
  std::size_t operator()(const QuantumIdentifier& g) const;
};

template <>
struct std::formatter<Quantum::NoSpaceString> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Quantum::NoSpaceString& q,
                              FmtContext& ctx) const {
    return tags.format(ctx, q.str());
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
    tags.format(ctx, q.isot);

    for (auto& v : q.state) tags.format(ctx, " "sv, v.first, " "sv, v.second);

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

    for (auto& v : q.state) tags.format(ctx, " "sv, v.first, " "sv, v.second);

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
