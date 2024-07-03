#pragma once

#include <matpack.h>

#include <format>
#include <limits>
#include <vector>

#include "array.h"
#include "configtypes.h"
#include "enums.h"
#include "lbl_lineshape_model.h"
#include "lbl_zeeman.h"
#include "quantum_numbers.h"

namespace lbl {
struct line {
  //! Einstein A coefficient
  Numeric a{};

  //! Line center
  Numeric f0{};

  //! Lower level energy
  Numeric e0{};

  //! Upper level statistical weight
  Numeric gu{};

  //! Lower level statistical weight
  Numeric gl{};

  //! Zeeman model
  zeeman::model z{};

  //! Line shape model
  line_shape::model ls{};

  //! Quantum numbers of this line
  QuantumNumberLocalState qn{};

  /*! Line strength in LTE divided by frequency-factor

  WARNING: 
  To agree with databases line strength, you must scale
  the output of this by f * (1 - exp(-hf/kT)) (c^2 / 8pi)

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric s(Numeric T, Numeric Q) const;

  [[nodiscard]] constexpr Numeric nlte_k(Numeric ru, Numeric rl) const {
    return (rl * gu / gl - ru) * a / Math::pow3(f0);
  }

  [[nodiscard]] constexpr Numeric dnlte_k_drl() const {
    return gu / gl * a / Math::pow3(f0);
  }

  [[nodiscard]] constexpr Numeric dnlte_k_dru() const {
    return -a / Math::pow3(f0);
  }

  [[nodiscard]] constexpr Numeric nlte_e(Numeric ru) const { return ru * a; }

  [[nodiscard]] constexpr Numeric dnlte_e_dru() const { return a; }

  [[nodiscard]] static constexpr Numeric dnlte_e_drl() { return 0; }

  /*! Derivative of s(T, Q) wrt to this->e0

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_de0(Numeric T, Numeric Q) const;

  //! The ratio of ds_de0 / s
  [[nodiscard]] constexpr Numeric ds_de0_s_ratio(Numeric T) const {
    return -1 / (Constant::k * T);
  }

  /*! Derivative of s(T, Q) wrt to this->f0

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_df0(Numeric T, Numeric Q) const;

  //! The ratio of ds_df0 / s
  [[nodiscard]] constexpr Numeric ds_df0_s_ratio() const { return -3 / f0; }

  /*! Derivative of s(T, Q) wrt to this->a

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_da(Numeric T, Numeric Q) const;

  /*! Derivative of s(T, Q) wrt to input t

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @param[in] dQ_dt Partition function derivative at temperature wrt t [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_dT(Numeric T, Numeric Q, Numeric dQ_dt) const;

  friend std::ostream& operator<<(std::ostream& os, const line& x);

  friend std::istream& operator>>(std::istream& is, line& x);
};

struct band_data {
  std::vector<line> lines{};

  LineByLineLineshape lineshape{LineByLineLineshape::VP_LTE};

  LineByLineCutoffType cutoff{LineByLineCutoffType::None};

  Numeric cutoff_value{std::numeric_limits<Numeric>::infinity()};

  [[nodiscard]] auto&& back() { return lines.back(); }
  [[nodiscard]] auto&& back() const { return lines.back(); }
  [[nodiscard]] auto&& front() { return lines.front(); }
  [[nodiscard]] auto&& front() const { return lines.front(); }
  [[nodiscard]] auto size() const { return lines.size(); }
  [[nodiscard]] auto begin() { return lines.begin(); }
  [[nodiscard]] auto begin() const { return lines.begin(); }
  [[nodiscard]] auto cbegin() const { return lines.cbegin(); }
  [[nodiscard]] auto end() { return lines.end(); }
  [[nodiscard]] auto end() const { return lines.end(); }
  [[nodiscard]] auto cend() const { return lines.cend(); }
  template <typename T>
  void push_back(T&& l) {
    lines.push_back(std::forward<T>(l));
  }
  template <typename... Ts>
  lbl::line& emplace_back(Ts&&... l) {
    return lines.emplace_back(std::forward<Ts>(l)...);
  }

  [[nodiscard]] constexpr Numeric get_cutoff_frequency() const {
    using enum LineByLineCutoffType;
    switch (cutoff) {
      case None:
        return std::numeric_limits<Numeric>::infinity();
      case ByLine:
        return cutoff_value;
    }
    return -1;
  }

  void sort(LineByLineVariable v = LineByLineVariable::f0);

  //! Gets all the lines between (f0-get_cutoff_frequency(), f1+get_cutoff_frequency())
  [[nodiscard]] std::pair<Size, std::span<const line>> active_lines(
      Numeric f0, Numeric f1) const;

  [[nodiscard]] Rational max(QuantumNumberType) const;

  friend std::ostream& operator<<(std::ostream& os, const band_data& x);
};

struct band {
  QuantumIdentifier key{"Ar-8"};
  band_data data{};

  friend std::ostream& operator<<(std::ostream& os, const band&);
};

struct line_pos {
  Size line;
  Size spec{std::numeric_limits<Size>::max()};
  Size iz{std::numeric_limits<Size>::max()};
};

//! The key to finding any absorption line
struct line_key {
  //! The band the line belongs to
  QuantumIdentifier band;

  //! The line count within the band
  Size line{std::numeric_limits<Size>::max()};

  //! The species index if (ls_var is not FINAL)
  Size spec{std::numeric_limits<Size>::max()};

  /* The variable to be used for the line shape derivative

  If ls_var is FINAL, then the var variable is used for the line
  parameter.  ls_var and var are not both allowed to be FINAL.
  */
  LineShapeModelVariable ls_var{static_cast<LineShapeModelVariable>(-1)};

  //! The line shape coefficient if ls_var is not FINAL
  LineShapeModelCoefficient ls_coeff{
      static_cast<LineShapeModelCoefficient>(-1)};

  /* The line parameter to be used for the line shape derivative
  
  If var is FINAL, then the ls_var variable is used for the line shape
  parameter.  ls_var and var are not both allowed to be FINAL.
  */
  LineByLineVariable var{static_cast<LineByLineVariable>(-1)};

  [[nodiscard]] auto operator<=>(const line_key&) const = default;

  friend std::ostream& operator<<(std::ostream& os, const line_key& x);
};

std::ostream& operator<<(std::ostream& os, const std::vector<line>& x);

std::ostream& operator<<(std::ostream& os, const std::vector<band>& x);
}  // namespace lbl

//! Support hashing of line keys
template <>
struct std::hash<lbl::line_key> {
  Size operator()(const lbl::line_key& x) const {
    return (std::hash<QuantumIdentifier>{}(x.band) << 32) ^
           std::hash<Size>{}(x.line) ^ std::hash<Size>{}(x.spec);
  }
};

using LblLineKey = lbl::line_key;

using AbsorptionBand = lbl::band;

//! A list of multiple bands
using ArrayOfAbsorptionBand = std::vector<lbl::band>;

template <>
struct std::formatter<lbl::line> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line& v, FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');

    const std::string_view sep = tags.sep();
    std::format_to(ctx.out(),
                   "{}{}{}{}{}{}{}{}{}",
                   v.a,
                   sep,
                   v.f0,
                   sep,
                   v.e0,
                   sep,
                   v.gu,
                   sep,
                   v.gl);

    if (not tags.short_str) {
      std::formatter<lbl::zeeman::model> z{};
      std::formatter<lbl::line_shape::model> ls{};
      std::formatter<QuantumNumberLocalState> qn{};

      make_compat(z, ls, qn);

      std::format_to(ctx.out(), "{}", sep);
      z.format(v.z, ctx);
      std::format_to(ctx.out(), "{}", sep);
      ls.format(v.ls, ctx);
      std::format_to(ctx.out(), "{}", sep);
      qn.format(v.qn, ctx);
    }

    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<lbl::band_data> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::band_data& v, FmtContext& ctx) const {
    std::formatter<std::vector<lbl::line>> lines{};
    std::formatter<LineByLineLineshape> lineshape{};
    std::formatter<LineByLineCutoffType> cutoff{};
    std::formatter<Numeric> cutoff_value{};

    make_compat(lines, lineshape, cutoff, cutoff_value);
    const auto sep = tags.sep();
    lineshape.format(v.lineshape, ctx);
    cutoff.format(v.cutoff, ctx);
    std::format_to(ctx.out(), "{}", sep);
    cutoff_value.format(v.cutoff_value, ctx);
    std::format_to(ctx.out(), "\n");
    return lines.format(v.lines, ctx);
  }
};

template <>
struct std::formatter<AbsorptionBand> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const AbsorptionBand& v, FmtContext& ctx) const {
    std::formatter<QuantumIdentifier> key;
    std::formatter<lbl::band_data> data;

    make_compat(key, data);

    key.format(v.key, ctx);
    std::format_to(ctx.out(), "\n");
    return data.format(v.data, ctx);
  }
};

template <>
struct std::formatter<lbl::line_key> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line_key& v, FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    std::formatter<QuantumIdentifier> band{};
    make_compat(band);
    band.format(v.band, ctx);
    std::format_to(ctx.out(), "{}{}{}{}", sep, v.line, sep, v.spec);

    return ctx.out();
  }
};
