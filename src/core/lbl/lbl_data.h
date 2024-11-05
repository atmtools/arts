#pragma once

#include <array.h>
#include <configtypes.h>
#include <enumsLineByLineCutoffType.h>
#include <enumsLineByLineLineshape.h>
#include <enumsLineByLineVariable.h>
#include <enumsLineShapeModelCoefficient.h>
#include <enumsLineShapeModelVariable.h>
#include <enumsQuantumNumberType.h>
#include <matpack.h>
#include <quantum_numbers.h>

#include <format>
#include <limits>
#include <unordered_set>
#include <vector>
#include <unordered_set>

#include "lbl_lineshape_model.h"
#include "lbl_zeeman.h"


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

  /** Compute the HITRAN linestrength for this line
   * 
   * @param hitran_s The HITRAN line strength
   * @param isot The isotope to use - required to get the correct partition function
   * @param T0 The reference temperature.  Defaults to HITRAN 296.0.
   * @return Numeric Hitran equivalent linestrength
   */
  [[nodiscard]] Numeric hitran_a(const Numeric hitran_s,
                                 const SpeciesIsotope& isot,
                                 const Numeric T0 = 296.0) const;

  /** The HITRAN equivalent line strength
   * 
   * @param isot The isotope to use
   * @param T0 The reference temperature.  Defaults to HITRAN 296.0.
   * @return Numeric The HITRAN equivalent line strength (including isotopoic ratio)
   */
  [[nodiscard]] Numeric hitran_s(const SpeciesIsotope& isot,
                                 const Numeric T0 = 296.0) const;

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
      case None:   return std::numeric_limits<Numeric>::infinity();
      case ByLine: return cutoff_value;
    }
    return -1;
  }

  void sort(LineByLineVariable v = LineByLineVariable::f0);

  //! Gets all the lines between (f0-get_cutoff_frequency(), f1+get_cutoff_frequency())
  [[nodiscard]] std::pair<Size, std::span<const line>> active_lines(
      Numeric f0, Numeric f1) const;

  [[nodiscard]] Rational max(QuantumNumberType) const;

  //! Returns true if the line is new for the band_data (based on quantum numbers)
  bool merge(const line& linedata);

  friend std::ostream& operator<<(std::ostream& os, const band_data& x);
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

  //! The species index if (ls_var is not invalid)
  Size spec{std::numeric_limits<Size>::max()};

  /* The variable to be used for the line shape derivative

  If ls_var is invalid, then the var variable is used for the line
  parameter.  ls_var and var are not both allowed to be invalid.
  */
  LineShapeModelVariable ls_var{static_cast<LineShapeModelVariable>(-1)};

  //! The line shape coefficient if ls_var is not invalid
  LineShapeModelCoefficient ls_coeff{
      static_cast<LineShapeModelCoefficient>(-1)};

  /* The line parameter to be used for the line shape derivative
  
  If var is invalid, then the ls_var variable is used for the line shape
  parameter.  ls_var and var are not both allowed to be invalid.
  */
  LineByLineVariable var{static_cast<LineByLineVariable>(-1)};

  [[nodiscard]] auto operator<=>(const line_key&) const = default;

  friend std::ostream& operator<<(std::ostream& os, const line_key& x);

  [[nodiscard]] Numeric& get_value(
      std::unordered_map<QuantumIdentifier, lbl::band_data>&) const;
  [[nodiscard]] const Numeric& get_value(
      const std::unordered_map<QuantumIdentifier, lbl::band_data>&) const;
};

std::ostream& operator<<(std::ostream& os, const std::vector<line>& x);

std::ostream& operator<<(
    std::ostream& os,
    const std::unordered_map<QuantumIdentifier, band_data>& x);

std::unordered_set<SpeciesEnum> species_in_bands(
    const std::unordered_map<QuantumIdentifier, band_data>& bands);
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

using AbsorptionBand = lbl::band_data;

//! A list of multiple bands
using AbsorptionBands = std::unordered_map<QuantumIdentifier, AbsorptionBand>;

template <>
struct std::formatter<lbl::line> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line& v, FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.a, sep, v.f0, sep, v.e0, sep, v.gu, sep, v.gl);
    if (not tags.short_str) tags.format(ctx, sep, v.z, sep, v.ls, sep, v.qn);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<lbl::band_data> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::band_data& v, FmtContext& ctx) const {
    const auto sep = tags.sep();

    tags.format(ctx, v.lineshape, sep, v.cutoff, sep, v.cutoff_value);
    if (not tags.short_str) tags.format(ctx, sep, v.lines);

    return ctx.out();
  }
};

template <>
struct std::formatter<lbl::line_key> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line_key& v, FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    return tags.format(ctx, v.band, sep, v.line, sep, v.spec);
  }
};
