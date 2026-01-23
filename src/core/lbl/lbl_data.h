#pragma once

#include <array.h>
#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <configtypes.h>
#include <enumsLineByLineCutoffType.h>
#include <enumsLineByLineLineshape.h>
#include <enumsLineByLineVariable.h>
#include <enumsLineShapeModelCoefficient.h>
#include <enumsLineShapeModelVariable.h>
#include <enumsQuantumNumberType.h>
#include <matpack.h>
#include <quantum.h>
#include <xml.h>

#include <format>
#include <limits>
#include <unordered_set>
#include <vector>

#include "enumsSpeciesEnum.h"
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
  QuantumState qn{};

  /*! Line strength in LTE divided by frequency-factor

  WARNING: 
  To agree with databases line strength, you must scale
  the output of this by f * (1 - exp(-hf/kT)) (c^2 / 8pi)

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric s(Numeric T, Numeric Q) const {
    return a * gu * std::exp(-e0 / (Constant::k * T)) / (Math::pow3(f0) * Q);
  }

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
  [[nodiscard]]
  Numeric ds_de0(Numeric T, Numeric Q) const {
    return -a * gu * std::exp(-e0 / (Constant::k * T)) /
           (Math::pow3(f0) * Constant::k * T * Q);
  }

  //! The ratio of ds_de0 / s
  [[nodiscard]] constexpr Numeric ds_de0_s_ratio(Numeric T) const {
    return -1 / (Constant::k * T);
  }

  /*! Derivative of s(T, Q) wrt to this->f0

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]]
  Numeric ds_df0(Numeric T, Numeric Q) const {
    return -3 * a * gu * std::exp(-e0 / (Constant::k * T)) /
           (Math::pow4(f0) * Q);
  }

  //! The ratio of ds_df0 / s
  [[nodiscard]] constexpr Numeric ds_df0_s_ratio() const { return -3 / f0; }

  /*! Derivative of s(T, Q) wrt to this->a

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]]
  Numeric ds_da(Numeric T, Numeric Q) const {
    return gu * std::exp(-e0 / (Constant::k * T)) / (Math::pow3(f0) * Q);
  }

  /*! Derivative of s(T, Q) wrt to input t

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @param[in] dQ_dt Partition function derivative at temperature wrt t [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_dT(Numeric T, Numeric Q, Numeric dQ_dT) const {
    return a * gu * (e0 * Q - Constant::k * Math::pow2(T) * dQ_dT) *
           std::exp(-e0 / (Constant::k * T)) /
           (Math::pow3(f0) * Constant::k * Math::pow2(T) * Math::pow2(Q));
  }

  /** Compute the HITRAN Einstein Coefficient for this line
   * 
   * @param hitran_s The HITRAN line strength
   * @param isot The isotope to use - required to get the correct partition function
   * @param T0 The reference temperature.  Defaults to HITRAN 296.0.
   * @return Numeric Hitran equivalent Einstein Coefficient
   */
  [[nodiscard]] Numeric hitran_a(const Numeric hitran_s,
                                 const SpeciesIsotope& isot,
                                 const Numeric T0 = 296.0) const;

  /** Compute the Einstein Coefficient for this line
   * 
   * @param s The line strength
   * @param isot The isotope to use - required to get the correct partition function
   * @param T0 The reference temperature.
   * @return Numeric equivalent Einstein Coefficient
   */
  [[nodiscard]] Numeric compute_a(const Numeric s,
                                  const SpeciesIsotope& isot,
                                  const Numeric T0) const;

  /** The HITRAN equivalent line strength
   * 
   * @param isot The isotope to use
   * @param T0 The reference temperature.  Defaults to HITRAN 296.0.
   * @return Numeric The HITRAN equivalent line strength (including isotopoic ratio)
   */
  [[nodiscard]] Numeric hitran_s(const SpeciesIsotope& isot,
                                 const Numeric T0 = 296.0) const;

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
};

struct line_pos {
  Size line;
  Size iz{std::numeric_limits<Size>::max()};
};

//! The key to finding any absorption line
struct line_key {
  //! The band the line belongs to
  QuantumIdentifier band;

  //! The line count within the band
  Size line{std::numeric_limits<Size>::max()};

  //! The species (if ls_var is not unused)
  SpeciesEnum spec{SpeciesEnum::unused};

  /* The variable to be used for the line shape derivative

  If ls_var is unused, then the var variable is used for the line
  parameter.  ls_var and var are not both allowed to be unused.
  */
  LineShapeModelVariable ls_var{LineShapeModelVariable::unused};

  //! The line shape coefficient if ls_var is not unused
  LineShapeModelCoefficient ls_coeff{LineShapeModelCoefficient::unused};

  /* The line parameter to be used for the line shape derivative
  
  If var is unused, then the ls_var variable is used for the line shape
  parameter.  ls_var and var are not both allowed to be unused.
  */
  LineByLineVariable var{LineByLineVariable::unused};

  [[nodiscard]] auto operator<=>(const line_key&) const = default;

  [[nodiscard]] Numeric& get_value(
      std::unordered_map<QuantumIdentifier, lbl::band_data>&) const;
  [[nodiscard]] const Numeric& get_value(
      const std::unordered_map<QuantumIdentifier, lbl::band_data>&) const;
};

/** Returns all species in the band, including those that are broadening species
 * 
 * @param bands The bands to search
 * @return std::unordered_set<SpeciesEnum> 
 */
std::unordered_set<SpeciesEnum> species_in_bands(
    const std::unordered_map<QuantumIdentifier, band_data>& bands);

/** Wraps keep_hitran_s for band_data per species, to remove all lines that are not in the keep map
 * 
 * @param bands The bands to use
 * @param keep A map of species to minimum hitran_s values to keep.  Missing species keep all their lines.
 * @param T0 The reference temperature.  Defaults to 296.0.
 */
void keep_hitran_s(std::unordered_map<QuantumIdentifier, band_data>& bands,
                   const std::unordered_map<SpeciesEnum, Numeric>& keep,
                   const Numeric T0 = 296);

/** Compute what lines should be kept.  Meant to be used in conjunction with keep_hitran_s.
 * 
 * The same percentile of lines are kept for all species
 * 
 * @param bands The bands to use
 * @param approx_percentile The percentile to keep [0, 100]
 * @param T0 The reference temperature.  Defaults to 296.0.
 * @return A map of species to minimum hitran_s values to keep
 */
std::unordered_map<SpeciesEnum, Numeric> percentile_hitran_s(
    const std::unordered_map<QuantumIdentifier, band_data>& bands,
    const Numeric approx_percentile,
    const Numeric T0 = 296);

/** Compute what lines should be kept.  Meant to be used in conjunction with keep_hitran_s.
 * 
 * Only species in the approx_percentile map are affected.  Otherwise acts like the pure index version.
 * 
 * @param bands The bands to use
 * @param approx_percentile The percentile to keep species: [0, 100]
 * @param T0 The reference temperature.  Defaults to 296.0.
 * @return A map of species to minimum hitran_s values to keep
 */
std::unordered_map<SpeciesEnum, Numeric> percentile_hitran_s(
    const std::unordered_map<QuantumIdentifier, band_data>& bands,
    const std::unordered_map<SpeciesEnum, Numeric>& approx_percentile,
    const Numeric T0 = 296);

Size count_lines(const std::unordered_map<QuantumIdentifier, lbl::band_data>&);
}  // namespace lbl

//! Support hashing of line keys
template <>
struct std::hash<lbl::line_key> {
  static std::size_t operator()(const lbl::line_key& x) {
    std::size_t seed = 0;

    boost::hash_combine(seed, std::hash<QuantumIdentifier>{}(x.band));
    boost::hash_combine(seed, std::hash<Size>{}(x.line));
    boost::hash_combine(seed, std::hash<SpeciesEnum>{}(x.spec));

    return seed;
  }
};

using LblLineKey = lbl::line_key;

using AbsorptionLine = lbl::line;

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

  [[nodiscard]] std::string to_string(const lbl::line& v) const;

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line& v, FmtContext& ctx) const {
    return tags.format(ctx, to_string(v));
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

  [[nodiscard]] std::string to_string(const lbl::band_data& v) const;

  template <class FmtContext>
  FmtContext::iterator format(const lbl::band_data& v, FmtContext& ctx) const {
    return tags.format(ctx, to_string(v));
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

template <>
struct xml_io_stream<AbsorptionBand> {
  static constexpr std::string_view type_name = "AbsorptionBand"sv;

  static void write(std::ostream& os,
                    const AbsorptionBand& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   AbsorptionBand& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<LblLineKey> {
  static constexpr std::string_view type_name = "LblLineKey"sv;

  static void write(std::ostream& os,
                    const LblLineKey& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, LblLineKey& x, bifstream* pbifs = nullptr);
};

template <>
std::optional<std::string> to_helper_string<AbsorptionBands>(
    const AbsorptionBands&);
