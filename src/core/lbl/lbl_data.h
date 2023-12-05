#pragma once

#include <matpack.h>

#include <limits>
#include <vector>

#include "configtypes.h"
#include "lbl_lineshape_model.h"
#include "lbl_zeeman.h"
#include "quantum_numbers.h"

namespace lbl {
ENUMCLASS(variable, char, f0, e0, a)

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
  [[nodiscard]] Numeric s(Numeric T, Numeric Q) const noexcept;

  [[nodiscard]] constexpr Numeric nlte_k(Numeric ru, Numeric rl) const noexcept {
    return (rl * gu / gl - ru) * a / Math::pow3(f0);
  }

  [[nodiscard]] constexpr Numeric dnlte_k_drl() const noexcept {
    return gu / gl * a / Math::pow3(f0);
  }

  [[nodiscard]] constexpr Numeric dnlte_k_dru() const noexcept {
    return - a / Math::pow3(f0);
  }

  [[nodiscard]] constexpr Numeric nlte_e(Numeric ru) const noexcept {
    return ru * a;
  }

  [[nodiscard]] constexpr Numeric dnlte_e_dru() const noexcept {
    return a;
  }

  [[nodiscard]] static constexpr Numeric dnlte_e_drl() noexcept {
    return 0;
  }

  /*! Derivative of s(T, Q) wrt to this->e0

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_de0(Numeric T, Numeric Q) const noexcept;

  //! The ratio of ds_de0 / s
  [[nodiscard]] constexpr Numeric ds_de0_s_ratio(Numeric T) const noexcept {
    return -1 / (Constant::k * T);
  }

  /*! Derivative of s(T, Q) wrt to this->f0

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_df0(Numeric T, Numeric Q) const noexcept;

  //! The ratio of ds_df0 / s
  [[nodiscard]] constexpr Numeric ds_df0_s_ratio() const noexcept {
    return -3 / f0;
  }

  /*! Derivative of s(T, Q) wrt to this->a

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_da(Numeric T, Numeric Q) const noexcept;

  /*! Derivative of s(T, Q) wrt to input t

  @param[in] T Temperature [K]
  @param[in] Q Partition function at temperature [-]
  @param[in] dQ_dt Partition function derivative at temperature wrt t [-]
  @return Line strength in LTE divided by frequency [per m^2]
  */
  [[nodiscard]] Numeric ds_dT(Numeric T,
                              Numeric Q,
                              Numeric dQ_dt) const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const line& x);

  friend std::istream& operator>>(std::istream& is, line& x);
};

ENUMCLASS(CutoffType, char, None, ByLine)

ENUMCLASS(Lineshape, char, VP_LTE, VP_LINE_NLTE, VP_ECS_MAKAROV)

struct band_data {
  std::vector<line> lines{};

  Lineshape lineshape{Lineshape::VP_LTE};

  CutoffType cutoff{CutoffType::None};

  Numeric cutoff_value{std::numeric_limits<Numeric>::infinity()};

  [[nodiscard]] auto&& back() noexcept { return lines.back(); }
  [[nodiscard]] auto&& back() const noexcept { return lines.back(); }
  [[nodiscard]] auto&& front() noexcept { return lines.front(); }
  [[nodiscard]] auto&& front() const noexcept { return lines.front(); }
  [[nodiscard]] auto size() const noexcept { return lines.size(); }
  [[nodiscard]] auto begin() noexcept { return lines.begin(); }
  [[nodiscard]] auto begin() const noexcept { return lines.begin(); }
  [[nodiscard]] auto cbegin() const noexcept { return lines.cbegin(); }
  [[nodiscard]] auto end() noexcept { return lines.end(); }
  [[nodiscard]] auto end() const noexcept { return lines.end(); }
  [[nodiscard]] auto cend() const noexcept { return lines.cend(); }
  template <typename T>
  void push_back(T&& l) noexcept {
    lines.push_back(std::forward<T>(l));
  }
  template <typename... Ts>
  lbl::line& emplace_back(Ts&&... l) noexcept {
    return lines.emplace_back(std::forward<Ts>(l)...);
  }

  [[nodiscard]] constexpr Numeric get_cutoff_frequency() const noexcept {
    using enum CutoffType;
    switch (cutoff) {
      case None:
        return std::numeric_limits<Numeric>::infinity();
      case ByLine:
        return cutoff_value;
      case FINAL:;  // Leave last
    }
    return -1;
  }

  void sort(variable v = variable::f0);

  //! Gets all the lines between (f0-get_cutoff_frequency(), f1+get_cutoff_frequency())
  [[nodiscard]] std::pair<Size, std::span<const line>> active_lines(
      Numeric f0, Numeric f1) const;

[[nodiscard]] Rational max(QuantumNumberType) const;

  friend std::ostream& operator<<(std::ostream& os, const band_data& x);
};

struct band {
  QuantumIdentifier key{};
  band_data data{};

  friend std::ostream& operator<<(std::ostream& os, const band&);
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
  line_shape::variable ls_var{line_shape::variable::FINAL};

  //! The line shape coefficient if ls_var is not FINAL
  temperature::coefficient ls_coeff{temperature::coefficient::FINAL};

  /* The line parameter to be used for the line shape derivative
  
  If var is FINAL, then the ls_var variable is used for the line shape
  parameter.  ls_var and var are not both allowed to be FINAL.
  */
  variable var{variable::FINAL};

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
using AbsorptionBands = std::vector<lbl::band>;
