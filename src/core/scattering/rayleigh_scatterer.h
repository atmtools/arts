#ifndef RAYLEIGH_SCATTERER_H_
#define RAYLEIGH_SCATTERER_H_

#include <atm.h>
#include <configtypes.h>
#include <enumsRayleighType.h>
#include <rtepack.h>

#include "general_tro_spectral.h"

namespace scattering {

/** Rayleigh scattering species selected by a model tag.
 */
struct RayleighScatterer {
  RayleighType type{RayleighType::EarthAir};
  Numeric diameter{};  // particle diameter [m]

  constexpr RayleighScatterer() = default;
  RayleighScatterer(RayleighType type, Numeric diameter = 0.0);

  /** Scattering and absorption coefficients for one frequency.
   *
   * Returns {scat_coeff, abs_coeff} [1/m].
   * Dispatches on `type`.
   */
  [[nodiscard]] std::pair<Numeric, Numeric> cross_sections(
      Numeric freq_hz, const AtmPoint& atm_point) const;

  /** TRO spectral bulk scattering properties (for DISORT etc.).
   *
   * Rayleigh phase function: P_11(θ) = (3/4)(1 + cos²θ).
   * Legendre expansion: β_0 = 1, β_2 = 1/10, others = 0.
   */
  [[nodiscard]] ScatteringTroSpectralVector
  get_bulk_scattering_properties_tro_spectral(const AtmPoint& atm_point,
                                              const Vector& f_grid,
                                              Index l) const;

  /** TRO gridded bulk scattering properties.
   *
   * Evaluates the Rayleigh phase matrix at each zenith angle
   * in the grid.  All 6 compact TRO elements are filled.
   */
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded>
  get_bulk_scattering_properties_tro_gridded(
      const AtmPoint& atm_point,
      const Vector& f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  /** Radar single-scattering data at exact backscatter (180°).
   *
   * Fills per-frequency backscatter phase matrix and extinction
   * propagation matrix.  These are bulk (nd-weighted) quantities.
   */
  void get_radar_single_scat(MuelmatVector& Z_back,
                             PropmatVector& K_ext,
                             const AtmPoint& atm_point,
                             const Vector& f_grid) const;

  [[nodiscard]] Numeric get_diameter() const { return diameter; }
  void set_diameter(Numeric d) { diameter = d; }

  friend std::ostream& operator<<(std::ostream& os, const RayleighScatterer& s);
};
}  // namespace scattering

template <>
struct std::formatter<scattering::RayleighScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const scattering::RayleighScatterer& v,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "RayleighScatterer"sv);
    }

    return tags.format(ctx, v.type, ", diameter=", v.diameter);
  }
};

template <>
struct xml_io_stream<scattering::RayleighScatterer> {
  static constexpr std::string_view type_name = "RayleighScatterer";

  static void write(std::ostream&,
                    const scattering::RayleighScatterer&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);

  static void read(std::istream&,
                   scattering::RayleighScatterer&,
                   bifstream* = nullptr);
};

#endif  // RAYLEIGH_SCATTERER_H_
