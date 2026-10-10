#ifndef HENYEY_GREENSTEIN_H_
#define HENYEY_GREENSTEIN_H_

#include <configtypes.h>
#include <operators.h>

#include <tuple>

#include "atm.h"
#include "bulk_scattering_properties.h"
#include "general_tro_spectral.h"
#include "properties.h"

using ExtSSACallback = CustomOperator<std::pair<Numeric, Numeric>, Numeric, const AtmPoint&>;

namespace scattering {
struct ExtinctionSSALookup {
  ScatteringSpeciesProperty extinction_field{};
  ScatteringSpeciesProperty ssa_field{};
  ExtinctionSSALookup(ScatteringSpeciesProperty extinction_field_, ScatteringSpeciesProperty ssa_field_);

  std::pair<Numeric, Numeric> operator()(Numeric, const AtmPoint& atm_point);
};

/** A Henyey-Greenstein scatterer with asymmetry parameter g
 *
 * The scattering matrix is F11 = F22 = p, F33 = F44 = p f and F12 = F34 = 0,
 * with p the Henyey-Greenstein phase function and
 * f = (3 cos(Theta) - cos^3(Theta)) / 2.  It produces no polarization, and it
 * is regular at forward and backward scattering, as a physical scattering
 * matrix is: F22 + F33 = p (1 + cos(Theta))^2 (2 - cos(Theta)) / 2 and
 * F22 - F33 = p (1 - cos(Theta))^2 (2 + cos(Theta)) / 2, so F33 = F22 at
 * 0 deg, F33 = -F22 and F44 = F11 - 2 F22 at 180 deg, and its expansion in
 * generalized spherical functions converges as fast as the Legendre series
 * of p.
 */
struct HenyeyGreensteinScatterer {
  ExtSSACallback ext_ssa_callback{};

  Numeric g = 0.0;

  HenyeyGreensteinScatterer() = default;
  HenyeyGreensteinScatterer(ExtSSACallback ext_ssa_callback, const Numeric& g);
  HenyeyGreensteinScatterer(ScatteringSpeciesProperty extinction_field,
                            ScatteringSpeciesProperty ssa_field,
                            const Numeric&            g_);

  HenyeyGreensteinScatterer(const HenyeyGreensteinScatterer&)                = default;
  HenyeyGreensteinScatterer(HenyeyGreensteinScatterer&&) noexcept            = default;
  HenyeyGreensteinScatterer& operator=(const HenyeyGreensteinScatterer&)     = default;
  HenyeyGreensteinScatterer& operator=(HenyeyGreensteinScatterer&&) noexcept = default;

  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded>
  get_bulk_scattering_properties_tro_gridded(const AtmPoint&,
                                             const Vector&                    f_grid,
                                             std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded>
  get_bulk_scattering_properties_tro_gridded_derivative(const AtmPoint&,
                                                        const Vector&,
                                                        std::shared_ptr<ZenithAngleGrid>,
                                                        const AtmKeyVal&) const;

  [[nodiscard]] ScatteringTroSpectralVector get_bulk_scattering_properties_tro_spectral(const AtmPoint&,
                                                                                        const Vector& f_grid,
                                                                                        Index         l) const;

  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
  get_bulk_scattering_properties_aro_gridded(const AtmPoint&,
                                             const Vector&                                f_grid,
                                             const Vector&                                za_inc_grid,
                                             const Vector&                                delta_aa_grid,
                                             std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;

  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
  get_bulk_scattering_properties_aro_gridded_derivative(const AtmPoint&,
                                                        const Vector&,
                                                        const Vector&,
                                                        const Vector&,
                                                        std::shared_ptr<scattering::ZenithAngleGrid>,
                                                        const AtmKeyVal&) const;

  /** The SHT of degree and order of the laboratory-frame phase matrix over the scattering directions, per
   *  incidence zenith angle [deg]: the closed form sampled exactly on the SHT's own grid and transformed */
  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral>
  get_bulk_scattering_properties_aro_spectral(
      const AtmPoint&, const Vector& f_grid, const Vector& za_inc_grid, Index degree, Index order) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame phase matrix, exactly
   *
   * See tro_lab_frame_fourier_modes; the phase matrix is the closed form at
   * every scattering angle.  The extinction matrix and absorption vector are
   * at the incidence zenith angles.
   *
   * @param za_inc_grid The incidence zenith angles [deg]
   * @param za_scat_grid The scattering zenith angles [deg], ascending
   * @param max_mode The highest mode, >= 0
   */
  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>
  get_bulk_scattering_properties_aro_fourier(const AtmPoint&,
                                             const Vector& f_grid,
                                             const Vector& za_inc_grid,
                                             const Vector& za_scat_grid,
                                             Index         max_mode) const;

  [[nodiscard]] Numeric get_g() const { return g; };
  void                  set_g(const Numeric& g_) { g = g_; };

  friend std::ostream& operator<<(std::ostream& os, const HenyeyGreensteinScatterer& scatterer);
};
}  // namespace scattering

template <> struct std::formatter<scattering::HenyeyGreensteinScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const scattering::HenyeyGreensteinScatterer& v, FmtContext& ctx) const {
    if (tags.names) { return tags.format(ctx, "HenyeyGreensteinScatterer"sv); }

    return tags.format(ctx, v.get_g());
  }
};

template <> struct xml_io_stream<scattering::HenyeyGreensteinScatterer> {
  static constexpr std::string_view type_name = "HenyeyGreensteinScatterer";

  static void write(std::ostream&,
                    const scattering::HenyeyGreensteinScatterer&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);

  static void read(std::istream&, scattering::HenyeyGreensteinScatterer&, bifstream* = nullptr);
};

#endif  // HENYEY_GREENSTEIN_H_
