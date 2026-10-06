#pragma once

#include "bulk_scattering_properties.h"
#include "general_tro_spectral.h"
#include "particle_habit.h"
#include "psd.h"

namespace scattering {

using PSD = std::variant<MonodispersePSD,
                         MGDSingleMoment,
                         MGDMass,
                         MGDTwoMoment,
                         DelanoeEtAl14,
                         FieldEtAl07,
                         McFarquharHeymsfield97,
                         BinnedPSD>;

/*** A scattering habit
 *
 * A scattering habit combines a particle habit with an additional PSD
 * and thus defines a mapping between atmospheric scattering species properties
 * and corresponding bulk skattering properties.
 */
class ScatteringHabit {
 public:
  ScatteringHabit() = default;
  ScatteringHabit(const ParticleHabit& particle_habit_,
                  const PSD&           psd_,
                  Numeric              mass_size_rel_a_ = -1.0,
                  Numeric              mass_size_rel_b_ = -1.0);

  BulkScatteringPropertiesTROGridded get_bulk_scattering_properties_tro_gridded(const AtmPoint&,
                                                                                const Vector& f_grid,
                                                                                const Numeric f_tol = 1e-3) const;

  BulkScatteringProperties<Format::TRO, Representation::Gridded> get_bulk_scattering_properties_tro_gridded(
      const AtmPoint& point, const Vector& f_grid, std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  BulkScatteringProperties<Format::TRO, Representation::Gridded> get_bulk_scattering_properties_tro_gridded_derivative(
      const AtmPoint&, const Vector&, std::shared_ptr<ZenithAngleGrid>, const AtmKeyVal&) const;

  /** The bulk Legendre series to degree at f_grid, summed over the particles
   *
   * Every particle must hold a TRO Legendre series of at least that degree.
   * The series are interpolated linearly in temperature and frequency.  Each
   * coefficient is the scattering-plane Mueller matrix on Y_l0 (see
   * tro_legendre.h).
   */
  ScatteringTroSpectralVector get_bulk_scattering_properties_tro_spectral(const AtmPoint&,
                                                                          const Vector& f_grid,
                                                                          const Index   degree) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame bulk phase matrix at the zenith angles [deg]
   *
   * See ssd_to_aro_fourier: TRO Legendre series and SHT ARO data give them
   * exactly at any scattering zenith angle, gridded ARO data on their own
   * zenith grids.  Gridded TRO data must be converted to a Legendre series
   * first.
   */
  BulkScatteringProperties<Format::ARO, Representation::Fourier> get_bulk_scattering_properties_aro_fourier(
      const AtmPoint&, const Vector& f_grid, const Vector& za_inc_grid, const Vector& za_scat_grid, Index max_mode)
      const;

  BulkScatteringProperties<Format::ARO, Representation::Gridded> get_bulk_scattering_properties_aro_gridded(
      const AtmPoint&                  point,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  BulkScatteringProperties<Format::ARO, Representation::Gridded> get_bulk_scattering_properties_aro_gridded_derivative(
      const AtmPoint&, const Vector&, const Vector&, const Vector&, std::shared_ptr<ZenithAngleGrid>, const AtmKeyVal&)
      const;

  //  BulkScatteringProperties<Format::TRO, Representation::Gridded>
  //  get_bulk_scattering_properties_tro_spectral(
  //      const AtmPoint&,
  //      const Vector& f_grid,
  //      Index l) const;

 private:
  //! The number density [m^-3] each particle stands for at point (see scattering::number_densities)
  Vector  number_densities(const AtmPoint& point) const;
  PSDData number_densities_with_derivatives(const AtmPoint& point) const;

  ParticleHabit particle_habit;
  Numeric       mass_size_rel_a, mass_size_rel_b;
  PSD           psd;
};

}  // namespace scattering

template <> struct std::formatter<scattering::ScatteringHabit> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const scattering::ScatteringHabit&, FmtContext& ctx) const {
    if (tags.names) { return tags.format(ctx, "ScatteringHabit"sv); }

    return tags.format(ctx);
  }
};

template <> struct xml_io_stream<scattering::ScatteringHabit> {
  static constexpr std::string_view type_name = "ScatteringHabit";

  static void write(std::ostream&, const scattering::ScatteringHabit&, bofstream* = nullptr, std::string_view = ""sv);

  static void read(std::istream&, scattering::ScatteringHabit&, bifstream* = nullptr);
};
