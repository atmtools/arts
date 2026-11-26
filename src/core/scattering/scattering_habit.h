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
 * and thus defines a mapping between scattering species properties
 * and corresponding bulk skattering properties.
 *
 * The explicit methods take the temperature [K] and the scattering species
 * properties that the PSD reads.  The derivatives take their derivatives by
 * the target as well.  The template wrappers take the point of a medium and
 * its spectral point, and read the temperature and the properties from the
 * point.  The habit has no spectral properties, so they do not read the
 * spectral point.
 */
class ScatteringHabit {
 public:
  ScatteringHabit() = default;
  ScatteringHabit(const ParticleHabit& particle_habit_,
                  const PSD&           psd_,
                  Numeric              mass_size_rel_a_ = -1.0,
                  Numeric              mass_size_rel_b_ = -1.0);

  //! The bulk properties on the particles' own scattering angles
  [[nodiscard]] BulkScatteringPropertiesTROGridded tro_gridded(Numeric                     temperature,
                                                               const ScatteringSpeciesMap& ssprops,
                                                               const Vector&               f_grid) const;

  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(
      Numeric                          temperature,
      const ScatteringSpeciesMap&      ssprops,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  //! The derivative, from the derivatives of the temperature and of the scattering species properties
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_derivative(
      Numeric                          temperature,
      const ScatteringSpeciesMap&      ssprops,
      Numeric                          d_temperature,
      const ScatteringSpeciesMap&      d_ssprops,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  /** The bulk Legendre series to degree at f_grid, summed over the particles
   *
   * Every particle must hold a TRO Legendre series of at least that degree.
   * The series are interpolated linearly in temperature and frequency.  Each
   * coefficient is the scattering-plane Mueller matrix on Y_l0 (see
   * tro_legendre.h).
   */
  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(Numeric                     temperature,
                                                         const ScatteringSpeciesMap& ssprops,
                                                         const Vector&               f_grid,
                                                         const Index                 degree) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame bulk phase matrix at the zenith angles [deg]
   *
   * See ssd_to_aro_fourier: TRO Legendre series and SHT ARO data give them
   * exactly at any scattering zenith angle, gridded ARO data on their own
   * zenith grids.  Gridded TRO data must be converted to a Legendre series
   * first.
   */
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(
      Numeric                     temperature,
      const ScatteringSpeciesMap& ssprops,
      const Vector&               f_grid,
      const Vector&               za_inc_grid,
      const Vector&               za_scat_grid,
      Index                       max_mode) const;

  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
      Numeric                          temperature,
      const ScatteringSpeciesMap&      ssprops,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  //! The derivative, from the derivatives of the temperature and of the scattering species properties
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_derivative(
      Numeric                          temperature,
      const ScatteringSpeciesMap&      ssprops,
      Numeric                          d_temperature,
      const ScatteringSpeciesMap&      d_ssprops,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringPropertiesTROGridded tro_gridded(const Point& point,
                                                               const SpectralPoint&,
                                                               const Vector& f_grid) const {
    return tro_gridded(point.temperature, scattering_species_properties(point), f_grid);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
    return tro_gridded(point.temperature, scattering_species_properties(point), f_grid, std::move(za_scat_grid));
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_derivative(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const medium_key<Point>&         target) const {
    return tro_gridded_derivative(point.temperature,
                                  scattering_species_properties(point),
                                  temperature_derivative(target),
                                  scattering_species_property_derivatives(target),
                                  f_grid,
                                  std::move(za_scat_grid));
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(const Point& point,
                                                         const SpectralPoint&,
                                                         const Vector& f_grid,
                                                         const Index   degree) const {
    return tro_spectral(point.temperature, scattering_species_properties(point), f_grid, degree);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(
      const Point& point,
      const SpectralPoint&,
      const Vector& f_grid,
      const Vector& za_inc_grid,
      const Vector& za_scat_grid,
      Index         max_mode) const {
    return aro_fourier(
        point.temperature, scattering_species_properties(point), f_grid, za_inc_grid, za_scat_grid, max_mode);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
    return aro_gridded(point.temperature,
                       scattering_species_properties(point),
                       f_grid,
                       za_inc_grid,
                       delta_aa_grid,
                       std::move(za_scat_grid));
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_derivative(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const medium_key<Point>&         target) const {
    return aro_gridded_derivative(point.temperature,
                                  scattering_species_properties(point),
                                  temperature_derivative(target),
                                  scattering_species_property_derivatives(target),
                                  f_grid,
                                  za_inc_grid,
                                  delta_aa_grid,
                                  std::move(za_scat_grid));
  }

 private:
  //! The number density [m^-3] each particle stands for (see scattering::number_densities)
  [[nodiscard]] Vector  number_densities(Numeric temperature, const ScatteringSpeciesMap& ssprops) const;
  [[nodiscard]] PSDData number_densities_with_derivatives(Numeric                     temperature,
                                                          const ScatteringSpeciesMap& ssprops) const;

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
