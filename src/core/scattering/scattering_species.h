#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <matpack.h>

#include <cmath>
#include <format>
#include <variant>

#include "bulk_scattering_properties.h"
#include "callback_scattering_species.h"
#include "gas_scattering.h"
#include "general_tro_spectral.h"
#include "henyey_greenstein.h"
#include "particle_habit.h"
#include "properties.h"
#include "psd.h"
#include "scattering_habit.h"

namespace scattering {

struct ScatteringDataSpec {};

using Species = std::variant<GasScatterer, HenyeyGreensteinScatterer, CallbackScatteringSpecies, ScatteringHabit>;

}  // namespace scattering

using ScatteringSpecies = scattering::Species;

template <scattering::Format format, scattering::Representation repr> using BulkScatteringProperties =
    scattering::BulkScatteringProperties<format, repr>;

/** Array of scattering species
  *
  * The array of scattering species holds quantities that scatter and provides
  * convenience functions to calculate their bulk scattering properties at the
  * point of a medium and its spectral point, AtmPoint with AtmSpectralPoint or
  * SubsurfacePoint with SubsurfaceSpectralPoint.
  */
struct ArrayOfScatteringSpecies {
  std::vector<scattering::Species> species;

  void add(const scattering::Species& species);
  void prepare_scattering_data(scattering::ScatteringDataSpec);

  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded> tro_gridded(
      const Point&                                 point,
      const SpectralPoint&                         spectral_point,
      const Vector&                                f_grid,
      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;

  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<scattering::Format::TRO, scattering::Representation::Gridded>
  tro_gridded_derivative(const Point&                                 point,
                         const SpectralPoint&                         spectral_point,
                         const Vector&                                f_grid,
                         std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
                         const scattering::medium_key<Point>&         target) const;

  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(const Point&         point,
                                                         const SpectralPoint& spectral_point,
                                                         const Vector&        f_grid,
                                                         Index                degree) const;

  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded> aro_gridded(
      const Point&                                 point,
      const SpectralPoint&                         spectral_point,
      const Vector&                                f_grid,
      const Vector&                                za_inc_grid,
      const Vector&                                delta_aa_grid,
      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;

  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
  aro_gridded_derivative(const Point&                                 point,
                         const SpectralPoint&                         spectral_point,
                         const Vector&                                f_grid,
                         const Vector&                                za_inc_grid,
                         const Vector&                                delta_aa_grid,
                         std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
                         const scattering::medium_key<Point>&         target) const;

  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral> aro_spectral(
      const Point&         point,
      const SpectralPoint& spectral_point,
      const Vector&        f_grid,
      const Vector&        za_inc_grid,
      Index                degree,
      Index                order) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame bulk phase matrix (see
   * PhaseMatrixData<ARO, Fourier>), at exactly the given zenith angles [deg], summed over the species
   */
  template <scattering::medium_point Point, scattering::spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier> aro_fourier(
      const Point&         point,
      const SpectralPoint& spectral_point,
      const Vector&        f_grid,
      const Vector&        za_inc_grid,
      const Vector&        za_scat_grid,
      Index                max_mode) const;
};

template <> struct std::formatter<ArrayOfScatteringSpecies> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const ArrayOfScatteringSpecies&, FmtContext& ctx) const {
    return ctx.out();
  }
};

using HenyeyGreensteinScatterer = scattering::HenyeyGreensteinScatterer;
using GasScatterer              = scattering::GasScatterer;
using ParticleHabit             = scattering::ParticleHabit;
using ScatteringHabit           = scattering::ScatteringHabit;
using PSD                       = scattering::PSD;

//! Naming the variant
template <> struct xml_io_stream_name<ScatteringSpecies> {
  static constexpr std::string_view name = "ScatteringSpecies"sv;
};

//! Required overload since it inherits
template <> struct xml_io_stream_name<ArrayOfScatteringSpecies> {
  static constexpr std::string_view name = "ArrayOfScatteringSpecies"sv;
};

template <> struct xml_io_stream_aggregate<ArrayOfScatteringSpecies> {
  static constexpr bool value = true;
};
