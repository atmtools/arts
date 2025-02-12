#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <matpack.h>

#include <cmath>
#include <format>
#include <variant>

#include "bulk_scattering_properties.h"
#include "henyey_greenstein.h"
#include "particle_habit.h"
#include "properties.h"
#include "general_tro_spectral.h"
#include "scattering_habit.h"

namespace scattering {


struct ScatteringDataSpec {};

using Species = std::variant<HenyeyGreensteinScatterer, ScatteringGeneralSpectralTRO>;

}  // namespace scattering

using ScatteringSpecies = scattering::Species;

template <scattering::Format format, scattering::Representation repr>
using BulkScatteringProperties =
    scattering::BulkScatteringProperties<format, repr>;

/** Array of scattering species
  *
  * The array of scattering species holds atmospheric quantities that scatter
  * and provides convenience functions to calculate their bulk scattering
  * properties.
  */
struct ArrayOfScatteringSpecies : public std::vector<scattering::Species> {
  void add(const scattering::Species& species) { push_back(species); }
  void prepare_scattering_data(scattering::ScatteringDataSpec) {}

  [[nodiscard]] BulkScatteringProperties<scattering::Format::TRO,
                                         scattering::Representation::Gridded>
  get_bulk_scattering_properties_tro_gridded(
      const AtmPoint& atm_point,
      const Vector& f_grid,
      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;

  [[nodiscard]] ScatteringTroSpectralVector
  get_bulk_scattering_properties_tro_spectral(const AtmPoint& atm_point,
                                              const Vector& f_grid,
                                              Index degree) const;

  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO,
                                         scattering::Representation::Gridded>
  get_bulk_scattering_properties_aro_gridded(
      const AtmPoint& atm_point,
      const Vector& f_grid,
      const Vector& za_inc_grid,
      const Vector& delta_aa_grid,
      std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const;

  [[nodiscard]] BulkScatteringProperties<scattering::Format::ARO,
                                         scattering::Representation::Spectral>
  get_bulk_scattering_properties_aro_spectral(const AtmPoint& atm_point,
                                              const Vector& f_grid,
                                              const Vector& za_inc_grid,
                                              Index degree,
                                              Index order) const;
};

inline std::ostream& operator<<(std::ostream& os,
                                const ArrayOfScatteringSpecies& /*species*/) {
  os << "An array of scattering species." << '\n';
  return os;
}

template <>
struct std::formatter<ArrayOfScatteringSpecies> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ArrayOfScatteringSpecies&,
                              FmtContext& ctx) const {
    return ctx.out();
  }
};

using HenyeyGreensteinScatterer = scattering::HenyeyGreensteinScatterer;
using ParticleHabit             = scattering::ParticleHabit;
using ScatteringHabit           = scattering::ScatteringHabit;
using PSD                       = scattering::PSD;

std::ostream& operator<<(
    std::ostream& os,
    const std::variant<HenyeyGreensteinScatterer,
                       scattering::ScatteringHabit>& /*species*/);
