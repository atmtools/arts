#pragma once

#include "psd.h"
#include "particle_habit.h"
#include "bulk_scattering_properties.h"
#include "general_tro_spectral.h"

namespace scattering {

using PSD = std::variant<MGDSingleMoment, BinnedPSD>;

/*** A scattering habit
 *
 * A scattering habit combines a particle habit with an additional PSD
 * and thus defines a mapping between atmospheric scattering species properties
 * and corresponding bulk skattering properties.
 */
class ScatteringHabit {
 public:

  ScatteringHabit() {};
  ScatteringHabit(const ParticleHabit &particle_habit_,
                  const PSD &psd_,
                  Numeric mass_size_rel_a_ = -1.0,
                  Numeric mass_size_rel_b_ = -1.0)
    : particle_habit(particle_habit_), mass_size_rel_a(mass_size_rel_a_), mass_size_rel_b(mass_size_rel_b_), psd(psd_) {
    if ((mass_size_rel_a < 0.0) || (mass_size_rel_b < 0.0)) {
      auto size_param = std::visit([](auto const& psd){return psd.get_size_parameter();}, psd);
      auto [sizes, mass_size_rel_a, mass_size_rel_b] = particle_habit.get_size_mass_info(size_param);
    }
  }

  BulkScatteringPropertiesTROGridded
  get_bulk_scattering_properties_tro_gridded(
      const AtmPoint&,
      const Vector& f_grid,
      const Numeric f_tol = 1e-3) const;

  ScatteringTroSpectralVector
  get_bulk_scattering_properties_tro_spectral(
      const AtmPoint&,
      const Vector& f_grid,
      const Numeric f_tol = 1e-3) const;

//  BulkScatteringProperties<Format::TRO, Representation::Gridded>
//  get_bulk_scattering_properties_tro_spectral(
//      const AtmPoint&,
//      const Vector& f_grid,
//      Index l) const;

 private:
  ParticleHabit particle_habit;
  Numeric mass_size_rel_a, mass_size_rel_b;
  PSD psd;
};


}

template <>
struct std::formatter<scattering::ScatteringHabit> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const scattering::ScatteringHabit& v,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "ScatteringHabit"sv);
    }

    return tags.format(ctx);
  }
};
