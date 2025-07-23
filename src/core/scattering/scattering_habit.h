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

  ScatteringHabit() = default;
  ScatteringHabit(const ParticleHabit &particle_habit_,
                  const PSD &psd_,
                  Numeric mass_size_rel_a_ = -1.0,
                  Numeric mass_size_rel_b_ = -1.0);

  BulkScatteringPropertiesTROGridded
  get_bulk_scattering_properties_tro_gridded(
      const AtmPoint&,
      const Vector& f_grid,
      const Numeric f_tol = 1e-3) const;

  ScatteringTroSpectralVector
  get_bulk_scattering_properties_tro_spectral(
      const AtmPoint&,
      const Vector& f_grid,
      const Index degree [[maybe_unused]]) const;

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
  FmtContext::iterator format(const scattering::ScatteringHabit&,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "ScatteringHabit"sv);
    }

    return tags.format(ctx);
  }
};

template <>
struct xml_io_stream<scattering::ScatteringHabit> {
  static constexpr std::string_view type_name = "ScatteringHabit";

  static void write(std::ostream&,
                    const scattering::ScatteringHabit&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);

  static void read(std::istream&,
                   scattering::ScatteringHabit&,
                   bifstream* = nullptr);
};
