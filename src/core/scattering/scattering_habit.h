#pragma once

#include "psd.h"
#include "particle_habit.h"
#include "bulk_scattering_properties.h"

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
  ScatteringHabit(ParticleHabit particle_habit_,
                  Numeric scat_species_a_,
                  Numeric scat_species_b_,
                  PSD psd_)
    : particle_habit(particle_habit_), scat_species_a(scat_species_a_), scat_species_b(scat_species_b_), psd(psd_) {}

  BulkScatteringPropertiesTROGridded
  get_bulk_scattering_properties_tro_gridded(
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
  Numeric scat_species_a, scat_species_b;
  PSD psd;
};


}
