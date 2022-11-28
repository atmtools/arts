/* Copyright (C) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmer.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, rite to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   bulk_particle_habit.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-15

  \brief Defines the BulkParticleHabit class, which combines a particle habit
  with a PND.

  The BulkParticleHabit class combines a particle habit, i.e. a collection of
  scattering particles, with a PND. This combination allows mapping bulk
  properties to scattering data and this is what is required to calculate the
  input for the RT routines.

  The BulkParticleHabit implements the generic interface for scattering defined
  by the ScatteringSpecies class.
*/
#pragma once

#include <iostream>
#include <memory>

#include "agenda_class.h"
#include "optproperties.h"
#include "scattering/maths.h"
#include "scattering.h"
#include "scattering/particle_habit.h"

using scattering::to_arts;

////////////////////////////////////////////////////////////////////////////////
// ScatteringParticle and ArrayOfScatteringParticle
////////////////////////////////////////////////////////////////////////////////

typedef scattering::Particle ScatteringParticle;
typedef Array<ScatteringParticle> ArrayOfScatteringParticle;

std::ostream &operator<<(std::ostream &out, const ScatteringParticle &);
std::ostream &operator<<(std::ostream &out, const ArrayOfScatteringParticle &);

////////////////////////////////////////////////////////////////////////////////
// Conversion from legacy data.
////////////////////////////////////////////////////////////////////////////////

/** Convert from legacy format to new scattering format.
 *
 * @param legacy_data The scattering data in legacy ARTS format to convert
 * to the new format.
 * @return SingleScatteringData object containing the single scattering data.
 */
namespace detail {
scattering::SingleScatteringData from_legacy_format(
    const SingleScatteringData &legacy_data);
}

////////////////////////////////////////////////////////////////////////////////
// BulkParticleHabit
////////////////////////////////////////////////////////////////////////////////
/** BulkParticleHabit
 *
 * A scattering habit describes an atmospheric scattering species through  an
 * explicit particle model of scattering properties at given particle sizes and
 * a parametrized PSD that describes the distribution of those sizes throughout
 * the atmosphere.
 *
 * This class implements the abstract interface for scattering species (defined
 * by ScatteringSpeciesImpl) and instances of this class are used as a elements
 * of the *scattering_species* WSM containing the atmosphere's scattering
 * species.
 *
 */
class BulkParticleHabit : public ScatteringSpeciesImpl {

    // Extracts agenda input from particle bulkprop field.
    Matrix get_agenda_input(Matrix pbp_field,
                            ArrayOfString pbf_names) const;

    // Get names of agenda input for which to calculate Jacobians.
    ArrayOfString get_dpnd_data_dx_names(ArrayOfRetrievalQuantity jacobian_quantities,
                                         bool jacobian_do) const;

 public:
  BulkParticleHabit();
  /** Create BulkParticleHabit.
   * @param name The name of the scattering species.
   * @param scat_data The ensemble scattering data describing the particles
   * which make up the habit.
   * @param meta_data The particle meta data corresponding to the particles in
   * scat_data
   * @pnd_agenda The agenda to use to calculate the number densities of each
   * particle.
   * @pnd_agenda_input The names of the properties from pbp_field to use as input
   * for the particle size distributions.
   */
  BulkParticleHabit(const String &name,
                  const ArrayOfSingleScatteringData &scat_data,
                  const ArrayOfScatteringMetaData &meta_data,
                  const Agenda &pnd_agenda,
                  const ArrayOfString &pnd_agenda_input);

  /** Create BulkParticleHabit.
   * @param name The name of the scattering species.
   * @pnd_agenda The agenda to use to calculate the number densities of each
   * particle.
   * @param pnd_agenda_input The names of the properties from pbp_field to use as input
   * for the particle size distributions.
   * @param particle_habit Particle model object containing the scattering data
   * describing the particles in the habit.
   */
  BulkParticleHabit(const String name,
                  const Agenda &pnd_agenda,
                  ArrayOfString pnd_agenda_input,
                  std::shared_ptr<scattering::ParticleHabit> particle_habit);

  /** Create BulkParticleHabit.
   * @param name The name of the scattering species.
   * @param scat_data The ensemble scattering data describing the particles
   * which make up the habit.
   * @param meta_data The particle meta data corresponding to the particles in
   * scat_data
   * @param index_start Start index of the books in particle_bulkprop_field
   * @param index_end Index pointing behind the first book after 'start_index'
   */
  BulkParticleHabit(const String &name,
                  const ArrayOfSingleScatteringData &scat_data,
                  Index index_start,
                  Index index_end);

  /** Create a BulkParticleHabit with a fixed particle number density.
   * @param name The name of the scattering species.
   * @param index_start Start index of the books in particle_bulkprop_field
   * @param index_end Index pointing behind the first book after 'start_index'
   * that doesn't contain any of the habits PND values.
   * @particle_habit Particle model object containing the scattering data
   * describing the particles in the habit.
   */
  BulkParticleHabit(const String name,
                  Index index_start,
                  Index index_end,
                  std::shared_ptr<scattering::ParticleHabit> particle_habit);

  BulkParticleHabit(const BulkParticleHabit &) = default;
  BulkParticleHabit &operator=(const BulkParticleHabit &) = default;
  BulkParticleHabit &operator=(BulkParticleHabit &&) = default;

  ~BulkParticleHabit();

  /// The masses of the particles in the habit.
  Vector get_particle_mass() const {
    return to_arts(particle_habit_->get_mass());
  }

  /// The volume equivalent diameters of the particles in the habit.
  Vector get_particle_d_eq() const {
    return to_arts(particle_habit_->get_d_eq());
  }

  /// The maximum diameters of the particles in the habit.
  Vector get_particle_d_max() const {
    return to_arts(particle_habit_->get_d_max());
  }

  Vector get_particle_area() const {
    return to_arts(particle_habit_->get_d_max());
  }

  /// Extra meta data and conver to ARTS legacy format.
  ArrayOfScatteringMetaData get_meta_data() const {
      auto extract_meta = [](const scattering::Particle &particle){
          return ScatteringMetaData{
              particle.get_name(),
              particle.get_source(),
              particle.get_refractive_index(),
              particle.get_mass(),
              particle.get_d_max(),
              particle.get_d_eq(),
              particle.get_d_aero()
              };
      };
      ArrayOfScatteringMetaData result{};
      std::transform(
          particle_habit_->get_particles().begin(),
          particle_habit_->get_particles().end(),
          std::back_inserter(result),
          extract_meta
          );
      return result;
  }

  /** Calculate bulk scattering properties for 1D atmosphere.
   *
   * @param pbp_field Matrix view containing the particle bulk properties for
   * all layers in the atmosphere.
   * @param pbp_names The names corresponding to the rows in pbp_field.
   * @param temperature Vector containing the temperatures in the atmosphere.
   * @param jacobian_quantities The quantities for which to compute the jacobian.
   * @param jacobian_do Whether or not to calculate the Jacobian.
   */
  BulkScatteringProperties calculate_bulk_properties(Workspace &ws,
                                                     ConstMatrixView pbp_field,
                                                     const ArrayOfString &pbp_names,
                                                     ConstVectorView temperature,
                                                     const ArrayOfRetrievalQuantity &jacobian_quantities,
                                                     bool jacobian_do) const;


  std::shared_ptr<ScatteringSpeciesImpl> prepare_scattering_data(ScatteringPropertiesSpec specs) const;

  std::pair<Matrix, Tensor3> get_absorption_and_extinction(
      Workspace &ws,
      ConstMatrixView pbp_field,
      const ArrayOfString &pbp_names,
      Numeric frequency,
      ConstVectorView temperature,
      Numeric lon_inc,
      Numeric lat_inc,
      const ArrayOfRetrievalQuantity &jacobian_quantities,
      bool jacobian_do) const;

  std::pair<Vector, Matrix> sample_incoming_direction(
      Workspace& ws,
      ConstMatrixView pbp_field,
      const ArrayOfString &pbp_names,
      Numeric frequency,
      ConstVectorView temperature,
      ConstVectorView los_scat_rev,
      Rng& rng,
      Numeric scat_coeff_tot) const;

  friend std::ostream &operator<<(std::ostream &out, const BulkParticleHabit &);


 private:
  String name_;
  std::shared_ptr<Agenda> pnd_agenda_;
  ArrayOfString pnd_agenda_input_;
  std::shared_ptr<scattering::ParticleHabit> particle_habit_;
  Index index_start_ = 0;
  Index index_end_ = 0;
};
