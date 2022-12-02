/* Copyright (C) 2020 Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>

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
  \file   scattering_species.h
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-18

  \brief This file defines the ScatteringSpecies and ArrayOfScatteringSpecies
         classes.

 The ScatteringSpecies class provides a generic interface for all objects
 representing scattering data in ARTS. Together with the
 ArrayOfScatteringSpecies class it implements the high-level interface through
 which all other ARTS components should interact with scattering data.
*/
#pragma once

#include <numbers>

#include "array.h"
//#include "eigen.h"
#include "optproperties.h"
#include "jacobian.h"
#include "matpackI.h"
#include "rng.h"
#include "scattering.h"


////////////////////////////////////////////////////////////////////////////////
// ScatteringSpecies
////////////////////////////////////////////////////////////////////////////////
/** A generic scattering species.
 *
 * The ScatteringSpecies class provides a generic high-level interface for
 * the representation of scattering data in ARTS.
 *
 * A ScatteringSpecies object is a light-weight container that can hold any
 * type of scattering data. The ScatteringSpecies object manages the
 * lifetime of the scattering data. However, copy semantics are shallow,
 * which means that a copied ScatteringSpecies object shares the underlying
 * data with its source.
 */
class ScatteringSpecies {
 public:

  /// Create empty scattering species.
  ScatteringSpecies() {}

  /// Create container for existing scattering data pointer.
  ScatteringSpecies(std::shared_ptr<ScatteringSpeciesImpl> scattering_habit)
      : impl_(scattering_habit) {}

  /** Prepare scattering data to calculate expected properties.
   *
   * This function prepares scattering data for the calculation of specific
   * properties. It returns a new scattering species object which was
   * converted to allow faster calculation of the scattering properties.
   *
   * @param specs ScatteringPropertiesSpec object describing the expected
   * scattering properties.
   * @return A new scattering species object which was prepared to compute
   * the expected scattering properties.
   */
  ScatteringSpecies prepare_scattering_data(ScatteringPropertiesSpec specs) const {
    return ScatteringSpecies(impl_->prepare_scattering_data(specs));
  }

  /** Calculate bulk scattering properties.
   *
   * @param ws The ARTS workspace to used for calculation of the bulk
   * properties.
   * @param pbp_field The particle bulk properties field flattened into a matrix
   * with rows corresponding to the different atmospheric layers.
   * @param pbp_names The names of the particle bulk properties corresponding
   * to the columns in pbp_field
   * @param temperature The atmospheric temperatures flattened into a vector.
   * @jacobian_quantities The quantities for which to compute the jacobian.
   * @jacobian_do The JacobianDo flag.
   * @return A BulkScatteringProperties object containing the bulk scattering
   * properties.
   */
  BulkScatteringProperties calculate_bulk_properties(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      const Vector temperature,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
    return impl_->calculate_bulk_properties(ws,
                                            pbp_field,
                                            pbp_names,
                                            temperature,
                                            jacobian_quantities,
                                            jacobian_do);
  }

  std::pair<Matrix, Tensor3> get_absorption_and_extinction(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      Numeric frequency,
      const Vector temperature,
      Numeric lon_inc,
      Numeric lat_inc,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
      return impl_->get_absorption_and_extinction(ws,
                                                  pbp_field,
                                                  pbp_names,
                                                  frequency,
                                                  temperature,
                                                  lon_inc,
                                                  lat_inc,
                                                  jacobian_quantities,
                                                  jacobian_do);
  }

  std::pair<Vector, Matrix> sample_incoming_direction(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      Numeric frequency,
      const Vector temperature,
      ConstVectorView los_scat_rev,
      Rng& rng,
      Numeric scat_coeff_tot) const {
      return impl_->sample_incoming_direction(
          ws, pbp_field, pbp_names,
          frequency, temperature, los_scat_rev,
          rng, scat_coeff_tot
          );
  }

  bool operator==(const ScatteringSpecies &other) const {
      if (impl_) {
          return impl_.get() == other.impl_.get();
      }
      return false;
  }

  friend std::ostream& operator<<(std::ostream& out, const ScatteringSpecies&);

 private:
  std::shared_ptr<ScatteringSpeciesImpl> impl_ = nullptr;
};

////////////////////////////////////////////////////////////////////////////////
// ArrayOfScatteringSpecies
////////////////////////////////////////////////////////////////////////////////

/** Array of ScatteringSpecies
 *
 * This is simply a container class for multiple scattering species. It
 * essentially extends the functions required to calculate bulk scattering
 * properties to an array containing multiple ScatteringSpecies objects.
 */
class ArrayOfScatteringSpecies final : public Array<ScatteringSpecies> {
 public:

  ArrayOfScatteringSpecies() noexcept = default;
  ArrayOfScatteringSpecies(size_t n) : Array<ScatteringSpecies>(n) {}
  ArrayOfScatteringSpecies(Index n, const ScatteringSpecies value)
      : Array<ScatteringSpecies>(n, value) {}
  ArrayOfScatteringSpecies(const Array<ScatteringSpecies> &arr)
      : Array<ScatteringSpecies>(arr) {}
  ArrayOfScatteringSpecies(Array<ScatteringSpecies>&& arr)
      : Array<ScatteringSpecies>(std::move(arr)) {}

  /** Prepare scattering data to calculate expected properties.
   *
   * This function prepares scattering data for the calculation of specific
   * properties. It returns a new scattering species object which was
   * converted to allow faster calculation of the scattering properties.
   *
   * @param specs ScatteringPropertiesSpec object describing the expected
   * scattering properties.
   * @return A new scattering species object which was prepared to compute
   * the expected scattering properties.
   */
  ArrayOfScatteringSpecies prepare_scattering_data(
      ScatteringPropertiesSpec specs) const {
    ArrayOfScatteringSpecies result(size());
    for (size_t i = 0; i < size(); ++i) {
      result[i] = this->operator[](i).prepare_scattering_data(specs);
    }
    result.prepared_ = true;
    result.phase_function_norm_ = specs.phase_function_norm;
    return result;
  }

  /** Calculate bulk scattering properties.
   *
   *
   * @param ws The ARTS workspace to used for calculation of the bulk
   * properties.
   * @param pbp_field The particle bulk properties field flattened into a matrix
   * with rows corresponding to the different atmospheric layers.
   * @param pbp_names The names of the particle bulk properties corresponding
   * to the columns in pbp_field
   * @param temperature The atmospheric temperatures flattened into a vector.
   * @jacobian_quantities The quantities for which to compute the jacobian.
   * @jacobian_do The JacobianDo flag.
   * @return A BulkScatteringProperties object containing the bulk scattering
   * properties.
   */
  BulkScatteringProperties calculate_bulk_properties(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      const Vector temperature,
      const ArrayOfRetrievalQuantity& jacobian_quantities,
      bool jacobian_do) const {
    if (!prepared_) {
      std::runtime_error(
          "The scattering species must be prepared using "
          "'prepare_scattering' data before the bulk "
          " properties can be computed.");
    }
    auto result =
        this->operator[](0).calculate_bulk_properties(ws,
                                                      pbp_field,
                                                      pbp_names,
                                                      temperature,
                                                      jacobian_quantities,
                                                      jacobian_do);
    for (size_t i = 1; i < this->size(); ++i) {
      result +=
          this->operator[](i).calculate_bulk_properties(ws,
                                                        pbp_field,
                                                        pbp_names,
                                                        temperature,
                                                        jacobian_quantities,
                                                        jacobian_do);
    }
    result.normalize(phase_function_norm_);
    return result;
  }

  std::pair<Matrix, Tensor3> get_absorption_and_extinction(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      Numeric frequency,
      const Vector temperature,
      Numeric lon_inc,
      Numeric lat_inc,
      const ArrayOfRetrievalQuantity &jacobian_quantities,
      bool jacobian_do) const
  {
      if (size() == 0) {
          Matrix abs(pbp_field.nrows(), 4);
          abs = 0.0;
          Tensor3 ext(pbp_field.nrows(), 4, 3);
          ext = 0.0;
          return std::make_pair(abs, ext);
      }
      auto result = this->operator[](0).get_absorption_and_extinction(
          ws, pbp_field, pbp_names, frequency,
          temperature, lon_inc, lat_inc, jacobian_quantities, jacobian_do
          );
      for (size_t i = 1; i < this->size(); ++i) {
          auto result_i = this->operator[](i).get_absorption_and_extinction(
              ws, pbp_field, pbp_names,
              frequency, temperature, lon_inc, lat_inc,
              jacobian_quantities, jacobian_do
              );
          std::get<0>(result) += std::get<0>(result_i);
          std::get<1>(result) += std::get<1>(result_i);
      }
      return result;
  }

  /** Sample incoming scattering direction.
   *
   * Sample incoming direction and phase matrix at a given atmospheric point.
   *
   * @param ws The current workspace. Required for potential Agenda evaluation.
   * @param pbp_field: The particle bulkproperties field.
   * @param pbp_names: Names of the quantities in the PBP field.
   * @param frequency: The current frequency.
   * @param temperature: Vector containing atmospheric temperatures.
   * @param los_scat_reverse: The reversed outgoing direction of the
   * scattering event.
   * @param lat_inc: Latitude components of the reverse incoming direction in
   * degrees.
   * @param rng: Random generator object to use to generate random numbers.
   * @param scat_coeff_tot: The total scattering coefficient at the current
   * location.
   *
   * @return A pair containing the reverse scattering direction and the
   * corresponding phase matrix.
   *
   */
  std::pair<Vector, Matrix> sample_incoming_direction(
      Workspace& ws,
      const MatrixView pbp_field,
      const ArrayOfString pbp_names,
      Numeric frequency,
      const Vector temperature,
      ConstVectorView los_scat_rev,
      Rng& rng,
      Numeric scat_coeff_tot) const {

      auto lat_inc = Conversion::deg2rad(los_scat_rev[0]);
      auto lon_inc = Conversion::deg2rad(los_scat_rev[1]);

      if (size() == 0) {
          throw std::runtime_error(
              "There are no scattering species, so there's no way to sample"
              " a scattering"
              );
      } else if (size() == 1) {
          return this->operator[](0).sample_incoming_direction(
              ws,
              pbp_field,
              pbp_names,
              frequency,
              temperature,
              los_scat_rev,
              rng,
              scat_coeff_tot);
      }

      double r = rng.draw() * scat_coeff_tot;
      double scat_coeff_acc = 0.0;

      // If there are several scatterers we MC integrate over them.
      // Note that this requires scaling of the scattering matrix.
      for (size_t i = 0; i < size(); ++i) {
          auto s = this->operator[](0);
          double abs, ext;
          auto abs_ext = s.get_absorption_and_extinction(
              ws,
              pbp_field,
              pbp_names,
              frequency,
              temperature,
              lon_inc,
              lat_inc,
              {},
              false);
          abs = std::get<0>(abs_ext)(0, 0);
          ext = std::get<1>(abs_ext)(0, 0, 0);
          scat_coeff_acc += ext - abs;
          if (r < scat_coeff_acc) {
              auto result = s.sample_incoming_direction(
                  ws,
                  pbp_field,
                  pbp_names,
                  frequency,
                  temperature,
                  los_scat_rev,
                  rng,
                  scat_coeff_tot);
              std::get<1>(result) *= (scat_coeff_tot / (ext - abs));
              return result;
          }
      }
  }

private:

  bool prepared_ = false;
  Numeric phase_function_norm_ = 4.0 * pi_v<Numeric>;

};

void Append(  // WS Generic Output:
    ArrayOfScatteringSpecies& out,
    const String& out_name,
    const ArrayOfScatteringSpecies& in,
    const String& direction,
    const String& in_name,
    const String& direction_name,
    const Verbosity& verbosity);

void Append(  // WS Generic Output:
    ArrayOfScatteringSpecies& out,
    const String& out_name,
    const ScatteringSpecies& in,
    const String& direction,
    const String& in_name,
    const String& direction_name,
    const Verbosity& verbosity);

void Select(  // WS Generic Output:
    ArrayOfScatteringSpecies& needles,
    // WS Generic Input:
    const ArrayOfScatteringSpecies& haystack,
    const ArrayOfIndex& needleind,
    const Verbosity& verbosity);
