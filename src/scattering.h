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
#include "optproperties.h"
#include "jacobian.h"
#include "matpackI.h"
#include "rng.h"
#include "scattering/single_scattering_data.h"

using std::numbers::pi_v;

class Workspace;

enum class Format {Gridded, Spectral};
enum class ReferenceFrame {ScatteringPlane, Lab};

/** Specification of required scattering properties.
 *
 * This struct holds requirements describing the scattering
 * data required by different scattering solvers.
 */
struct ScatteringPropertiesSpec {
    ScatteringPropertiesSpec(const Vector& f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             Index l_max_,
                             Index m_max_=0,
                             Numeric phase_function_norm_=1.0);

    ScatteringPropertiesSpec(scattering::math::ConstVectorPtr<Numeric> f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             scattering::math::ConstVectorPtr<Numeric> lon_scat_,
                             scattering::ConstLatitudeGridPtr<Numeric> lat_scat_,
                             Numeric phase_function_norm=1.0);
    ScatteringPropertiesSpec(const Vector& f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             const Vector& lon_scat_,
                             const Vector& lat_scat_,
                             Numeric phase_function_norm_=1.0)
        : ScatteringPropertiesSpec(std::make_shared<scattering::math::Vector<Numeric>>(scattering::to_eigen(f_grid)),
                               frame,
                               n_stokes_,
                                   std::make_shared<scattering::math::Vector<Numeric>>(scattering::to_eigen(lon_scat_)),
                                   std::make_shared<scattering::IrregularLatitudeGrid<Numeric>>(scattering::to_eigen(lat_scat_)),
                               phase_function_norm_) {}
    ScatteringPropertiesSpec(scattering::math::ConstVectorPtr<Numeric> f_grid,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             scattering::math::ConstVectorPtr<Numeric> lat_inc_,
                             scattering::math::ConstVectorPtr<Numeric> lon_scat_,
                             scattering::ConstLatitudeGridPtr<Numeric> lat_scat_,
                             Numeric phase_function_norm_=1.0);
    ScatteringPropertiesSpec(const Vector& f_grid_,
                             ReferenceFrame frame_,
                             Index n_stokes_,
                             const Vector& lat_inc_,
                             const Vector& lon_scat_,
                             const Vector& lat_scat_,
                             Numeric phase_function_norm_=1.0)
        : ScatteringPropertiesSpec(std::make_shared<scattering::math::Vector<Numeric>>(scattering::to_eigen(f_grid_)),
                               frame_,
                               n_stokes_,
                                   std::make_shared<scattering::math::Vector<Numeric>>(scattering::to_eigen(lat_inc_)),
                                   std::make_shared<scattering::math::Vector<Numeric>>(scattering::to_eigen(lon_scat_)),
                                   std::make_shared<scattering::IrregularLatitudeGrid<Numeric>>(scattering::to_eigen(lat_scat_)),
                               phase_function_norm_) {}

    Format format;
    ReferenceFrame frame;
    Index n_stokes = 0;
    Index l_max = 0;
    Index m_max = 0;
    Numeric phase_function_norm = 1.0;
    scattering::math::ConstVectorPtr<Numeric> lon_inc = nullptr;
    scattering::math::ConstVectorPtr<Numeric> lat_inc = nullptr;
    scattering::math::ConstVectorPtr<Numeric> lon_scat = nullptr;
    scattering::ConstLatitudeGridPtr<Numeric> lat_scat = nullptr;
    scattering::math::ConstVectorPtr<Numeric> f_grid = nullptr;
    Index n_angs_frame_conversion = 32;
};


/** BulkScatteringProperties
 *
 * The class represents the bulk scattering properties of the atmosphere. It is the
 * result of calculating the scattering properties all points along a propagation
 * path. The scattering data for each point in the atmospheric grid is represent by a
 * SingleScatteringData object. The data for the whole grid is stored as a
 * 1-dimensional array.
 */
class BulkScatteringProperties {
 public:
  /** Create BulkScatteringProperties
     *
     * @param Array of scattering::SingleScatteringData objects containing
     * the bulk scattering properties at all points along the propagation
     * path.
     */
  BulkScatteringProperties(Array<scattering::SingleScatteringData> data)
      : data_(data),
        n_freqs_(data[0].get_f_grid().size()),
        stokes_dim_(data[0].get_stokes_dim()) {}

  /** Extracts extinction coefficients from bulk properties.
   *
   * @return Tensor4 with dimensions [n_points, n_freqs, n_lon_inc, n_lat_inc]
   * containing the extinction coefficient for all points along the propagation
   * path, frequencies in f_grid, and incoming azimuth and zenith angles.
   */
  Tensor4 get_extinction_coeff() const;

  /** Extracts extinction matrices from bulk properties.
   *
   * @return Tensor6 with dimensions
   * [n_points, n_freqs, n_lon_inc, n_lat_inc, n_stokes, n_stokes]
   * containing the extinction matrices for all points along the propagation
   * path, frequencies in f_grid, and incoming azimuth and zenith angles.
   */
  Tensor6 get_extinction_matrix() const;
  Tensor6 get_extinction_matrix(Index stokes_dim) const;

  /** Extracts absorption coefficients from bulk properties.
   *
   * @return Tensor5 with dimensions [n_points, n_freqs, n_lon_inc]
   * containing the absorption coefficient for all points along the
   * propagation path, frequencies in f_grid, and incoming azimuth
   * and zenith angles.
   */
  Tensor4 get_absorption_coeff() const;

  /** Extract absorption vector from bulk properties.
   *
   * @return Tensor5 with dimensions [n_points, n_freqs, n_lon_inc, stokes]
   * containing the absorption vector for all points along the propagation
   * path, frequencies in f_grid, incoming azimuth and zenith angles and the
   * stokes components.
   */
  Tensor5 get_absorption_vector() const;
  Tensor5 get_absorption_vector(Index stokes_dim) const;

  /** Get spectral coefficients.
   *
   * Return spectral components of SHT-transformed scattering matrix computed
   * using orthonormal basis functions.
   *
   * @return Tensor5 containing Legendre coefficients along
   * columns, incoming scattering angles along rows and pages,
   * propagation-path points along paragraphs and frequencies 
   * along books.
   */
  Tensor5 get_spectral_coeffs() const;

  /** Get Legendre coefficients of phase function.
   * @return Tensor5 containing Legendre coefficients along
   * columns, incoming scattering angles along rows and pages,
   * propagation-path points along paragraphs and frequencies 
   * along books.
   */
  Tensor5 get_legendre_coeffs() const;

  /** Return phase function.
   * @return Tensor3 containing the phase functions elements along
   * column, the atmospheric layers along rows and the frequencies along
   * the pages.
   */
  Tensor3 get_phase_function() const;

  /** Return full phase matrix.
   *
   * Returns the full phase matrix for all points along the propagation paths.
   * @return Tensor7 containing the full phase matrix with frequencies along
   * the first axis, propagation path points along the second, incoming zenith,
   * scattering azimuth and scattering zenith along axis 3, 4, 5, respectively.
   * The last two axes correspond to the rows and columns of the phase matrix.
   */
  Tensor7 get_phase_matrix() const;
  /** Return full phase matrix for a given stokes dimension.
   *
   * Returns the full phase matrix for all points along the propagation paths.
   * @return Tensor7 containing the full phase matrix with frequencies along
   * the first axis, propagation path points along the second, incoming zenith,
   * scattering azimuth and scattering zenith along axis 3, 4, 5, respectively.
   * The last two axes correspond to the rows and columns of the phase matrix.
   */
  Tensor7 get_phase_matrix(Index stokes_dim) const;

  BulkScatteringProperties &operator+=(const BulkScatteringProperties &other) {
    for (Index i = 0; i < data_.size(); ++i) {
      data_[i] += other.data_[i];
    }
    return *this;
  }

  /** Normalize phase matrix data.
   *
   * Normalizes the scattering-angle integral of phase matrix to the
   * given value.
   *
   * @param norm The desired integral of the phase matrix.
   */
  void normalize(Numeric norm) {
    for (auto &d : data_) {
      d.normalize(norm);
    }
  }

  /** Downsample phase matrix along scattering azimuth angles.
   *
   * This function downsamples the phase matrix to a coarser grid by
   * averaging it along the scattering azimuth grid.
   *
   * @param lon_scat: The coarser grid to which to downsample the
   * phase matrix.
   */
  void downsample_lon_scat(const Vector &lon_scat) {
      auto lon_scat_ptr = std::make_shared<scattering::math::Vector<Numeric>>(scattering::to_eigen(lon_scat));
    for (Index i = 0; i < data_.size(); ++i) {
      data_[i] = data_[i].downsample_lon_scat(lon_scat_ptr);
    }
  }

 private:
  Array<scattering::SingleScatteringData> data_;
  Index n_freqs_;
  Index stokes_dim_;
};

////////////////////////////////////////////////////////////////////////////////
// Abstract interface
////////////////////////////////////////////////////////////////////////////////

/** Abstract interface for scattering species representations.
 *
 * This abstract class defines the generic interface that all 
 *
 */
class ScatteringSpeciesImpl {
public:
    virtual ~ScatteringSpeciesImpl(){};
    ScatteringSpeciesImpl(){};
    ScatteringSpeciesImpl(const ScatteringSpeciesImpl &) = default;

    virtual std::shared_ptr<ScatteringSpeciesImpl> prepare_scattering_data(
        ScatteringPropertiesSpec specs) const = 0;
    virtual BulkScatteringProperties calculate_bulk_properties(
        Workspace &ws,
        ConstMatrixView pbp_field,
        const ArrayOfString &pbf_names,
        ConstVectorView temperature,
        const ArrayOfRetrievalQuantity &jacobian_quantities,
        bool jacobian_do) const = 0;

    virtual std::pair<Matrix, Tensor3> get_absorption_and_extinction(
        Workspace &ws,
        ConstMatrixView pbp_field,
        const ArrayOfString &pbf_names,
        Numeric frequency,
        ConstVectorView temperature,
        Numeric lon_inc,
        Numeric lat_inc,
        const ArrayOfRetrievalQuantity &jacobian_quantities,
        bool jacobian_do
        ) const = 0;

    virtual std::pair<Vector, Matrix> sample_incoming_direction(
        Workspace &ws,
        ConstMatrixView pbp_field,
        const ArrayOfString &pbf_names,
        Numeric frequency,
        ConstVectorView temperature,
        ConstVectorView los_scat_rev,
        Rng& rng,
        Numeric scat_coeff_tot
        ) const = 0;
};

