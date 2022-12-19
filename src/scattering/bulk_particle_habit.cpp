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
  \file   scattering_habit.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-15

  \brief  Implementation of scattering_habit.h
*/
#include <iomanip>
#include "rte.h"
#include "auto_md.h"
#include "arts_conversions.h"
#include "scattering/bulk_particle_habit.h"

const String SCATSPECIES_MAINTAG = "Scattering species";

namespace detail {

using scattering::EigenTensor;

/** Extract backward scattering coefficient from phase matrix.
 *
 * @param phase_matrix Phase matrix in scattering format.
 * @return Eigen tensor containing the backscattering coefficient in
 *     scattering-compatible format.
 */
EigenTensor<7> extract_backward_scattering_coeff(
    const EigenTensor<7> &phase_matrix) {
  auto dimensions = phase_matrix.dimensions();
  auto backward_scattering_coeff =
      phase_matrix.chip<4>(0).chip<4>(dimensions[5] - 1).chip<4>(0);
  dimensions[4] = 1;
  dimensions[5] = 1;
  dimensions[6] = 1;
  return backward_scattering_coeff.reshape(dimensions);
}

/** Extract forward scattering coefficient from phase matrix.
 *
 * @param phase_matrix Phase matrix in scattering format.
 * @return Eigen tensor containing the forward scattering coefficient in
 *     scattering-compatible format.
 */
EigenTensor<7> extract_forward_scattering_coeff(
    const EigenTensor<7> &phase_matrix) {
  auto dimensions = phase_matrix.dimensions();
  auto backward_scattering_coeff =
      phase_matrix.chip<4>(0).chip<4>(0).chip<4>(0);
  dimensions[4] = 1;
  dimensions[5] = 1;
  dimensions[6] = 1;
  return backward_scattering_coeff.reshape(dimensions);
}

/** Convert ARTS to scattering format.
 *
 * This function converts ARTS phase matrix data to scattering format.
 *
 * @param tensor ARTS phase matrix data.
 * @return Eigen tensor containing the phase matrix in scattering-compatible format.
 */
EigenTensor<7> from_legacy_format(const Tensor7 &tensor) {
    EigenTensor<7> tensor_eigen = scattering::to_eigen(tensor);
    std::array<Eigen::Index, 7> shuffle_dimensions = {0, 1, 5, 4, 3, 2, 6};
    EigenTensor<7> result = tensor_eigen.shuffle(shuffle_dimensions);
    return tensor_eigen.shuffle(shuffle_dimensions);
}

/** Convert ARTS to scattering format.
 *
 * This function converts a ARTS extinction matrix or absorption
 * vector data to scattering format
 *
 * @param tensor The input data to convert.
 * @return Eigen tensor containing the data in scattering-compatible format.
 */
EigenTensor<7> from_legacy_format(const Tensor5 &tensor) {
    EigenTensor<5> tensor_eigen = scattering::to_eigen(tensor);
  auto result = scattering::math::unsqueeze<4, 5>(tensor_eigen);
  std::array<Eigen::Index, 7> shuffle_dimensions = {0, 1, 3, 2, 4, 5, 6};
  return result.shuffle(shuffle_dimensions);
}

/** Convert ARTS legacy scattering data to new scattering format.
 *
 * @param arts_data SingleScatteringData object containing the scattering
 * data in legacy format.
 */
scattering::SingleScatteringData from_legacy_format(
    const SingleScatteringData &arts_data) {
    scattering::math::Vector<Numeric> f_grid = scattering::to_eigen(arts_data.f_grid);
    scattering::math::Vector<Numeric> t_grid = scattering::to_eigen(arts_data.T_grid);
    scattering::math::Vector<Numeric> lon_inc = scattering::math::Vector<Numeric>::Constant(1, 0.0);
    scattering::math::Vector<Numeric> lat_inc = scattering::math::Vector<Numeric>::Constant(1, 0.0);
    scattering::math::Vector<Numeric> lon_scat = scattering::to_eigen(arts_data.aa_grid) * M_PI / 180.0;
  if (lon_scat.size() == 0) {
      lon_scat = lon_inc;
  }
  scattering::math::Vector<Numeric> lat_scat = scattering::to_eigen(arts_data.za_grid);
  lat_scat *= (M_PI / 180.0);
  if (arts_data.ptype > PTYPE_TOTAL_RND) {
      lat_inc = lat_scat;
  }
  auto phase_matrix = from_legacy_format(arts_data.pha_mat_data);
  auto extinction_matrix = from_legacy_format(arts_data.ext_mat_data);
  auto absorption_vector = from_legacy_format(arts_data.abs_vec_data);
  auto backward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);
  auto forward_scattering_coeff =
      extract_backward_scattering_coeff(phase_matrix);

  auto dimensions = phase_matrix.dimensions();

  assert(phase_matrix.dimension(0) == f_grid.size());
  assert(phase_matrix.dimension(1) == t_grid.size());
  assert(phase_matrix.dimension(4) == lon_scat.size());
  assert(phase_matrix.dimension(5) == lat_scat.size());

  return scattering::SingleScatteringData(
      f_grid, t_grid,
      lon_inc, lat_inc,
      lon_scat, lat_scat,
      phase_matrix,
      extinction_matrix,
      absorption_vector,
      backward_scattering_coeff,
      forward_scattering_coeff);
}

}  // namespace detail


////////////////////////////////////////////////////////////////////////////////
// BulkParticleHabit
////////////////////////////////////////////////////////////////////////////////


BulkParticleHabit::BulkParticleHabit() {}
BulkParticleHabit::~BulkParticleHabit() {}

BulkParticleHabit::BulkParticleHabit(
    const String &name,
    const ArrayOfSingleScatteringData &arts_scat_data,
    const ArrayOfScatteringMetaData &meta_data,
    const Agenda &pnd_agenda,
    const ArrayOfString &pnd_agenda_input)
    : name_(name),
      pnd_agenda_(std::make_shared<Agenda>(pnd_agenda)),
      pnd_agenda_input_(pnd_agenda_input)
{
  Index n = arts_scat_data.size();
  Index m = meta_data.size();
  std::vector<scattering::Particle> particles{};
  scattering::math::Vector<Numeric> d_eq(n), d_max(n), mass(n);

  for (size_t i = 0; i < n; ++i) {

    auto scattering_data = detail::from_legacy_format(arts_scat_data[i]);
    auto meta = meta_data[i];
    particles.emplace_back(meta.mass,
                           meta.diameter_volume_equ,
                           meta.diameter_max,
                           scattering_data);
  }
  particle_habit_ = std::make_shared<scattering::ParticleHabit>(particles);
}

BulkParticleHabit::BulkParticleHabit(
    const String &name,
    const ArrayOfSingleScatteringData &arts_scat_data,
    Index index_start,
    Index index_end)
    : name_(name),
      index_start_(index_start),
      index_end_(index_end)
{
    Index n = arts_scat_data.size();
    std::vector<scattering::Particle> particles{};
    scattering::math::Vector<Numeric> d_eq(n), d_max(n), mass(n);

    for (size_t i = 0; i < n; ++i) {

        auto scattering_data = detail::from_legacy_format(arts_scat_data[i]);
        particles.emplace_back(0.0, 0.0, 0.0, scattering_data);
    }
    particle_habit_ = std::make_shared<scattering::ParticleHabit>(particles);
}

BulkParticleHabit::BulkParticleHabit(
    const String name,
    const Agenda &pnd_agenda,
    ArrayOfString pnd_agenda_input,
    std::shared_ptr<scattering::ParticleHabit> particle_model)
    : name_(name),
      pnd_agenda_(std::make_shared<Agenda>(pnd_agenda)),
      pnd_agenda_input_(pnd_agenda_input),
      particle_habit_(particle_model) {}

BulkParticleHabit::BulkParticleHabit(
    const String name,
    Index index_start,
    Index index_end,
    std::shared_ptr<scattering::ParticleHabit> particle_model)
    : name_(name),
      index_start_(index_start),
      index_end_(index_end),
      particle_habit_(particle_model) {}

ArrayOfString BulkParticleHabit::get_dpnd_data_dx_names(
    ArrayOfRetrievalQuantity jacobian_quantities, bool jacobian_do) const {
  ArrayOfString dpnd_data_dx_names = {};
  //if (jacobian_do) {
  //  for (auto &jq : jacobian_quantities) {
  //    if (jq.MainTag() == SCATSPECIES_MAINTAG) {
  //      if (jq.Subtag() == name_) {
  //        dpnd_data_dx_names.push_back(jq.SubSubtag());
  //      }
  //    }
  //  }
  //}
  return dpnd_data_dx_names;
}

Matrix BulkParticleHabit::get_agenda_input(Matrix pbp_field,
                                         ArrayOfString pbf_names) const {
    assert(pbp_field.nrows() == pbf_names.size());

  Index n_inputs = pnd_agenda_input_.size();
  ArrayOfIndex row_indices(n_inputs);

  for (Index i = 0; i < n_inputs; ++i) {
    row_indices[i] = find_first(pbf_names, pnd_agenda_input_[i]);
    if (row_indices[i] < 0) {
      ostringstream os;
      os << "Particle habit requires input " << pnd_agenda_input_[i]
         << "\",\nbut this quantity "
         << "could not found in *particle_bulkprop_names*.\n"
         << "(Note that temperature must be written as \"Temperature\")";
      throw runtime_error(os.str());
    }
  }

  Matrix agenda_input(pbp_field.ncols(), n_inputs);
  for (size_t i = 0; i < n_inputs; ++i) {
      agenda_input(joker, i) = pbp_field(row_indices[i], joker);
  }

  return agenda_input;
}

std::shared_ptr<ScatteringSpeciesImpl> BulkParticleHabit::prepare_scattering_data(ScatteringPropertiesSpec specs) const {

    scattering::ParticleHabit formatted = particle_habit_->set_stokes_dim(specs.n_stokes);


    if (specs.frame == ReferenceFrame::Lab) {
        Index n = 37 * 2 - 1; //specs.lat_scat.nelem();
        scattering::math::Vector<Numeric> lon_scat(n);
        Numeric dx = 2.0 * M_PI / static_cast<Numeric>(n - 1);
        lon_scat[0] = 0.0;
        for (Index i = 1; i < n; ++i) {
            lon_scat[i] = lon_scat[i - 1] + dx;
        }

        formatted = formatted.to_lab_frame(specs.lat_inc,
                                           std::make_shared<scattering::math::Vector<Numeric>>(lon_scat),
                                           specs.lat_scat,
                                           specs.n_stokes);
    }

    if (specs.format == Format::Spectral) {
        // Data must be regridded to ensure conformity with grids expected
        // by SHT transform.
        formatted = formatted.regrid().to_spectral(specs.l_max, specs.m_max);
    } else {
        formatted = formatted.downsample_scattering_angles(specs.lon_scat,
                                                           specs.lat_scat);
        formatted = formatted.to_gridded(specs.lon_inc,
                                         specs.lat_inc,
                                         specs.lon_scat,
                                         specs.lat_scat);
    }

    auto new_model = std::make_shared<scattering::ParticleHabit>(formatted.interpolate_frequency(specs.f_grid));
    auto result = std::make_shared<BulkParticleHabit>(name_, *pnd_agenda_, pnd_agenda_input_, new_model);
    return result;
}

BulkScatteringProperties BulkParticleHabit::calculate_bulk_properties(
    Workspace &ws,
    ConstMatrixView pbp_field,
    const ArrayOfString &pbf_names,
    ConstVectorView temperature,
    const ArrayOfRetrievalQuantity &jacobian_quantities,
    bool jacobian_do) const {

  // If pnd_agenda_ is not null, PND field is calculated using the provided
  // PSD.
  if (pnd_agenda_) {
    Matrix pnd_data;
    Tensor3 dpnd_data_dx;
    Matrix input = get_agenda_input(pbp_field, pbf_names);
    ArrayOfString dpnd_data_dx_names =
        get_dpnd_data_dx_names(jacobian_quantities, jacobian_do);

    pnd_agendaExecute(ws,
                      pnd_data,
                      dpnd_data_dx,
                      Array<ArrayOfScatteringMetaData>{get_meta_data()},
                      temperature,
                      input,
                      pnd_agenda_input_,
                      dpnd_data_dx_names,
                      *pnd_agenda_);

    auto n_levels = pnd_data.nrows();
    Array<scattering::SingleScatteringData> bulk_properties(n_levels);
    for (Index i = 0; i < n_levels; ++i) {
        scattering::math::Vector<Numeric> number_densities = scattering::to_eigen(pnd_data(i, joker));
      bulk_properties[i] = particle_habit_->calculate_bulk_properties(
          temperature[i], number_densities);
    }
    return BulkScatteringProperties(bulk_properties);
  } else {
      if (pbp_field.ncols() < particle_habit_->size()) {
          throw std::runtime_error(
              "If a scattering habit has no PSD agenda, the pbp-field must "
              "have at least as many columns as the particle model has "
              "elements."
              );
      }
      auto col_start = index_start_;
      ConstMatrixView pnd_data = pbp_field(joker, Range(col_start, pbp_field.ncols()));
      auto n_levels = pnd_data.ncols();
      Array<scattering::SingleScatteringData> bulk_properties(n_levels);
      for (Index i = 0; i < n_levels; ++i) {
          scattering::math::Vector<Numeric> number_densities = scattering::to_eigen(pnd_data(i, joker));
          bulk_properties[i] = particle_habit_->calculate_bulk_properties(
              temperature[i], number_densities);
      }
      return BulkScatteringProperties(bulk_properties);
  }
}

std::pair<Matrix, Tensor3> BulkParticleHabit::get_absorption_and_extinction(
    Workspace &ws,
    ConstMatrixView pbp_field,
    const ArrayOfString &pbf_names,
    Numeric frequency,
    ConstVectorView temperature,
    Numeric lon_inc,
    Numeric lat_inc,
    const ArrayOfRetrievalQuantity &jacobian_quantities,
    bool jacobian_do) const {

  // If pnd_agenda_ is not null, PND field is calculated using the provided
  // PSD.
  if (pnd_agenda_.get()) {
    Matrix pnd_data;
    Tensor3 dpnd_data_dx;
    Matrix input = get_agenda_input(pbp_field, pbf_names);
    ArrayOfString dpnd_data_dx_names =
        get_dpnd_data_dx_names(jacobian_quantities, jacobian_do);

    pnd_agendaExecute(ws,
                      pnd_data,
                      dpnd_data_dx,
                      Array<ArrayOfScatteringMetaData>{get_meta_data()},
                      temperature,
                      input,
                      pnd_agenda_input_,
                      dpnd_data_dx_names,
                      *pnd_agenda_);

    auto n_levels = pnd_data.nrows();
    auto stokes_dim = particle_habit_->get_stokes_dim();
    Matrix absorption(n_levels, stokes_dim);
    Tensor3 extinction(n_levels, stokes_dim, stokes_dim);
    for (Index i = 0; i < n_levels; ++i) {
        scattering::math::Vector<Numeric> number_densities = scattering::to_eigen(pnd_data(i, joker));
      auto abs = particle_habit_->get_absorption_vector(
          frequency,
          temperature[i],
          lon_inc,
          lat_inc,
          number_densities,
          stokes_dim);
      absorption(i, joker) = to_arts(abs);
      auto ext = particle_habit_->get_extinction_matrix(
          frequency,
          temperature[i],
          lon_inc,
          lat_inc,
          number_densities,
          stokes_dim);
      extinction(i, joker, joker) = to_arts(ext);
    }
    return std::make_pair(absorption, extinction);
  } else {
      if (pbp_field.ncols() < particle_habit_->size()) {
          throw std::runtime_error(
              "If a scattering habit has no PSD agenda, the pbp-field must "
              "have at least as many columns as the particle model has "
              "elements."
              );
      }
      auto col_start = index_start_;
      ConstMatrixView pnd_data = pbp_field(joker, Range(col_start, pbp_field.ncols()));
      auto n_levels = pnd_data.nrows();
      auto stokes_dim = particle_habit_->get_stokes_dim();
      Matrix absorption(n_levels, stokes_dim);
      Tensor3 extinction(n_levels, stokes_dim, stokes_dim);
      std::cout << "N_LEVELS :: " << n_levels << std::endl;
      for (Index i = 0; i < n_levels; ++i) {
          scattering::math::Vector<Numeric> number_densities = scattering::to_eigen(pnd_data(i, joker));
          auto abs = particle_habit_->get_absorption_vector(
              frequency,
              temperature[i],
              lon_inc,
              lat_inc,
              number_densities,
              stokes_dim);
          absorption(i, joker) = to_arts(abs);
          auto ext = particle_habit_->get_extinction_matrix(
                  frequency,
                  temperature[i],
                  lon_inc,
                  lat_inc,
                  number_densities,
                  stokes_dim);
          extinction(i, joker, joker) = to_arts(ext);
      }
      return std::make_pair(absorption, extinction);
  }
}


std::pair<Vector, Matrix> BulkParticleHabit::sample_incoming_direction(
    Workspace &ws,
    ConstMatrixView pbp_field,
    const ArrayOfString &pbp_names,
    Numeric frequency,
    ConstVectorView temperature,
    ConstVectorView los_scat_rev,
    Rng& rng,
    Numeric scat_coeff_tot) const {

  // Determine PND.
  Matrix pnd_data;
  if (pnd_agenda_.get()) {
    Tensor3 dpnd_data_dx;
    Matrix input = get_agenda_input(pbp_field, pbp_names);
    ArrayOfString dpnd_data_dx_names =
        get_dpnd_data_dx_names({}, false);
    pnd_agendaExecute(ws,
                      pnd_data,
                      dpnd_data_dx,
                      Array<ArrayOfScatteringMetaData>{get_meta_data()},
                      temperature,
                      input,
                      pnd_agenda_input_,
                      dpnd_data_dx_names,
                      *pnd_agenda_);
  } else {
     auto col_start = index_start_;
     pnd_data = pbp_field(joker, Range(col_start, pbp_field.ncols()));
  }

  // Prepare input for sampling.
  Vector los_scat{3};
  mirror_los(los_scat, los_scat_rev, 3);
  VectorView pnd = pnd_data(0, joker);

  auto stokes_dim = particle_habit_->get_stokes_dim();

  bool found = false;
  Numeric lon_inc = 0.0;
  Numeric lat_inc = 0.0;
  Numeric lon_scat = Conversion::deg2rad(los_scat[1]);
  Numeric lat_scat = Conversion::deg2rad(los_scat[0]);

  scattering::math::Vector<Numeric> number_densities = scattering::to_eigen(pnd_data(0, joker));
  Numeric phase_function_max = particle_habit_->get_phase_function_maximum_inc(
      frequency,
      temperature[0],
      lon_scat,
      lat_scat,
      number_densities
      );

  Numeric Z_00 = 0.0;
  Matrix Z;
  std::cout << "ANGS SCAT ::" << lon_scat << " / " << lat_scat << std::endl;
  std::cout << "PND :: " << number_densities << std::endl;

  while (!found) {
    lon_inc = 2.0 * M_PI * rng.draw() - M_PI;
    lat_inc = acos(-1.0 + 2.0 * rng.draw());

    Z_00 = particle_habit_
        ->get_phase_matrix(frequency,
                          temperature[0],
                          lon_inc,
                          lat_inc,
                          lon_scat,
                          lat_scat,
                          number_densities,
                          1)(0, 0);
    double r = rng.draw();
    std::cout << Z_00 << " / " << r  << " / " << phase_function_max << std::endl;
    if (r <= Z_00 / phase_function_max) {
      found = true;
      Z = to_arts(
          particle_habit_->get_phase_matrix(
              frequency,
              temperature[0],
              lon_inc,
              lat_inc,
              lon_scat,
              lat_scat,
              number_densities,
              stokes_dim)
          );
    }
  }

  Vector los_inc(2), los_inc_rev(2);
  los_inc[0] = Conversion::rad2deg(lat_inc);
  los_inc[1] = Conversion::rad2deg(lon_inc);
  mirror_los(los_inc_rev, los_inc, 3);
  Z *= scat_coeff_tot / Z(0, 0);

  return std::make_pair(los_inc_rev, Z);
}

std::ostream &operator<<(std::ostream &output, const BulkParticleHabit &habit) {
  output << "Scattering habit: " << habit.name_ << std::endl;
  output << "\t Particle d_eq: " << habit.get_particle_d_eq() << std::endl;
  output << "\t Particle d_max: " << habit.get_particle_d_max() << std::endl;
  output << std::endl;
  return output;
}
