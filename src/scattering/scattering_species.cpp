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
  \file   scattering.cc
  \author Simon Pfreundschuh <simon.pfreundschuh@chalmers.se>
  \date   2020-09-25

  \brief  Implementation of scattering.h
*/
#include <memory>

#include "scattering/maths.h"
#include "scattering/scattering_species.h"

using scattering::math::tensor_index;

////////////////////////////////////////////////////////////////////////////////
// ScatteringPropertiesSpec
////////////////////////////////////////////////////////////////////////////////

ScatteringPropertiesSpec::ScatteringPropertiesSpec(const Vector& f_grid_,
                                                   ReferenceFrame frame_,
                                                   Index stokes_dim,
                                                   Index l_max_,
                                                   Index m_max_,
                                                   Numeric phase_function_norm_)
    : format(Format::Spectral),
      n_stokes(stokes_dim),
      l_max(l_max_),
      m_max(m_max_),
      phase_function_norm(phase_function_norm_),
      f_grid(std::make_shared<EigenVector>(to_eigen(f_grid_))) {}

ScatteringPropertiesSpec::ScatteringPropertiesSpec(EigenVectorPtr f_grid_,
                                                   ReferenceFrame frame_,
                                                   Index stokes_dim,
                                                   EigenVectorPtr lon_scat_,
                                                   scattering::LatitudeGridPtr<Numeric> lat_scat_,
                                                   Numeric phase_function_norm_)
    : format(Format::Gridded),
      n_stokes(stokes_dim),
      phase_function_norm(phase_function_norm_),
      lon_inc(std::make_shared<EigenVector>(1, 1)),
      lat_inc(std::make_shared<EigenVector>(1, 1)),
      lon_scat(lon_scat_),
      lat_scat(lat_scat_),
      f_grid(f_grid_) {
  (*lon_inc)[0] = M_PI;
  (*lat_inc)[0] = 0.5 * M_PI;
}

ScatteringPropertiesSpec::ScatteringPropertiesSpec(EigenVectorPtr f_grid_,
                                                   ReferenceFrame frame_,
                                                   Index stokes_dim,
                                                   EigenVectorPtr lat_inc_,
                                                   EigenVectorPtr lon_scat_,
                                                   scattering::LatitudeGridPtr<Numeric> lat_scat_,
                                                   Numeric phase_function_norm_)
    : format(Format::Gridded),
      frame(frame_),
      n_stokes(stokes_dim),
      phase_function_norm(phase_function_norm_),
      lon_inc(std::make_shared<EigenVector>(1, 1)),
      lat_inc(lat_inc_),
      lon_scat(lon_scat_),
      lat_scat(lat_scat_),
      f_grid(f_grid_) {
    (*lon_inc)[0] = M_PI;
    std::cout << f_grid << " / " << lon_inc << " / " << lat_inc << " / " << lon_scat << " / " << lat_scat << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
// BulkScatteringProperties
////////////////////////////////////////////////////////////////////////////////

Tensor4 BulkScatteringProperties::get_extinction_coeff() const {
  auto n_layers = data_.size();
  auto n_freqs = data_[0].get_n_freqs();
  auto n_lon_inc = data_[0].get_n_lon_inc();
  auto n_lat_inc = data_[0].get_n_lat_inc();
  Tensor4 result(n_freqs, n_layers, n_lon_inc, n_lat_inc);
  for (Index i = 0; i < n_layers; ++i) {
    // This is a rank-5 tensor, but degenerate along the axis corresponding
    // to temperature, so we can copy it directly into the output.
    auto exintction = data_[i].get_extinction_coeff();
    for (Index j = 0; j < n_freqs; ++j) {
      auto exintction_f = tensor_index<1>(exintction, {j});
      auto input_start = exintction_f.data();
      auto input_end = exintction_f.data() + exintction_f.size();
      auto result_view = result(j, i, joker, joker);
      std::copy(input_start, input_end, result_view.get_c_array());
    }
  }
  return result;
}

Tensor6 BulkScatteringProperties::get_extinction_matrix(
    Index stokes_dim) const {
  auto n_layers = data_.size();
  auto n_freqs = data_[0].get_n_freqs();
  auto n_lon_inc = data_[0].get_n_lon_inc();
  auto n_lat_inc = data_[0].get_n_lat_inc();
  Tensor6 result(
      n_freqs, n_layers, n_lon_inc, n_lat_inc, stokes_dim, stokes_dim, 0.0);
  for (Index i = 0; i < n_layers; ++i) {
    // This is a rank-8 tensor, but degenerate along the axis corresponding
    // to temperature as well as the scattering angles, so we can copy it
    // directly into the output.
    auto extinction = data_[i].get_extinction_matrix(stokes_dim);
    for (Index j = 0; j < n_freqs; ++j) {
      auto extinction_f = tensor_index<1>(extinction, {j});
      auto input_start = extinction_f.data();
      auto input_end = extinction_f.data() + extinction_f.size();
      auto result_view = result(j, i, joker, joker, joker, joker);
      std::copy(input_start, input_end, result_view.get_c_array());
    }
  }
  return result;
}

Tensor6 BulkScatteringProperties::get_extinction_matrix() const {
  return get_extinction_matrix(stokes_dim_);
}

Tensor4 BulkScatteringProperties::get_absorption_coeff() const {
  auto n_layers = data_.size();
  auto n_freqs = data_[0].get_n_freqs();
  auto n_lon_inc = data_[0].get_n_lon_inc();
  auto n_lat_inc = data_[0].get_n_lat_inc();
  Tensor4 result(n_freqs, n_layers, n_lon_inc, n_lat_inc);
  for (Index i = 0; i < n_layers; ++i) {
    // This is a rank-5 tensor, but degenerate along the axis corresponding
    // to temperature, so we can copy it directly into the output.
    auto absorption = data_[i].get_absorption_coeff();
    for (Index j = 0; j < n_freqs; ++j) {
      auto absorption_f = tensor_index<1>(absorption, {j});
      auto input_start = absorption_f.data();
      auto input_end = absorption_f.data() + absorption_f.size();
      auto result_view = result(j, i, joker, joker);
      std::copy(input_start, input_end, result_view.get_c_array());
    }
  }
  return result;
}

Tensor5 BulkScatteringProperties::get_absorption_vector(
    Index stokes_dim) const {
  auto n_layers = data_.size();
  auto n_freqs = data_[0].get_n_freqs();
  auto n_lon_inc = data_[0].get_n_lon_inc();
  auto n_lat_inc = data_[0].get_n_lat_inc();
  Tensor5 result(n_freqs, n_layers, n_lon_inc, n_lat_inc, stokes_dim);
  for (Index i = 0; i < n_layers; ++i) {
    // This is a rank-5 tensor, but degenerate along the axis corresponding
    // to temperature, so we can copy it directly into the output.
    auto absorption = data_[i].get_absorption_vector(stokes_dim);
    for (Index j = 0; j < n_freqs; ++j) {
      auto absorption_f = tensor_index<1>(absorption, {j});
      auto input_start = absorption_f.data();
      auto input_end = absorption_f.data() + absorption_f.size();
      auto result_view = result(j, i, joker, joker, joker);
      std::copy(input_start, input_end, result_view.get_c_array());
    }
  }
  return result;
}

Tensor5 BulkScatteringProperties::get_absorption_vector() const {
  return get_absorption_vector(stokes_dim_);
}

Tensor5 BulkScatteringProperties::get_spectral_coeffs() const {
  auto n_layers = data_.size();
  auto n_coeffs = data_[0].get_phase_function_spectral().dimension(4);
  auto n_freqs = data_[0].get_n_freqs();
  auto n_lon_inc = data_[0].get_n_lon_inc();
  auto n_lat_inc = data_[0].get_n_lat_inc();
  Tensor5 result(n_freqs_, n_layers, n_lon_inc, n_lat_inc, n_coeffs);
  for (Index i = 0; i < n_layers; ++i) {
    EigenTensor<5> coeffs = data_[i].get_phase_function_spectral().real();
    for (Index j = 0; j < n_freqs; ++j) {
      auto coeffs_f = tensor_index<1>(coeffs, {j});
      auto input_start = coeffs_f.data();
      auto input_end = coeffs_f.data() + coeffs_f.size();
      EigenTensor<5> phase_matrix =
          data_[i].get_phase_function_spectral().real();
      auto result_view = result(j, i, joker, joker, joker);
      std::copy(input_start, input_end, result_view.get_c_array());
    }
  }
  return result;
}

Tensor5 BulkScatteringProperties::get_legendre_coeffs() const {
  auto n_coeffs = data_[0].get_phase_function_spectral().dimension(4);
  auto result = get_spectral_coeffs();
  // Need to go from SHT coefficients to Legendre coefficients.
  for (Index i = 0; i < n_coeffs; ++i) {
    result(joker, joker, joker, joker, i) *= sqrt(4.0 * M_PI / (2.0 * i + 1.0));
  }
  return result;
}

Tensor7 BulkScatteringProperties::get_phase_matrix(Index stokes_dim) const {
  auto n_layers = data_.size();
  auto n_lat_inc = data_[0].get_n_lat_inc();
  auto n_lon_scat = data_[0].get_n_lon_scat();
  auto n_lat_scat = data_[0].get_n_lat_scat();
  Tensor7 result(n_freqs_,
                 n_layers,
                 n_lat_inc,
                 n_lon_scat,
                 n_lat_scat,
                 stokes_dim,
                 stokes_dim,
                 0.0);
  for (Index i = 0; i < n_layers; ++i) {
    EigenTensor<8> z = data_[i].get_phase_matrix(stokes_dim);
    for (Index j = 0; j < n_freqs_; ++j) {
      for (Index k = 0; k < n_lat_inc; ++k) {
        for (Index l = 0; l < n_lon_scat; ++l) {
          for (Index m = 0; m < n_lat_scat; ++m) {
            for (Index n = 0; n < stokes_dim; ++n) {
              for (Index o = 0; o < stokes_dim; ++o) {
                result(j, i, k, l, m, n, o) =
                    z.coeffRef(j, 0, 0, k, l, m, n, o);
              }
            }
          }
        }
      }
    }
  }
  return result;
}

Tensor7 BulkScatteringProperties::get_phase_matrix() const {
  return get_phase_matrix(stokes_dim_);
}

Tensor3 BulkScatteringProperties::get_phase_function() const {
  auto n = data_.size();
  Index n_spectral_coeffs_ = data_[0].get_phase_function().dimension(5);
  Tensor3 result(n_freqs_, n, n_spectral_coeffs_);
  for (Index i = 0; i < n; ++i) {
    EigenTensor<6> phase_matrix = data_[i].get_phase_function();
    result(joker, i, joker) = MatrixView(phase_matrix.data(),
                                         Range(0, n_freqs_, n_spectral_coeffs_),
                                         Range(0, n_spectral_coeffs_, 1));
  }
  return result;
}
