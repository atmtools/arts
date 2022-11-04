#include "scattering/single_scattering_data.h"

namespace scattering {

//
// Gridded
//

SingleScatteringData::SingleScatteringData(
    math::VectorPtr<double> f_grid,
    math::VectorPtr<double> t_grid,
    math::VectorPtr<double> lon_inc,
    math::VectorPtr<double> lat_inc,
    math::VectorPtr<double> lon_scat,
    std::shared_ptr<LatitudeGrid<double>> lat_scat,
    math::TensorPtr<double, 7> phase_matrix,
    math::TensorPtr<double, 7> extinction_matrix,
    math::TensorPtr<double, 7> absorption_vector,
    math::TensorPtr<double, 7> backward_scattering_coeff,
    math::TensorPtr<double, 7> forward_scattering_coeff)
    : data_(new SingleScatteringDataGridded<double>(f_grid,
                                                    t_grid,
                                                    lon_inc,
                                                    lat_inc,
                                                    lon_scat,
                                                    lat_scat,
                                                    phase_matrix,
                                                    extinction_matrix,
                                                    absorption_vector,
                                                    backward_scattering_coeff,
                                                    forward_scattering_coeff)) {}

SingleScatteringData::SingleScatteringData(
    scattering::math::Vector<double> f_grid,
    scattering::math::Vector<double> t_grid,
    scattering::math::Vector<double> lon_inc,
    scattering::math::Vector<double> lat_inc,
    scattering::math::Vector<double> lon_scat,
    scattering::math::Vector<double> lat_scat,
    ParticleType type)
    : SingleScatteringData(
          std::make_shared<math::Vector<double>>(f_grid),
          std::make_shared<math::Vector<double>>(t_grid),
          std::make_shared<math::Vector<double>>(lon_inc),
          std::make_shared<math::Vector<double>>(lat_inc),
          std::make_shared<math::Vector<double>>(lon_scat),
          std::make_shared<IrregularLatitudeGrid<double>>(lat_scat),
          std::make_shared<math::Tensor<double, 7>>(
              std::array<Index, 7>{f_grid.size(),
                                   t_grid.size(),
                                   lon_inc.size(),
                                   lat_inc.size(),
                                   lon_scat.size(),
                                   lat_scat.size(),
                                   detail::get_n_phase_matrix_elements(type)}),
          std::make_shared<math::Tensor<double, 7>>(std::array<Index, 7>{
              f_grid.size(),
              t_grid.size(),
              lon_inc.size(),
              lat_inc.size(),
              1,
              1,
              detail::get_n_extinction_matrix_elements(type)}),
          std::make_shared<math::Tensor<double, 7>>(std::array<Index, 7>{
              f_grid.size(),
              t_grid.size(),
              lon_inc.size(),
              lat_inc.size(),
              1,
              1,
              detail::get_n_absorption_vector_elements(type)}),
          std::make_shared<math::Tensor<double, 7>>(
              std::array<Index, 7>{f_grid.size(),
                                   t_grid.size(),
                                   lon_inc.size(),
                                   lat_inc.size(),
                                   1,
                                   1,
                                   1}),
          std::make_shared<math::Tensor<double, 7>>(
              std::array<Index, 7>{f_grid.size(),
                                   t_grid.size(),
                                   lon_inc.size(),
                                   lat_inc.size(),
                                   1,
                                   1,
                                   1})) {}

//
// Spectral
//

SingleScatteringData::SingleScatteringData(
    math::VectorPtr<double> f_grid,
    math::VectorPtr<double> t_grid,
    math::VectorPtr<double> lon_inc,
    math::VectorPtr<double> lat_inc,
    std::shared_ptr<sht::SHT> sht_scat,
    math::TensorPtr<std::complex<double>, 6> phase_matrix,
    math::TensorPtr<std::complex<double>, 6> extinction_matrix,
    math::TensorPtr<std::complex<double>, 6> absorption_vector,
    math::TensorPtr<std::complex<double>, 6> backward_scattering_coeff,
    math::TensorPtr<std::complex<double>, 6> forward_scattering_coeff)
    : data_(new SingleScatteringDataSpectral<double>(f_grid,
                                                     t_grid,
                                                     lon_inc,
                                                     lat_inc,
                                                     sht_scat,
                                                     phase_matrix,
                                                     extinction_matrix,
                                                     absorption_vector,
                                                     backward_scattering_coeff,
                                                     forward_scattering_coeff)) {}

SingleScatteringData::SingleScatteringData(
    math::Vector<double> f_grid,
    math::Vector<double> t_grid,
    math::Vector<double> lon_inc,
    math::Vector<double> lat_inc,
    sht::SHT sht_scat,
    math::Tensor<std::complex<double>, 6> phase_matrix,
    math::Tensor<std::complex<double>, 6> extinction_matrix,
    math::Tensor<std::complex<double>, 6> absorption_vector,
    math::Tensor<std::complex<double>, 6> backward_scattering_coeff,
    math::Tensor<std::complex<double>, 6> forward_scattering_coeff)
    : data_(new SingleScatteringDataSpectral<double>(
                std::make_shared<math::Vector<double>>(f_grid),
                std::make_shared<math::Vector<double>>(t_grid),
                std::make_shared<math::Vector<double>>(lon_inc),
                std::make_shared<math::Vector<double>>(lat_inc),
                std::make_shared<sht::SHT>(sht_scat),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(phase_matrix),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(extinction_matrix),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(absorption_vector),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(backward_scattering_coeff),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(forward_scattering_coeff))) {}


SingleScatteringData::SingleScatteringData(
    math::Vector<double> f_grid,
    math::Vector<double> t_grid,
    math::Vector<double> lon_inc,
    math::Vector<double> lat_inc,
    math::Tensor<std::complex<double>, 6> phase_matrix,
    math::Tensor<std::complex<double>, 6> extinction_matrix,
    math::Tensor<std::complex<double>, 6> absorption_vector,
    math::Tensor<std::complex<double>, 6> backward_scattering_coeff,
    math::Tensor<std::complex<double>, 6> forward_scattering_coeff)
    : data_(new SingleScatteringDataSpectral<double>(
                std::make_shared<math::Vector<double>>(f_grid),
                std::make_shared<math::Vector<double>>(t_grid),
                std::make_shared<math::Vector<double>>(lon_inc),
                std::make_shared<math::Vector<double>>(lat_inc),
                std::make_shared<sht::SHT>(
                    sht::SHT(sht::SHT::calc_l_max(phase_matrix.dimension(4)))
                    ),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(phase_matrix),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(extinction_matrix),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(absorption_vector),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(backward_scattering_coeff),
                std::make_shared<math::Tensor<std::complex<double>, 6>>(forward_scattering_coeff))) {}


SingleScatteringData::SingleScatteringData(
    scattering::math::Vector<double> f_grid,
    scattering::math::Vector<double> t_grid,
    scattering::math::Vector<double> lon_inc,
    scattering::math::Vector<double> lat_inc,
    Index l_max,
    ParticleType type)
    : SingleScatteringData(
          std::make_shared<math::Vector<double>>(f_grid),
          std::make_shared<math::Vector<double>>(t_grid),
          std::make_shared<math::Vector<double>>(lon_inc),
          std::make_shared<math::Vector<double>>(lat_inc),
          std::make_shared<sht::SHT>(l_max,
                                     l_max,
                                     2 * l_max + 2,
                                     2 * l_max + 2),
          std::make_shared<math::Tensor<std::complex<double>, 6>>(
              math::zeros<std::complex<double>>(
                  f_grid.size(),
                  t_grid.size(),
                  lon_inc.size(),
                  lat_inc.size(),
                  sht::SHT::calc_n_spectral_coeffs(l_max, l_max),
                  detail::get_n_phase_matrix_elements(type))),
          std::make_shared<math::Tensor<std::complex<double>, 6>>(
              math::zeros<std::complex<double>>(
                  f_grid.size(),
                  t_grid.size(),
                  lon_inc.size(),
                  lat_inc.size(),
                  1,
                  detail::get_n_extinction_matrix_elements(type))),
          std::make_shared<math::Tensor<std::complex<double>, 6>>(
              math::zeros<std::complex<double>>(
                  f_grid.size(),
                  t_grid.size(),
                  lon_inc.size(),
                  lat_inc.size(),
                  1,
                  detail::get_n_absorption_vector_elements(type))),
          std::make_shared<math::Tensor<std::complex<double>, 6>>(
              math::zeros<std::complex<double>>(
              f_grid.size(),
                                   t_grid.size(),
                                   lon_inc.size(),
                                   lat_inc.size(),
                                   1,
                      1)),
          std::make_shared<math::Tensor<std::complex<double>, 6>>(
          math::zeros<std::complex<double>>(
              f_grid.size(),
                                   t_grid.size(),
                                   lon_inc.size(),
                                   lat_inc.size(),
                                   1,
              1))) {}

}
