/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   gas_scattering.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-03-15

  \brief  Implementation file for functions related to gas scattering.
*/

#include "gas_scattering.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <physics_funcs.h>

#include <cmath>

#include "sht.h"

namespace {
/** The Rayleigh scattering matrix at scattering angle theta [rad], with F11 integrating to 4 pi over the sphere */
rtepack::compact_planar_muelmat rayleigh_scattering_matrix(Numeric theta_rad, Numeric depolarization_factor) {
  using Math::pow2;

  const Numeric delta       = (1.0 - depolarization_factor) / (1.0 + 0.5 * depolarization_factor);
  const Numeric delta_prime = (1.0 - 2.0 * depolarization_factor) / (1.0 - depolarization_factor);
  const Numeric cos_theta   = cos(theta_rad);

  return {0.75 * delta * (1.0 + pow2(cos_theta)) + 1.0 - delta,  // F11
          -0.75 * delta * pow2(sin(theta_rad)),                  // F12
          0.75 * delta * (1.0 + pow2(cos_theta)),                // F22
          1.5 * delta * cos_theta,                               // F33
          0.0,                                                   // F34
          1.5 * delta * delta_prime * cos_theta};                // F44
}
}  // namespace

Vector calc_rayleighPhaMat(const Numeric& theta_rad, const Numeric& depolarization_factor) {
  using Constant::pi;

  ARTS_USER_ERROR_IF(theta_rad != std::clamp<Numeric>(theta_rad, 0.0, pi),
                     "Error in calc_rayleighPhaMat: Scattering angle must be in the range [0, pi], is {}",
                     theta_rad);
  ARTS_USER_ERROR_IF(depolarization_factor < 0.0 or depolarization_factor >= 1.0,
                     "The depolarization factor must be in [0, 1), is {}",
                     depolarization_factor);

  const auto f = rayleigh_scattering_matrix(theta_rad, depolarization_factor);
  return Vector{f[0], f[1], f[2], f[3], f[4], f[5]};
}

namespace scattering {
namespace {

Numeric air_simple_cross_section(Numeric frequency) {
  static constexpr std::array coefficients{3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};

  const Numeric wavelength_um = Conversion::freq2wavelen(frequency) * 1e6;
  Numeric       sum           = 0.0;
  Numeric       inverse_power = 1.0;
  for (const auto coefficient : coefficients) {
    sum           += coefficient * inverse_power;
    inverse_power /= Math::pow2(wavelength_um);
  }
  return 1e-32 * sum / Math::pow4(wavelength_um);
}

/** The scattering matrix at scattering angle theta [rad], with F11 integrating to 4 pi over the sphere */
rtepack::compact_planar_muelmat normalized_phase_matrix(const GasScatteringPhaseMatrix& phase_matrix,
                                                        Numeric                         scattering_angle) {
  return std::visit(
      [scattering_angle](const auto& phase) -> rtepack::compact_planar_muelmat {
        using Phase = std::remove_cvref_t<decltype(phase)>;
        if constexpr (std::is_same_v<Phase, IsotropicGasScattering>) {
          return {1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
        } else {
          return rayleigh_scattering_matrix(scattering_angle, phase.depolarization_factor);
        }
      },
      phase_matrix);
}

Numeric scattering_coefficient(const GasScatteringCoefficient& coefficient,
                               Numeric                         frequency,
                               const AtmPoint&                 atm_point) {
  return std::visit([&](const auto& model) { return model(frequency, atm_point); }, coefficient);
}

/** The relative derivative of the scattering coefficient, which is proportional to number density */
Numeric relative_derivative(const AtmPoint& atm_point, const AtmKeyVal& target) {
  if (target == AtmKeyVal{AtmKey::p}) return 1.0 / atm_point.pressure;
  if (target == AtmKeyVal{AtmKey::t}) return -1.0 / atm_point.temperature;
  return 0.0;
}

}  // namespace

Numeric AirSimpleGasScattering::operator()(Numeric frequency, const AtmPoint& atm_point) const {
  return air_simple_cross_section(frequency) * number_density(atm_point.pressure, atm_point.temperature);
}

ConstantGasScattering::ConstantGasScattering(Numeric cross_section_)
    : cross_section(cross_section_){ARTS_USER_ERROR_IF(
          cross_section < 0.0, "The gas-scattering cross-section must be non-negative, is {}", cross_section)}

      Numeric ConstantGasScattering::operator()(Numeric, const AtmPoint& atm_point) const {
  return cross_section * number_density(atm_point.pressure, atm_point.temperature);
}

RayleighGasScattering::RayleighGasScattering(Numeric depolarization_factor_)
    : depolarization_factor(depolarization_factor_){
          ARTS_USER_ERROR_IF(depolarization_factor < 0.0 or depolarization_factor >= 1.0,
                             "The depolarization factor must be in [0, 1), is {}",
                             depolarization_factor)}

      GasScatterer::GasScatterer(GasScatteringCoefficient coefficient_, GasScatteringPhaseMatrix phase_matrix_)
    : coefficient(std::move(coefficient_)), phase_matrix(std::move(phase_matrix_)) {}

BulkScatteringProperties<Format::TRO, Representation::Gridded> GasScatterer::get_bulk_scattering_properties_tro_gridded(
    const AtmPoint& atm_point, const Vector& f_grid, std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  ARTS_USER_ERROR_IF(not za_scat_grid, "A scattering zenith-angle grid is required")

  auto                                                           t_grid_ptr = std::make_shared<Vector>(Vector{0.0});
  auto                                                           f_grid_ptr = std::make_shared<Vector>(f_grid);
  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> phase{t_grid_ptr, f_grid_ptr, za_scat_grid};
  ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> extinction{t_grid_ptr, f_grid_ptr};
  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> absorption{t_grid_ptr, f_grid_ptr};

  const auto angles = grid_vector(*za_scat_grid);
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    const Numeric coefficient_per_m = scattering_coefficient(coefficient, f_grid[iv], atm_point);
    extinction[0, iv, 0]            = coefficient_per_m;
    for (Size ia = 0; ia < angles.size(); ++ia) {
      const auto normalized = normalized_phase_matrix(phase_matrix, Conversion::deg2rad(angles[ia]));
      for (Index is = 0; is < 6; ++is) {
        phase[0, iv, ia, is] = coefficient_per_m * normalized[is] / (4.0 * Constant::pi);
      }
    }
  }

  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = std::move(extinction),
          .absorption_vector = std::move(absorption)};
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
GasScatterer::get_bulk_scattering_properties_tro_gridded_derivative(const AtmPoint&                  atm_point,
                                                                    const Vector&                    f_grid,
                                                                    std::shared_ptr<ZenithAngleGrid> za_scat_grid,
                                                                    const AtmKeyVal&                 target) const {
  auto out  = get_bulk_scattering_properties_tro_gridded(atm_point, f_grid, std::move(za_scat_grid));
  out      *= relative_derivative(atm_point, target);
  return out;
}

ScatteringTroSpectralVector GasScatterer::get_bulk_scattering_properties_tro_spectral(const AtmPoint& atm_point,
                                                                                      const Vector&   f_grid,
                                                                                      Index           degree) const {
  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be non-negative, is {}", degree)

  SpecmatMatrix phase(f_grid.size(), degree + 1, Specmat{0.0});
  PropmatVector extinction(f_grid.size());
  StokvecVector absorption(f_grid.size());

  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    const Numeric k    = scattering_coefficient(coefficient, f_grid[iv], atm_point);
    const Numeric amp  = 0.5 * Constant::inv_sqrt_pi * k;
    extinction[iv].A() = k;

    std::visit(
        [&](const auto& model) {
          using Model = std::remove_cvref_t<decltype(model)>;
          // Each coefficient is the scattering-plane Mueller matrix of the degree
          const auto set = [&](Index l, const rtepack::compact_planar_muelmat& F) {
            phase[iv, l] = Specmat{F.expand().data};
          };
          if constexpr (std::is_same_v<Model, IsotropicGasScattering>) {
            set(0, {amp, 0.0, amp, amp, 0.0, amp});
          } else {
            const Numeric delta       = (1.0 - model.depolarization_factor) / (1.0 + 0.5 * model.depolarization_factor);
            const Numeric delta_prime = (1.0 - 2.0 * model.depolarization_factor) / (1.0 - model.depolarization_factor);

            // F11 = 1 + delta/2 P2, F12 = delta/2(P2 - 1),
            // F22 = delta(1 + P2/2), F33 = 3 delta/2 P1, F44 = 3 delta delta'/2 P1
            set(0, {amp, -0.5 * delta * amp, delta * amp, 0.0, 0.0, 0.0});
            if (degree >= 1)
              set(1,
                  {0.0,
                   0.0,
                   0.0,
                   0.5 * std::numbers::sqrt3 * delta * amp,
                   0.0,
                   0.5 * std::numbers::sqrt3 * delta * delta_prime * amp});
            if (degree >= 2) {
              const Numeric p2 = 0.5 * delta * amp / std::sqrt(5.0);
              set(2, {p2, p2, p2, 0.0, 0.0, 0.0});
            }
          }
        },
        phase_matrix);
  }

  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = std::move(extinction),
          .absorption_vector = std::move(absorption)};
}

BulkScatteringProperties<Format::ARO, Representation::Gridded> GasScatterer::get_bulk_scattering_properties_aro_gridded(
    const AtmPoint&                  atm_point,
    const Vector&                    f_grid,
    const Vector&                    za_inc_grid,
    const Vector&                    delta_aa_grid,
    std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  ARTS_USER_ERROR_IF(not za_scat_grid, "A scattering zenith-angle grid is required")

  auto t_grid_ptr = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);
  auto za_inc_ptr = std::make_shared<Vector>(za_inc_grid);
  ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> extinction{t_grid_ptr, f_grid_ptr};
  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> absorption{t_grid_ptr, f_grid_ptr};
  Vector                                                              scale(f_grid.size());
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    extinction[0, iv, 0] = scattering_coefficient(coefficient, f_grid[iv], atm_point);
    scale[iv]            = extinction[0, iv, 0] / (4.0 * Constant::pi);
  }

  // The lab-frame phase matrix at the exact scattering angle of every direction pair
  const auto closed_form = [&](Numeric theta, matpack::data_t<Numeric, 3>& scattering_matrix) {
    const auto normalized = normalized_phase_matrix(phase_matrix, theta);
    for (Size iv = 0; iv < f_grid.size(); ++iv) {
      for (Index is = 0; is < 6; ++is) scattering_matrix[0, iv, is] = scale[iv] * normalized[is];
    }
  };
  auto phase = tro_lab_frame<Numeric>(t_grid_ptr,
                                      f_grid_ptr,
                                      za_inc_ptr,
                                      std::make_shared<Vector>(delta_aa_grid),
                                      std::move(za_scat_grid),
                                      closed_form);

  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = extinction.to_lab_frame(za_inc_ptr),
          .absorption_vector = absorption.to_lab_frame(za_inc_ptr)};
}

BulkScatteringProperties<Format::ARO, Representation::Gridded>
GasScatterer::get_bulk_scattering_properties_aro_gridded_derivative(const AtmPoint&                  atm_point,
                                                                    const Vector&                    f_grid,
                                                                    const Vector&                    za_inc_grid,
                                                                    const Vector&                    delta_aa_grid,
                                                                    std::shared_ptr<ZenithAngleGrid> za_scat_grid,
                                                                    const AtmKeyVal&                 target) const {
  auto out = get_bulk_scattering_properties_aro_gridded(
      atm_point, f_grid, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
  out *= relative_derivative(atm_point, target);
  return out;
}

BulkScatteringProperties<Format::ARO, Representation::Spectral>
GasScatterer::get_bulk_scattering_properties_aro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, const Vector& za_inc_grid, Index degree, Index order) const {
  auto sht = sht::provider.get_instance_lm(degree, order);
  return get_bulk_scattering_properties_aro_gridded(atm_point,
                                                    f_grid,
                                                    za_inc_grid,
                                                    *sht->get_aa_grid_ptr(),
                                                    std::make_shared<ZenithAngleGrid>(sht->get_zenith_angle_grid()))
      .to_spectral(degree, order);
}

BulkScatteringProperties<Format::ARO, Representation::Fourier> GasScatterer::get_bulk_scattering_properties_aro_fourier(
    const AtmPoint& atm_point,
    const Vector&   f_grid,
    const Vector&   za_inc_grid,
    const Vector&   za_scat_grid,
    Index           max_mode) const {
  auto t_grid_ptr = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);
  auto za_inc_ptr = std::make_shared<Vector>(za_inc_grid);
  ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> extinction{t_grid_ptr, f_grid_ptr};
  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> absorption{t_grid_ptr, f_grid_ptr};
  Vector                                                              scale(f_grid.size());
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    extinction[0, iv, 0] = scattering_coefficient(coefficient, f_grid[iv], atm_point);
    scale[iv]            = extinction[0, iv, 0] / (4.0 * Constant::pi);
  }

  const auto closed_form = [&](Numeric theta, matpack::data_t<Numeric, 3>& scattering_matrix) {
    const auto normalized = normalized_phase_matrix(phase_matrix, theta);
    for (Size iv = 0; iv < f_grid.size(); ++iv) {
      for (Index is = 0; is < 6; ++is) scattering_matrix[0, iv, is] = scale[iv] * normalized[is];
    }
  };
  // F11 / scale integrates to 4 pi, so the phase integral is the scattering coefficient
  Matrix integral(1, f_grid.size());
  for (Size iv = 0; iv < f_grid.size(); ++iv) integral[0, iv] = extinction[0, iv, 0];
  auto phase = tro_lab_frame_fourier_modes<Numeric>(
      t_grid_ptr,
      f_grid_ptr,
      za_inc_ptr,
      std::make_shared<const ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat_grid)),
      max_mode,
      integral,
      closed_form);

  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = extinction.to_lab_frame(za_inc_ptr).to_fourier(),
          .absorption_vector = absorption.to_lab_frame(za_inc_ptr).to_fourier()};
}

}  // namespace scattering

namespace {
template <typename T>
void write_empty_gas_scattering(std::ostream& os, std::string_view type_name, std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

template <typename T> void read_empty_gas_scattering(std::istream& is, std::string_view type_name) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
}  // namespace

void xml_io_stream<scattering::AirSimpleGasScattering>::write(std::ostream& os,
                                                              const scattering::AirSimpleGasScattering&,
                                                              bofstream*,
                                                              std::string_view name) {
  write_empty_gas_scattering<scattering::AirSimpleGasScattering>(os, type_name, name);
}
void xml_io_stream<scattering::AirSimpleGasScattering>::read(std::istream& is,
                                                             scattering::AirSimpleGasScattering&,
                                                             bifstream*) {
  read_empty_gas_scattering<scattering::AirSimpleGasScattering>(is, type_name);
}

void xml_io_stream<scattering::IsotropicGasScattering>::write(std::ostream& os,
                                                              const scattering::IsotropicGasScattering&,
                                                              bofstream*,
                                                              std::string_view name) {
  write_empty_gas_scattering<scattering::IsotropicGasScattering>(os, type_name, name);
}
void xml_io_stream<scattering::IsotropicGasScattering>::read(std::istream& is,
                                                             scattering::IsotropicGasScattering&,
                                                             bifstream*) {
  read_empty_gas_scattering<scattering::IsotropicGasScattering>(is, type_name);
}

void xml_io_stream<scattering::ConstantGasScattering>::write(std::ostream&                            os,
                                                             const scattering::ConstantGasScattering& x,
                                                             bofstream*                               pbofs,
                                                             std::string_view                         name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, x.cross_section, pbofs);
  tag.write_to_end_stream(os);
}
void xml_io_stream<scattering::ConstantGasScattering>::read(std::istream&                      is,
                                                            scattering::ConstantGasScattering& x,
                                                            bifstream*                         pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  xml_read_from_stream(is, x.cross_section, pbifs);
  tag.read_from_stream(is);
  tag.check_end_name(type_name);
  ARTS_USER_ERROR_IF(x.cross_section < 0.0, "The gas-scattering cross-section must be non-negative")
}

void xml_io_stream<scattering::RayleighGasScattering>::write(std::ostream&                            os,
                                                             const scattering::RayleighGasScattering& x,
                                                             bofstream*                               pbofs,
                                                             std::string_view                         name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, x.depolarization_factor, pbofs);
  tag.write_to_end_stream(os);
}
void xml_io_stream<scattering::RayleighGasScattering>::read(std::istream&                      is,
                                                            scattering::RayleighGasScattering& x,
                                                            bifstream*                         pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  xml_read_from_stream(is, x.depolarization_factor, pbifs);
  tag.read_from_stream(is);
  tag.check_end_name(type_name);
  ARTS_USER_ERROR_IF(x.depolarization_factor < 0.0 or x.depolarization_factor >= 1.0,
                     "The depolarization factor must be in [0, 1)")
}

void xml_io_stream<scattering::GasScatterer>::write(std::ostream&                   os,
                                                    const scattering::GasScatterer& x,
                                                    bofstream*                      pbofs,
                                                    std::string_view                name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, x.coefficient, pbofs);
  xml_write_to_stream(os, x.phase_matrix, pbofs);
  tag.write_to_end_stream(os);
}
void xml_io_stream<scattering::GasScatterer>::read(std::istream& is, scattering::GasScatterer& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  xml_read_from_stream(is, x.coefficient, pbifs);
  xml_read_from_stream(is, x.phase_matrix, pbifs);
  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
