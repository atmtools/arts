#include "radar.h"

#include <arts_constants.h>
#include <atm_path.h>
#include <physics_funcs.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>

namespace radar {
Numeric liebe93_k2(const Numeric frequency, const Numeric temperature) {
  const Complex permittivity = permittivity_water_liebe93(frequency, temperature);
  return std::norm((permittivity - 1.0) / (permittivity + 2.0));
}

Numeric ze_factor(const Numeric frequency, const Numeric ze_tref, const Numeric k2) {
  const Numeric dielectric = k2 < 0.0 ? liebe93_k2(frequency, ze_tref) : k2;
  if (dielectric <= 0.0) arts::user_error("Radar k2 must be positive, or negative to request Liebe-93; got {}", k2);
  const Numeric wavelength = Constant::speed_of_light / frequency;
  return 4e18 * std::pow(wavelength, 4) / (std::pow(Constant::pi, 4) * dielectric);
}

Numeric unit_factor(const RadarUnit unit, const Numeric frequency, const Numeric ze_tref, const Numeric k2) {
  switch (unit) {
    case RadarUnit::unit: return 1.0;
    case RadarUnit::Ze:
    case RadarUnit::dBZe: return ze_factor(frequency, ze_tref, k2);
  }
  std::unreachable();
}

RadarRangeMode resolve(const RadarRangeMode mode, const Numeric largest_limit) {
  if (mode != RadarRangeMode::Legacy) return mode;
  return largest_limit > 1.0 ? RadarRangeMode::Altitude : RadarRangeMode::RoundTripTime;
}

PathRange path_range(const ArrayOfPropagationPathPoint& path, const Vector3& radar_pos, const Vector2& ellipsoid) {
  const Size np = path.size();
  PathRange  out{
       .layer = np > 1 ? path::distance(path, ellipsoid) : Vector{}, .distance = Vector(np), .time = Vector(np)};
  if (np == 0) return out;

  out.distance[0] = path::distance(radar_pos, path.front().pos, ellipsoid);
  out.time[0]     = 2.0 * out.distance[0] / Constant::speed_of_light;
  for (Size ip = 1; ip < np; ++ip) {
    out.distance[ip] = out.distance[ip - 1] + out.layer[ip - 1];
    out.time[ip] =
        out.time[ip - 1] + out.layer[ip - 1] * (path[ip - 1].ngroup + path[ip].ngroup) / Constant::speed_of_light;
  }
  return out;
}

Vector range_coordinate(const PathRange& range, const ArrayOfPropagationPathPoint& path, const RadarRangeMode mode) {
  switch (mode) {
    case RadarRangeMode::Altitude: {
      Vector out(path.size());
      for (Size ip = 0; ip < path.size(); ++ip) out[ip] = path[ip].pos[0];
      return out;
    }
    case RadarRangeMode::Distance:      return range.distance;
    case RadarRangeMode::RoundTripTime: return range.time;
    case RadarRangeMode::Legacy:        arts::user_error("The Legacy radar range mode must be resolved before use");
  }
  std::unreachable();
}

void add_gate_samples(std::vector<rtepack::radar_sample>& samples,
                      const ConstVectorView&              coordinate,
                      const Numeric                       lo,
                      const Numeric                       hi,
                      const Size                          gate) {
  const Numeric width = hi - lo;
  for (Size i = 1; i < coordinate.size(); i++) {
    const Numeric c0 = coordinate[i - 1];
    const Numeric c1 = coordinate[i];
    const Numeric a  = std::max(lo, std::min(c0, c1));
    const Numeric b  = std::min(hi, std::max(c0, c1));
    if (not(a < b)) continue;
    samples.push_back(
        {.layer = i, .x = (std::midpoint(a, b) - c0) / (c1 - c0), .weight = (b - a) / width, .gate = gate});
  }
}

Size sort_samples(std::vector<rtepack::radar_sample>& samples) {
  stdr::stable_sort(samples, {}, &rtepack::radar_sample::layer);
  return samples.empty() ? 0 : samples.back().layer + 1;
}

Vector gate_coverage(const std::span<const rtepack::radar_sample>& samples, const Size ngate) {
  Vector out(ngate, 0.0);
  for (const auto& sample : samples) out[sample.gate] += sample.weight;
  return out;
}

Vector gate_sum(const std::span<const rtepack::radar_sample>& samples, const ConstVectorView& level, const Size ngate) {
  Vector out(ngate, 0.0);
  for (const auto& [layer, x, weight, gate] : samples)
    out[gate] += weight * ((1.0 - x) * level[layer - 1] + x * level[layer]);
  return out;
}

namespace {
Numeric normalized_delta_aa(Numeric aa) {
  while (aa < -180.0) aa += 360.0;
  while (aa > 180.0) aa -= 360.0;
  return aa;
}

bool at_pole(const Numeric za) { return std::abs(std::sin(za * Constant::pi / 180.0)) < 1e-12; }
}  // namespace

Optics optics(const ArrayOfPropagationPathPoint& path,
              const Size                         np,
              const AtmField&                    atm_field,
              const AscendingGrid&               freq_grid,
              const ArrayOfScatteringSpecies&    scattering_species,
              const JacobianTargets&             jac_targets,
              const Numeric                      pext_scaling,
              const GasModel&                    gas) {
  const Size nf = freq_grid.size();
  const Size nq = jac_targets.target_count();

  if (np > path.size()) arts::user_error("Cannot evaluate {} radar levels of a path with {} points", np, path.size());

  Optics out;
  out.K_out.resize(np);
  out.K_ret.resize(np);
  out.dK_out.resize(np);
  out.dK_ret.resize(np);
  out.Z.resize(np, MuelmatVector(nf, Muelmat{}));
  out.dZ.resize(np, MuelmatMatrix(nq, nf, Muelmat{}));
  out.gas_extinction.resize(nf, np);
  out.particle_extinction.resize(nf, np);
  if (np == 0) return out;

  const ArrayOfAtmPoint atm_path =
      forward_atm_path(ArrayOfPropagationPathPoint(path.begin(), path.begin() + np), atm_field);

  for (Size ip = 0; ip < np; ++ip) {
    // Propagation-path LOS is the propagation direction toward the receiver;
    // its mirror is the transmitted propagation direction.  ARO phase
    // matrices take the latter as incident and the former as scattered.
    const PropagationPathPoint& receive_point  = path[ip];
    PropagationPathPoint        transmit_point = path[ip];
    transmit_point.los                         = path::mirror(receive_point.los);

    PropmatVector gas_out, gas_ret;
    PropmatMatrix gas_out_jac, gas_ret_jac;
    gas(gas_out, gas_out_jac, transmit_point, atm_path[ip]);
    gas(gas_ret, gas_ret_jac, receive_point, atm_path[ip]);
    if (gas_out.size() != nf or gas_ret.size() != nf or static_cast<Size>(gas_out_jac.nrows()) != nq or
        static_cast<Size>(gas_ret_jac.nrows()) != nq or static_cast<Size>(gas_out_jac.ncols()) != nf or
        static_cast<Size>(gas_ret_jac.ncols()) != nf)
      arts::user_error(
          "The gas model must give {} propagation matrices and a [{}, {}] Jacobian, got {}, {:B,}, {}, and {:B,}",
          nf,
          nq,
          nf,
          gas_out.size(),
          gas_out_jac.shape(),
          gas_ret.size(),
          gas_ret_jac.shape());

    // Azimuth is undefined at zenith and nadir.  Use its canonical value at
    // either pole to avoid a 0*NaN rotation without changing the scattering
    // angle or the physical polarization basis.
    const Vector2& transmit_los = transmit_point.los;
    const Vector2& receive_los  = receive_point.los;
    const Vector   delta_aa{at_pole(transmit_los[0]) or at_pole(receive_los[0])
                                ? 0.0
                                : normalized_delta_aa(receive_los[1] - transmit_los[1])};
    const Vector   transmit_za{transmit_los[0]};
    const Vector   receive_za{receive_los[0]};
    const auto     transmit_za_grid =
        std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(transmit_za));
    const auto receive_za_grid =
        std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(receive_za));

    const auto bulk_out = scattering_species.get_bulk_scattering_properties_aro_gridded(
        atm_path[ip], freq_grid, transmit_za, delta_aa, receive_za_grid);
    const auto bulk_ret = scattering_species.get_bulk_scattering_properties_aro_gridded(
        atm_path[ip], freq_grid, receive_za, Vector{0.0}, transmit_za_grid);
    if (not bulk_out.phase_matrix.has_value())
      arts::user_error("Radar single scattering requires phase matrices from all scattering species");

    const auto& phase = *bulk_out.phase_matrix;
    out.K_out[ip]     = gas_out;
    out.K_ret[ip]     = gas_ret;
    for (Size iv = 0; iv < nf; ++iv) {
      const Propmat particle_out       = bulk_out.extinction_matrix[0, iv, 0].expand();
      out.K_out[ip][iv]               += pext_scaling * particle_out;
      out.K_ret[ip][iv]               += pext_scaling * bulk_ret.extinction_matrix[0, iv, 0].expand();
      out.Z[ip][iv]                    = phase[0, iv, 0, 0, 0];
      out.gas_extinction[iv, ip]       = gas_out[iv].A();
      out.particle_extinction[iv, ip]  = pext_scaling * particle_out.A();
    }

    out.dK_out[ip] = std::move(gas_out_jac);
    out.dK_ret[ip] = std::move(gas_ret_jac);
    for (const auto& target : jac_targets.atm) {
      const auto dbulk_out = scattering_species.get_bulk_scattering_properties_aro_gridded_derivative(
          atm_path[ip], freq_grid, transmit_za, delta_aa, receive_za_grid, target.type);
      const auto dbulk_ret = scattering_species.get_bulk_scattering_properties_aro_gridded_derivative(
          atm_path[ip], freq_grid, receive_za, Vector{0.0}, transmit_za_grid, target.type);
      if (not dbulk_out.phase_matrix) arts::user_error("Radar derivatives require scattering phase matrices");

      const auto& dphase = *dbulk_out.phase_matrix;
      for (Size iv = 0; iv < nf; ++iv) {
        out.dK_out[ip][target.target_pos, iv] += pext_scaling * dbulk_out.extinction_matrix[0, iv, 0].expand();
        out.dK_ret[ip][target.target_pos, iv] += pext_scaling * dbulk_ret.extinction_matrix[0, iv, 0].expand();
        out.dZ[ip][target.target_pos, iv]      = dphase[0, iv, 0, 0, 0];
      }
    }
  }

  return out;
}

void single_scattering(StokvecMatrix&                                I,
                       StokvecTensor4&                               dI,
                       const Stokvec&                                transmitted,
                       const Optics&                                 optics,
                       const ConstVectorView&                        layer,
                       const std::span<const rtepack::radar_sample>& samples,
                       const Size                                    ngate,
                       const TransmittanceOption                     option) {
  const Size np = optics.size();
  const Size nf = optics.gas_extinction.nrows();
  const Size nq = np == 0 ? 0 : optics.dK_out.front().nrows();

  if (np > layer.size() + 1)
    arts::user_error("The radar optics have {} levels but there are only {} layers", np, layer.size());

  I.resize(ngate, nf);
  dI.resize(ngate, np, nq, nf);
  rtepack::radar_single_scattering(I,
                                   dI,
                                   transmitted,
                                   optics.K_out,
                                   optics.dK_out,
                                   optics.K_ret,
                                   optics.dK_ret,
                                   optics.Z,
                                   optics.dZ,
                                   layer[Range(0, np == 0 ? 0 : np - 1)],
                                   samples,
                                   option);
}
}  // namespace radar
