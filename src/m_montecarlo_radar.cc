#include <arts_constants.h>
#include <atm_path.h>
#include <mc_antenna.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
Numeric normalized_delta_aa(Numeric aa) {
  while (aa < -180.0) aa += 360.0;
  while (aa > 180.0) aa -= 360.0;
  return aa;
}

Numeric ze_factor(const Numeric frequency, const Numeric k2) {
  const Numeric wavelength = Constant::speed_of_light / frequency;
  return 4e18 * std::pow(wavelength, 4) / (std::pow(Constant::pi, 4) * k2);
}

struct OpticalPoint {
  Propmat outgoing_extinction{};
  Propmat return_extinction{};
  Muelmat backscatter{0.0};
};

OpticalPoint optical_point(const Workspace&                ws,
                           const Numeric                   frequency,
                           const PropagationPathPoint&     ray_point,
                           const AtmPoint&                 atm_point,
                           const ArrayOfScatteringSpecies& scattering_species,
                           const Agenda&                   spectral_propmat_agenda) {
  const AscendingGrid freq_grid{frequency};

  PropmatVector gas_return;
  StokvecVector src_return;
  PropmatMatrix gas_return_jac;
  StokvecMatrix src_return_jac;
  spectral_propmat_agendaExecute(ws,
                                 gas_return,
                                 src_return,
                                 gas_return_jac,
                                 src_return_jac,
                                 freq_grid,
                                 Vector3{0.0, 0.0, 0.0},
                                 JacobianTargets{},
                                 {},
                                 ray_point,
                                 atm_point,
                                 spectral_propmat_agenda);

  PropagationPathPoint outgoing_point = ray_point;
  outgoing_point.los                  = path::mirror(ray_point.los);
  PropmatVector gas_outgoing;
  StokvecVector src_outgoing;
  PropmatMatrix gas_outgoing_jac;
  StokvecMatrix src_outgoing_jac;
  spectral_propmat_agendaExecute(ws,
                                 gas_outgoing,
                                 src_outgoing,
                                 gas_outgoing_jac,
                                 src_outgoing_jac,
                                 freq_grid,
                                 Vector3{0.0, 0.0, 0.0},
                                 JacobianTargets{},
                                 {},
                                 outgoing_point,
                                 atm_point,
                                 spectral_propmat_agenda);

  const Vector2 transmit_los = outgoing_point.los;
  const Vector2 receive_los  = ray_point.los;
  const auto    at_pole      = [](const Numeric za) { return std::abs(std::sin(za * Constant::pi / 180.0)) < 1e-12; };
  const Vector  delta_aa{at_pole(transmit_los[0]) or at_pole(receive_los[0])
                             ? 0.0
                             : normalized_delta_aa(receive_los[1] - transmit_los[1])};
  auto          receive_za =
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{receive_los[0]}));
  const auto bulk_outgoing = scattering_species.get_bulk_scattering_properties_aro_gridded(
      atm_point, freq_grid, Vector{transmit_los[0]}, delta_aa, receive_za);
  auto transmit_za =
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{transmit_los[0]}));
  const auto bulk_return = scattering_species.get_bulk_scattering_properties_aro_gridded(
      atm_point, freq_grid, Vector{receive_los[0]}, Vector{0.0}, transmit_za);

  if (not bulk_outgoing.phase_matrix.has_value())
    arts::user_error("MCRadar requires phase matrices from all scattering species");

  const auto& phase = *bulk_outgoing.phase_matrix;

  OpticalPoint out;
  out.outgoing_extinction  = gas_outgoing[0];
  out.outgoing_extinction += bulk_outgoing.extinction_matrix[0, 0, 0].expand();
  out.return_extinction    = gas_return[0];
  out.return_extinction   += bulk_return.extinction_matrix[0, 0, 0].expand();
  out.backscatter          = phase[0, 0, 0, 0, 0];
  return out;
}

}  // namespace

void MCRadar(const Workspace&                ws,
             StokvecVector&                  radar_signal,
             StokvecVector&                  radar_error,
             const AtmField&                 atm_field,
             const SurfaceField&             surf_field,
             const ArrayOfScatteringSpecies& scattering_species,
             const MCAntenna&                mc_antenna,
             const Agenda&                   ray_path_observer_agenda,
             const Agenda&                   spectral_propmat_agenda,
             const Numeric&                  frequency,
             const Vector3&                  sensor_pos,
             const Vector2&                  sensor_los,
             const Stokvec&                  mc_y_tx,
             const AscendingGrid&            range_bins,
             const Index&                    mc_seed,
             const Index&                    mc_max_iter,
             const Index&                    mc_max_scatorder,
             const Numeric&                  k2,
             const String&                   unit) try {
  ARTS_TIME_REPORT

  if (frequency <= 0.0) arts::user_error("frequency must be positive");
  if (range_bins.size() < 2) arts::user_error("range_bins must contain at least two edges");
  if (range_bins.front() < 0.0) arts::user_error("range_bins cannot contain negative distances");
  if (mc_max_iter <= 0) arts::user_error("mc_max_iter must be positive");
  if (mc_max_scatorder != 1)
    arts::user_error("This first ARTS3 MCRadar implementation supports mc_max_scatorder=1; got {}", mc_max_scatorder);
  if (scattering_species.species.empty()) arts::user_error("MCRadar requires at least one scattering species");
  if (unit != "1" and unit != "Ze") arts::user_error(R"(unit must be either "1" or "Ze")");
  if (k2 <= 0.0) arts::user_error("k2 must be positive");

  const Index nbins = static_cast<Index>(range_bins.size()) - 1;
  radar_signal.resize(nbins);
  radar_error.resize(nbins);
  std::fill(radar_signal.begin(), radar_signal.end(), Stokvec{0.0});
  std::fill(radar_error.begin(), radar_error.end(), Stokvec{0.0});
  StokvecVector squared_sum(nbins, Stokvec{0.0});

  RandomNumberGenerator<> rng(mc_seed);
  // The general RNG avoids reusing seeds process-wide.  A forward model with
  // an explicit seed instead promises reproducibility across repeated calls.
  rng.force_seed(mc_seed);
  const Matrix33 ant_to_enu = rotmat_enu(sensor_los);
  Matrix33       enu_to_ant;
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j) enu_to_ant[i, j] = ant_to_enu[j, i];

  for (Index iter = 0; iter < mc_max_iter; ++iter) {
    StokvecVector photon(nbins, Stokvec{0.0});
    const auto [sampled_los, ray_rotation] = mc_antenna.draw_los(rng, ant_to_enu, sensor_los);
    const Numeric antenna_weight           = mc_antenna.return_los(ray_rotation, enu_to_ant);
    const Muelmat tx_rotation              = rotmat_stokes(1.0, 1.0, ant_to_enu, ray_rotation);
    const Muelmat rx_rotation              = rotmat_stokes(-1.0, 1.0, ray_rotation, ant_to_enu);
    const Stokvec transmitted              = tx_rotation * mc_y_tx;

    ArrayOfPropagationPathPoint ray_path;
    ray_path_observer_agendaExecute(ws, ray_path, sensor_pos, sampled_los, ray_path_observer_agenda);
    if (ray_path.size() < 2) continue;

    const ArrayOfAtmPoint atm_path  = forward_atm_path(ray_path, atm_field);
    const Vector          distances = path::distance(ray_path, surf_field.ellipsoid);
    Muelmat               outbound  = Muelmat::id();
    Muelmat               returning = Muelmat::id();
    Numeric               travelled = 0.0;
    OpticalPoint          previous =
        optical_point(ws, frequency, ray_path.front(), atm_path.front(), scattering_species, spectral_propmat_agenda);

    for (Size ip = 1; ip < ray_path.size(); ++ip) {
      const Numeric ds = distances[ip - 1];
      if (ds <= 0.0) continue;

      const OpticalPoint current =
          optical_point(ws, frequency, ray_path[ip], atm_path[ip], scattering_species, spectral_propmat_agenda);
      const Muelmat outgoing_segment = rtepack::tran(previous.outgoing_extinction, current.outgoing_extinction, ds)();
      const Muelmat return_segment   = rtepack::tran(current.return_extinction, previous.return_extinction, ds)();
      // Integrate only the overlap with each gate.  Interpolate optical
      // properties to its midpoint, including the extinction along the
      // partial segment between the previous path point and that midpoint.
      const auto  first_edge = stdr::upper_bound(range_bins, travelled);
      const Index first_bin  = std::max<Index>(0, first_edge - range_bins.begin() - 1);
      for (Index ibin = first_bin; ibin < nbins and range_bins[ibin] < travelled + ds; ++ibin) {
        const Numeric begin = std::max(travelled, range_bins[ibin]);
        const Numeric end   = std::min(travelled + ds, range_bins[ibin + 1]);
        if (end <= begin) continue;
        const Numeric midpoint = 0.5 * (begin + end) - travelled;
        const Numeric fraction = midpoint / ds;
        const Propmat outgoing_extinction =
            (1.0 - fraction) * previous.outgoing_extinction + fraction * current.outgoing_extinction;
        const Propmat return_extinction =
            (1.0 - fraction) * previous.return_extinction + fraction * current.return_extinction;
        const Muelmat outgoing_midpoint =
            rtepack::tran(previous.outgoing_extinction, outgoing_extinction, midpoint)() * outbound;
        const Muelmat return_midpoint =
            returning * rtepack::tran(return_extinction, previous.return_extinction, midpoint)();
        const Muelmat backscatter = (1.0 - fraction) * previous.backscatter + fraction * current.backscatter;
        const Stokvec contribution =
            antenna_weight * (end - begin) *
            (rx_rotation * rtepack::radar_return(return_midpoint, backscatter, outgoing_midpoint, transmitted));
        photon[ibin] += contribution;
      }

      travelled += ds;
      if (travelled >= range_bins.back()) break;
      outbound  = outgoing_segment * outbound;
      returning = returning * return_segment;
      previous  = current;
    }

    for (Index b = 0; b < nbins; ++b) {
      const Numeric width = range_bins[b + 1] - range_bins[b];
      for (Index s = 0; s < 4; ++s) {
        const Numeric value  = photon[b][s] / width;
        radar_signal[b][s]  += value;
        squared_sum[b][s]   += value * value;
      }
    }
  }

  const Numeric factor = unit == "Ze" ? ze_factor(frequency, k2) / (2.0 * Constant::pi) : 1.0;
  for (Index b = 0; b < nbins; ++b) {
    for (Index s = 0; s < 4; ++s) {
      const Numeric mean     = radar_signal[b][s] / static_cast<Numeric>(mc_max_iter);
      const Numeric variance = std::max(0.0, squared_sum[b][s] / static_cast<Numeric>(mc_max_iter) - mean * mean);
      radar_signal[b][s]     = factor * mean;
      radar_error[b][s]      = factor * std::sqrt(variance / static_cast<Numeric>(mc_max_iter));
    }
  }
}
ARTS_METHOD_ERROR_CATCH
