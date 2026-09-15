#include <arts_constants.h>
#include <atm_path.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <sensor_meta_info.h>
#include <time_report.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <string_view>
#include <vector>

namespace {

Numeric liebe93_k2(const Numeric frequency, const Numeric temperature) {
  ARTS_USER_ERROR_IF(
      temperature < Constant::temperature_at_0c - 40.0 or temperature > Constant::temperature_at_0c + 100.0,
      "Automatic radar k2 requires {} K <= ze_tref <= {} K; got {} K",
      Constant::temperature_at_0c - 40.0,
      Constant::temperature_at_0c + 100.0,
      temperature)
  ARTS_USER_ERROR_IF(frequency < 10e9 or frequency > 1000e9,
                     "Automatic radar k2 uses Liebe-93 and requires 10 GHz <= frequency <= 1000 GHz; got {} Hz",
                     frequency)

  const Numeric theta = 1.0 - 300.0 / temperature;
  const Numeric e0    = 77.66 - 103.3 * theta;
  const Numeric e1    = 0.0671 * e0;
  const Numeric f1    = 20.2 + 146.0 * theta + 316.0 * theta * theta;
  const Numeric e2    = 3.52;
  const Numeric f2    = 39.8 * f1;
  const Complex ifGHz{0.0, frequency / 1e9};
  const Complex n  = std::sqrt(e2 + (e1 - e2) / (1.0 - ifGHz / f2) + (e0 - e1) / (1.0 - ifGHz / f1));
  const Complex n2 = n * n;
  return std::norm((n2 - 1.0) / (n2 + 2.0));
}

Numeric ze_factor(const Numeric frequency, const Numeric ze_tref, const Numeric k2) {
  const Numeric dielectric = k2 < 0.0 ? liebe93_k2(frequency, ze_tref) : k2;
  ARTS_USER_ERROR_IF(dielectric <= 0.0, "Radar k2 must be positive, or negative to request Liebe-93; got {}", k2)
  const Numeric wavelength = Constant::speed_of_light / frequency;
  return 4e18 * std::pow(wavelength, 4) / (std::pow(Constant::pi, 4) * dielectric);
}

struct OpticalProfile {
  std::vector<PropmatVector> outgoing_extinction;
  std::vector<PropmatVector> return_extinction;
  std::vector<PropmatVector> gas_extinction;
  std::vector<PropmatVector> particle_extinction;
  std::vector<MuelmatVector> backscatter;
  std::vector<PropmatMatrix> outgoing_extinction_jac;
  std::vector<PropmatMatrix> return_extinction_jac;
  std::vector<MuelmatMatrix> backscatter_jac;
};

Numeric normalized_delta_aa(Numeric aa) {
  while (aa < -180.0) aa += 360.0;
  while (aa > 180.0) aa -= 360.0;
  return aa;
}

OpticalProfile optical_profile(const Workspace&                   ws,
                               const AscendingGrid&               frequency,
                               const ArrayOfPropagationPathPoint& ray_path,
                               const AtmField&                    atm_field,
                               const ArrayOfScatteringSpecies&    scattering_species,
                               const Agenda&                      spectral_propmat_agenda,
                               const JacobianTargets&             jac_targets,
                               const Numeric                      pext_scaling) {
  const Size np = ray_path.size();
  const Size nf = frequency.size();
  const Size nq = jac_targets.target_count();

  OpticalProfile out;
  out.outgoing_extinction.resize(np, PropmatVector(nf));
  out.return_extinction.resize(np, PropmatVector(nf));
  out.gas_extinction.resize(np, PropmatVector(nf));
  out.particle_extinction.resize(np, PropmatVector(nf));
  out.backscatter.resize(np, MuelmatVector(nf, Muelmat{0.0}));
  out.outgoing_extinction_jac.resize(np, PropmatMatrix(nq, nf));
  out.return_extinction_jac.resize(np, PropmatMatrix(nq, nf));
  out.backscatter_jac.resize(np, MuelmatMatrix(nq, nf, Muelmat{0.0}));

  const ArrayOfAtmPoint atm_path = forward_atm_path(ray_path, atm_field);

  for (Size ip = 0; ip < np; ++ip) {
    PropmatVector gas_return;
    StokvecVector source_return;
    PropmatMatrix gas_return_jac;
    StokvecMatrix source_return_jac;
    spectral_propmat_agendaExecute(ws,
                                   gas_return,
                                   source_return,
                                   gas_return_jac,
                                   source_return_jac,
                                   frequency,
                                   Vector3{0.0, 0.0, 0.0},
                                   jac_targets,
                                   {},
                                   ray_path[ip],
                                   atm_path[ip],
                                   spectral_propmat_agenda);
    PropagationPathPoint outgoing_point = ray_path[ip];
    outgoing_point.los                  = path::mirror(outgoing_point.los);
    PropmatVector gas_outgoing;
    StokvecVector source_outgoing;
    PropmatMatrix gas_outgoing_jac;
    StokvecMatrix source_outgoing_jac;
    spectral_propmat_agendaExecute(ws,
                                   gas_outgoing,
                                   source_outgoing,
                                   gas_outgoing_jac,
                                   source_outgoing_jac,
                                   frequency,
                                   Vector3{0.0, 0.0, 0.0},
                                   jac_targets,
                                   {},
                                   outgoing_point,
                                   atm_path[ip],
                                   spectral_propmat_agenda);

    // Propagation-path LOS is the propagation direction toward the receiver;
    // its mirror is the transmitted propagation direction.  ARO phase
    // matrices take the latter as incident and the former as scattered.
    const Vector2 receive_los  = ray_path[ip].los;
    const Vector2 transmit_los = outgoing_point.los;
    // Azimuth is undefined at zenith and nadir.  Use its canonical value at
    // either pole to avoid a 0*NaN rotation without changing the scattering
    // angle or the physical polarization basis.
    const auto   at_pole = [](const Numeric za) { return std::abs(std::sin(za * Constant::pi / 180.0)) < 1e-12; };
    const Vector delta_aa{at_pole(transmit_los[0]) or at_pole(receive_los[0])
                              ? 0.0
                              : normalized_delta_aa(receive_los[1] - transmit_los[1])};
    auto         receive_za =
        std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{receive_los[0]}));
    const auto bulk_outgoing = scattering_species.get_bulk_scattering_properties_aro_gridded(
        atm_path[ip], frequency, Vector{transmit_los[0]}, delta_aa, receive_za);
    auto transmit_za =
        std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{transmit_los[0]}));
    const auto bulk_return = scattering_species.get_bulk_scattering_properties_aro_gridded(
        atm_path[ip], frequency, Vector{receive_los[0]}, Vector{0.0}, transmit_za);
    ARTS_USER_ERROR_IF(not bulk_outgoing.phase_matrix.has_value(),
                       "Radar single scattering requires phase matrices from all scattering species")

    const auto& phase                      = *bulk_outgoing.phase_matrix;
    const auto  phase_coeffs               = phase.get_const_coeff_vector_view();
    const auto  outgoing_extinction_coeffs = bulk_outgoing.extinction_matrix.get_const_coeff_vector_view();
    const auto  return_extinction_coeffs   = bulk_return.extinction_matrix.get_const_coeff_vector_view();
    for (Size iv = 0; iv < nf; ++iv) {
      const auto& outgoing             = outgoing_extinction_coeffs[0, iv, 0];
      const auto& returning            = return_extinction_coeffs[0, iv, 0];
      out.gas_extinction[ip][iv]       = gas_outgoing[iv];
      out.particle_extinction[ip][iv]  = Propmat{outgoing[0], outgoing[1], 0.0, 0.0, 0.0, 0.0, outgoing[2]};
      out.outgoing_extinction[ip][iv]  = gas_outgoing[iv];
      out.outgoing_extinction[ip][iv] += pext_scaling * out.particle_extinction[ip][iv];
      out.return_extinction[ip][iv]    = gas_return[iv];
      out.return_extinction[ip][iv] +=
          pext_scaling * Propmat{returning[0], returning[1], 0.0, 0.0, 0.0, 0.0, returning[2]};
      out.backscatter[ip][iv] = phase_coeffs[0, iv, 0, 0, 0];
    }
    for (const auto& target : jac_targets.atm) {
      auto receive_za_jac =
          std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{receive_los[0]}));
      const auto dbulk_outgoing = scattering_species.get_bulk_scattering_properties_aro_gridded_derivative(
          atm_path[ip], frequency, Vector{transmit_los[0]}, delta_aa, receive_za_jac, target.type);
      auto transmit_za_jac =
          std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{transmit_los[0]}));
      const auto dbulk_return = scattering_species.get_bulk_scattering_properties_aro_gridded_derivative(
          atm_path[ip], frequency, Vector{receive_los[0]}, Vector{0.0}, transmit_za_jac, target.type);
      ARTS_USER_ERROR_IF(not dbulk_outgoing.phase_matrix, "Radar derivatives require scattering phase matrices")
      const auto dphase_coeffs               = dbulk_outgoing.phase_matrix->get_const_coeff_vector_view();
      const auto doutgoing_extinction_coeffs = dbulk_outgoing.extinction_matrix.get_const_coeff_vector_view();
      const auto dreturn_extinction_coeffs   = dbulk_return.extinction_matrix.get_const_coeff_vector_view();
      for (Size iv = 0; iv < nf; ++iv) {
        out.outgoing_extinction_jac[ip][target.target_pos, iv] = gas_outgoing_jac[target.target_pos, iv];
        const auto& doutgoing                                  = doutgoing_extinction_coeffs[0, iv, 0];
        const auto& dreturning                                 = dreturn_extinction_coeffs[0, iv, 0];
        out.outgoing_extinction_jac[ip][target.target_pos, iv] +=
            pext_scaling * Propmat{doutgoing[0], doutgoing[1], 0.0, 0.0, 0.0, 0.0, doutgoing[2]};
        out.return_extinction_jac[ip][target.target_pos, iv] = gas_return_jac[target.target_pos, iv];
        out.return_extinction_jac[ip][target.target_pos, iv] +=
            pext_scaling * Propmat{dreturning[0], dreturning[1], 0.0, 0.0, 0.0, 0.0, dreturning[2]};
        out.backscatter_jac[ip][target.target_pos, iv] = dphase_coeffs[0, iv, 0, 0, 0];
      }
    }
  }

  return out;
}

enum class RadarRangeMode : char { Altitude, Distance, RoundTripTime };

RadarRangeMode range_mode(const String& mode, const Matrix& limits) {
  if (mode == "Altitude") return RadarRangeMode::Altitude;
  if (mode == "Distance") return RadarRangeMode::Distance;
  if (mode == "RoundTripTime") return RadarRangeMode::RoundTripTime;
  if (mode == "Legacy") {
    Numeric largest = -std::numeric_limits<Numeric>::infinity();
    for (Index i = 0; i < limits.nrows(); ++i)
      for (Index j = 0; j < limits.ncols(); ++j) largest = std::max(largest, limits[i, j]);
    return largest > 1.0 ? RadarRangeMode::Altitude : RadarRangeMode::RoundTripTime;
  }
  ARTS_USER_ERROR("Unknown radar range_mode \"{}\". Expected Altitude, Distance, RoundTripTime, or Legacy.", mode)
}

Vector range_coordinate(const ArrayOfPropagationPathPoint& ray_path,
                        const Vector3&                     observer,
                        const Vector2&                     ellipsoid,
                        const RadarRangeMode               mode) {
  const Size np = ray_path.size();
  Vector     out(np);
  if (mode == RadarRangeMode::Altitude) {
    for (Size ip = 0; ip < np; ++ip) out[ip] = ray_path[ip].pos[0];
    return out;
  }

  const Numeric initial = path::distance(observer, ray_path.front().pos, ellipsoid);
  if (mode == RadarRangeMode::Distance) {
    out[0]               = initial;
    const Vector segment = path::distance(ray_path, ellipsoid);
    for (Size ip = 1; ip < np; ++ip) out[ip] = out[ip - 1] + segment[ip - 1];
    return out;
  }

  out[0]               = 2.0 * initial / Constant::speed_of_light;
  const Vector segment = path::distance(ray_path, ellipsoid);
  for (Size ip = 1; ip < np; ++ip) {
    out[ip] =
        out[ip - 1] + segment[ip - 1] * (ray_path[ip - 1].ngroup + ray_path[ip].ngroup) / Constant::speed_of_light;
  }
  return out;
}

template <class T> std::optional<T> bin_average(const Vector& x, const std::vector<T>& y, Numeric lo, Numeric hi) {
  ARTS_USER_ERROR_IF(x.size() != y.size(), "Internal radar range/profile size mismatch")
  if (x.size() < 2) return std::nullopt;

  const bool    increasing = x.back() > x.front();
  const Numeric xmin       = std::min(x.front(), x.back());
  const Numeric xmax       = std::max(x.front(), x.back());
  lo                       = std::max(lo, xmin);
  hi                       = std::min(hi, xmax);
  if (not(lo < hi)) return std::nullopt;

  auto interpolate = [&](const Numeric q) {
    Size i = 0;
    if (increasing) {
      while (i + 1 < x.size() and x[i + 1] < q) ++i;
    } else {
      while (i + 1 < x.size() and x[i + 1] > q) ++i;
    }
    if (i + 1 >= x.size()) return y.back();
    const Numeric w = (q - x[i]) / (x[i + 1] - x[i]);
    return (1.0 - w) * y[i] + w * y[i + 1];
  };

  std::vector<Numeric> q{lo};
  q.reserve(x.size() + 2);
  for (auto v : x)
    if (v > lo and v < hi) q.push_back(v);
  q.push_back(hi);
  std::ranges::sort(q);

  T integral{0.0};
  T past = interpolate(q.front());
  for (Size i = 1; i < q.size(); ++i) {
    const T now  = interpolate(q[i]);
    integral    += 0.5 * (q[i] - q[i - 1]) * (past + now);
    past         = now;
  }
  return integral * (1.0 / (hi - lo));
}

struct PreparedSimulation {
  const AscendingGrid*        frequency{};
  const SensorPosLosVector*   poslos{};
  Size                        iposlos{};
  ArrayOfPropagationPathPoint ray_path{};
};

struct RadarResult {
  Vector measurement;
  Matrix jacobian;
  Matrix auxiliary;
};

RadarResult radar_forward(const Workspace&                       ws,
                          const ArrayOfSensorObsel&              sensor,
                          const Matrix&                          limits,
                          const std::vector<PreparedSimulation>& simulations,
                          const AtmField&                        atm_field,
                          const SurfaceField&                    surf_field,
                          const ArrayOfScatteringSpecies&        scattering_species,
                          const Agenda&                          spectral_propmat_agenda,
                          const Stokvec&                         transmitted,
                          const RadarRangeMode                   coordinate_mode,
                          const String&                          unit,
                          const Numeric                          ze_tref,
                          const Numeric                          k2,
                          const Numeric                          dbze_min,
                          const Numeric                          pext_scaling,
                          const ArrayOfString&                   aux_vars,
                          const JacobianTargets&                 jac_targets) {
  const Size        ny = sensor.size();
  const Size        nq = jac_targets.target_count();
  RadarResult       out{.measurement = Vector(ny, 0.0),
                        .jacobian    = Matrix(ny, jac_targets.x_size(), 0.0),
                        .auxiliary   = Matrix(aux_vars.size(), ny, 0.0)};
  std::vector<bool> seen(ny, false);
  Vector            scalar_weight_sum(ny, 0.0);

  for (const auto& simulation : simulations) {
    const auto& frequency = *simulation.frequency;
    const auto& poslos    = *simulation.poslos;
    const auto& ray_path  = simulation.ray_path;
    const Size  iposlos   = simulation.iposlos;
    if (ray_path.size() < 2) continue;

    const OpticalProfile optical = optical_profile(
        ws, frequency, ray_path, atm_field, scattering_species, spectral_propmat_agenda, jac_targets, pext_scaling);
    const Size np = ray_path.size();
    const Size nf = frequency.size();

    std::vector<StokvecVector>              attenuated(np, StokvecVector(nf));
    std::vector<StokvecVector>              unattenuated(np, StokvecVector(nf));
    std::vector<StokvecVector>              gas_ext(np, StokvecVector(nf));
    std::vector<StokvecVector>              particle_ext(np, StokvecVector(nf));
    std::vector<Muelmat>                    outgoing_transmission(nf, Muelmat::id());
    std::vector<Muelmat>                    return_transmission(nf, Muelmat::id());
    std::vector<MuelmatMatrix>              outgoing_transmission_jac(nf, MuelmatMatrix(np, nq, Muelmat{0.0}));
    std::vector<MuelmatMatrix>              return_transmission_jac(nf, MuelmatMatrix(np, nq, Muelmat{0.0}));
    std::vector<std::vector<StokvecVector>> attenuated_jac(
        np * nq, std::vector<StokvecVector>(np, StokvecVector(nf, Stokvec{0.0})));

    for (Size ip = 0; ip < np; ++ip) {
      if (ip > 0) {
        const Numeric ds = path::distance(ray_path[ip - 1].pos, ray_path[ip].pos, surf_field.ellipsoid);
        for (Size iv = 0; iv < nf; ++iv) {
          const Muelmat       previous_outgoing = outgoing_transmission[iv];
          const Muelmat       previous_return   = return_transmission[iv];
          const rtepack::tran outgoing_segment(
              optical.outgoing_extinction[ip - 1][iv], optical.outgoing_extinction[ip][iv], ds);
          const rtepack::tran return_segment(
              optical.return_extinction[ip][iv], optical.return_extinction[ip - 1][iv], ds);
          const Muelmat outgoing_segment_transmission = outgoing_segment();
          const Muelmat return_segment_transmission   = return_segment();
          for (Size jp = 0; jp < np; ++jp)
            for (Size iq = 0; iq < nq; ++iq) {
              outgoing_transmission_jac[iv][jp, iq] =
                  outgoing_segment_transmission * outgoing_transmission_jac[iv][jp, iq];
              return_transmission_jac[iv][jp, iq] = return_transmission_jac[iv][jp, iq] * return_segment_transmission;
            }
          for (Size iq = 0; iq < nq; ++iq) {
            outgoing_transmission_jac[iv][ip - 1, iq] +=
                outgoing_segment.deriv(outgoing_segment_transmission,
                                       optical.outgoing_extinction[ip - 1][iv],
                                       optical.outgoing_extinction[ip][iv],
                                       optical.outgoing_extinction_jac[ip - 1][iq, iv],
                                       ds,
                                       0.0) *
                previous_outgoing;
            outgoing_transmission_jac[iv][ip, iq] += outgoing_segment.deriv(outgoing_segment_transmission,
                                                                            optical.outgoing_extinction[ip - 1][iv],
                                                                            optical.outgoing_extinction[ip][iv],
                                                                            optical.outgoing_extinction_jac[ip][iq, iv],
                                                                            ds,
                                                                            0.0) *
                                                     previous_outgoing;
            return_transmission_jac[iv][ip - 1, iq] +=
                previous_return * return_segment.deriv(return_segment_transmission,
                                                       optical.return_extinction[ip][iv],
                                                       optical.return_extinction[ip - 1][iv],
                                                       optical.return_extinction_jac[ip - 1][iq, iv],
                                                       ds,
                                                       0.0);
            return_transmission_jac[iv][ip, iq] +=
                previous_return * return_segment.deriv(return_segment_transmission,
                                                       optical.return_extinction[ip][iv],
                                                       optical.return_extinction[ip - 1][iv],
                                                       optical.return_extinction_jac[ip][iq, iv],
                                                       ds,
                                                       0.0);
          }
          outgoing_transmission[iv] = outgoing_segment_transmission * previous_outgoing;
          return_transmission[iv]   = previous_return * return_segment_transmission;
        }
      }
      for (Size iv = 0; iv < nf; ++iv) {
        unattenuated[ip][iv] = optical.backscatter[ip][iv] * transmitted;
        attenuated[ip][iv]   = rtepack::radar_return(
            return_transmission[iv], optical.backscatter[ip][iv], outgoing_transmission[iv], transmitted);
        gas_ext[ip][iv]      = Stokvec{optical.gas_extinction[ip][iv].A(), 0.0, 0.0, 0.0};
        particle_ext[ip][iv] = Stokvec{pext_scaling * optical.particle_extinction[ip][iv].A(), 0.0, 0.0, 0.0};
        for (Size jp = 0; jp < np; ++jp) {
          for (Size iq = 0; iq < nq; ++iq) {
            const Muelmat dbackscatter = jp == ip ? optical.backscatter_jac[ip][iq, iv] : Muelmat{0.0};
            const Stokvec derivative   = rtepack::radar_return_derivative(return_transmission[iv],
                                                                          return_transmission_jac[iv][jp, iq],
                                                                          optical.backscatter[ip][iv],
                                                                          dbackscatter,
                                                                          outgoing_transmission[iv],
                                                                          outgoing_transmission_jac[iv][jp, iq],
                                                                          transmitted);
            attenuated_jac[jp * nq + iq][ip][iv] = derivative;
          }
        }
      }
    }

    const Vector  coordinate = range_coordinate(ray_path, poslos[iposlos].pos, surf_field.ellipsoid, coordinate_mode);
    const Numeric background = ray_path.back().los_type == PathPositionType::space
                                   ? 0.0
                                   : (ray_path.back().los_type == PathPositionType::surface ? 1.0 : 2.0);

    for (Size iy = 0; iy < ny; ++iy) {
      const auto& obsel = sensor[iy];
      if (not obsel.same_freqs(frequency)) continue;

      StokvecVector              bin_signal(nf, Stokvec{0.0});
      StokvecVector              bin_backscatter(nf, Stokvec{0.0});
      StokvecVector              bin_gas(nf, Stokvec{0.0});
      StokvecVector              bin_particle(nf, Stokvec{0.0});
      std::vector<StokvecVector> bin_jac(np * nq, StokvecVector(nf, Stokvec{0.0}));
      bool                       inside = false;
      for (Size iv = 0; iv < nf; ++iv) {
        std::vector<Stokvec> a(np), b(np), g(np), p(np);
        for (Size ip = 0; ip < np; ++ip) {
          a[ip] = attenuated[ip][iv];
          b[ip] = unattenuated[ip][iv];
          g[ip] = gas_ext[ip][iv];
          p[ip] = particle_ext[ip][iv];
        }
        if (auto avg = bin_average(coordinate, a, limits[iy, 0], limits[iy, 1])) {
          inside              = true;
          bin_signal[iv]      = *avg;
          bin_backscatter[iv] = *bin_average(coordinate, b, limits[iy, 0], limits[iy, 1]);
          bin_gas[iv]         = *bin_average(coordinate, g, limits[iy, 0], limits[iy, 1]);
          bin_particle[iv]    = *bin_average(coordinate, p, limits[iy, 0], limits[iy, 1]);

          const Numeric factor  = unit == "1" ? 1.0 : ze_factor(frequency[iv], ze_tref, k2);
          bin_signal[iv]       *= factor;
          bin_backscatter[iv]  *= factor;
          for (Size jp = 0; jp < np; ++jp) {
            for (Size iq = 0; iq < nq; ++iq) {
              std::vector<Stokvec> derivative(np);
              for (Size kp = 0; kp < np; ++kp) derivative[kp] = attenuated_jac[jp * nq + iq][kp][iv];
              bin_jac[jp * nq + iq][iv] = factor * *bin_average(coordinate, derivative, limits[iy, 0], limits[iy, 1]);
            }
          }
        }
      }
      if (not inside) {
        if (not seen[iy]) {
          out.measurement[iy] = std::numeric_limits<Numeric>::quiet_NaN();
          for (Size ia = 0; ia < aux_vars.size(); ++ia)
            out.auxiliary[ia, iy] = std::numeric_limits<Numeric>::quiet_NaN();
        }
        continue;
      }

      if (not seen[iy] or std::isnan(out.measurement[iy])) {
        out.measurement[iy] = 0.0;
        for (Size ia = 0; ia < aux_vars.size(); ++ia) out.auxiliary[ia, iy] = 0.0;
      }
      const Numeric signal  = obsel.sumup(bin_signal, iposlos);
      out.measurement[iy]  += signal;
      for (const auto& target : jac_targets.atm) {
        const auto& field = atm_field[target.type];
        for (Size jp = 0; jp < np; ++jp) {
          const Numeric local = obsel.sumup(bin_jac[jp * nq + target.target_pos], iposlos);
          for (const auto& [index, weight] : field.flat_weight(ray_path[jp].pos))
            out.jacobian[iy, target.x_start + index] += weight * local;
        }
      }

      const Numeric scalar_weight  = obsel.sumup(Stokvec{1.0, 0.0, 0.0, 0.0}, iposlos);
      scalar_weight_sum[iy]       += scalar_weight;
      for (Size ia = 0; ia < aux_vars.size(); ++ia) {
        if (aux_vars[ia] == "Radiative background") {
          out.auxiliary[ia, iy] += scalar_weight * background;
        } else if (aux_vars[ia] == "Backscattering") {
          out.auxiliary[ia, iy] += obsel.sumup(bin_backscatter, iposlos);
        } else if (aux_vars[ia] == "Abs species extinction") {
          out.auxiliary[ia, iy] += obsel.sumup(bin_gas, iposlos);
        } else if (aux_vars[ia] == "Particle extinction") {
          out.auxiliary[ia, iy] += obsel.sumup(bin_particle, iposlos);
        }
      }
      seen[iy] = true;
    }
  }

  for (Size ia = 0; ia < aux_vars.size(); ++ia) {
    if (aux_vars[ia] == "Backscattering") continue;
    for (Size iy = 0; iy < ny; ++iy)
      if (seen[iy] and scalar_weight_sum[iy] != 0.0) out.auxiliary[ia, iy] /= scalar_weight_sum[iy];
  }

  if (unit == "dBZe") {
    const Numeric ze_min = std::pow(10.0, dbze_min / 10.0);
    for (Size iy = 0; iy < ny; ++iy) {
      if (std::isnan(out.measurement[iy])) continue;
      if (out.measurement[iy] <= ze_min) {
        out.measurement[iy]     = dbze_min;
        out.jacobian[iy, joker] = 0.0;
      } else {
        const Numeric factor     = 10.0 / (std::log(10.0) * out.measurement[iy]);
        out.jacobian[iy, joker] *= factor;
        out.measurement[iy]      = 10.0 * std::log10(out.measurement[iy]);
      }
    }
    for (Size ia = 0; ia < aux_vars.size(); ++ia) {
      if (aux_vars[ia] != "Backscattering") continue;
      for (Size iy = 0; iy < ny; ++iy) {
        auto& value = out.auxiliary[ia, iy];
        if (not std::isnan(value)) value = value <= ze_min ? dbze_min : 10.0 * std::log10(value);
      }
    }
  }

  return out;
}

}  // namespace

void measurement_sensorAddSimpleRadar(ArrayOfSensorObsel&    measurement_sensor,
                                      ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                      Matrix&                radar_range_limits,
                                      const AscendingGrid&   freq_grid,
                                      const Vector3&         pos,
                                      const Vector2&         los,
                                      const Stokvec&         pol,
                                      const AscendingGrid&   range_bins) try {
  ARTS_TIME_REPORT
  ARTS_USER_ERROR_IF(freq_grid.empty(), "Radar frequency grid cannot be empty")
  ARTS_USER_ERROR_IF(range_bins.size() < 2, "Radar range_bins needs at least two edges")

  const Size old_n = measurement_sensor.size();
  ARTS_USER_ERROR_IF(radar_range_limits.nrows() != 0 and
                         (radar_range_limits.ncols() != 2 or static_cast<Size>(radar_range_limits.nrows()) != old_n),
                     "radar_range_limits must be empty or have shape [{}, 2], got {:B,}",
                     old_n,
                     radar_range_limits.shape())

  const Size nbins = range_bins.size() - 1;
  const Size nnew  = freq_grid.size() * nbins;
  Matrix     new_limits(old_n + nnew, 2);
  if (old_n) new_limits[Range(0, old_n), joker] = radar_range_limits;

  auto poslos = std::make_shared<const SensorPosLosVector>(SensorPosLosVector{SensorPosLos{.pos = pos, .los = los}});
  measurement_sensor.reserve(old_n + nnew);
  measurement_sensor_meta.reserve(measurement_sensor_meta.size() + freq_grid.size());

  Size row = old_n;
  for (Size iv = 0; iv < freq_grid.size(); ++iv) {
    auto                        one_frequency = std::make_shared<const AscendingGrid>(AscendingGrid{freq_grid[iv]});
    sensor::SparseStokvecMatrix weights(1, 1);
    weights[0, 0] = pol;

    for (Size ib = 0; ib < nbins; ++ib) {
      measurement_sensor.emplace_back(one_frequency, poslos, weights);
      new_limits[row, 0] = range_bins[ib];
      new_limits[row, 1] = range_bins[ib + 1];
      ++row;
    }

    SortedGriddedField1 meta;
    meta.data_name     = std::format("radar-{}-Hz", freq_grid[iv]);
    meta.gridname<0>() = "range";
    Vector centers(nbins);
    for (Size ib = 0; ib < nbins; ++ib) centers[ib] = std::midpoint(range_bins[ib], range_bins[ib + 1]);
    meta.grid<0>() = AscendingGrid{std::move(centers)};
    meta.data.resize(nbins);
    meta.data = 0.0;
    measurement_sensor_meta.push_back(SensorMetaInfo{.data = std::move(meta)});
  }
  radar_range_limits = std::move(new_limits);
}
ARTS_METHOD_ERROR_CATCH

void measurement_vecFromRadarSingleScattering(const Workspace&                ws,
                                              Vector&                         measurement_vec,
                                              Matrix&                         measurement_jac,
                                              Matrix&                         radar_aux,
                                              const ArrayOfSensorObsel&       measurement_sensor,
                                              const Matrix&                   radar_range_limits,
                                              const JacobianTargets&          jac_targets,
                                              const AtmField&                 atm_field,
                                              const SurfaceField&             surf_field,
                                              const ArrayOfScatteringSpecies& scattering_species,
                                              const Agenda&                   ray_path_observer_agenda,
                                              const Agenda&                   spectral_propmat_agenda,
                                              const Stokvec&                  transmitted_stokes,
                                              const String&                   range_mode_string,
                                              const String&                   unit,
                                              const Numeric&                  ze_tref,
                                              const Numeric&                  k2,
                                              const Numeric&                  dbze_min,
                                              const Numeric&                  pext_scaling,
                                              const ArrayOfString&            aux_vars) try {
  ARTS_TIME_REPORT
  const Size ny = measurement_sensor.size();
  ARTS_USER_ERROR_IF(ny == 0, "measurement_sensor cannot be empty for a radar calculation")
  ARTS_USER_ERROR_IF(radar_range_limits.nrows() != static_cast<Index>(ny) or radar_range_limits.ncols() != 2,
                     "radar_range_limits must have shape [{}, 2], got {:B,}",
                     ny,
                     radar_range_limits.shape())
  for (Size iy = 0; iy < ny; ++iy) {
    measurement_sensor[iy].check();
    ARTS_USER_ERROR_IF(not(radar_range_limits[iy, 0] < radar_range_limits[iy, 1]),
                       "Radar range limits in row {} are not strictly increasing: [{}, {}]",
                       iy,
                       radar_range_limits[iy, 0],
                       radar_range_limits[iy, 1])
  }
  ARTS_USER_ERROR_IF(transmitted_stokes.I() != 1.0,
                     "The transmitted Stokes I component must equal one, got {}",
                     transmitted_stokes.I())
  ARTS_USER_ERROR_IF(unit != "1" and unit != "Ze" and unit != "dBZe",
                     "Radar unit must be \"1\", \"Ze\", or \"dBZe\"; got \"{}\"",
                     unit)
  ARTS_USER_ERROR_IF(
      pext_scaling < 0.0 or pext_scaling > 2.0, "pext_scaling must be between 0 and 2; got {}", pext_scaling)
  ARTS_USER_ERROR_IF(scattering_species.species.empty(), "Radar single scattering requires scattering species")
  ARTS_USER_ERROR_IF(not jac_targets.surf.empty() or not jac_targets.subsurf.empty() or not jac_targets.line.empty() or
                         not jac_targets.sensor.empty() or not jac_targets.error.empty(),
                     "The radar forward model currently supports atmospheric Jacobian targets only")
  for (const auto& name : aux_vars) {
    ARTS_USER_ERROR_IF(name != "Radiative background" and name != "Backscattering" and
                           name != "Abs species extinction" and name != "Particle extinction",
                       "Unknown radar auxiliary variable \"{}\"",
                       name)
  }

  const RadarRangeMode            coordinate_mode = range_mode(range_mode_string, radar_range_limits);
  const SensorSimulations         collected       = collect_simulations(measurement_sensor);
  std::vector<PreparedSimulation> simulations;
  simulations.reserve(collected.size());
  for (const auto& simulation : collected) {
    const auto& [frequency, poslos, iposlos] = simulation;
    PreparedSimulation prepared{.frequency = &frequency, .poslos = &poslos, .iposlos = iposlos};
    ray_path_observer_agendaExecute(
        ws, prepared.ray_path, poslos[iposlos].pos, poslos[iposlos].los, ray_path_observer_agenda);
    simulations.push_back(std::move(prepared));
  }

  const auto calculate = [&](const AtmField& field) {
    return radar_forward(ws,
                         measurement_sensor,
                         radar_range_limits,
                         simulations,
                         field,
                         surf_field,
                         scattering_species,
                         spectral_propmat_agenda,
                         transmitted_stokes,
                         coordinate_mode,
                         unit,
                         ze_tref,
                         k2,
                         dbze_min,
                         pext_scaling,
                         aux_vars,
                         jac_targets);
  };

  RadarResult base = calculate(atm_field);
  measurement_vec  = std::move(base.measurement);
  measurement_jac  = std::move(base.jacobian);
  radar_aux        = std::move(base.auxiliary);

  const Size nx = jac_targets.x_size();
  if (nx == 0) return;

  Vector state(nx, 0.0);
  for (const auto& target : jac_targets.atm) target.update_state(state, atm_field);
  for (const auto& target : jac_targets.atm) target.update_jac(measurement_jac, state, atm_field);
}
ARTS_METHOD_ERROR_CATCH

void model_state_vecFromRadarOnionPeeling(const Workspace&                ws,
                                          Vector&                         model_state_vec,
                                          Vector&                         measurement_vec_fit,
                                          Matrix&                         measurement_jac,
                                          AtmField&                       atm_field,
                                          const Vector&                   measurement_vec,
                                          const ArrayOfSensorObsel&       measurement_sensor,
                                          const Matrix&                   radar_range_limits,
                                          const JacobianTargets&          jac_targets,
                                          const SurfaceField&             surf_field,
                                          const ArrayOfScatteringSpecies& scattering_species,
                                          const Agenda&                   ray_path_observer_agenda,
                                          const Agenda&                   spectral_propmat_agenda,
                                          const Stokvec&                  transmitted_stokes,
                                          const String&                   range_mode_string,
                                          const String&                   unit,
                                          const Numeric&                  ze_tref,
                                          const Numeric&                  k2,
                                          const Numeric&                  dbze_min,
                                          const Numeric&                  pext_scaling,
                                          const Numeric&                  measurement_noise_floor,
                                          const Numeric&                  state_min,
                                          const Numeric&                  state_max,
                                          const Numeric&                  max_step,
                                          const Numeric&                  tolerance,
                                          const Index&                    max_iterations,
                                          const Index&                    max_sweeps) try {
  ARTS_TIME_REPORT
  const Size ny = measurement_sensor.size();
  const Size nx = jac_targets.x_size();
  ARTS_USER_ERROR_IF(radar_range_limits.nrows() != static_cast<Index>(ny) or radar_range_limits.ncols() != 2,
                     "radar_range_limits must have shape [{}, 2], got {:B,}",
                     ny,
                     radar_range_limits.shape())
  ARTS_USER_ERROR_IF(nx == 0 or jac_targets.atm.empty(),
                     "Radar onion peeling requires at least one finalized atmospheric Jacobian target")
  ARTS_USER_ERROR_IF(not jac_targets.surf.empty() or not jac_targets.subsurf.empty() or not jac_targets.line.empty() or
                         not jac_targets.sensor.empty() or not jac_targets.error.empty(),
                     "Radar onion peeling currently changes atmospheric state coordinates only")
  ARTS_USER_ERROR_IF(state_min > state_max, "state_min ({}) exceeds state_max ({})", state_min, state_max)
  ARTS_USER_ERROR_IF(max_step <= 0.0, "max_step must be positive, got {}", max_step)
  ARTS_USER_ERROR_IF(tolerance <= 0.0, "tolerance must be positive, got {}", tolerance)
  ARTS_USER_ERROR_IF(max_iterations < 1, "max_iterations must be positive, got {}", max_iterations)
  ARTS_USER_ERROR_IF(max_sweeps < 1, "max_sweeps must be positive, got {}", max_sweeps)

  // The state is deliberately obtained and applied through the normal target
  // transformations.  The peeled profile is therefore an ordinary ARTS model
  // state and can be passed directly to OEM (including transformed targets).
  model_state_vec.resize(nx);
  model_state_vec = 0.0;
  for (const auto& target : jac_targets.atm) target.update_state(model_state_vec, atm_field);

  const auto calculate = [&]() {
    Matrix radar_aux;
    measurement_vecFromRadarSingleScattering(ws,
                                             measurement_vec_fit,
                                             measurement_jac,
                                             radar_aux,
                                             measurement_sensor,
                                             radar_range_limits,
                                             jac_targets,
                                             atm_field,
                                             surf_field,
                                             scattering_species,
                                             ray_path_observer_agenda,
                                             spectral_propmat_agenda,
                                             transmitted_stokes,
                                             range_mode_string,
                                             unit,
                                             ze_tref,
                                             k2,
                                             dbze_min,
                                             pext_scaling,
                                             {});
  };

  calculate();

  // Onion order is distance from the sensor.  Distance and round-trip time
  // already increase away from it.  For altitude gates, use the altitude of
  // the observation geometry so both downward- and upward-looking radars are
  // ordered correctly.  This ordering is local to every sensor geometry;
  // interleaving independent profiles is harmless because their Jacobian
  // blocks have zero cross-sensitivity.
  std::vector<Size> order(ny);
  stdr::iota(order, Size{0});
  const RadarRangeMode coordinate_mode = range_mode(range_mode_string, radar_range_limits);
  const auto           distance_key    = [&](const Size iy) {
    const Numeric center = std::midpoint(radar_range_limits[iy, 0], radar_range_limits[iy, 1]);
    if (coordinate_mode != RadarRangeMode::Altitude) return center;
    const auto& poslos = measurement_sensor[iy].poslos_grid();
    ARTS_USER_ERROR_IF(poslos.empty(), "Radar observation {} has no position/LOS geometry", iy)
    return std::abs(center - poslos.front().pos[0]);
  };
  stdr::stable_sort(order, [&](const Size a, const Size b) { return distance_key(a) < distance_key(b); });

  std::vector<bool>                  peeled(nx, false);
  std::vector<std::pair<Size, Size>> gate_state;
  for (const Size iy : order) {
    if (not std::isfinite(measurement_vec[iy]) or measurement_vec[iy] <= measurement_noise_floor) continue;

    // Associate this range gate with its strongest not-yet-retrieved local
    // atmospheric state sensitivity.  This uses the same mapped Jacobian as
    // OEM and avoids imposing a legacy pressure-grid or two-species layout.
    Size    ix_best = nx;
    Numeric best    = 0.0;
    for (Size ix = 0; ix < nx; ++ix) {
      const Numeric sensitivity = std::abs(measurement_jac[iy, ix]);
      if (not peeled[ix] and std::isfinite(sensitivity) and sensitivity > best) {
        best    = sensitivity;
        ix_best = ix;
      }
    }
    if (ix_best == nx) continue;

    bool converged = false;
    for (Index iteration = 0; iteration < max_iterations; ++iteration) {
      const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
      if (std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]))) {
        converged = true;
        break;
      }

      const Numeric derivative = measurement_jac[iy, ix_best];
      ARTS_USER_ERROR_IF(not std::isfinite(derivative) or derivative == 0.0,
                         "Radar onion peeling encountered a zero/non-finite derivative for gate {} and state {}",
                         iy,
                         ix_best)
      Numeric step             = residual / derivative;
      step                     = std::clamp(step, -max_step, max_step);
      model_state_vec[ix_best] = std::clamp(model_state_vec[ix_best] + step, state_min, state_max);
      for (const auto& target : jac_targets.atm) target.update_model(atm_field, model_state_vec);
      calculate();
    }
    if (not converged) {
      const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
      converged = std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]));
    }
    ARTS_USER_ERROR_IF(not converged,
                       "Radar onion peeling did not converge for gate {} (observed {}, fitted {}) after {} iterations",
                       iy,
                       measurement_vec[iy],
                       measurement_vec_fit[iy],
                       max_iterations)
    peeled[ix_best] = true;
    gate_state.emplace_back(iy, ix_best);
  }

  // GeodeticField3 interpolation means adjacent model-state coordinates can
  // weakly affect a gate already visited.  Repeated outward sweeps are the
  // continuous-grid counterpart of ARTS2's layer-by-layer peeling and recover
  // mutual consistency without replacing the analytical Jacobian by a table.
  auto fit_is_converged = [&]() {
    return std::ranges::all_of(gate_state, [&](const auto& pair) {
      const Size    iy       = pair.first;
      const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
      return std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]));
    });
  };
  bool all_converged = fit_is_converged();
  for (Index sweep = 1; sweep < max_sweeps and not all_converged; ++sweep) {
    for (const auto [iy, ix] : gate_state) {
      for (Index iteration = 0; iteration < max_iterations; ++iteration) {
        const Numeric residual = measurement_vec[iy] - measurement_vec_fit[iy];
        if (std::abs(residual) <= tolerance * std::max(Numeric{1e-30}, std::abs(measurement_vec[iy]))) break;
        const Numeric derivative = measurement_jac[iy, ix];
        ARTS_USER_ERROR_IF(not std::isfinite(derivative) or derivative == 0.0,
                           "Radar onion peeling encountered a zero/non-finite derivative for gate {} and state {}",
                           iy,
                           ix)
        Numeric step        = std::clamp(residual / derivative, -max_step, max_step);
        model_state_vec[ix] = std::clamp(model_state_vec[ix] + step, state_min, state_max);
        for (const auto& target : jac_targets.atm) target.update_model(atm_field, model_state_vec);
        calculate();
      }
    }
    all_converged = fit_is_converged();
  }
  ARTS_USER_ERROR_IF(
      not all_converged, "Radar onion peeling did not produce a mutually consistent fit after {} sweeps", max_sweeps)

  // Ensure the returned field, fit, and Jacobian all represent precisely the
  // returned state, ready for a subsequent OEM calculation.
  for (const auto& target : jac_targets.atm) target.update_model(atm_field, model_state_vec);
  calculate();
}
ARTS_METHOD_ERROR_CATCH
