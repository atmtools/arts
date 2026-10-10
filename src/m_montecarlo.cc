#include <arts_constants.h>
#include <atm_path.h>
#include <lagrange_interp.h>
#include <mc_antenna.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <scattering_species.h>
#include <workspace.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>

namespace {

struct PolarizedOptics {
  Propmat extinction{};
  Stokvec absorption{};
  Stokvec emission{};
};

PolarizedOptics polarized_optics(const Workspace&                ws,
                                 Numeric                         frequency,
                                 const PropagationPathPoint&     point,
                                 const AtmPoint&                 atm,
                                 const ArrayOfScatteringSpecies& scattering_species,
                                 const Agenda&                   spectral_propmat_agenda,
                                 const JacobianTargets&          jac_targets,
                                 Index                           local_jac_target,
                                 Numeric                         local_jac_step) {
  const AscendingGrid freq_grid{frequency};
  PropmatVector       gas;
  StokvecVector       src;
  PropmatMatrix       gas_jac;
  StokvecMatrix       src_jac;
  spectral_propmat_agendaExecute(ws,
                                 gas,
                                 src,
                                 gas_jac,
                                 src_jac,
                                 freq_grid,
                                 Vector3{0.0, 0.0, 0.0},
                                 jac_targets,
                                 {},
                                 point,
                                 atm,
                                 spectral_propmat_agenda);

  if (local_jac_target >= 0) {
    ARTS_USER_ERROR_IF(static_cast<Size>(local_jac_target) >= jac_targets.target_count(),
                       "Local Monte Carlo Jacobian target {} is out of bounds",
                       local_jac_target)
    gas[0] += local_jac_step * gas_jac[local_jac_target, 0];
    src[0] += local_jac_step * src_jac[local_jac_target, 0];
  }

  Propmat particle_ext{};
  Stokvec particle_abs{};
  if (not scattering_species.species.empty()) {
    // Path points already store the direction of propagation.  In particular,
    // observer agendas mirror their external LOS when constructing the path.
    const Numeric za = point.los[0];
    auto za_grid     = std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{za}));
    const auto bulk  = scattering_species.get_bulk_scattering_properties_aro_gridded(
        atm, freq_grid, Vector{za}, Vector{0.0}, std::move(za_grid));
    particle_ext.A() = bulk.extinction_matrix[0, 0, 0, 0];
    particle_ext.B() = bulk.extinction_matrix[0, 0, 0, 1];
    particle_ext.W() = bulk.extinction_matrix[0, 0, 0, 2];
    particle_abs.I() = bulk.absorption_vector[0, 0, 0, 0];
    particle_abs.Q() = bulk.absorption_vector[0, 0, 0, 1];
  }

  const Stokvec absorption = Stokvec{gas[0].A(), gas[0].B(), gas[0].C(), gas[0].D()} + particle_abs;

  return {.extinction = gas[0] + particle_ext,
          .absorption = absorption,
          .emission   = rtepack::level_emission(absorption, src[0], frequency, atm.temperature)};
}

//! Linear interpolation of a cyclic coordinate along the shorter arc from ya to yb, cycled to the cycler's range
template <lagrange_interp::cyclic cycler>
Numeric cyclic_interp(Numeric ya, Numeric yb, const lagrange_interp::lag_t<1>& lag) {
  ya                    = cycler::cycle(ya);
  const Numeric shorter = cycler::cycle(cycler::cycle(yb) - ya + cycler::midpoint()) - cycler::midpoint();
  return cycler::cycle(lagrange_interp::interp(Vector2{ya, ya + shorter}, lag));
}

//! The path point at the fraction x in [0, 1] of the segment from a to b
PropagationPathPoint interpolate(const PropagationPathPoint& a, const PropagationPathPoint& b, Numeric x) {
  static constexpr std::array<Numeric, 2> segment{0.0, 1.0};
  const lagrange_interp::lag_t<1>         lag(segment, x, lagrange_interp::ascending_grid_t{});
  const auto linear = [&lag](Numeric ya, Numeric yb) { return lagrange_interp::interp(Vector2{ya, yb}, lag); };

  PropagationPathPoint p = a;
  p.pos_type             = PathPositionType::atm;
  p.los_type             = PathPositionType::atm;
  p.pos[0]               = linear(a.pos[0], b.pos[0]);
  p.pos[1]               = linear(a.pos[1], b.pos[1]);
  p.pos[2]               = cyclic_interp<lagrange_interp::loncross>(a.pos[2], b.pos[2], lag);
  // The interpolated zenith angle lies between the end points; the clamp removes the round-off of the weights,
  // which can otherwise leave [0, 180] deg
  p.los[0] = std::clamp(linear(a.los[0], b.los[0]), std::min(a.los[0], b.los[0]), std::max(a.los[0], b.los[0]));
  p.los[1] = cyclic_interp<lagrange_interp::azicross>(a.los[1], b.los[1], lag);
  p.nreal  = linear(a.nreal, b.nreal);
  p.ngroup = linear(a.ngroup, b.ngroup);
  return p;
}

Stokvec boundary_radiance(const Workspace&            ws,
                          Numeric                     frequency,
                          const PropagationPathPoint& point,
                          const SurfaceField&         surf_field,
                          const SubsurfaceField&      subsurf_field,
                          const Agenda&               spectral_rad_space_agenda,
                          const Agenda&               spectral_rad_surface_agenda) {
  StokvecVector       value;
  StokvecMatrix       jac;
  const AscendingGrid grid{frequency};
  if (point.los_type == PathPositionType::space) {
    spectral_rad_space_agendaExecute(ws, value, jac, grid, JacobianTargets{}, point, spectral_rad_space_agenda);
  } else if (point.los_type == PathPositionType::surface) {
    spectral_rad_surface_agendaExecute(
        ws, value, jac, grid, JacobianTargets{}, point, surf_field, subsurf_field, spectral_rad_surface_agenda);
  } else {
    ARTS_USER_ERROR("Monte Carlo path ended without reaching space or the surface")
  }
  return value[0];
}

Vector2 isotropic_los(RandomNumberGenerator<>& rng) {
  auto uniform = rng.get<std::uniform_real_distribution>(0.0, 1.0);
  return {Conversion::rad2deg(std::acos(1.0 - 2.0 * uniform())), 360.0 * uniform() - 180.0};
}

struct Collision {
  bool                 happened{false};
  PropagationPathPoint point{};
  AtmPoint             atm{};
  PolarizedOptics      optics{};
  PropagationPathPoint end{};
  Muelmat              conditional_transport{Muelmat::id()};
};

Collision next_collision(const Workspace&                ws,
                         RandomNumberGenerator<>&        rng,
                         Numeric                         frequency,
                         const Vector3&                  pos,
                         const Vector2&                  los,
                         const AtmField&                 atm_field,
                         const SurfaceField&             surf_field,
                         const ArrayOfScatteringSpecies& scattering_species,
                         const Agenda&                   ray_path_observer_agenda,
                         const Agenda&                   spectral_propmat_agenda,
                         const JacobianTargets&          jac_targets,
                         Index                           local_jac_target,
                         Numeric                         local_jac_step) {
  ArrayOfPropagationPathPoint path;
  ray_path_observer_agendaExecute(ws, path, pos, los, ray_path_observer_agenda);
  ARTS_USER_ERROR_IF(path.empty(), "The ray-path agenda returned an empty path")
  if (path.size() == 1) return {.end = path.back()};

  const auto    atm         = forward_atm_path(path, atm_field);
  const auto    distances   = path::distance(path, surf_field.ellipsoid);
  auto          uniform     = rng.get<std::uniform_real_distribution>(0.0, 1.0);
  const Numeric target      = -std::log(std::max(uniform(), std::numeric_limits<Numeric>::min()));
  Numeric       tau         = 0.0;
  Muelmat       conditional = Muelmat::id();
  auto          previous    = polarized_optics(ws,
                                               frequency,
                                               path[0],
                                               atm[0],
                                               scattering_species,
                                               spectral_propmat_agenda,
                                               jac_targets,
                                               local_jac_target,
                                               local_jac_step);

  for (Size i = 1; i < path.size(); ++i) {
    const auto    current = polarized_optics(ws,
                                             frequency,
                                             path[i],
                                             atm[i],
                                             scattering_species,
                                             spectral_propmat_agenda,
                                             jac_targets,
                                             local_jac_target,
                                             local_jac_step);
    const Numeric dtau    = 0.5 * (previous.extinction.A() + current.extinction.A()) * distances[i - 1];
    if (dtau > 0.0 and tau + dtau >= target) {
      const Numeric x           = std::clamp((target - tau) / dtau, 0.0, 1.0);
      auto          point       = interpolate(path[i - 1], path[i], x);
      auto          atm_point   = atm_field.at(point.pos);
      auto          optics      = polarized_optics(ws,
                                                   frequency,
                                                   point,
                                                   atm_point,
                                                   scattering_species,
                                                   spectral_propmat_agenda,
                                                   jac_targets,
                                                   local_jac_target,
                                                   local_jac_step);
      const Numeric ds          = x * distances[i - 1];
      const Numeric partial_tau = 0.5 * (previous.extinction.A() + optics.extinction.A()) * ds;
      const Muelmat partial     = rtepack::tran(previous.extinction, optics.extinction, ds)();
      conditional               = conditional * (std::exp(partial_tau) * partial);
      return {.happened              = true,
              .point                 = point,
              .atm                   = std::move(atm_point),
              .optics                = optics,
              .end                   = path.back(),
              .conditional_transport = conditional};
    }
    if (distances[i - 1] > 0.0) {
      const Muelmat segment = rtepack::tran(previous.extinction, current.extinction, distances[i - 1])();
      conditional           = conditional * (std::exp(std::max(0.0, dtau)) * segment);
    }
    tau      += std::max(0.0, dtau);
    previous  = current;
  }
  return {.end = path.back(), .conditional_transport = conditional};
}

Muelmat phase_weight(const ArrayOfScatteringSpecies& species,
                     const AtmPoint&                 atm,
                     Numeric                         frequency,
                     const Vector2&                  old_los,
                     const Vector2&                  new_los,
                     Numeric                         scattering_coefficient) {
  // old_los and new_los use the external observer-LOS convention.  Convert
  // both to propagation directions before evaluating the phase matrix.
  const Vector2 outgoing = path::mirror(old_los);
  const Vector2 incoming = path::mirror(new_los);
  Numeric       delta_aa = outgoing[1] - incoming[1];
  while (delta_aa < -180.0) delta_aa += 360.0;
  while (delta_aa > 180.0) delta_aa -= 360.0;
  auto za = std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{outgoing[0]}));
  const auto bulk = species.get_bulk_scattering_properties_aro_gridded(
      atm, AscendingGrid{frequency}, Vector{incoming[0]}, Vector{delta_aa}, za);
  ARTS_USER_ERROR_IF(not bulk.phase_matrix, "MCGeneral failed to transform the phase matrix to the laboratory frame")
  const auto flat = (*bulk.phase_matrix)[0, 0, 0, 0, 0, joker];
  Muelmat    out{0.0};
  for (Index i = 0; i < 4; ++i)
    for (Index j = 0; j < 4; ++j) out[i, j] = flat[4 * i + j];
  return (4.0 * Constant::pi / scattering_coefficient) * out;
}

}  // namespace

void MCGeneral(const Workspace&                ws,
               Stokvec&                        mc_spectral_rad,
               Stokvec&                        mc_error,
               Index&                          mc_iteration_count,
               const AtmField&                 atm_field,
               const SurfaceField&             surf_field,
               const SubsurfaceField&          subsurf_field,
               const ArrayOfScatteringSpecies& scattering_species,
               const MCAntenna&                mc_antenna,
               const Agenda&                   ray_path_observer_agenda,
               const Agenda&                   spectral_propmat_agenda,
               const Agenda&                   spectral_rad_space_agenda,
               const Agenda&                   spectral_rad_surface_agenda,
               const Numeric&                  frequency,
               const Vector3&                  sensor_pos,
               const Vector2&                  sensor_los,
               const Index&                    mc_seed,
               const Index&                    mc_min_iter,
               const Index&                    mc_max_iter,
               const Index&                    mc_max_scatorder,
               const Numeric&                  mc_std_err,
               const Numeric&                  mc_max_time,
               const JacobianTargets&          local_jac_targets,
               const Index&                    local_jac_target,
               const Numeric&                  local_jac_step) try {
  ARTS_TIME_REPORT
  ARTS_USER_ERROR_IF(frequency <= 0.0, "frequency must be positive")
  ARTS_USER_ERROR_IF(mc_min_iter < 1 or mc_max_iter < mc_min_iter, "Require 1 <= mc_min_iter <= mc_max_iter")
  ARTS_USER_ERROR_IF(mc_max_scatorder < 0, "mc_max_scatorder cannot be negative")
  ARTS_USER_ERROR_IF(mc_std_err < 0.0 or mc_max_time < 0.0, "Stopping thresholds cannot be negative")

  RandomNumberGenerator<> rng(mc_seed);
  rng.force_seed(mc_seed);
  const Matrix33 ant_to_enu = rotmat_enu(sensor_los);
  Stokvec        sum{0.0}, sumsq{0.0};
  const auto     started = std::chrono::steady_clock::now();
  mc_iteration_count     = 0;

  for (; mc_iteration_count < mc_max_iter; ++mc_iteration_count) {
    const auto [initial_los, ray_rotation] = mc_antenna.draw_los(rng, ant_to_enu, sensor_los);
    Vector3 pos                            = sensor_pos;
    Vector2 los                            = initial_los;
    Muelmat weight                         = rotmat_stokes(-1.0, -1.0, ray_rotation, ant_to_enu);
    Stokvec sample{0.0};

    for (Index order = 0;;) {
      const auto collision = next_collision(ws,
                                            rng,
                                            frequency,
                                            pos,
                                            los,
                                            atm_field,
                                            surf_field,
                                            scattering_species,
                                            ray_path_observer_agenda,
                                            spectral_propmat_agenda,
                                            local_jac_targets,
                                            local_jac_target,
                                            local_jac_step);
      weight               = weight * collision.conditional_transport;
      if (not collision.happened) {
        sample = weight * boundary_radiance(ws,
                                            frequency,
                                            collision.end,
                                            surf_field,
                                            subsurf_field,
                                            spectral_rad_space_agenda,
                                            spectral_rad_surface_agenda);
        break;
      }

      auto          uniform    = rng.get<std::uniform_real_distribution>(0.0, 1.0);
      const Numeric extinction = collision.optics.extinction.A();
      const Numeric scattering = std::max(0.0, extinction - collision.optics.absorption.I());
      if (scattering <= 0.0 or uniform() >= scattering / extinction) {
        const Numeric absorption = collision.optics.absorption.I();
        ARTS_USER_ERROR_IF(absorption <= 0.0, "Invalid non-positive absorption collision probability")
        sample = weight * (collision.optics.emission / absorption);
        break;
      }
      if (++order > mc_max_scatorder) break;

      const Vector2 new_los = isotropic_los(rng);
      weight = weight * phase_weight(scattering_species, collision.atm, frequency, los, new_los, scattering);
      pos    = collision.point.pos;
      los    = new_los;
    }

    sum += sample;
    for (Index s = 0; s < 4; ++s) sumsq[s] += sample[s] * sample[s];
    const Index n = mc_iteration_count + 1;
    if (n < mc_min_iter) continue;
    const Numeric nn = static_cast<Numeric>(n);

    Numeric largest_error = 0.0;
    for (Index s = 0; s < 4; ++s) {
      const Numeric mean = sum[s] / nn;
      largest_error      = std::max(largest_error, std::sqrt(std::max(0.0, sumsq[s] / nn - mean * mean) / nn));
    }
    const Numeric elapsed = std::chrono::duration<Numeric>(std::chrono::steady_clock::now() - started).count();
    if ((mc_std_err > 0.0 and largest_error <= mc_std_err) or (mc_max_time > 0.0 and elapsed >= mc_max_time)) {
      ++mc_iteration_count;
      break;
    }
  }

  ARTS_USER_ERROR_IF(mc_iteration_count == 0, "MCGeneral produced no samples")
  const Numeric n = static_cast<Numeric>(mc_iteration_count);
  for (Index s = 0; s < 4; ++s) {
    mc_spectral_rad[s]     = sum[s] / n;
    const Numeric variance = std::max(0.0, sumsq[s] / n - mc_spectral_rad[s] * mc_spectral_rad[s]);
    mc_error[s]            = std::sqrt(variance / n);
  }
}
ARTS_METHOD_ERROR_CATCH

void MCGeneral(const Workspace&                ws,
               Stokvec&                        mc_spectral_rad,
               Stokvec&                        mc_error,
               Index&                          mc_iteration_count,
               const AtmField&                 atm_field,
               const SurfaceField&             surf_field,
               const SubsurfaceField&          subsurf_field,
               const ArrayOfScatteringSpecies& scattering_species,
               const MCAntenna&                mc_antenna,
               const Agenda&                   ray_path_observer_agenda,
               const Agenda&                   spectral_propmat_agenda,
               const Agenda&                   spectral_rad_space_agenda,
               const Agenda&                   spectral_rad_surface_agenda,
               const Numeric&                  frequency,
               const Vector3&                  sensor_pos,
               const Vector2&                  sensor_los,
               const Index&                    mc_seed,
               const Index&                    mc_min_iter,
               const Index&                    mc_max_iter,
               const Index&                    mc_max_scatorder,
               const Numeric&                  mc_std_err,
               const Numeric&                  mc_max_time) {
  MCGeneral(ws,
            mc_spectral_rad,
            mc_error,
            mc_iteration_count,
            atm_field,
            surf_field,
            subsurf_field,
            scattering_species,
            mc_antenna,
            ray_path_observer_agenda,
            spectral_propmat_agenda,
            spectral_rad_space_agenda,
            spectral_rad_surface_agenda,
            frequency,
            sensor_pos,
            sensor_los,
            mc_seed,
            mc_min_iter,
            mc_max_iter,
            mc_max_scatorder,
            mc_std_err,
            mc_max_time,
            JacobianTargets{},
            -1,
            0.0);
}

void spectral_radMonteCarlo(const Workspace&                ws,
                            StokvecVector&                  spectral_rad,
                            StokvecMatrix&                  spectral_rad_jac,
                            ArrayOfPropagationPathPoint&    ray_path,
                            const AscendingGrid&            freq_grid,
                            const JacobianTargets&          jac_targets,
                            const Vector3&                  obs_pos,
                            const Vector2&                  obs_los,
                            const AtmField&                 atm_field,
                            const SurfaceField&             surf_field,
                            const SubsurfaceField&          subsurf_field,
                            const ArrayOfScatteringSpecies& scat_species,
                            const Agenda&                   ray_path_observer_agenda,
                            const Agenda&                   spectral_propmat_agenda,
                            const Agenda&                   spectral_rad_space_agenda,
                            const Agenda&                   spectral_rad_surface_agenda,
                            const Index&                    mc_seed,
                            const Index&                    mc_min_iter,
                            const Index&                    mc_max_iter,
                            const Index&                    mc_max_scatorder,
                            const Numeric&                  mc_std_err,
                            const Numeric&                  mc_max_time) try {
  ARTS_TIME_REPORT

  ray_path_observer_agendaExecute(ws, ray_path, obs_pos, obs_los, ray_path_observer_agenda);
  ARTS_USER_ERROR_IF(ray_path.empty(), "The ray-path agenda returned an empty path")

  spectral_rad.resize(freq_grid.size());
  spectral_rad_jac.resize(jac_targets.x_size(), freq_grid.size());
  spectral_rad_jac = 0.0;

  MCAntenna pencil;
  pencil.set_pencil_beam();

  const auto run = [&](Stokvec&               out,
                       const Numeric          frequency,
                       const Index            seed,
                       const AtmField&        atm,
                       const SurfaceField&    surf,
                       const SubsurfaceField& subsurf,
                       const Index            local_target = -1,
                       const Numeric          local_step   = 0.0) {
    Stokvec error;
    Index   iterations;
    MCGeneral(ws,
              out,
              error,
              iterations,
              atm,
              surf,
              subsurf,
              scat_species,
              pencil,
              ray_path_observer_agenda,
              spectral_propmat_agenda,
              spectral_rad_space_agenda,
              spectral_rad_surface_agenda,
              frequency,
              obs_pos,
              obs_los,
              seed,
              mc_min_iter,
              mc_max_iter,
              mc_max_scatorder,
              mc_std_err,
              mc_max_time,
              jac_targets,
              local_target,
              local_step);
  };

  for (Size i = 0; i < freq_grid.size(); ++i) {
    run(spectral_rad[i], freq_grid[i], mc_seed + static_cast<Index>(i), atm_field, surf_field, subsurf_field);
  }

  const auto finite_difference = [&]<typename Target, typename Field>(const Target& target, const Field& field) {
    ARTS_USER_ERROR_IF(target.d == 0.0, "Monte Carlo Jacobian perturbation must be nonzero for target {}", target.type)

    Vector state(jac_targets.x_size(), 0.0);
    target.update_state(state, field);

    for (Size j = 0; j < target.x_size; ++j) {
      Field  perturbed = field;
      Vector perturbed_state{state};
      perturbed_state[target.x_start + j] += target.d;
      target.update_model(perturbed, perturbed_state);

      for (Size i = 0; i < freq_grid.size(); ++i) {
        Stokvec value;
        if constexpr (std::same_as<Field, AtmField>) {
          run(value, freq_grid[i], mc_seed + static_cast<Index>(i), perturbed, surf_field, subsurf_field);
        } else if constexpr (std::same_as<Field, SurfaceField>) {
          run(value, freq_grid[i], mc_seed + static_cast<Index>(i), atm_field, perturbed, subsurf_field);
        } else {
          run(value, freq_grid[i], mc_seed + static_cast<Index>(i), atm_field, surf_field, perturbed);
        }
        spectral_rad_jac[target.x_start + j, i] = (value - spectral_rad[i]) / target.d;
      }
    }
  };

  for (const auto& target : jac_targets.atm) finite_difference(target, atm_field);
  for (const auto& target : jac_targets.surf) finite_difference(target, surf_field);
  for (const auto& target : jac_targets.subsurf) finite_difference(target, subsurf_field);

  for (const auto& target : jac_targets.line) {
    ARTS_USER_ERROR_IF(target.x_size != 1,
                       "Monte Carlo spectroscopic target {} has {} state elements; expected one",
                       target.type,
                       target.x_size)
    ARTS_USER_ERROR_IF(target.d == 0.0, "Monte Carlo Jacobian perturbation must be nonzero for target {}", target.type)

    for (Size i = 0; i < freq_grid.size(); ++i) {
      Stokvec value;
      run(value,
          freq_grid[i],
          mc_seed + static_cast<Index>(i),
          atm_field,
          surf_field,
          subsurf_field,
          static_cast<Index>(target.target_pos),
          target.d);
      spectral_rad_jac[target.x_start, i] = (value - spectral_rad[i]) / target.d;
    }
  }
}
ARTS_METHOD_ERROR_CATCH
