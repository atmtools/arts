#include <atm.h>
#include <new_jacobian.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <ppath_struct.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

void dradEmpty(StokvecMatrix &drad,
               const Vector &f_grid,
               const JacobianTargets &jacobian_targets) {
  drad.resize(jacobian_targets.x_size(), f_grid.nelem());
  drad = Stokvec{0.0, 0.0, 0.0, 0.0};
}

void dradFromBackground(StokvecMatrix &drad,
                        const StokvecMatrix &background_drad,
                        const MuelmatVector &background_transmittance) {
  ARTS_USER_ERROR_IF(
      background_drad.ncols() != background_transmittance.nelem(),
      "background_drad must have same number of rows as the "
      "size of jacobian_targets")

  //! The radiance derivative shape is the background shape
  drad.resize(background_drad.shape());

  //! Set the background radiance derivative as that which is seen after "this" swath
  for (Index i = 0; i < drad.nrows(); i++) {
    std::transform(background_transmittance.begin(),
                   background_transmittance.end(),
                   background_drad[i].begin(),
                   drad[i].begin(),
                   std::multiplies<>());
  }
}

void dradAddPathPropagation(StokvecMatrix &drad,
                            const ArrayOfStokvecMatrix &ppvar_drad,
                            const JacobianTargets &jacobian_targets,
                            const AtmField &atm_field,
                            const Ppath &ppath) {
  const auto np = ppvar_drad.size();
  const auto nj = drad.nrows();
  const auto nf = drad.ncols();

  ARTS_USER_ERROR_IF(
      static_cast<Size>(drad.nrows()) != jacobian_targets.x_size(),
      "Bad size of drad, it's inner dimension should match the size of jacobian_targets")

  ARTS_USER_ERROR_IF(
      static_cast<Size>(ppath.pos.nrows()) != np,
      "ppath.pos must have same number of rows as the size of ppvar_drad")
  ARTS_USER_ERROR_IF(ppath.pos.ncols() != 3, "ppath.pos must have 3 columns")

  for (auto &dr : ppvar_drad) {
    ARTS_USER_ERROR_IF(
        dr.nrows() != nj,
        "ppvar_drad elements must have same number of rows as the size of "
        "jacobian_targets")
    ARTS_USER_ERROR_IF(
        dr.ncols() != nf,
        "ppvar_drad elements must have same number of columns as the size of "
        "cumulative_transmission")
  }

  //! Checks that the jacobian_targets can be used and throws if not
  jacobian_targets.throwing_check(nj);

  //! The altitude, latitude and longitude vectors must be copied because of how atm_field works
  const auto [alt, lat, lon] = [&]() {
    const auto n = ppath.pos.nrows();
    std::array<Vector, 3> g;
    for (auto &v : g) v.resize(n);
    for (Index i = 0; i < n; i++) {
      g[0][i] = ppath.pos[i][0];
      g[1][i] = ppath.pos[i][1];
      g[2][i] = ppath.pos[i][2];
    }
    return g;
  }();

  //! The derivative part from the atmosphere
  for (auto &atm_block : jacobian_targets.atm()) {
    ARTS_USER_ERROR_IF(not atm_field.contains(atm_block.type),
                       "No ",
                       atm_block.type,
                       " in atm_field but in jacobian_targets")
    const auto &data = atm_field[atm_block.type];
    const auto weights = data.flat_weights(alt, lat, lon);
    ARTS_ASSERT(weights.size() == np)

    for (Size j = 0; j < np; j++) {
      for (auto &w : weights[j]) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.x_start;
          ARTS_ASSERT(i < static_cast<Size>(nj))
          std::transform(
              ppvar_drad[j][atm_block.target_pos].begin(),
              ppvar_drad[j][atm_block.target_pos].end(),
              drad[i].begin(),
              drad[i].begin(),
              [x = w.second](auto &a, auto &b) { return x * a + b; });
        }
      }
    }
  }
}

void dradAddPathPropagation2(StokvecMatrix &drad,
                             const ArrayOfStokvecMatrix &ppvar_drad,
                             const JacobianTargets &jacobian_targets,
                             const AtmField &atm_field,
                             const ArrayOfPropagationPathPoint &rad_path) {
  const auto np = ppvar_drad.size();
  const auto nj = drad.nrows();
  const auto nf = drad.ncols();

  ARTS_USER_ERROR_IF(
      static_cast<Size>(drad.nrows()) != jacobian_targets.x_size(),
      "Bad size of drad, it's inner dimension should match the size of jacobian_targets")

  ARTS_USER_ERROR_IF(
      rad_path.size() != np,
      "ppath.pos must have same number of rows as the size of ppvar_drad")

  for (auto &dr : ppvar_drad) {
    ARTS_USER_ERROR_IF(
        dr.nrows() != nj,
        "ppvar_drad elements must have same number of rows as the size of "
        "jacobian_targets")
    ARTS_USER_ERROR_IF(
        dr.ncols() != nf,
        "ppvar_drad elements must have same number of columns as the size of "
        "cumulative_transmission")
  }

  //! Checks that the jacobian_targets can be used and throws if not
  jacobian_targets.throwing_check(nj);

  //! The altitude, latitude and longitude vectors must be copied because of how atm_field works
  const auto [alt, lat, lon] = [&]() {
    std::array<Vector, 3> g{Vector(np), Vector(np), Vector(np)};
    for (Size i = 0; i < np; i++) {
      g[0][i] = rad_path[i].pos[0];
      g[1][i] = rad_path[i].pos[1];
      g[2][i] = rad_path[i].pos[2];
    }
    return g;
  }();

  //! The derivative part from the atmosphere
  for (auto &atm_block : jacobian_targets.atm()) {
    ARTS_USER_ERROR_IF(not atm_field.contains(atm_block.type),
                       "No ",
                       atm_block.type,
                       " in atm_field but in jacobian_targets")
    const auto &data = atm_field[atm_block.type];
    const auto weights = data.flat_weights(alt, lat, lon);
    ARTS_ASSERT(weights.size() == np)

    for (Size j = 0; j < np; j++) {
      for (auto &w : weights[j]) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.x_start;
          ARTS_ASSERT(i < static_cast<Size>(nj))
          std::transform(
              ppvar_drad[j][atm_block.target_pos].begin(),
              ppvar_drad[j][atm_block.target_pos].end(),
              drad[i].begin(),
              drad[i].begin(),
              [x = w.second](auto &a, auto &b) { return x * a + b; });
        }
      }
    }
  }
}

void radFromPathPropagation(StokvecVector &rad,
                            const ArrayOfStokvecVector &ppvar_rad) {
  ARTS_USER_ERROR_IF(ppvar_rad.empty(), "Empty ppvar_rad")
  rad = ppvar_rad.front();
}

void radStandardEmission(const Workspace &ws,
                         StokvecVector &rad,
                         StokvecMatrix &drad,
                         const Vector &f_grid,
                         const JacobianTargets &jacobian_targets,
                         const AtmField &atm_field,
                         const Ppath &ppath,
                         const Agenda &space_radiation_agenda,
                         const Agenda &surface_radiation_agenda,
                         const Agenda &stop_distance_radiation_agenda,
                         const Agenda &propmat_clearsky_agenda,
                         const Numeric &rte_alonglos_v,
                         const Index &hse_derivative) try {
  background_radFromPath(ws,
                         rad,
                         drad,
                         f_grid,
                         jacobian_targets,
                         ppath,
                         space_radiation_agenda,
                         surface_radiation_agenda,
                         stop_distance_radiation_agenda);

  ArrayOfAtmPoint ppvar_atm;
  ppvar_atmFromPath(ppvar_atm, ppath, atm_field);

  ArrayOfVector ppvar_f;
  ppvar_fFromPath(ppvar_f, f_grid, ppath, ppvar_atm, rte_alonglos_v);

  ArrayOfPropmatVector ppvar_propmat;
  ArrayOfStokvecVector ppvar_nlte;
  ArrayOfPropmatMatrix ppvar_dpropmat;
  ArrayOfStokvecMatrix ppvar_dnlte;
  ppvar_propmatCalc(ws,
                    ppvar_propmat,
                    ppvar_nlte,
                    ppvar_dpropmat,
                    ppvar_dnlte,
                    propmat_clearsky_agenda,
                    jacobian_targets,
                    ppvar_f,
                    ppath,
                    ppvar_atm);

  ArrayOfStokvecVector ppvar_src;
  ArrayOfStokvecMatrix ppvar_dsrc;
  ppvar_srcFromPropmat(ppvar_src,
                       ppvar_dsrc,
                       ppvar_propmat,
                       ppvar_nlte,
                       ppvar_dpropmat,
                       ppvar_dnlte,
                       ppvar_f,
                       ppvar_atm,
                       jacobian_targets);

  ArrayOfMuelmatVector ppvar_tramat;
  ArrayOfArrayOfMuelmatMatrix ppvar_dtramat;
  Vector ppvar_distance;
  ArrayOfArrayOfVector ppvar_ddistance;
  ppvar_tramatCalc(ppvar_tramat,
                   ppvar_dtramat,
                   ppvar_distance,
                   ppvar_ddistance,
                   ppvar_propmat,
                   ppvar_dpropmat,
                   ppath,
                   ppvar_atm,
                   jacobian_targets,
                   hse_derivative);

  ArrayOfMuelmatVector ppvar_cumtramat;
  ppvar_cumtramatForward(ppvar_cumtramat, ppvar_tramat);

  ArrayOfStokvecVector ppvar_rad;
  ArrayOfStokvecMatrix ppvar_drad;
  ppvar_radCalcEmission(ppvar_rad,
                        ppvar_drad,
                        rad,
                        ppvar_src,
                        ppvar_dsrc,
                        ppvar_tramat,
                        ppvar_cumtramat,
                        ppvar_dtramat);

  MuelmatVector background_transmittance;
  background_transmittanceFromPathPropagationBack(background_transmittance,
                                                  ppvar_cumtramat);

  radFromPathPropagation(rad, ppvar_rad);

  dradFromBackground(drad, drad, background_transmittance);

  dradAddPathPropagation(drad, ppvar_drad, jacobian_targets, atm_field, ppath);
}
ARTS_METHOD_ERROR_CATCH

void radStandardEmission2(const Workspace &ws,
                          StokvecVector &rad,
                          StokvecMatrix &drad,
                          const Vector &f_grid,
                          const JacobianTargets &jacobian_targets,
                          const AtmField &atm_field,
                          const SurfaceField &surface_field,
                          const ArrayOfPropagationPathPoint &rad_path,
                          const Agenda &space_radiation_agenda,
                          const Agenda &surface_radiation_agenda,
                          const Agenda &propmat_clearsky_agenda,
                          const Numeric &rte_alonglos_v,
                          const Index &hse_derivative) try {
  background_radFromPath2(ws,
                          rad,
                          drad,
                          f_grid,
                          jacobian_targets,
                          rad_path,
                          space_radiation_agenda,
                          surface_radiation_agenda);

  ArrayOfAtmPoint ppvar_atm;
  ppvar_atmFromPath2(ppvar_atm, rad_path, atm_field);

  ArrayOfVector ppvar_f;
  ppvar_fFromPath2(ppvar_f, f_grid, rad_path, ppvar_atm, rte_alonglos_v);

  ArrayOfPropmatVector ppvar_propmat;
  ArrayOfStokvecVector ppvar_nlte;
  ArrayOfPropmatMatrix ppvar_dpropmat;
  ArrayOfStokvecMatrix ppvar_dnlte;
  ppvar_propmatCalc2(ws,
                     ppvar_propmat,
                     ppvar_nlte,
                     ppvar_dpropmat,
                     ppvar_dnlte,
                     propmat_clearsky_agenda,
                     jacobian_targets,
                     ppvar_f,
                     rad_path,
                     ppvar_atm);

  ArrayOfStokvecVector ppvar_src;
  ArrayOfStokvecMatrix ppvar_dsrc;
  ppvar_srcFromPropmat(ppvar_src,
                       ppvar_dsrc,
                       ppvar_propmat,
                       ppvar_nlte,
                       ppvar_dpropmat,
                       ppvar_dnlte,
                       ppvar_f,
                       ppvar_atm,
                       jacobian_targets);

  ArrayOfMuelmatVector ppvar_tramat;
  ArrayOfArrayOfMuelmatMatrix ppvar_dtramat;
  Vector ppvar_distance;
  ArrayOfArrayOfVector ppvar_ddistance;
  ppvar_tramatCalc2(ppvar_tramat,
                    ppvar_dtramat,
                    ppvar_distance,
                    ppvar_ddistance,
                    ppvar_propmat,
                    ppvar_dpropmat,
                    rad_path,
                    ppvar_atm,
                    surface_field,
                    jacobian_targets,
                    hse_derivative);

  ArrayOfMuelmatVector ppvar_cumtramat;
  ppvar_cumtramatForward(ppvar_cumtramat, ppvar_tramat);

  ArrayOfStokvecVector ppvar_rad;
  ArrayOfStokvecMatrix ppvar_drad;
  ppvar_radCalcEmission(ppvar_rad,
                        ppvar_drad,
                        rad,
                        ppvar_src,
                        ppvar_dsrc,
                        ppvar_tramat,
                        ppvar_cumtramat,
                        ppvar_dtramat);

  MuelmatVector background_transmittance;
  background_transmittanceFromPathPropagationBack(background_transmittance,
                                                  ppvar_cumtramat);

  radFromPathPropagation(rad, ppvar_rad);

  dradFromBackground(drad, drad, background_transmittance);

  dradAddPathPropagation2(
      drad, ppvar_drad, jacobian_targets, atm_field, rad_path);
}
ARTS_METHOD_ERROR_CATCH
