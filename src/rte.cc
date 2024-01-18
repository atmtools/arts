/**
  @file   rte.cc
  @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  @date   2002-05-29

  \brief  Functions to solve radiative transfer tasks.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "rte.h"

#include <workspace.h>

#include <algorithm>
#include <cmath>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "arts_conversions.h"
#include "atm.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "jacobian.h"
#include "math_funcs.h"
#include "matpack_concepts.h"
#include "new_jacobian.h"
#include "path_point.h"
#include "physics_funcs.h"
#include "refraction.h"
#include "rtepack.h"
#include "special_interp.h"
#include "species_tags.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric TEMP_0_C=Constant::temperature_at_0c;

/*===========================================================================
  === The functions in alphabetical order
  ===========================================================================*/

void adapt_stepwise_partial_derivatives(
    PropmatMatrix& dK_dx,
    StokvecMatrix& dS_dx,
    const JacobianTargets& jacobian_targets,
    const ConstVectorView& ppath_f_grid,
    const ConstVectorView& ppath_line_of_sight) {
  // All relevant quantities are extracted first
  DEBUG_ONLY(const Size nq = jacobian_targets.target_count();)
  DEBUG_ONLY(const Index nv = ppath_f_grid.size();)
  ARTS_ASSERT(nq == static_cast<Size>(dK_dx.nrows()))
  ARTS_ASSERT(nq == static_cast<Size>(dS_dx.nrows()))
  ARTS_ASSERT(nv == dK_dx.ncols())
  ARTS_ASSERT(nv == dS_dx.ncols())

  for (auto w : {Atm::Key::wind_u, Atm::Key::wind_v, Atm::Key::wind_w}) {
    if (auto wind_pair = jacobian_targets.find<Jacobian::AtmTarget>(w);
        wind_pair.first) {
      const auto i = wind_pair.second->target_pos;
      ARTS_ASSERT(i < nq)

      const Vector scale =
          get_stepwise_f_partials(ppath_line_of_sight, ppath_f_grid, w);

      std::transform(scale.begin(),
                     scale.end(),
                     dK_dx[i].begin(),
                     dK_dx[i].begin(),
                     std::multiplies<>());
      std::transform(scale.begin(),
                     scale.end(),
                     dS_dx[i].begin(),
                     dS_dx[i].begin(),
                     std::multiplies<>());
    }
  }
}

void apply_iy_unit(MatrixView iy,
                   const String& iy_unit,
                   const ConstVectorView& f_grid,
                   const Numeric& n,
                   const ArrayOfIndex& i_pol) {
  // The code is largely identical between the two apply_iy_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Size ns = iy.ncols();

  ARTS_ASSERT(f_grid.size() == nf);
  ARTS_ASSERT(i_pol.size() == ns);

  if (iy_unit == "1") {
    if (n != 1) {
      iy *= (n * n);
    }
  }

  else if (iy_unit == "RJBT") {
    for (Index iv = 0; iv < nf; iv++) {
      const Numeric scfac = invrayjean(1, f_grid[iv]);
      for (Size is = 0; is < ns; is++) {
        if (i_pol[is] < 5)  // Stokes components
        {
          iy(iv, is) *= scfac;
        } else  // Measuement single pols
        {
          iy(iv, is) *= 2 * scfac;
        }
      }
    }
  }

  else if (iy_unit == "PlanckBT") {
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = ns - 1; is >= 0; is--)  // Order must here be reversed
      {
        if (i_pol[is] == 1) {
          iy(iv, is) = invplanck(iy(iv, is), f_grid[iv]);
        } else if (i_pol[is] < 5) {
          ARTS_ASSERT(i_pol[0] == 1);
          iy(iv, is) = invplanck(0.5 * (iy(iv, 0) + iy(iv, is)), f_grid[iv]) -
                       invplanck(0.5 * (iy(iv, 0) - iy(iv, is)), f_grid[iv]);
        } else {
          iy(iv, is) = invplanck(2 * iy(iv, is), f_grid[iv]);
        }
      }
    }
  }

  else if (iy_unit == "W/(m^2 m sr)") {
    for (Index iv = 0; iv < nf; iv++) {
      const Numeric scfac = n * n * f_grid[iv] * (f_grid[iv] / SPEED_OF_LIGHT);
      for (Size is = 0; is < ns; is++) {
        iy(iv, is) *= scfac;
      }
    }
  }

  else if (iy_unit == "W/(m^2 m-1 sr)") {
    iy *= (n * n * SPEED_OF_LIGHT);
  }

  else {
    ARTS_USER_ERROR (
      "Unknown option: iy_unit = \"", iy_unit, "\"\n"
      "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
      "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\"")
  }
}

void apply_iy_unit2(Tensor3View J,
                    const ConstMatrixView& iy,
                    const String& iy_unit,
                    const ConstVectorView& f_grid,
                    const Numeric& n,
                    const ArrayOfIndex& i_pol) {
  // The code is largely identical between the two apply_iy_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();
  const Index np = J.npages();

  ARTS_ASSERT(J.nrows() == nf);
  ARTS_ASSERT(J.ncols() == ns);
  ARTS_ASSERT(f_grid.size() == nf);
  ARTS_ASSERT(i_pol.size() == static_cast<Size>(ns));

  if (iy_unit == "1") {
    if (n != 1) {
      J *= (n * n);
    }
  }

  else if (iy_unit == "RJBT") {
    for (Index iv = 0; iv < nf; iv++) {
      const Numeric scfac = invrayjean(1, f_grid[iv]);
      for (Index is = 0; is < ns; is++) {
        if (i_pol[is] < 5)  // Stokes componenets
        {
          for (Index ip = 0; ip < np; ip++) {
            J(ip, iv, is) *= scfac;
          }
        } else  // Measuement single pols
        {
          for (Index ip = 0; ip < np; ip++) {
            J(ip, iv, is) *= 2 * scfac;
          }
        }
      }
    }
  }

  else if (iy_unit == "PlanckBT") {
    for (Index iv = 0; iv < f_grid.size(); iv++) {
      for (Index is = ns - 1; is >= 0; is--) {
        Numeric scfac = 1;
        if (i_pol[is] == 1) {
          scfac = dinvplanckdI(iy(iv, is), f_grid[iv]);
        } else if (i_pol[is] < 5) {
          ARTS_ASSERT(i_pol[0] == 1);
          scfac = dinvplanckdI(0.5 * (iy(iv, 0) + iy(iv, is)), f_grid[iv]) +
                  dinvplanckdI(0.5 * (iy(iv, 0) - iy(iv, is)), f_grid[iv]);
        } else {
          scfac = dinvplanckdI(2 * iy(iv, is), f_grid[iv]);
        }
        //
        for (Index ip = 0; ip < np; ip++) {
          J(ip, iv, is) *= scfac;
        }
      }
    }
  }

  else if (iy_unit == "W/(m^2 m sr)") {
    for (Index iv = 0; iv < nf; iv++) {
      const Numeric scfac = n * n * f_grid[iv] * (f_grid[iv] / SPEED_OF_LIGHT);
      for (Index ip = 0; ip < np; ip++) {
        for (Index is = 0; is < ns; is++) {
          J(ip, iv, is) *= scfac;
        }
      }
    }
  }

  else if (iy_unit == "W/(m^2 m-1 sr)") {
    J *= (n * n * SPEED_OF_LIGHT);
  }

  else {
    ARTS_USER_ERROR (
      "Unknown option: iy_unit = \"", iy_unit, "\"\n"
      "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
      "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\"")
  }
}

Numeric dotprod_with_los(const ConstVectorView& los,
                         const Numeric& u,
                         const Numeric& v,
                         const Numeric& w) {
  // Strength of field
  const Numeric f = sqrt(u * u + v * v + w * w);

  // Zenith and azimuth angle for field (in radians)
  const Numeric za_f = acos(w / f);
  const Numeric aa_f = atan2(u, v);

  // Zenith and azimuth angle for photon direction (in radians)
  Vector los_p;
  mirror_los(los_p, los);
  const Numeric za_p = Conversion::deg2rad(1) * los_p[0];
  const Numeric aa_p = Conversion::deg2rad(1) * los_p[1];

  return f * (cos(za_f) * cos(za_p) + sin(za_f) * sin(za_p) * cos(aa_f - aa_p));
}

void get_stepwise_blackbody_radiation(VectorView B, VectorView dB_dT,
                                      const ConstVectorView &ppath_f_grid,
                                      const Numeric &ppath_temperature,
                                      const bool &do_temperature_derivative) {
  std::transform(ppath_f_grid.begin(), ppath_f_grid.end(), B.begin(),
                 [T = ppath_temperature](auto &&f) { return planck(f, T); });

  if (do_temperature_derivative)
    std::transform(
        ppath_f_grid.begin(), ppath_f_grid.end(), dB_dT.begin(),
        [T = ppath_temperature](auto &&f) { return dplanck_dt(f, T); });
}

void get_stepwise_blackbody_radiation(
    Vector &B, Matrix &dB, const Vector &ppath_f_grid,
    const Numeric &ppath_temperature,
    const ArrayOfRetrievalQuantity &jacobian_quantities,
    const bool j_analytical_do) {
  B.resize(ppath_f_grid.size());
  std::transform(ppath_f_grid.begin(), ppath_f_grid.end(), B.begin(),
                 [T = ppath_temperature](auto &&f) { return planck(f, T); });

  if (j_analytical_do) {
    dB.resize(jacobian_quantities.size(), ppath_f_grid.size());
    for (Size i = 0; i < jacobian_quantities.size(); ++i) {
      if (jacobian_quantities[i] == Jacobian::Atm::Temperature) {
        std::transform(
            ppath_f_grid.begin(), ppath_f_grid.end(), dB[i].begin(),
            [T = ppath_temperature](auto &&f) { return dplanck_dt(f, T); });
      } else {
        dB[i] = 0.0;
      }
    }
  }
}

void get_stepwise_clearsky_propmat(
    const Workspace& ws,
    PropmatVector& K,
    StokvecVector& S,
    PropmatMatrix& dK_dx,
    StokvecMatrix& dS_dx,
    const Agenda& propmat_clearsky_agenda,
    const JacobianTargets& jacobian_targets,
    const Vector& ppath_f_grid,
    const PropagationPathPoint& path_point,
    const AtmPoint& atm_point) {
  static const ArrayOfSpeciesTag select_abs_species{};
  static const ArrayOfRetrievalQuantity jacobian_quantities_empty{};

  // Perform the propagation matrix computations
  propmat_clearsky_agendaExecute(
      ws,
      K,
      S,
      dK_dx,
      dS_dx,
      jacobian_targets,
      select_abs_species,
      ppath_f_grid,
      path_point,
      atm_point,
      propmat_clearsky_agenda);

  const Vector sensor_like_los{path::mirror(path_point.los)};
  adapt_stepwise_partial_derivatives(
      dK_dx, dS_dx, jacobian_targets, ppath_f_grid, sensor_like_los);
}

Vector get_stepwise_f_partials(const ConstVectorView& line_of_sight,
                               const ConstVectorView& f_grid,
                               const Atm::Key wind_type) {
  // Doppler relevant velocity
  Numeric dv_doppler_dx = 0.0;
  
  Vector deriv(f_grid);
  
  switch (wind_type) {
    case Atm::Key::wind_u:
      dv_doppler_dx =
          (dotprod_with_los(line_of_sight, 1, 0, 0));
      break;
    case Atm::Key::wind_v:
      dv_doppler_dx =
          (dotprod_with_los(line_of_sight, 0, 1, 0));
      break;
    case Atm::Key::wind_w:
      dv_doppler_dx =
          (dotprod_with_los(line_of_sight, 0, 0, 1));
      break;
    default:
      ARTS_ASSERT(false, "Not allowed to call this function without a wind parameter as wind_type");
      break;
  }
  
  deriv *= - dv_doppler_dx / Constant::c;
  return deriv;
}

void get_stepwise_scattersky_propmat(
    StokvecVector& ap,
    PropmatVector& Kp,
    StokvecMatrix& dap_dx,
    PropmatMatrix& dKp_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstMatrixView& ppath_1p_pnd,  // the ppath_pnd at this ppath point
    const ArrayOfMatrix&
        ppath_dpnd_dx,  // the full ppath_dpnd_dx, ie all ppath points
    const Index ppath_1p_id,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstVectorView& ppath_line_of_sight,
    const ConstVectorView& ppath_temperature,
    const bool& jacobian_do) {
  const Size nf = Kp.size();

  ArrayOfArrayOfSingleScatteringData scat_data_mono;

  // Direction of outgoing scattered radiation (which is reversed to
  // LOS). Only used for extracting scattering properties.
  Vector dir;
  mirror_los(dir, ppath_line_of_sight);
  Matrix dir_array(1, 2, 0.);
  dir_array(0, joker) = dir;

  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor5 ext_mat_bulk;
  Tensor4 abs_vec_bulk;
  Index ptype_bulk;

  // get per-scat-elem data here. and fold with pnd.
  // keep per-scat-elem data to fold with dpnd_dx further down in
  // analyt-jac-loop.
  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      Vector{ppath_temperature},
                      dir_array,
                      -1);

  opt_prop_ScatSpecBulk(ext_mat_ssbulk,
                        abs_vec_ssbulk,
                        ptype_ssbulk,
                        ext_mat_Nse,
                        abs_vec_Nse,
                        ptypes_Nse,
                        ppath_1p_pnd,
                        t_ok);
  opt_prop_Bulk(ext_mat_bulk,
                abs_vec_bulk,
                ptype_bulk,
                ext_mat_ssbulk,
                abs_vec_ssbulk,
                ptype_ssbulk);

  for (Size iv = 0; iv < nf; iv++) {
    ap[iv] = rtepack::to_stokvec(abs_vec_bulk(iv, 0, 0, joker));
    Kp[iv] = rtepack::to_propmat(ext_mat_bulk(iv, 0, 0, joker, joker));
  }

  if (jacobian_do)
    FOR_ANALYTICAL_JACOBIANS_DO(
        if (ppath_dpnd_dx[iq].empty()) {
          dap_dx[iq] = 0.0;
          dKp_dx[iq] = 0.0;
        } else {
          // check, whether we have any non-zero ppath_dpnd_dx in this
          // pnd-affecting x? might speed up things a little bit.
          opt_prop_ScatSpecBulk(ext_mat_ssbulk,
                                abs_vec_ssbulk,
                                ptype_ssbulk,
                                ext_mat_Nse,
                                abs_vec_Nse,
                                ptypes_Nse,
                                ppath_dpnd_dx[iq](joker, Range(ppath_1p_id, 1)),
                                t_ok);
          opt_prop_Bulk(ext_mat_bulk,
                        abs_vec_bulk,
                        ptype_bulk,
                        ext_mat_ssbulk,
                        abs_vec_ssbulk,
                        ptype_ssbulk);
          for (Size iv = 0; iv < nf; iv++) {
              dap_dx[iq][iv] = rtepack::to_stokvec(abs_vec_bulk(iv, 0, 0, joker));
              dKp_dx[iq][iv] = rtepack::to_propmat(ext_mat_bulk(iv, 0, 0, joker, joker));
          }
        })
}

void get_stepwise_scattersky_source(
    StokvecVector& Sp,
    StokvecMatrix& dSp_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstVectorView& ppath_1p_pnd,  // the ppath_pnd at this ppath point
    const ArrayOfMatrix&
        ppath_dpnd_dx,  // the full ppath_dpnd_dx, ie all ppath points
    const Index ppath_1p_id,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstTensor7View& cloudbox_field,
    const ConstVectorView& za_grid,
    const ConstVectorView& aa_grid,
    const ConstMatrixView& ppath_line_of_sight,
    const GridPos& ppath_pressure,
    const Vector& temperature,
    const bool& jacobian_do,
    const Index& t_interp_order) {
  ARTS_USER_ERROR ("This function handles so far only 1D atmospheres.");

  const Index nf = Sp.size();
  const Index ne = ppath_1p_pnd.size();
  ARTS_ASSERT(TotalNumberOfElements(scat_data) == ne);
  const Index nza = za_grid.size();
  const Index naa = aa_grid.size();
  const Index nq = jacobian_do ? jacobian_quantities.size() : 0;

  // interpolate incident field to this ppath point (no need to do this
  // separately per scatelem)
  GridPos gp_p;
  gridpos_copy(gp_p, ppath_pressure);
  Vector itw_p(2);
  interpweights(itw_p, gp_p);
  Tensor3 inc_field(nf, nza, 4, 0.);
  for (Index iv = 0; iv < nf; iv++) {
    for (Index iza = 0; iza < nza; iza++) {
      for (Index i = 0; i < 4; i++) {
        inc_field(iv, iza, i) =
            interp(itw_p, cloudbox_field(iv, joker, 0, 0, iza, 0, i), gp_p);
      }
    }
  }

  // create matrix of incident directions (flat representation of the
  // za_grid * aa_grid matrix)
  Matrix idir(nza * naa, 2);
  Index ia = 0;
  for (Index iza = 0; iza < nza; iza++) {
    for (Index iaa = 0; iaa < naa; iaa++) {
      idir(ia, 0) = za_grid[iza];
      idir(ia, 1) = aa_grid[iaa];
      ia++;
    }
  }

  // setting prop (aka scattered) direction
  Matrix pdir(1, 2);
  //if( ppath_line_of_sight.ncols()==2 )
  //  pdir(0,joker) = ppath_line_of_sight;
  //else // 1D case only (currently the only handled case). azimuth not defined.
  {
    pdir(0, 0) = ppath_line_of_sight(0, 0);
    pdir(0, 1) = 0.;
  }

  // some further variables needed for pha_mat extraction
  Index nf_ssd = scat_data[0][0].pha_mat_data.nlibraries();
  Index duplicate_freqs = ((nf == nf_ssd) ? 0 : 1);
  Tensor6 pha_mat_1se(nf_ssd, 1, 1, nza * naa, 4, 4);
  Vector t_ok(1);
  Index ptype;
  Tensor3 scat_source_1se(ne, nf, 4, 0.);

  Index ise_flat = 0;
  for (Size i_ss = 0; i_ss < scat_data.size(); i_ss++) {
    for (Size i_se = 0; i_se < scat_data[i_ss].size(); i_se++) {
      // determine whether we have some valid pnd for this
      // scatelem (in pnd or dpnd)
      Index val_pnd = 0;
      if (ppath_1p_pnd[ise_flat] != 0) {
        val_pnd = 1;
      } else if (jacobian_do) {
        for (Index iq = 0; (!val_pnd) && (iq < nq); iq++) {
          if ((not(jacobian_quantities[iq] == Jacobian::Type::Sensor) and not(jacobian_quantities[iq] == Jacobian::Special::SurfaceString)) &&
              !ppath_dpnd_dx[iq].empty() &&
              ppath_dpnd_dx[iq](ise_flat, ppath_1p_id) != 0) {
            val_pnd = 1;
          }
        }
      }

      if (val_pnd) {
        pha_mat_1ScatElem(pha_mat_1se,
                          ptype,
                          t_ok,
                          scat_data[i_ss][i_se],
                          temperature,
                          pdir,
                          idir,
                          0,
                          t_interp_order);
        ARTS_USER_ERROR_IF (t_ok[0] == 0,
            "Interpolation error for (flat-array) scattering "
            "element #", ise_flat, "\n"
            "at location/temperature point #", ppath_1p_id, "\n")

        Index this_iv = 0;
        for (Index iv = 0; iv < nf; iv++) {
          if (!duplicate_freqs) {
            this_iv = iv;
          }
          Tensor3 product_fields(nza, naa, 4, 0.);

          ia = 0;
          for (Index iza = 0; iza < nza; iza++) {
            for (Index iaa = 0; iaa < naa; iaa++) {
              for (Index i = 0; i < 4; i++) {
                for (Index j = 0; j < 4; j++) {
                  product_fields(iza, iaa, i) +=
                      pha_mat_1se(this_iv, 0, 0, ia, i, j) *
                      inc_field(iv, iza, j);
                }
              }
              ia++;
            }
          }

          for (Index i = 0; i < 4; i++) {
            scat_source_1se(ise_flat, iv, i) = AngIntegrate_trapezoid(
                product_fields(joker, joker, i), za_grid, aa_grid);
          }
        }  // for iv
      }    // if val_pnd

      ise_flat++;

    }  // for i_se
  }    // for i_ss

  for (Index iv = 0; iv < nf; iv++) {
    Vector scat_source(4, 0.);
    for (ise_flat = 0; ise_flat < ne; ise_flat++) {
      for (Index i = 0; i < 4; i++) {
        scat_source[i] +=
            scat_source_1se(ise_flat, iv, i) * ppath_1p_pnd[ise_flat];
      }
    }

    Sp[iv] = rtepack::to_stokvec(scat_source);

    if (jacobian_do) {
      FOR_ANALYTICAL_JACOBIANS_DO(
          if (ppath_dpnd_dx[iq].empty()) { dSp_dx[iq] = 0.0; } else {
            scat_source = 0.;
            for (ise_flat = 0; ise_flat < ne; ise_flat++) {
              for (Index i = 0; i < 4; i++) {
                scat_source[i] += scat_source_1se(ise_flat, iv, i) *
                                  ppath_dpnd_dx[iq](ise_flat, ppath_1p_id);
                dSp_dx[iq][iv] = rtepack::to_stokvec(scat_source);
              }
            }
          })
    }
  }  // for iv
}

void iy_transmittance_mult(Tensor3& iy_trans_total,
                          const ConstTensor3View& iy_trans_old,
                          const ConstTensor3View& iy_trans_new) {
  const Index nf = iy_trans_old.npages();
  const Index ns = iy_trans_old.ncols();

  ARTS_ASSERT(ns == iy_trans_old.nrows());
  ARTS_ASSERT(nf == iy_trans_new.npages());
  ARTS_ASSERT(ns == iy_trans_new.nrows());
  ARTS_ASSERT(ns == iy_trans_new.ncols());

  iy_trans_total.resize(nf, ns, ns);

  for (Index iv = 0; iv < nf; iv++) {
    mult(iy_trans_total(iv, joker, joker),
         iy_trans_old(iv, joker, joker),
         iy_trans_new(iv, joker, joker));
  }
}

void iy_transmittance_mult(Matrix& iy_new,
                          const ConstTensor3View& iy_trans,
                          const ConstMatrixView& iy_old) {
  const Index nf = iy_trans.npages();
  const Index ns = iy_trans.ncols();

  ARTS_ASSERT(ns == iy_trans.nrows());
  ARTS_ASSERT(nf == iy_old.nrows());
  ARTS_ASSERT(ns == iy_old.ncols());

  iy_new.resize(nf, ns);

  for (Index iv = 0; iv < nf; iv++) {
    mult(iy_new(iv, joker), iy_trans(iv, joker, joker), iy_old(iv, joker));
  }
}

void mirror_los(Vector& los_mirrored,
                const ConstVectorView& los) {
  los_mirrored.resize(2);
  //
  los_mirrored[0] = 180 - los[0];
  los_mirrored[1] = los[1] + 180;
  if (los_mirrored[1] > 180) {
    los_mirrored[1] -= 360;
  }
}

void mueller_modif2stokes(Matrix &Cs) {
  //
  Cs.resize(4, 4);
  Cs(0, 0) = 1;
  Cs(0, 1) = Cs(1, 0) = 1;
  Cs(1, 1) = -1;
  Cs(0, 2) = Cs(1, 2) = Cs(2, 0) = Cs(2, 1) = 0;
  Cs(2, 2) = 1;
  Cs(0, 3) = Cs(1, 3) = Cs(2, 3) = Cs(3, 0) = Cs(3, 1) = Cs(3, 2) = 0;
  Cs(3, 3) = 1;
}

void mueller_rotation(Matrix& L,
                      const Numeric& rotangle) {
  //
  L.resize(4, 4);
  L(0, 0) = 1;
    const Numeric alpha = 2 * Conversion::deg2rad(1) * rotangle;
    const Numeric c2 = cos(alpha);
    L(0,1) = L(1,0) = 0;
    L(1,1) = c2;
      const Numeric s2 = sin(alpha);
      L(0,2) = L(2,0) = 0;
      L(1,2) = s2;
      L(2,1) = -s2;      
      L(2,2) = c2;
        L(0,3) = L(1,3) = L(2,3) = L(3,0) = L(3,1) = L(3,2) = 0;
        L(3,3) = 1;   
}

void mueller_stokes2modif(Matrix& Cm) {
  //
  Cm.resize(4, 4);
  Cm(0,0) = 0.5;
    Cm(0,1) = Cm(1,0) = 0.5;
    Cm(1,1) = -0.5;
      Cm(0,2) = Cm(1,2) = Cm(2,0) = Cm(2,1) = 0;
      Cm(2,2) = 1;
        Cm(0,3) = Cm(1,3) = Cm(2,3) = Cm(3,0) = Cm(3,1) = Cm(3,2) = 0;
        Cm(3,3) = 1;   
}

void ze_cfac(Vector& fac,
             const Vector& f_grid,
             const Numeric& ze_tref,
             const Numeric& k2) {
  const Index nf = f_grid.size();

  ARTS_ASSERT(fac.size() == nf);

  // Refractive index for water (if needed)
  Matrix complex_n(0, 0);
  if (k2 <= 0) {
    complex_n_water_liebe93(complex_n, f_grid, ze_tref);
  }

  // Common conversion factor
  static constexpr Numeric a = 4e18 / Math::pow4(Constant::pi);

  for (Index iv = 0; iv < nf; iv++) {
    // Reference dielectric factor.
    Numeric K2;
    if (k2 >= 0) {
      K2 = k2;
    } else {
      Complex n(complex_n(iv, 0), complex_n(iv, 1));
      Complex n2 = n * n;
      Complex K = (n2 - Numeric(1.0)) / (n2 + Numeric(2.0));
      Numeric absK = abs(K);
      K2 = absK * absK;
    }

    // Wavelength
    Numeric la = SPEED_OF_LIGHT / f_grid[iv];

    fac[iv] = a * la * la * la * la / K2;
  }
}
