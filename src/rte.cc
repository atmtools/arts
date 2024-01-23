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
#include "path_point.h"
#include "physics_funcs.h"
#include "refraction.h"
#include "rtepack.h"
#include "special_interp.h"
#include "species_tags.h"

inline constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;
inline constexpr Numeric TEMP_0_C = Constant::temperature_at_0c;

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
    ARTS_USER_ERROR("Unknown option: iy_unit = \"",
                    iy_unit,
                    "\"\n"
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
    ARTS_USER_ERROR("Unknown option: iy_unit = \"",
                    iy_unit,
                    "\"\n"
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

void get_stepwise_blackbody_radiation(VectorView B,
                                      VectorView dB_dT,
                                      const ConstVectorView& ppath_f_grid,
                                      const Numeric& ppath_temperature,
                                      const bool& do_temperature_derivative) {
  std::transform(ppath_f_grid.begin(),
                 ppath_f_grid.end(),
                 B.begin(),
                 [T = ppath_temperature](auto&& f) { return planck(f, T); });

  if (do_temperature_derivative)
    std::transform(
        ppath_f_grid.begin(),
        ppath_f_grid.end(),
        dB_dT.begin(),
        [T = ppath_temperature](auto&& f) { return dplanck_dt(f, T); });
}

void get_stepwise_clearsky_propmat(const Workspace& ws,
                                   PropmatVector& K,
                                   StokvecVector& S,
                                   PropmatMatrix& dK_dx,
                                   StokvecMatrix& dS_dx,
                                   const Agenda& propagation_matrix_agenda,
                                   const JacobianTargets& jacobian_targets,
                                   const Vector& ppath_f_grid,
                                   const PropagationPathPoint& path_point,
                                   const AtmPoint& atm_point) {
  static const ArrayOfSpeciesTag select_abs_species{};

  // Perform the propagation matrix computations
  propagation_matrix_agendaExecute(ws,
                                   K,
                                   S,
                                   dK_dx,
                                   dS_dx,
                                   jacobian_targets,
                                   select_abs_species,
                                   ppath_f_grid,
                                   path_point,
                                   atm_point,
                                   propagation_matrix_agenda);

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
      dv_doppler_dx = (dotprod_with_los(line_of_sight, 1, 0, 0));
      break;
    case Atm::Key::wind_v:
      dv_doppler_dx = (dotprod_with_los(line_of_sight, 0, 1, 0));
      break;
    case Atm::Key::wind_w:
      dv_doppler_dx = (dotprod_with_los(line_of_sight, 0, 0, 1));
      break;
    default:
      ARTS_ASSERT(
          false,
          "Not allowed to call this function without a wind parameter as wind_type");
      break;
  }

  deriv *= -dv_doppler_dx / Constant::c;
  return deriv;
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

void mirror_los(Vector& los_mirrored, const ConstVectorView& los) {
  los_mirrored.resize(2);
  //
  los_mirrored[0] = 180 - los[0];
  los_mirrored[1] = los[1] + 180;
  if (los_mirrored[1] > 180) {
    los_mirrored[1] -= 360;
  }
}

void mueller_modif2stokes(Matrix& Cs) {
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

void mueller_rotation(Matrix& L, const Numeric& rotangle) {
  //
  L.resize(4, 4);
  L(0, 0) = 1;
  const Numeric alpha = 2 * Conversion::deg2rad(1) * rotangle;
  const Numeric c2 = cos(alpha);
  L(0, 1) = L(1, 0) = 0;
  L(1, 1) = c2;
  const Numeric s2 = sin(alpha);
  L(0, 2) = L(2, 0) = 0;
  L(1, 2) = s2;
  L(2, 1) = -s2;
  L(2, 2) = c2;
  L(0, 3) = L(1, 3) = L(2, 3) = L(3, 0) = L(3, 1) = L(3, 2) = 0;
  L(3, 3) = 1;
}

void mueller_stokes2modif(Matrix& Cm) {
  //
  Cm.resize(4, 4);
  Cm(0, 0) = 0.5;
  Cm(0, 1) = Cm(1, 0) = 0.5;
  Cm(1, 1) = -0.5;
  Cm(0, 2) = Cm(1, 2) = Cm(2, 0) = Cm(2, 1) = 0;
  Cm(2, 2) = 1;
  Cm(0, 3) = Cm(1, 3) = Cm(2, 3) = Cm(3, 0) = Cm(3, 1) = Cm(3, 2) = 0;
  Cm(3, 3) = 1;
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
