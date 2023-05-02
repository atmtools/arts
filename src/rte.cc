/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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
#include "arts_constexpr_math.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic_OLD.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_concepts.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath_OLD.h"
#include "refraction.h"
#include "special_interp.h"
#include "species_tags.h"
#include <cmath>
#include <stdexcept>

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric TEMP_0_C=Constant::temperature_at_0c;

/*===========================================================================
  === The functions in alphabetical order
  ===========================================================================*/

void adapt_stepwise_partial_derivatives(
    ArrayOfPropagationMatrix& dK_dx,
    ArrayOfStokesVector& dS_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstVectorView& ppath_f_grid,
    const ConstVectorView& ppath_line_of_sight,
    const Index& lte,
    const Index& atmosphere_dim,
    const bool& jacobian_do) {
  if (not jacobian_do) return;

  // All relevant quantities are extracted first
  const Index nq = jacobian_quantities.nelem();

  // Computational temporary vector
  Vector a;

  for (Index i = 0; i < nq; i++) {
    if (jacobian_quantities[i] == Jacobian::Type::Sensor or jacobian_quantities[i] == Jacobian::Special::SurfaceString) continue;

    if (jacobian_quantities[i].is_wind()) {
      const auto scale = get_stepwise_f_partials(ppath_line_of_sight, ppath_f_grid, jacobian_quantities[i].Target().atm, atmosphere_dim);
      dK_dx[i] *= scale;
      if (not lte) dS_dx[i] *= scale;
    }
  }
}

void adjust_los(VectorView los, const Index& atmosphere_dim) {
  if (atmosphere_dim == 1) {
    if (los[0] < 0) {
      los[0] = -los[0];
    } else if (los[0] > 180) {
      los[0] = 360 - los[0];
    }
  } else if (atmosphere_dim == 2) {
    if (los[0] < -180) {
      los[0] = los[0] + 360;
    } else if (los[0] > 180) {
      los[0] = los[0] - 360;
    }
  } else {
    // If any of the angles out-of-bounds, use cart2zaaa to resolve
    if (abs(los[0] - 90) > 90 || abs(los[1]) > 180) {
      Numeric dx, dy, dz;
      zaaa2cart(dx, dy, dz, los[0], los[1]);
      cart2zaaa(los[0], los[1], dx, dy, dz);
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
  const Index ns = iy.ncols();

  ARTS_ASSERT(f_grid.nelem() == nf);
  ARTS_ASSERT(i_pol.nelem() == ns);

  if (iy_unit == "1") {
    if (n != 1) {
      iy *= (n * n);
    }
  }

  else if (iy_unit == "RJBT") {
    for (Index iv = 0; iv < nf; iv++) {
      const Numeric scfac = invrayjean(1, f_grid[iv]);
      for (Index is = 0; is < ns; is++) {
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
      for (Index is = 0; is < ns; is++) {
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
  ARTS_ASSERT(f_grid.nelem() == nf);
  ARTS_ASSERT(i_pol.nelem() == ns);

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
    for (Index iv = 0; iv < f_grid.nelem(); iv++) {
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

void bending_angle1d(Numeric& alpha, const Ppath& ppath) {
  Numeric theta;
  if (ppath.dim < 3) {
    theta = abs(ppath.start_pos[1] - ppath.end_pos[1]);
  } else {
    theta = sphdist(ppath.start_pos[1],
                    ppath.start_pos[2],
                    ppath.end_pos[1],
                    ppath.end_pos[2]);
  }

  // Eq 17 in Kursinski et al., TAO, 2000:
  alpha = ppath.start_los[0] - ppath.end_los[0] + theta;

  // This as
  // phi_r = 180 - ppath.end_los[0]
  // phi_t = ppath.start_los[0]
}

/** Just to avoid duplicatuion of code in *defocusing_general*.
   
    rte_los is mainly an input, but is also returned "adjusted" (with zenith
    and azimuth angles inside defined ranges) 
 
    @param[out]   pos                 Position of ppath at optical distance lo0
    @param[in,out]   rte_los          Direction for transmitted signal 
                                      (disturbed from nominal value)
    @param[out]   rte_pos             Position of transmitter.
    @param[out]   background          Raditaive background of ppath.
    @param[in]    lo0                 Optical path length between transmitter 
                                      and receiver.
    @param[in]    ppath_step_agenda   As the WSV with the same name.
    @param[in]    atmosphere_dim      As the WSV with the same name.
    @param[in]    p_grid              As the WSV with the same name.
    @param[in]    lat_grid            As the WSV with the same name.
    @param[in]    lon_grid            As the WSV with the same name.
    @param[in]    z_field             As the WSV with the same name.
    @param[in]    f_grid              As the WSV with the same name.
    @param[in]    refellipsoid        As the WSV with the same name.
    @param[in]    z_surface           As the WSV with the same name.
    @param[in]    verbosity           As the WSV with the same name.

    @author Patrick Eriksson 
    @date   2012-04-11
 */
void defocusing_general_sub(Workspace& ws,
                            Vector& pos,
                            Vector& rte_los,
                            Index& background,
                            const Vector& rte_pos,
                            const Numeric& lo0,
                            const Agenda& ppath_step_agenda,
                            const Numeric& ppath_lmax,
                            const Numeric& ppath_lraytrace,
                            const Index& atmosphere_dim,
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid,
                            const Tensor3& z_field,
                            const Vector& f_grid,
                            const Vector& refellipsoid,
                            const Matrix& z_surface,
                            const Verbosity& verbosity) {
  // Special treatment of 1D around zenith/nadir
  // (zenith angles outside [0,180] are changed by *adjust_los*)
  bool invert_lat = false;
  if (atmosphere_dim == 1 && (rte_los[0] < 0 || rte_los[0] > 180)) {
    invert_lat = true;
  }

  // Handle cases where angles have moved out-of-bounds due to disturbance
  adjust_los(rte_los, atmosphere_dim);

  // Calculate the ppath for disturbed rte_los
  Ppath ppx;
  //
  ppath_calc(ws,
             ppx,
             ppath_step_agenda,
             atmosphere_dim,
             p_grid,
             lat_grid,
             lon_grid,
             z_field,
             f_grid,
             refellipsoid,
             z_surface,
             0,
             ArrayOfIndex(0),
             rte_pos,
             rte_los,
             ppath_lmax,
             ppath_lraytrace,
             false,
             verbosity);
  //
  background = ppath_what_background(ppx);

  // Calcualte cumulative optical path for ppx
  Vector lox(ppx.np);
  Index ilast = ppx.np - 1;
  lox[0] = ppx.end_lstep;
  for (Index i = 1; i <= ilast; i++) {
    lox[i] =
        lox[i - 1] + ppx.lstep[i - 1] * (ppx.nreal[i - 1] + ppx.nreal[i]) / 2.0;
  }

  pos.resize(max(Index(2), atmosphere_dim));

  // Reciever at a longer distance (most likely out in space):
  if (lox[ilast] < lo0) {
    const Numeric dl = lo0 - lox[ilast];
    if (atmosphere_dim < 3) {
      Numeric x, z, dx, dz;
      poslos2cart(
          x, z, dx, dz, ppx.r[ilast], ppx.pos(ilast, 1), ppx.los(ilast, 0));
      cart2pol(pos[0],
               pos[1],
               x + dl * dx,
               z + dl * dz,
               ppx.pos(ilast, 1),
               ppx.los(ilast, 0));
    } else {
      Numeric x, y, z, dx, dy, dz;
      poslos2cart(x,
                  y,
                  z,
                  dx,
                  dy,
                  dz,
                  ppx.r[ilast],
                  ppx.pos(ilast, 1),
                  ppx.pos(ilast, 2),
                  ppx.los(ilast, 0),
                  ppx.los(ilast, 1));
      cart2sph(pos[0],
               pos[1],
               pos[2],
               x + dl * dx,
               y + dl * dy,
               z + dl * dz,
               ppx.pos(ilast, 1),
               ppx.pos(ilast, 2),
               ppx.los(ilast, 0),
               ppx.los(ilast, 1));
    }
  }

  // Interpolate to lo0
  else {
    GridPos gp;
    Vector itw(2);
    gridpos(gp, lox, lo0);
    interpweights(itw, gp);
    //
    pos[0] = interp(itw, ppx.r, gp);
    pos[1] = interp(itw, ppx.pos(joker, 1), gp);
    if (atmosphere_dim == 3) {
      pos[2] = interp(itw, ppx.pos(joker, 2), gp);
    }
  }

  if (invert_lat) {
    pos[1] = -pos[1];
  }
}

void defocusing_general(Workspace& ws,
                        Numeric& dlf,
                        const Agenda& ppath_step_agenda,
                        const Index& atmosphere_dim,
                        const Vector& p_grid,
                        const Vector& lat_grid,
                        const Vector& lon_grid,
                        const Tensor3& z_field,
                        const Vector& f_grid,
                        const Vector& refellipsoid,
                        const Matrix& z_surface,
                        const Ppath& ppath,
                        const Numeric& ppath_lmax,
                        const Numeric& ppath_lraytrace,
                        const Numeric& dza,
                        const Verbosity& verbosity) {
  // Optical and physical path between transmitter and reciver
  Numeric lo = ppath.start_lstep + ppath.end_lstep;
  Numeric lp = lo;
  for (Index i = 0; i <= ppath.np - 2; i++) {
    lp += ppath.lstep[i];
    lo += ppath.lstep[i] * (ppath.nreal[i] + ppath.nreal[i + 1]) / 2.0;
  }
  // Extract rte_pos and rte_los
  const Vector rte_pos{ppath.start_pos[Range(0, atmosphere_dim)]};
  //
  Vector rte_los0(max(Index(1), atmosphere_dim - 1)), rte_los;
  mirror_los(rte_los, ppath.start_los, atmosphere_dim);
  rte_los0 = rte_los[Range(0, max(Index(1), atmosphere_dim - 1))];

  // A new ppath with positive zenith angle off-set
  //
  Vector pos1;
  Index backg1;
  //
  rte_los = rte_los0;
  rte_los[0] += dza;
  //
  defocusing_general_sub(ws,
                         pos1,
                         rte_los,
                         backg1,
                         rte_pos,
                         lo,
                         ppath_step_agenda,
                         ppath_lmax,
                         ppath_lraytrace,
                         atmosphere_dim,
                         p_grid,
                         lat_grid,
                         lon_grid,
                         z_field,
                         f_grid,
                         refellipsoid,
                         z_surface,
                         verbosity);

  // Same thing with negative zenit angle off-set
  Vector pos2;
  Index backg2;
  //
  rte_los = rte_los0;  // Use rte_los0 as rte_los can have been "adjusted"
  rte_los[0] -= dza;
  //
  defocusing_general_sub(ws,
                         pos2,
                         rte_los,
                         backg2,
                         rte_pos,
                         lo,
                         ppath_step_agenda,
                         ppath_lmax,
                         ppath_lraytrace,
                         atmosphere_dim,
                         p_grid,
                         lat_grid,
                         lon_grid,
                         z_field,
                         f_grid,
                         refellipsoid,
                         z_surface,
                         verbosity);

  // Calculate distance between pos1 and 2, and derive the loss factor
  // All appears OK:
  if (backg1 == backg2) {
    Numeric l12;
    if (atmosphere_dim < 3) {
      distance2D(l12, pos1[0], pos1[1], pos2[0], pos2[1]);
    } else {
      distance3D(l12, pos1[0], pos1[1], pos1[2], pos2[0], pos2[1], pos2[2]);
    }
    //
    dlf = lp * 2 * Conversion::deg2rad(1) * dza / l12;
  }
  // If different backgrounds, then only use the second calculation
  else {
    Numeric l12;
    if (atmosphere_dim == 1) {
      const Numeric r = refellipsoid[0];
      distance2D(l12, r + ppath.end_pos[0], 0, pos2[0], pos2[1]);
    } else if (atmosphere_dim == 2) {
      const Numeric r = refell2r(refellipsoid, ppath.end_pos[1]);
      distance2D(l12, r + ppath.end_pos[0], ppath.end_pos[1], pos2[0], pos2[1]);
    } else {
      const Numeric r = refell2r(refellipsoid, ppath.end_pos[1]);
      distance3D(l12,
                 r + ppath.end_pos[0],
                 ppath.end_pos[1],
                 ppath.end_pos[2],
                 pos2[0],
                 pos2[1],
                 pos2[2]);
    }
    //
    dlf = lp * Conversion::deg2rad(1) * dza / l12;
  }
}

void defocusing_sat2sat(Workspace& ws,
                        Numeric& dlf,
                        const Agenda& ppath_step_agenda,
                        const Index& atmosphere_dim,
                        const Vector& p_grid,
                        const Vector& lat_grid,
                        const Vector& lon_grid,
                        const Tensor3& z_field,
                        const Vector& f_grid,
                        const Vector& refellipsoid,
                        const Matrix& z_surface,
                        const Ppath& ppath,
                        const Numeric& ppath_lmax,
                        const Numeric& ppath_lraytrace,
                        const Numeric& dza,
                        const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF (ppath.end_los[0] < 90 || ppath.start_los[0] > 90,
        "The function *defocusing_sat2sat* can only be used "
        "for limb sounding geometry.");

  // Index of tangent point
  Index it;
  find_tanpoint(it, ppath);
  ARTS_ASSERT(it >= 0);

  // Length between tangent point and transmitter/reciver
  Numeric lt = ppath.start_lstep, lr = ppath.end_lstep;
  for (Index i = it; i <= ppath.np - 2; i++) {
    lt += ppath.lstep[i];
  }
  for (Index i = 0; i < it; i++) {
    lr += ppath.lstep[i];
  }

  // Bending angle and impact parameter for centre ray
  Numeric alpha0, a0;
  bending_angle1d(alpha0, ppath);
  alpha0 *= Conversion::deg2rad(1);
  a0 = ppath.constant;

  // Azimuth loss term (Eq 18.5 in Kursinski et al.)
  const Numeric lf = lr * lt / (lr + lt);
  const Numeric alt = 1 / (1 - alpha0 * lf / refellipsoid[0]);

  // Calculate two new ppaths to get dalpha/da
  Numeric alpha1, a1, alpha2, a2, dada;
  Ppath ppt;
  Vector rte_pos{ppath.end_pos[Range(0, atmosphere_dim)]};
  Vector rte_los{ppath.end_los};
  //
  rte_los[0] -= dza;
  adjust_los(rte_los, atmosphere_dim);
  ppath_calc(ws,
             ppt,
             ppath_step_agenda,
             atmosphere_dim,
             p_grid,
             lat_grid,
             lon_grid,
             z_field,
             f_grid,
             refellipsoid,
             z_surface,
             0,
             ArrayOfIndex(0),
             rte_pos,
             rte_los,
             ppath_lmax,
             ppath_lraytrace,
             false,
             verbosity);
  bending_angle1d(alpha2, ppt);
  alpha2 *= Conversion::deg2rad(1);
  a2 = ppt.constant;
  //
  rte_los[0] += 2 * dza;
  adjust_los(rte_los, atmosphere_dim);
  ppath_calc(ws,
             ppt,
             ppath_step_agenda,
             atmosphere_dim,
             p_grid,
             lat_grid,
             lon_grid,
             z_field,
             f_grid,
             refellipsoid,
             z_surface,
             0,
             ArrayOfIndex(0),
             rte_pos,
             rte_los,
             ppath_lmax,
             ppath_lraytrace,
             false,
             verbosity);
  // This path can hit the surface. And we need to check if ppt is OK.
  // (remember this function only deals with sat-to-sat links and OK
  // background here is be space)
  // Otherwise use the centre ray as the second one.
  if (ppath_what_background(ppt) == 1) {
    bending_angle1d(alpha1, ppt);
    alpha1 *= Conversion::deg2rad(1);
    a1 = ppt.constant;
    dada = (alpha2 - alpha1) / (a2 - a1);
  } else {
    dada = (alpha2 - alpha0) / (a2 - a0);
  }

  // Zenith loss term (Eq 18 in Kursinski et al.)
  const Numeric zlt = 1 / (1 - dada * lf);

  // Total defocusing loss
  dlf = zlt * alt;
}

Numeric dotprod_with_los(const ConstVectorView& los,
                         const Numeric& u,
                         const Numeric& v,
                         const Numeric& w,
                         const Index& atmosphere_dim) {
  // Strength of field
  const Numeric f = sqrt(u * u + v * v + w * w);

  // Zenith and azimuth angle for field (in radians)
  const Numeric za_f = acos(w / f);
  const Numeric aa_f = atan2(u, v);

  // Zenith and azimuth angle for photon direction (in radians)
  Vector los_p;
  mirror_los(los_p, los, atmosphere_dim);
  const Numeric za_p = Conversion::deg2rad(1) * los_p[0];
  const Numeric aa_p = Conversion::deg2rad(1) * los_p[1];

  return f * (cos(za_f) * cos(za_p) + sin(za_f) * sin(za_p) * cos(aa_f - aa_p));
}

void get_iy(Workspace& ws,
            Matrix& iy,
            const Index& cloudbox_on,
            const Vector& f_grid,
            const AtmField& atm_field,
            const Vector& rte_pos,
            const Vector& rte_los,
            const Vector& rte_pos2,
            const String& iy_unit,
            const Agenda& iy_main_agenda) {
  ArrayOfTensor3 diy_dx;
  ArrayOfMatrix iy_aux;
  Ppath ppath;
  Vector geo_pos;
  Tensor3 iy_transmittance(0, 0, 0);
  const Index iy_agenda_call1 = 1;
  const ArrayOfString iy_aux_vars(0);
  const Index iy_id = 0;
  const Index jacobian_do = 0;

  iy_main_agendaExecute(ws,
                        iy,
                        iy_aux,
                        ppath,
                        diy_dx,
                        geo_pos,
                        iy_agenda_call1,
                        iy_transmittance,
                        iy_aux_vars,
                        iy_id,
                        iy_unit,
                        cloudbox_on,
                        jacobian_do,
                        f_grid,
                        atm_field,
                        rte_pos,
                        rte_los,
                        rte_pos2,
                        iy_main_agenda);
}

void get_iy_of_background(Workspace& ws,
                          Matrix& iy,
                          ArrayOfTensor3& diy_dx,
                          const Tensor3& iy_transmittance,
                          const Index& iy_id,
                          const Index& jacobian_do,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Ppath& ppath,
                          const Vector& rte_pos2,
                          const AtmField& atm_field,
                          const Index& cloudbox_on,
                          const Index& stokes_dim,
                          const Vector& f_grid,
                          const String& iy_unit,
                          const Tensor3& surface_props_data,
                          const Agenda& iy_main_agenda,
                          const Agenda& iy_space_agenda,
                          const Agenda& iy_surface_agenda,
                          const Agenda& iy_cloudbox_agenda,
                          const Index& iy_agenda_call1,
                          const Verbosity& verbosity) {
  CREATE_OUT3;

  // Some sizes
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;

  // Set rtp_pos and rtp_los to match the last point in ppath.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector rtp_pos, rtp_los;
  rtp_pos.resize(3);
  rtp_pos = ppath.pos(np - 1, Range(0, 3));
  rtp_los.resize(ppath.los.ncols());
  rtp_los = ppath.los(np - 1, joker);

  out3 << "Radiative background: " << ppath.background << "\n";

  // Handle the different background cases
  //
  String agenda_name;
  //
  switch (ppath_what_background(ppath)) {
    case 1:  //--- Space ----------------------------------------------------
    {
      agenda_name = "iy_space_agenda";
      chk_not_empty(agenda_name, iy_space_agenda);
      iy_space_agendaExecute(ws, iy, f_grid, rtp_pos, rtp_los, iy_space_agenda);
    } break;

    case 2:  //--- The surface -----------------------------------------------
    {
      agenda_name = "iy_surface_agenda";
      chk_not_empty(agenda_name, iy_surface_agenda);
      //
      const Index los_id = iy_id % (Index)1000;
      Index iy_id_new = iy_id + (Index)9 * los_id;
      //
      // Surface jacobian stuff:
      ArrayOfString dsurface_names(0);
      if (jacobian_do && iy_agenda_call1) {
        for (Index i = 0; i < jacobian_quantities.nelem(); i++) {
          if (jacobian_quantities[i] == Jacobian::Special::SurfaceString) {
            dsurface_names.push_back(jacobian_quantities[i].Subtag());
          }
        }
      }
      ArrayOfTensor4 dsurface_rmatrix_dx(dsurface_names.nelem());
      ArrayOfMatrix dsurface_emission_dx(dsurface_names.nelem());
      //
      iy_surface_agendaExecute(ws,
                               iy,
                               diy_dx,
                               dsurface_rmatrix_dx,
                               dsurface_emission_dx,
                               iy_unit,
                               iy_transmittance,
                               iy_id_new,
                               cloudbox_on,
                               jacobian_do,
                               iy_main_agenda,
                               f_grid,
                               atm_field,
                               rtp_pos,
                               rtp_los,
                               rte_pos2,
                               surface_props_data,
                               dsurface_names,
                               iy_surface_agenda);
    } break;

    case 3:  //--- Cloudbox boundary or interior ------------------------------
    case 4: {
      agenda_name = "iy_cloudbox_agenda";
      chk_not_empty(agenda_name, iy_cloudbox_agenda);
      iy_cloudbox_agendaExecute(
          ws, iy, f_grid, rtp_pos, rtp_los, iy_cloudbox_agenda);
    } break;

    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      ARTS_ASSERT(false);
  }

  ARTS_USER_ERROR_IF (iy.ncols() != stokes_dim || iy.nrows() != nf,
      "The size of *iy* returned from *", agenda_name, "* is\n"
      "not correct:\n"
      "  expected size = [", nf, ",", stokes_dim, "]\n"
      "  size of iy    = [", iy.nrows(), ",", iy.ncols(), "]\n")
}

void get_ppath_cloudvars(ArrayOfIndex& clear2cloudy,
                         Matrix& ppath_pnd,
                         ArrayOfMatrix& ppath_dpnd_dx,
                         const Ppath& ppath,
                         const Index& atmosphere_dim,
                         const ArrayOfIndex& cloudbox_limits,
                         const Tensor4& pnd_field,
                         const ArrayOfTensor4& dpnd_field_dx) {
  const Index np = ppath.np;

  // Pnd along the ppath
  ppath_pnd.resize(pnd_field.nbooks(), np);
  ppath_pnd = 0;
  ppath_dpnd_dx.resize(dpnd_field_dx.nelem());

  bool any_dpnd = false;
  for (Index iq = 0; iq < dpnd_field_dx.nelem(); iq++) {
    if (dpnd_field_dx[iq].empty()) {
      ppath_dpnd_dx[iq].resize(0, 0);
    } else {
      any_dpnd = true;
      ppath_dpnd_dx[iq].resize(pnd_field.nbooks(), np);
    }
  }

  // A variable that can map from ppath to particle containers.
  // If outside cloudbox or all (d)pnd=0, this variable holds -1.
  clear2cloudy.resize(np);

  // Determine ppath_pnd and ppath_dpnd_dx
  Index nin = 0;
  for (Index ip = 0; ip < np; ip++)  // PPath point
  {
    Matrix itw(1, Index(pow(2.0, Numeric(atmosphere_dim))));

    ArrayOfGridPos gpc_p(1), gpc_lat(1), gpc_lon(1);
    GridPos gp_lat, gp_lon;
    if (atmosphere_dim >= 2) {
      gridpos_copy(gp_lat, ppath.gp_lat[ip]);
    }
    if (atmosphere_dim == 3) {
      gridpos_copy(gp_lon, ppath.gp_lon[ip]);
    }

    if (is_gp_inside_cloudbox(ppath.gp_p[ip],
                              gp_lat,
                              gp_lon,
                              cloudbox_limits,
                              true,
                              atmosphere_dim)) {
      interp_cloudfield_gp2itw(itw(0, joker),
                               gpc_p[0],
                               gpc_lat[0],
                               gpc_lon[0],
                               ppath.gp_p[ip],
                               gp_lat,
                               gp_lon,
                               atmosphere_dim,
                               cloudbox_limits);
      for (Index i = 0; i < pnd_field.nbooks(); i++) {
        interp_atmfield_by_itw(ExhaustiveVectorView{ppath_pnd(i, ip)},
                               atmosphere_dim,
                               pnd_field(i, joker, joker, joker),
                               gpc_p,
                               gpc_lat,
                               gpc_lon,
                               itw);
      }
      bool any_ppath_dpnd = false;
      if (any_dpnd) {
        for (Index iq = 0; iq < dpnd_field_dx.nelem();
             iq++)  // Jacobian parameter
        {
          if (!dpnd_field_dx[iq].empty()) {
            for (Index i = 0; i < pnd_field.nbooks();
                 i++)  // Scattering element
            {
              interp_atmfield_by_itw(ExhaustiveVectorView{ppath_dpnd_dx[iq](i, ip)},
                                     atmosphere_dim,
                                     dpnd_field_dx[iq](i, joker, joker, joker),
                                     gpc_p,
                                     gpc_lat,
                                     gpc_lon,
                                     itw);
            }
            if (max(ppath_dpnd_dx[iq](joker, ip)) > 0. ||
                min(ppath_dpnd_dx[iq](joker, ip)) < 0.)
              any_ppath_dpnd = true;
          }
        }
      }
      if (max(ppath_pnd(joker, ip)) > 0. || min(ppath_pnd(joker, ip)) < 0. ||
          any_ppath_dpnd) {
        clear2cloudy[ip] = nin;
        nin++;
      } else {
        clear2cloudy[ip] = -1;
      }
    } else {
      clear2cloudy[ip] = -1;
    }
  }
}

void get_ppath_f(Matrix& ppath_f,
                 const Ppath& ppath,
                 const ConstVectorView& f_grid,
                 const Index& atmosphere_dim,
                 const Numeric& rte_alonglos_v,
                 const ConstMatrixView& ppath_wind) {
  // Sizes
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;

  ppath_f.resize(nf, np);

  // Doppler relevant velocity
  //
  for (Index ip = 0; ip < np; ip++) {
    // Start by adding rte_alonglos_v (most likely sensor effects)
    Numeric v_doppler = rte_alonglos_v;

    // Include wind
    if (ppath_wind(1, ip) != 0 || ppath_wind(0, ip) != 0 ||
        ppath_wind(2, ip) != 0) {
      // The dot product below is valid for the photon direction. Winds
      // along this direction gives a positive contribution.
      v_doppler += dotprod_with_los(ppath.los(ip, joker),
                                    ppath_wind(0, ip),
                                    ppath_wind(1, ip),
                                    ppath_wind(2, ip),
                                    atmosphere_dim);
    }

    // Determine frequency grid
    if (v_doppler == 0) {
      ppath_f(joker, ip) = f_grid;
    } else {
      // Positive v_doppler means that sensor measures lower rest
      // frequencies
      const Numeric a = 1 - v_doppler / SPEED_OF_LIGHT;
      for (Index iv = 0; iv < nf; iv++) {
        ppath_f(iv, ip) = a * f_grid[iv];
      }
    }
  }
}

Range get_rowindex_for_mblock(const Sparse& sensor_response,
                              const Index& mblock_index) {
  const Index n1y = sensor_response.nrows();
  return Range(n1y * mblock_index, n1y);
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

void get_stepwise_clearsky_propmat(
    Workspace &ws, PropagationMatrix &K, StokesVector &S,
    ArrayOfPropagationMatrix &dK_dx, ArrayOfStokesVector &dS_dx,
    const Agenda &propmat_clearsky_agenda,
    const ArrayOfRetrievalQuantity &jacobian_quantities,
    const Vector &ppath_f_grid, const Vector &ppath_line_of_sight,
    const AtmPoint &atm_point, const bool jacobian_do) {
  static const ArrayOfSpeciesTag select_abs_species{};
  static const ArrayOfRetrievalQuantity jacobian_quantities_empty{};

  // Perform the propagation matrix computations
  propmat_clearsky_agendaExecute(
      ws, K, S, dK_dx, dS_dx,
      jacobian_do ? jacobian_quantities : jacobian_quantities_empty,
      select_abs_species, ppath_f_grid, ppath_line_of_sight, atm_point,
      propmat_clearsky_agenda);

  adapt_stepwise_partial_derivatives(dK_dx, dS_dx, jacobian_quantities,
                                     ppath_f_grid, ppath_line_of_sight,
                                     S.allZeroes(), 3, jacobian_do);
}

Vector get_stepwise_f_partials(const ConstVectorView& line_of_sight,
                               const ConstVectorView& f_grid,
                               const Jacobian::Atm wind_type,
                               const Index& atmosphere_dim) {
  // Doppler relevant velocity
  Numeric dv_doppler_dx = 0.0;
  
  Vector deriv(f_grid);
  
  switch (wind_type) {
    case Jacobian::Atm::WindMagnitude:
      dv_doppler_dx = 1.0;
      break;
    case Jacobian::Atm::WindU:
      dv_doppler_dx =
          (dotprod_with_los(line_of_sight, 1, 0, 0, atmosphere_dim));
      break;
    case Jacobian::Atm::WindV:
      dv_doppler_dx =
          (dotprod_with_los(line_of_sight, 0, 1, 0, atmosphere_dim));
      break;
    case Jacobian::Atm::WindW:
      dv_doppler_dx =
          (dotprod_with_los(line_of_sight, 0, 0, 1, atmosphere_dim));
      break;
    default:
      ARTS_ASSERT(false, "Not allowed to call this function without a wind parameter as wind_type");
      break;
  }
  
  deriv *= - dv_doppler_dx / Constant::c;
  return deriv;
}

void get_stepwise_scattersky_propmat(
    StokesVector& ap,
    PropagationMatrix& Kp,
    ArrayOfStokesVector& dap_dx,
    ArrayOfPropagationMatrix& dKp_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstMatrixView& ppath_1p_pnd,  // the ppath_pnd at this ppath point
    const ArrayOfMatrix&
        ppath_dpnd_dx,  // the full ppath_dpnd_dx, ie all ppath points
    const Index ppath_1p_id,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstVectorView& ppath_line_of_sight,
    const ConstVectorView& ppath_temperature,
    const Index& atmosphere_dim,
    const bool& jacobian_do) {
  const Index nf = Kp.NumberOfFrequencies(), stokes_dim = Kp.StokesDimensions();

  //StokesVector da_aux(nf, stokes_dim);
  //PropagationMatrix dK_aux(nf, stokes_dim);

  ArrayOfArrayOfSingleScatteringData scat_data_mono;

  // Direction of outgoing scattered radiation (which is reversed to
  // LOS). Only used for extracting scattering properties.
  Vector dir;
  mirror_los(dir, ppath_line_of_sight, atmosphere_dim);
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
                      stokes_dim,
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

  const Index nf_ssd = abs_vec_bulk.nbooks();  // number of freqs in extracted
                                               // optprops. if 1, we need to
                                               // duplicate the ext/abs output.

  for (Index iv = 0; iv < nf; iv++) {
    if (nf_ssd > 1) {
      ap.SetAtPosition(abs_vec_bulk(iv, 0, 0, joker), iv);
      Kp.SetAtPosition(ext_mat_bulk(iv, 0, 0, joker, joker), iv);
    } else {
      ap.SetAtPosition(abs_vec_bulk(0, 0, 0, joker), iv);
      Kp.SetAtPosition(ext_mat_bulk(0, 0, 0, joker, joker), iv);
    }
  }

  if (jacobian_do)
    FOR_ANALYTICAL_JACOBIANS_DO(
        if (ppath_dpnd_dx[iq].empty()) {
          dap_dx[iq].SetZero();
          dKp_dx[iq].SetZero();
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
          for (Index iv = 0; iv < nf; iv++) {
            if (nf_ssd > 1) {
              dap_dx[iq].SetAtPosition(abs_vec_bulk(iv, 0, 0, joker), iv);
              dKp_dx[iq].SetAtPosition(ext_mat_bulk(iv, 0, 0, joker, joker),
                                       iv);
            } else {
              dap_dx[iq].SetAtPosition(abs_vec_bulk(0, 0, 0, joker), iv);
              dKp_dx[iq].SetAtPosition(ext_mat_bulk(0, 0, 0, joker, joker), iv);
            }
          }
        })
}

void get_stepwise_scattersky_source(
    StokesVector& Sp,
    ArrayOfStokesVector& dSp_dx,
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
    const Index& atmosphere_dim,
    const bool& jacobian_do,
    const Index& t_interp_order) {
  ARTS_USER_ERROR_IF (atmosphere_dim != 1,
                      "This function handles so far only 1D atmospheres.");

  const Index nf = Sp.NumberOfFrequencies();
  const Index stokes_dim = Sp.StokesDimensions();
  const Index ne = ppath_1p_pnd.nelem();
  ARTS_ASSERT(TotalNumberOfElements(scat_data) == ne);
  const Index nza = za_grid.nelem();
  const Index naa = aa_grid.nelem();
  const Index nq = jacobian_do ? jacobian_quantities.nelem() : 0;

  // interpolate incident field to this ppath point (no need to do this
  // separately per scatelem)
  GridPos gp_p;
  gridpos_copy(gp_p, ppath_pressure);
  Vector itw_p(2);
  interpweights(itw_p, gp_p);
  Tensor3 inc_field(nf, nza, stokes_dim, 0.);
  for (Index iv = 0; iv < nf; iv++) {
    for (Index iza = 0; iza < nza; iza++) {
      for (Index i = 0; i < stokes_dim; i++) {
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
  Tensor6 pha_mat_1se(nf_ssd, 1, 1, nza * naa, stokes_dim, stokes_dim);
  Vector t_ok(1);
  Index ptype;
  Tensor3 scat_source_1se(ne, nf, stokes_dim, 0.);

  Index ise_flat = 0;
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
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
          Tensor3 product_fields(nza, naa, stokes_dim, 0.);

          ia = 0;
          for (Index iza = 0; iza < nza; iza++) {
            for (Index iaa = 0; iaa < naa; iaa++) {
              for (Index i = 0; i < stokes_dim; i++) {
                for (Index j = 0; j < stokes_dim; j++) {
                  product_fields(iza, iaa, i) +=
                      pha_mat_1se(this_iv, 0, 0, ia, i, j) *
                      inc_field(iv, iza, j);
                }
              }
              ia++;
            }
          }

          for (Index i = 0; i < stokes_dim; i++) {
            scat_source_1se(ise_flat, iv, i) = AngIntegrate_trapezoid(
                product_fields(joker, joker, i), za_grid, aa_grid);
          }
        }  // for iv
      }    // if val_pnd

      ise_flat++;

    }  // for i_se
  }    // for i_ss

  for (Index iv = 0; iv < nf; iv++) {
    Vector scat_source(stokes_dim, 0.);
    for (ise_flat = 0; ise_flat < ne; ise_flat++) {
      for (Index i = 0; i < stokes_dim; i++) {
        scat_source[i] +=
            scat_source_1se(ise_flat, iv, i) * ppath_1p_pnd[ise_flat];
      }
    }

    Sp.SetAtPosition(scat_source, iv);

    if (jacobian_do) {
      FOR_ANALYTICAL_JACOBIANS_DO(
          if (ppath_dpnd_dx[iq].empty()) { dSp_dx[iq].SetZero(); } else {
            scat_source = 0.;
            for (ise_flat = 0; ise_flat < ne; ise_flat++) {
              for (Index i = 0; i < stokes_dim; i++) {
                scat_source[i] += scat_source_1se(ise_flat, iv, i) *
                                  ppath_dpnd_dx[iq](ise_flat, ppath_1p_id);
                dSp_dx[iq].SetAtPosition(scat_source, iv);
              }
            }
          })
    }
  }  // for iv
}

void iyb_calc_body(bool& failed,
                   String& fail_msg,
                   ArrayOfArrayOfMatrix& iy_aux_array,
                   Workspace& ws,
                   Ppath& ppath,
                   Vector& iyb,
                   ArrayOfMatrix& diyb_dx,
                   Vector& geo_pos,
                   const Index& mblock_index,
                   const AtmField& atm_field,
                   const Index& cloudbox_on,
                   const Index& stokes_dim,
                   const Matrix& sensor_pos,
                   const Matrix& sensor_los,
                   const Matrix& transmitter_pos,
                   const Matrix& mblock_dlos,
                   const String& iy_unit,
                   const Agenda& iy_main_agenda,
                   const Index& j_analytical_do,
                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                   const ArrayOfArrayOfIndex& jacobian_indices,
                   const Vector& f_grid,
                   const ArrayOfString& iy_aux_vars,
                   const Index& ilos,
                   const Index& nf) {
constexpr Index atmosphere_dim = 3;

  // The try block here is necessary to correctly handle
  // exceptions inside the parallel region.
  try {
    //--- LOS of interest
    //
    Vector los(sensor_los.ncols());
    //
    los = sensor_los(mblock_index, joker);
    if (mblock_dlos.ncols() == 1) {
      los[0] += mblock_dlos(ilos, 0);
      adjust_los(los, atmosphere_dim);
    } else {
      add_za_aa(los[0],
                los[1],
                los[0],
                los[1],
                mblock_dlos(ilos, 0),
                mblock_dlos(ilos, 1));
    }

    //--- rtp_pos 1 and 2
    //
    Vector rtp_pos, rtp_pos2(0);
    //
    rtp_pos = sensor_pos(mblock_index, joker);
    if (!transmitter_pos.empty()) {
      rtp_pos2 = transmitter_pos(mblock_index, joker);
    }

    // Calculate iy and associated variables
    //
    Matrix iy;
    ArrayOfTensor3 diy_dx;
    Tensor3 iy_transmittance(0, 0, 0);
    const Index iy_agenda_call1 = 1;
    const Index iy_id =
        (Index)1e6 * (mblock_index + 1) + (Index)1e3 * (ilos + 1);
    //
    iy_main_agendaExecute(ws,
                          iy,
                          iy_aux_array[ilos],
                          ppath,
                          diy_dx,
                          geo_pos,
                          iy_agenda_call1,
                          iy_transmittance,
                          iy_aux_vars,
                          iy_id,
                          iy_unit,
                          cloudbox_on,
                          j_analytical_do,
                          f_grid,
                          atm_field,
                          rtp_pos,
                          los,
                          rtp_pos2,
                          iy_main_agenda);

    // Start row in iyb etc. for present LOS
    //
    const Index row0 = ilos * nf * stokes_dim;

    // Jacobian part
    //
    if (j_analytical_do) {
      FOR_ANALYTICAL_JACOBIANS_DO2(
          for (Index ip = 0;
               ip < jacobian_indices[iq][1] - jacobian_indices[iq][0] + 1;
               ip++) {
            for (Index is = 0; is < stokes_dim; is++) {
              diyb_dx[iq](Range(row0 + is, nf, stokes_dim), ip) =
                  diy_dx[iq](ip, joker, is);
            }
          })
    }

    // iy : copy to iyb
    for (Index is = 0; is < stokes_dim; is++) {
      iyb[Range(row0 + is, nf, stokes_dim)] = iy(joker, is);
    }

  }  // End try

  catch (const std::exception& e) {
#pragma omp critical(iyb_calc_fail)
    {
      fail_msg = e.what();
      failed = true;
    }
  }
}

void iyb_calc(Workspace& ws,
              Vector& iyb,
              ArrayOfVector& iyb_aux,
              ArrayOfMatrix& diyb_dx,
              Matrix& geo_pos_matrix,
              const Index& mblock_index,
              const AtmField& atm_field,
              const Index& cloudbox_on,
              const Index& stokes_dim,
              const Vector& f_grid,
              const Matrix& sensor_pos,
              const Matrix& sensor_los,
              const Matrix& transmitter_pos,
              const Matrix& mblock_dlos,
              const String& iy_unit,
              const Agenda& iy_main_agenda,
              const Index& j_analytical_do,
              const ArrayOfRetrievalQuantity& jacobian_quantities,
              const ArrayOfArrayOfIndex& jacobian_indices,
              const ArrayOfString& iy_aux_vars,
              const Verbosity& verbosity) {
  CREATE_OUT3;

  // Sizes
  const Index nf = f_grid.nelem();
  const Index nlos = mblock_dlos.nrows();
  const Index niyb = nf * nlos * stokes_dim;
  // Set up size of containers for data of 1 measurement block.
  // (can not be made below due to parallalisation)
  iyb.resize(niyb);
  //
  if (j_analytical_do) {
    diyb_dx.resize(jacobian_indices.nelem());
    FOR_ANALYTICAL_JACOBIANS_DO2(diyb_dx[iq].resize(
        niyb, jacobian_indices[iq][1] - jacobian_indices[iq][0] + 1);)
  } else {
    diyb_dx.resize(0);
  }
  // Assume empty geo_pos.
  geo_pos_matrix.resize(nlos, 5);
  geo_pos_matrix = NAN;

  // For iy_aux we don't know the number of quantities, and we have to store
  // all outout
  ArrayOfArrayOfMatrix iy_aux_array(nlos);

  String fail_msg;
  bool failed = false;
  if (nlos >= arts_omp_get_max_threads() || nlos * 10 >= nf) {
    out3 << "  Parallelizing los loop (" << nlos << " iterations, " << nf
         << " frequencies)\n";

    WorkspaceOmpParallelCopyGuard wss{ws};

    // Start of actual calculations
#pragma omp parallel for if (!arts_omp_in_parallel()) firstprivate(wss)
    for (Index ilos = 0; ilos < nlos; ilos++) {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      Ppath ppath;
      Vector geo_pos;
      iyb_calc_body(failed,
                    fail_msg,
                    iy_aux_array,
                    wss,
                    ppath,
                    iyb,
                    diyb_dx,
                    geo_pos,
                    mblock_index,
                    atm_field,
                    cloudbox_on,
                    stokes_dim,
                    sensor_pos,
                    sensor_los,
                    transmitter_pos,
                    mblock_dlos,
                    iy_unit,
                    iy_main_agenda,
                    j_analytical_do,
                    jacobian_quantities,
                    jacobian_indices,
                    f_grid,
                    iy_aux_vars,
                    ilos,
                    nf);

      if (geo_pos.nelem()) geo_pos_matrix(ilos, joker) = geo_pos;

      // Skip remaining iterations if an error occurred
      if (failed) continue;
    }
  } else {
    out3 << "  Not parallelizing los loop (" << nlos << " iterations, " << nf
         << " frequencies)\n";

    for (Index ilos = 0; ilos < nlos; ilos++) {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      Ppath ppath;
      Vector geo_pos;
      iyb_calc_body(failed,
                    fail_msg,
                    iy_aux_array,
                    ws,
                    ppath,
                    iyb,
                    diyb_dx,
                    geo_pos,
                    mblock_index,
                    atm_field,
                    cloudbox_on,
                    stokes_dim,
                    sensor_pos,
                    sensor_los,
                    transmitter_pos,
                    mblock_dlos,
                    iy_unit,
                    iy_main_agenda,
                    j_analytical_do,
                    jacobian_quantities,
                    jacobian_indices,
                    f_grid,
                    iy_aux_vars,
                    ilos,
                    nf);

      if (geo_pos.nelem()) geo_pos_matrix(ilos, joker) = geo_pos;

      // Skip remaining iterations if an error occurred
      if (failed) continue;
    }
  }

  ARTS_USER_ERROR_IF (failed,
                      "Run-time error in function: iyb_calc\n", fail_msg);

  // Compile iyb_aux
  //
  const Index nq = iy_aux_array[0].nelem();
  iyb_aux.resize(nq);
  //
  for (Index q = 0; q < nq; q++) {
    iyb_aux[q].resize(niyb);
    //
    for (Index ilos = 0; ilos < nlos; ilos++) {
      const Index row0 = ilos * nf * stokes_dim;
      for (Index iv = 0; iv < nf; iv++) {
        const Index row1 = row0 + iv * stokes_dim;
        const Index i1 = min(iv, iy_aux_array[ilos][q].nrows() - 1);
        for (Index is = 0; is < stokes_dim; is++) {
          Index i2 = min(is, iy_aux_array[ilos][q].ncols() - 1);
          iyb_aux[q][row1 + is] = iy_aux_array[ilos][q](i1, i2);
        }
      }
    }
  }
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
                const ConstVectorView& los,
                const Index& atmosphere_dim) {
  los_mirrored.resize(2);
  //
  if (atmosphere_dim == 1) {
    los_mirrored[0] = 180 - los[0];
    los_mirrored[1] = 180;
  } else if (atmosphere_dim == 2) {
    los_mirrored[0] = 180 - fabs(los[0]);
    if (los[0] >= 0) {
      los_mirrored[1] = 180;
    } else {
      los_mirrored[1] = 0;
    }
  } else if (atmosphere_dim == 3) {
    los_mirrored[0] = 180 - los[0];
    los_mirrored[1] = los[1] + 180;
    if (los_mirrored[1] > 180) {
      los_mirrored[1] -= 360;
    }
  }
}

void muellersparse_rotation(Sparse& H,
                            const Index& stokes_dim,
                            const Numeric& rotangle) {
  ARTS_ASSERT(stokes_dim > 1);
  ARTS_ASSERT(stokes_dim <= 4);
  ARTS_ASSERT(H.nrows() == stokes_dim);
  ARTS_ASSERT(H.ncols() == stokes_dim);
  ARTS_ASSERT(H(0, 1) == 0);
  ARTS_ASSERT(H(1, 0) == 0);
  //
  H.rw(0, 0) = 1;
  const Numeric a = Conversion::cosd(2 * rotangle);
  H.rw(1, 1) = a;
  if (stokes_dim > 2) {
    ARTS_ASSERT(H(2, 0) == 0);
    ARTS_ASSERT(H(0, 2) == 0);

    const Numeric b = Conversion::sind(2 * rotangle);
    H.rw(1, 2) = b;
    H.rw(2, 1) = -b;
    H.rw(2, 2) = a;
    if (stokes_dim > 3) {
      // More values should be checked, but to save time we just ARTS_ASSERT one
      ARTS_ASSERT(H(2, 3) == 0);
      H.rw(3, 3) = 1;
    }
  }
}

void mueller_modif2stokes(Matrix& Cs,
                          const Index& stokes_dim) {
  ARTS_ASSERT(stokes_dim >= 1);
  ARTS_ASSERT(stokes_dim <= 4);
  //
  Cs.resize(stokes_dim, stokes_dim);
  Cs(0,0) = 1;
  if (stokes_dim > 1 ) {
    Cs(0,1) = Cs(1,0) = 1;
    Cs(1,1) = -1;
    if (stokes_dim > 2 ) {
      Cs(0,2) = Cs(1,2) = Cs(2,0) = Cs(2,1) = 0;
      Cs(2,2) = 1;
      if (stokes_dim > 3 ) {
        Cs(0,3) = Cs(1,3) = Cs(2,3) = Cs(3,0) = Cs(3,1) = Cs(3,2) = 0;
        Cs(3,3) = 1;       
      }
    }
  }
}

void mueller_rotation(Matrix& L,
                      const Index& stokes_dim,
                      const Numeric& rotangle) {
  ARTS_ASSERT(stokes_dim >= 1);
  ARTS_ASSERT(stokes_dim <= 4);
  //
  L.resize(stokes_dim, stokes_dim);
  L(0, 0) = 1;
  if (stokes_dim > 1 ) {
    const Numeric alpha = 2 * Conversion::deg2rad(1) * rotangle;
    const Numeric c2 = cos(alpha);
    L(0,1) = L(1,0) = 0;
    L(1,1) = c2;
    if (stokes_dim > 2 ) {
      const Numeric s2 = sin(alpha);
      L(0,2) = L(2,0) = 0;
      L(1,2) = s2;
      L(2,1) = -s2;      
      L(2,2) = c2;
      if (stokes_dim > 3 ) {
        L(0,3) = L(1,3) = L(2,3) = L(3,0) = L(3,1) = L(3,2) = 0;
        L(3,3) = 1;       
      }
    }
  }
}

void mueller_stokes2modif(Matrix& Cm,
                          const Index& stokes_dim) {
  ARTS_ASSERT(stokes_dim >= 1);
  ARTS_ASSERT(stokes_dim <= 4);
  //
  Cm.resize(stokes_dim, stokes_dim);
  Cm(0,0) = 0.5;
  if (stokes_dim > 1 ) {
    Cm(0,1) = Cm(1,0) = 0.5;
    Cm(1,1) = -0.5;
    if (stokes_dim > 2 ) {
      Cm(0,2) = Cm(1,2) = Cm(2,0) = Cm(2,1) = 0;
      Cm(2,2) = 1;
      if (stokes_dim > 3 ) {
        Cm(0,3) = Cm(1,3) = Cm(2,3) = Cm(3,0) = Cm(3,1) = Cm(3,2) = 0;
        Cm(3,3) = 1;       
      }
    }
  }
}

void pos2true_latlon(Numeric& lat,
                     Numeric& lon,
                     const Index& atmosphere_dim,
                     const ConstVectorView& lat_grid,
                     const ConstVectorView& lat_true,
                     const ConstVectorView& lon_true,
                     const ConstVectorView& pos) {
  ARTS_ASSERT(pos.nelem() == atmosphere_dim);

  if (atmosphere_dim == 1) {
    ARTS_ASSERT(lat_true.nelem() == 1);
    ARTS_ASSERT(lon_true.nelem() == 1);
    //
    lat = lat_true[0];
    lon = lon_true[0];
  }

  else if (atmosphere_dim == 2) {
    ARTS_ASSERT(lat_true.nelem() == lat_grid.nelem());
    ARTS_ASSERT(lon_true.nelem() == lat_grid.nelem());
    GridPos gp;
    Vector itw(2);
    gridpos(gp, lat_grid, pos[1]);
    interpweights(itw, gp);
    lat = interp(itw, lat_true, gp);
    lon = interp(itw, lon_true, gp);
  }

  else {
    lat = pos[1];
    lon = pos[2];
  }
}

void rtmethods_jacobian_finalisation(
    Workspace& ws,
    ArrayOfTensor3& diy_dx,
    ArrayOfTensor3& diy_dpath,
    const Index& ns,
    const Index& nf,
    const Index& np,
    const Ppath& ppath,
    const ArrayOfAtmPoint& ppvar_atm,
    const Index& iy_agenda_call1,
    const Tensor3& iy_transmittance,
    const Agenda& water_p_eq_agenda,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfIndex& jac_species_i) {
  // Weight with iy_transmittance
  if (!iy_agenda_call1) {
    Matrix X, Y;
    //
    FOR_ANALYTICAL_JACOBIANS_DO(
        Y.resize(ns, diy_dpath[iq].npages());
        for (Index iv = 0; iv < nf; iv++) {
          X = transpose(diy_dpath[iq](joker, iv, joker));
          mult(Y, iy_transmittance(iv, joker, joker), X);
          diy_dpath[iq](joker, iv, joker) = transpose(Y);
        })
  }

  // Handle abs species retrieval units, both internally and impact on T-jacobian
  //
  Tensor3 water_p_eq(0, 0, 0);
  //
  // Conversion for abs species itself
  for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
    // Let x be VMR, and z the selected retrieval unit.
    // We have then that diy/dz = diy/dx * dx/dz
    //
    if (not(jacobian_quantities[iq] == Jacobian::Type::Sensor) and
        not(jacobian_quantities[iq] == Jacobian::Special::SurfaceString) and
        jac_species_i[iq] >= 0) {
      if (jacobian_quantities[iq].Mode() == "vmr") {
      }

      else if (jacobian_quantities[iq].Mode() == "rel") {
        // Here x = vmr*z
        for (Index ip = 0; ip < np; ip++) {
          diy_dpath[iq](ip, joker, joker) *= ppvar_atm[ip][abs_species[jac_species_i[iq]]];
        }
      }

      else if (jacobian_quantities[iq].Mode() == "nd") {
        // Here x = z/nd_tot
        for (Index ip = 0; ip < np; ip++) {
          diy_dpath[iq](ip, joker, joker) /=
              number_density(ppvar_atm[ip].pressure, ppvar_atm[ip].temperature);
        }
      }

      else if (jacobian_quantities[iq].Mode() == "rh") {
        // Here x = (p_sat/p) * z
        AtmField atm_field;
        Tensor3 t_data(ppvar_atm.nelem(), 1, 1);
        for (Index ip = 0; ip < np; ip++) {
          t_data(ip, 0, 0) = ppvar_atm[ip].temperature;
        }
        ARTS_ASSERT(false)
        water_p_eq_agendaExecute(ws, water_p_eq, atm_field, water_p_eq_agenda);
        for (Index ip = 0; ip < np; ip++) {
          diy_dpath[iq](ip, joker, joker) *= water_p_eq(ip, 0, 0) / ppvar_atm[ip].pressure;
        }
      }

      else if (jacobian_quantities[iq].Mode() == "q") {
        // Here we use the approximation of x = z/0.622
        diy_dpath[iq](joker, joker, joker) /= 0.622;
      }

      else {
        ARTS_ASSERT(0);
      }
    }
  }

  // Correction of temperature Jacobian
  for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
    // Let a be unit for abs species, and iy = f(T,a(T))
    // We have then that diy/dT = df/dT + df/da*da/dT
    // diy_dpath holds already df/dT. Remains is to add
    // df/da*da/dT for which abs species having da/dT != 0
    // This is only true for "nd" and "rh"
    //
    if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
      // Loop abs species, again
      for (Index ia = 0; ia < jacobian_quantities.nelem(); ia++) {
        if (jac_species_i[ia] >= 0) {
          if (jacobian_quantities[ia].Mode() == "nd") {
            for (Index ip = 0; ip < np; ip++) {
              Matrix ddterm{diy_dpath[ia](ip, joker, joker)};
              ddterm *= ppvar_atm[ip][abs_species[jac_species_i[ia]]] *
                        (number_density(ppvar_atm[ip].pressure, ppvar_atm[ip].temperature + 1) -
                         number_density(ppvar_atm[ip].pressure, ppvar_atm[ip].pressure));
              diy_dpath[iq](ip, joker, joker) += ddterm;
            }
          } else if (jacobian_quantities[ia].Mode() == "rh") {
            Tensor3 t_data(ppvar_atm.nelem(), 1, 1);
            for (Index ip = 0; ip < np; ip++) {
              t_data(ip, 0, 0) = ppvar_atm[ip].temperature;
            }
            // Calculate water sat. pressure if not already done
            if (water_p_eq.npages() == 0) {
              AtmField atm_field;
              ARTS_ASSERT(false)
              water_p_eq_agendaExecute(
                  ws, water_p_eq, atm_field, water_p_eq_agenda);
            }
            // Sat.pressure for +1K
            Tensor3 water_p_eq1K;
            t_data(joker, 0, 0) += 1;
              AtmField atm_field;
              ARTS_ASSERT(false)
            water_p_eq_agendaExecute(
                ws, water_p_eq1K, atm_field, water_p_eq_agenda);

            for (Index ip = 0; ip < np; ip++) {
              const Numeric p_eq = water_p_eq(ip, 0, 0);
              const Numeric p_eq1K = water_p_eq1K(ip, 0, 0);
              Matrix ddterm{diy_dpath[ia](ip, joker, joker)};
              ddterm *= ppvar_atm[ip][abs_species[jac_species_i[ia]]] *
                        (ppvar_atm[ip].pressure / pow(p_eq, 2.0)) * (p_eq1K - p_eq);
              diy_dpath[iq](ip, joker, joker) += ddterm;
            }
          }
        }
      }
    }
  }

  // Map to retrieval grids
  FOR_ANALYTICAL_JACOBIANS_DO(diy_from_path_to_rgrids(diy_dx[iq],
                                                      jacobian_quantities[iq],
                                                      diy_dpath[iq],
                                                      ppath);)
}

void rtmethods_unit_conversion(
    Matrix& iy,
    ArrayOfTensor3& diy_dx,
    Tensor3& ppvar_iy,
    const Vector& f_grid,
    const Ppath& ppath,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Index& j_analytical_do,
    const String& iy_unit) {
  // Determine refractive index to use for the n2 radiance law
  Numeric n = 1.0;  // First guess is that sensor is in space
  //
  if (ppath.end_lstep == 0)  // If true, sensor inside the atmosphere
  {
    n = ppath.nreal.back();
  }

  const Index ns = iy.ncols();

  // Polarisation index variable
  ArrayOfIndex i_pol(ns);
  for (Index is = 0; is < ns; is++) {
    i_pol[is] = is + 1;
  }

  // Jacobian part (must be converted to Tb before iy for PlanckBT)
  //
  if (j_analytical_do) {
    FOR_ANALYTICAL_JACOBIANS_DO2(
        apply_iy_unit2(diy_dx[iq], iy, iy_unit, f_grid, n, i_pol);)
  }

  // iy
  apply_iy_unit(iy, iy_unit, f_grid, n, i_pol);

  // ppvar_iy
  for (Index ip = 0; ip < ppath.np; ip++) {
    apply_iy_unit(
        ppvar_iy(joker, joker, ip), iy_unit, f_grid, ppath.nreal[ip], i_pol);
  }
}

void yCalc_mblock_loop_body(bool& failed,
                            String& fail_msg,
                            ArrayOfArrayOfVector& iyb_aux_array,
                            Workspace& ws,
                            Vector& y,
                            Vector& y_f,
                            ArrayOfIndex& y_pol,
                            Matrix& y_pos,
                            Matrix& y_los,
                            Matrix& y_geo,
                            Matrix& jacobian,
                            const AtmField& atm_field,
                            const Index& cloudbox_on,
                            const Index& stokes_dim,
                            const Vector& f_grid,
                            const Matrix& sensor_pos,
                            const Matrix& sensor_los,
                            const Matrix& transmitter_pos,
                            const Matrix& mblock_dlos,
                            const Sparse& sensor_response,
                            const Vector& sensor_response_f,
                            const ArrayOfIndex& sensor_response_pol,
                            const Matrix& sensor_response_dlos,
                            const String& iy_unit,
                            const Agenda& iy_main_agenda,
                            const Agenda& jacobian_agenda,
                            const Index& jacobian_do,
                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                            const ArrayOfArrayOfIndex& jacobian_indices,
                            const ArrayOfString& iy_aux_vars,
                            const Verbosity& verbosity,
                            const Index& mblock_index,
                            const Index& n1y,
                            const Index& j_analytical_do) {
  try {
    // Calculate monochromatic pencil beam data for 1 measurement block
    //
    Vector iyb, iyb_error, yb(n1y);
    ArrayOfMatrix diyb_dx;
    Matrix geo_pos_matrix;
    //
    iyb_calc(ws,
             iyb,
             iyb_aux_array[mblock_index],
             diyb_dx,
             geo_pos_matrix,
             mblock_index,
             atm_field,
             cloudbox_on,
             stokes_dim,
             f_grid,
             sensor_pos,
             sensor_los,
             transmitter_pos,
             mblock_dlos,
             iy_unit,
             iy_main_agenda,
             j_analytical_do,
             jacobian_quantities,
             jacobian_indices,
             iy_aux_vars,
             verbosity);

    // Apply sensor response matrix on iyb, and put into y
    //
    const Range rowind = get_rowindex_for_mblock(sensor_response, mblock_index);
    const Index row0 = rowind.offset;
    //
    mult(yb, sensor_response, iyb);
    //
    y[rowind] = yb;  // *yb* also used below, as input to jacobian_agenda

    // Fill information variables. And search for NaNs in *y*.
    //
    for (Index i = 0; i < n1y; i++) {
      const Index ii = row0 + i;
      ARTS_USER_ERROR_IF (std::isnan(y[ii]),
                          "One or several NaNs found in *y*.");
      y_f[ii] = sensor_response_f[i];
      y_pol[ii] = sensor_response_pol[i];
      y_pos(ii, joker) = sensor_pos(mblock_index, joker);
      y_los(ii, joker) = sensor_los(mblock_index, joker);
      y_los(ii, 0) += sensor_response_dlos(i, 0);
      if (sensor_response_dlos.ncols() > 1) {
        y_los(ii, 1) += sensor_response_dlos(i, 1);
      }
    }

    // Apply sensor response matrix on diyb_dx, and put into jacobian
    // (that is, analytical jacobian part)
    //
    if (j_analytical_do) {
      FOR_ANALYTICAL_JACOBIANS_DO2(
          mult(jacobian(rowind,
                        Range(jacobian_indices[iq][0],
                              jacobian_indices[iq][1] -
                                  jacobian_indices[iq][0] + 1)),
               sensor_response,
               diyb_dx[iq]);)
    }

    // Calculate remaining parts of *jacobian*
    //
    if (jacobian_do) {
      jacobian_agendaExecute(
          ws, jacobian, mblock_index, iyb, yb, jacobian_agenda);
    }

    // Handle geo-positioning
    if (!std::isnan(geo_pos_matrix(0, 0)))  // No data are flagged as NaN
    {
      // We set geo_pos based on the max value in sensor_response
      const Index nfs = f_grid.nelem() * stokes_dim;
      for (Index i = 0; i < n1y; i++) {
        Index jmax = -1;
        Numeric rmax = -99e99;
        for (Index j = 0; j < sensor_response.ncols(); j++) {
          if (sensor_response(i, j) > rmax) {
            rmax = sensor_response(i, j);
            jmax = j;
          }
        }
        const auto jhit = Index(floor(jmax / nfs));
        y_geo(row0 + i, joker) = geo_pos_matrix(jhit, joker);
      }
    }
  }

  catch (const std::exception& e) {
#pragma omp critical(yCalc_fail)
    {
      fail_msg = e.what();
      failed = true;
    }
  }
}

void ze_cfac(Vector& fac,
             const Vector& f_grid,
             const Numeric& ze_tref,
             const Numeric& k2) {
  const Index nf = f_grid.nelem();

  ARTS_ASSERT(fac.nelem() == nf);

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
