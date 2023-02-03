/* Copyright (C) 2005-2012 Cory Davis <cory@met.ed.ac.uk>
                            
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
 * @file   mc_antenna.cc
 * @author Cory Davis <cdavis@staffmail.ed.ac.uk>
 * @date   2005-12-01
 *
 * @brief  Monte Carlo Antenna implementation
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "mc_antenna.h"
#include <cfloat>
#include <sstream>
#include "arts_constants.h"
#include "arts_conversions.h"
#include "matpack_math.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

Numeric ran_gaussian(Rng& rng, const Numeric sigma) {
  Numeric x, y, r2;

  do {
    /* choose x,y in uniform square (-1,-1) to (+1,+1) */

    x = -1 + 2 * rng.draw();
    y = -1 + 2 * rng.draw();

    /* see if it is in the unit circle */
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt(-2.0 * log(r2) / r2);
}

Numeric ran_uniform(Rng& rng) {
  Numeric phi;

  phi = 360.0 * rng.draw() - 180.0;

  return phi;
}

void rotmat_enu(MatrixView R_ant2enu, ConstVectorView prop_los)

{
  Numeric cza, sza, caa, saa;

  cza = cos(prop_los[0] * DEG2RAD);
  sza = sin(prop_los[0] * DEG2RAD);
  caa = cos(prop_los[1] * DEG2RAD);
  saa = sin(prop_los[1] * DEG2RAD);

  R_ant2enu(0, 0) = -cza * saa;
  R_ant2enu(0, 1) = caa;
  R_ant2enu(0, 2) = sza * saa;

  R_ant2enu(1, 0) = -cza * caa;
  R_ant2enu(1, 1) = -saa;
  R_ant2enu(1, 2) = sza * caa;

  R_ant2enu(2, 0) = sza;
  R_ant2enu(2, 1) = 0.0;
  R_ant2enu(2, 2) = cza;
}

void rotmat_stokes(MatrixView R_pra,
                   const Index& stokes_dim,
                   const Numeric& f1_dir,
                   const Numeric& f2_dir,
                   ConstMatrixView R_f1,
                   ConstMatrixView R_f2)

{
  const Numeric flip = f1_dir * f2_dir;
  Numeric cos_pra1, sin_pra1, cos_pra2, sin_pra2;

  cos_pra1 = R_f1(joker, 0) * R_f2(joker, 0);
  sin_pra1 = f2_dir * (R_f1(joker, 0) * R_f2(joker, 1));
  sin_pra2 = f1_dir * (R_f1(joker, 1) * R_f2(joker, 0));
  cos_pra2 = f1_dir * f2_dir * (R_f1(joker, 1) * R_f2(joker, 1));

  R_pra = 0.0;
  R_pra(0, 0) = 1.0;
  if (stokes_dim > 1) {
    R_pra(1, 1) = 2 * cos_pra1 * cos_pra1 - 1.0;
    if (stokes_dim > 2) {
      R_pra(1, 2) = flip * 2 * cos_pra1 * sin_pra1;
      R_pra(2, 1) = 2 * cos_pra2 * sin_pra2;
      R_pra(2, 2) = flip * (2 * cos_pra2 * cos_pra2 - 1.0);
      if (stokes_dim > 3) {
        R_pra(3, 3) = flip * 1.0;
      }
    }
  }
}

void MCAntenna::set_pencil_beam() { atype = ANTENNA_TYPE_PENCIL_BEAM; }

void MCAntenna::set_gaussian(const Numeric& za_sigma, const Numeric& aa_sigma) {
  atype = ANTENNA_TYPE_GAUSSIAN;
  sigma_za = za_sigma;
  sigma_aa = aa_sigma;
}

void MCAntenna::set_gaussian_fwhm(const Numeric& za_fwhm,
                                  const Numeric& aa_fwhm) {
  atype = ANTENNA_TYPE_GAUSSIAN;
  sigma_za = za_fwhm / 2.3548;
  sigma_aa = aa_fwhm / 2.3548;
}

void MCAntenna::set_lookup(ConstVectorView za_grid_,
                           ConstVectorView aa_grid_,
                           ConstMatrixView G_lookup_) {
  atype = ANTENNA_TYPE_LOOKUP;
  za_grid = za_grid_;
  aa_grid = aa_grid_;
  G_lookup = G_lookup_;
}

void MCAntenna::return_los(Numeric& wgt,
                           ConstMatrixView R_return,
                           ConstMatrixView R_enu2ant) const {
  Numeric z, term_el, term_az;
  Numeric ant_el, ant_az;
  Vector k_vhk(3);

  switch (atype) {
    case ANTENNA_TYPE_PENCIL_BEAM:
      wgt = 1.0;
      break;

    case ANTENNA_TYPE_GAUSSIAN:

      mult(k_vhk, R_enu2ant, R_return(joker, 2));

      // Assume Gaussian is narrow enough that response is 0 beyond 90 degrees
      // Same assumption is made for drawing samples (draw_los)
      if (k_vhk[2] > 0) {
        ant_el = atan(k_vhk[0] / k_vhk[2]) * RAD2DEG;
        ant_az = atan(k_vhk[1] / k_vhk[2]) * RAD2DEG;
        term_el = ant_el / sigma_za;
        term_az = ant_az / sigma_aa;
        z = term_el * term_el + term_az * term_az;
        wgt = exp(-0.5 * z);
      } else {
        wgt = 0.0;
      }
      break;

    default:
      ARTS_USER_ERROR ("invalid Antenna type.")
  }
}

void MCAntenna::draw_los(VectorView sampled_rte_los,
                         MatrixView R_los,
                         Rng& rng,
                         ConstMatrixView R_ant2enu,
                         ConstVectorView bore_sight_los) const {
  Numeric ant_el, ant_az, ant_r;
  Numeric tel, taz;
  Vector k_vhk(3);

  switch (atype) {
    case ANTENNA_TYPE_PENCIL_BEAM:
      sampled_rte_los = bore_sight_los;
      R_los = R_ant2enu;
      break;

    case ANTENNA_TYPE_GAUSSIAN:

      ant_el = 91;
      ant_az = 91;

      // Assume Gaussian is narrow enough that response is 0 beyond 90 degrees
      // Same assumption is made for radar return samples (return_los)
      while (ant_el >= 90) {
        ant_el = ran_gaussian(rng, sigma_za);
      }
      while (ant_az >= 90) {
        ant_az = ran_gaussian(rng, sigma_aa);
      }

      // Propagation direction
      tel = tan(DEG2RAD * ant_el);
      taz = tan(DEG2RAD * ant_az);
      ant_r = sqrt(1 + tel * tel + taz * taz);
      k_vhk[0] = tel / ant_r;
      k_vhk[1] = taz / ant_r;
      k_vhk[2] = (Numeric)1.0 / ant_r;
      mult(R_los(joker, 2), R_ant2enu, k_vhk);

      sampled_rte_los[0] = acos(R_los(2, 2)) * RAD2DEG;

      // Horizontal polarization basis
      // If drawn los is at zenith or nadir, assume same azimuth as boresight
      if (((Numeric)1.0 - abs(R_los(2, 2))) < DBL_EPSILON) {
        // H is aligned with H of bs, use row not column because tranpose
        R_los(joker, 1) = R_ant2enu(1, joker);
        sampled_rte_los[1] = bore_sight_los[1];
      } else {
        const Vector uhat{0.0, 0.0, 1.0};
        Numeric magh;
        sampled_rte_los[1] = atan2(R_los(0, 2), R_los(1, 2)) * RAD2DEG;
        cross3(R_los(joker, 1), R_los(joker, 2), uhat);
        magh = sqrt(R_los(joker, 1) * R_los(joker, 1));
        R_los(joker, 1) /= magh;
      }

      // Vertical polarization basis
      cross3(R_los(joker, 0), R_los(joker, 1), R_los(joker, 2));

      break;

    default:
      ARTS_USER_ERROR ("invalid Antenna type.")
  }
}

ostream& operator<<(ostream& os, const MCAntenna&) {
  os << "MCAntenna: Output operator not implemented";
  return os;
}
