/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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
   USA.
*/
/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   star.cc
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Functions needed for radiative transfer with direct sources.
*/

#include "star.h"
#include "auto_md.h"
#include "agenda_class.h"
#include "propagationmatrix.h"
#include "geodetic.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions
  ===========================================================================*/

std::ostream& operator<<(std::ostream& os, const Star& star) {
  os << "Star: " << star.description;
  os << " Radius: " << star.radius << "m ";
  os << " Distance: " << star.distance << "m \n";
  os << " Latitude: " << star.latitude << "° \n";
  os << " Longitude: " << star.longitude << "° \n";
  os << " Spectrum [W/m2/Hz]: \n" << star.spectrum ;
  return os;
}

void get_scattered_starsource(Workspace& ws,
                              RadiationVector& scattered_starlight,
                              ArrayOfRadiationVector& dscattered_starlight,
                              const Vector& f_grid,
                              const Numeric& p,
                              const Numeric& T,
                              const Vector& vmr,
                              const Matrix& transmitted_starlight,
                              const Vector& in_los,
                              const Vector& out_los,
                              const Agenda& gas_scattering_agenda) {
  PropagationMatrix K_sca;
  TransmissionMatrix sca_mat;
  Vector sca_fct_dummy;

  // calculate gas scattering properties
  gas_scattering_agendaExecute(ws,
                               K_sca,
                               sca_mat,
                               sca_fct_dummy,
                               f_grid,
                               p,
                               T,
                               vmr,
                               in_los,
                               out_los,
                               0,
                               gas_scattering_agenda);

  //some basic quantities
  Index ns = transmitted_starlight.ncols();
  Index nf = f_grid.nelem();

  //allocate and resize
  Vector temp(ns);
  RadiationVector scattered_starlight_temp(1, ns);

  // Calculate the scattered radiation
  for (Index i_f = 0; i_f < nf; i_f++) {
    scattered_starlight_temp = transmitted_starlight(i_f, joker);
    scattered_starlight_temp.leftMul(sca_mat);

    for (Index j = 0; j < ns; j++) {
      scattered_starlight(i_f, j) =
          scattered_starlight_temp(0, j) * K_sca.Kjj(0, 0)[i_f];
    }
  }
}

void get_star_background(Matrix& iy,
                         const Star& star,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& refellipsoid) {

  //Calculate earth centric radial component of star_pos and rtp_pos.
  const Numeric R_star = star.distance;
  const Numeric R_rte = refell2r(refellipsoid, rtp_pos[1]) + rtp_pos[0];

  //Transform to cartesian coordinate system
  Numeric r_star_x, r_star_y, r_star_z;
  Numeric r_rte_x, r_rte_y, r_rte_z;
  Numeric r_los_x, r_los_y, r_los_z;

  // r_star
  sph2cart(r_star_x, r_star_y, r_star_z, R_star, star.latitude, star.longitude);

  // r_rte, r_los
  poslos2cart(r_rte_x,
              r_rte_y,
              r_rte_z,
              r_los_x,
              r_los_y,
              r_los_z,
              R_rte,
              rtp_pos[1],
              rtp_pos[2],
              rtp_los[0],
              rtp_los[1]);

  //Calculate vector of line of sight and unit vector pointing from
  //ppath point to the star.
  const Numeric r_ps_x = r_star_x - r_rte_x;
  const Numeric r_ps_y = r_star_y - r_rte_y;
  const Numeric r_ps_z = r_star_z - r_rte_z;

  //abs value of r_ps
  const Numeric r_ps =
      sqrt(r_ps_x * r_ps_x + r_ps_y * r_ps_y + r_ps_z * r_ps_z);

  //abs value of r_los
  const Numeric r_glos =
      sqrt(r_los_x * r_los_x + r_los_y * r_los_y + r_los_z * r_los_z);

  //Calculate angle beta between line of sight and the line between ppath point and the star
  //using scalar product
  const Numeric cos_beta =
      (r_ps_x * r_los_x + r_ps_y * r_los_y + r_ps_z * r_los_z) / (r_ps * r_glos);
  const Numeric beta = acos(cos_beta);

  // angular diameter of star
  const Numeric alpha = atan(star.radius / r_ps);

  //Check if we see the star. We see the star if the angle beta is smaller than
  // the angular diameter alpha of the star.
  if (beta <= alpha) {
    //Here we assume that the star radiates isotropically.
    Matrix star_radiance = star.spectrum;
    star_radiance /= PI;

    iy += star_radiance;
  }
}
