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
   USA. */

/*===========================================================================
  ===  File description
  ===========================================================================*/

#include "agenda_class.h"
#include "physics_funcs.h"
#include "constants.h"
#include "gas_scattering.h"
#include "optproperties.h"
#include "rte.h"
#include "transmissionmatrix.h"
#include <cmath>

using Constant::pi;
using Constant::boltzmann_constant;


/*!
  \file   m_gas_scattering.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-02-08

  \brief  Workspace functions related to gas scattering.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringOff(Workspace& ws,
                       Index& gas_scattering_do,
                       Agenda& gas_scattering_agenda,
                       const Verbosity&) {
  // set flag to False (default)
  gas_scattering_do = 0;

  gas_scattering_agenda = Agenda(ws);
  gas_scattering_agenda.set_name("gas_scattering_agenda");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringCoefXsecConst(PropagationMatrix& sca_coef,
                                 const Vector& f_grid,
                                 const Numeric& rtp_pressure,
                                 const Numeric& rtp_temperature,
                                 const Index& stokes_dim,
                                 const Numeric& ConstXsec,
                                 const Verbosity&) {
  // Some basic sizes
  const Index nf = f_grid.nelem();

  // Number density
  Numeric N = number_density(rtp_pressure, rtp_temperature);

  //Vector of constant cross sections
  Vector Xsec(nf, ConstXsec);

  // set coefficients
  PropagationMatrix sca_coef_temp(nf, stokes_dim);
  sca_coef_temp.SetZero();
  sca_coef_temp.Kjj() += Xsec;
  sca_coef_temp.Kjj() *= N;

  sca_coef=sca_coef_temp;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringCoefAirSimple(PropagationMatrix& sca_coef,
                                  const Vector& f_grid,
                                  const Numeric& rtp_pressure,
                                  const Numeric& rtp_temperature,
                                  const Index& stokes_dim,
                                  const Verbosity&) {

  Vector coefficients{3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};

  // Number density
  Numeric N = number_density(rtp_pressure, rtp_temperature);

  PropagationMatrix sca_coef_temp(f_grid.nelem(), stokes_dim);
  sca_coef_temp.SetZero();

  Matrix eye(stokes_dim, stokes_dim, 0.0);
  for (Index i = 0; i < stokes_dim; i++){
    eye(i,i) = 1;
  }

  for (Index f = 0; f < f_grid.nelem(); f++) {
    Matrix Xsec = eye;
    Numeric wavelen = Conversion::freq2wavelen(f_grid[f]) * 1e6;
    Numeric sigma;
    Numeric sum = 0;
    for (Index i = 0; i < 4; i++){
      sum += coefficients[i] * pow(wavelen, -2.0 * static_cast<double>(i));
    }
    sigma = pow(wavelen, -4) * sum * 1e-32;
    Xsec *= sigma;
    sca_coef_temp.SetAtPosition(Xsec, f);
  }

  sca_coef_temp.Kjj() *= N;
  sca_coef=sca_coef_temp;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringMatrixIsotropic(TransmissionMatrix& sca_mat,
                                   Vector& sca_fct_legendre,
                                   const Vector& in_los,
                                   const Vector& out_los,
                                   const Index& stokes_dim,
                                   const Index& gas_scattering_output_type,
                                   const Verbosity&) {
  //out
  if (gas_scattering_output_type) {
    sca_fct_legendre.resize(1);
    sca_fct_legendre = 1.;

  } else {
    if (in_los.nelem() > 0 && out_los.nelem() > 0) {
      TransmissionMatrix sca_mat_temp(1, stokes_dim);
      sca_mat_temp.setIdentity();
      sca_mat = sca_mat_temp;
    } else {
      // set the scattering matrics empty in case the in and out los are empty
      sca_mat = TransmissionMatrix();
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scatteringMatrixRayleigh(TransmissionMatrix& sca_mat,
                                  Vector& sca_fct_legendre,
                                  const Vector& in_los,
                                  const Vector& out_los,
                                  const Index& stokes_dim,
                                  const Index& gas_scattering_output_type,
                                  const Numeric& depolarization_factor,
                                  const Verbosity&) {

  ARTS_USER_ERROR_IF(in_los.nelem() != out_los.nelem(),
    "The length of the vectors of incoming and outgoing direction must be the same.")



    if (gas_scattering_output_type) {
    sca_fct_legendre.resize(3);
    sca_fct_legendre = {1, 0, 0.1};

  } else {

    //if in_los or out_los is empty then sca_mat is empty.
    if (in_los.nelem()>0 && out_los.nelem()>0){

      Index atmosphere_dim = 1;
      if (in_los.nelem()==2){
        atmosphere_dim=3;
      }

      //For the scattering calculation we need the propagation direction of incoming
      //and outgoing radiation. Therefore we have to convert the line of sights to
      //propagation directions.
      Vector in_prop;
      Vector out_prop;

      mirror_los(in_prop, in_los,atmosphere_dim);
      mirror_los(out_prop, out_los,atmosphere_dim);

      // calc_scatteringAngle() between in_los and out_los
      Numeric za_inc = in_prop[0];
      Numeric za_sca = out_prop[0];

      Numeric aa_inc = 0;
      Numeric aa_sca = 0;
      if (atmosphere_dim==3){
        aa_inc = in_prop[1];
        aa_sca = out_prop[1];
      }

      Numeric theta_rad = scat_angle(za_sca, aa_sca, za_inc, aa_inc);

      // Rayleigh phase matrix in scattering system
      Vector pha_mat_int = calc_rayleighPhaMat(theta_rad, stokes_dim);

      // transform the phase matrix
      Matrix pha_mat(stokes_dim, stokes_dim, 0.0);

      // account for depolarization factor
      Numeric delta =
          (1 - depolarization_factor) / (1 + 0.5 * depolarization_factor);
      Numeric delta_prime =
          (1 - 2 * depolarization_factor) / (1 - depolarization_factor);
      Vector depol(6, 0.0);

      switch (stokes_dim) {
        case 4:
          depol[5] = 1;
          pha_mat_int[5] *= delta_prime;
          [[fallthrough]];
        case 1:
          depol[0] = 1;
          depol *= 1 - delta;
          pha_mat_int += depol;
      }

      pha_mat_labCalc(pha_mat, pha_mat_int, out_los[0], out_los[1], in_los[0], in_los[1], theta_rad);

      TransmissionMatrix sca_mat_temp(pha_mat);

      sca_mat = sca_mat_temp;
    } else {
      // set the scattering matrics empty in case the in and out los are empty
      sca_mat = TransmissionMatrix();
    }
  }
}
