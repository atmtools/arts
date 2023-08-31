/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <workspace.h>
#include "physics_funcs.h"
#include "arts_conversions.h"
#include "gas_scattering.h"
#include "optproperties.h"
#include "rte.h"
#include "rtepack.h"
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
void gas_scatteringOff(Index& gas_scattering_do,
                       Agenda& gas_scattering_agenda) {
  // set flag to False (default)
  gas_scattering_do = 0;

  gas_scattering_agenda = Agenda("gas_scattering_agenda");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scattering_coefXsecConst(PropmatVector& gas_scattering_coef,
                                 const Vector& f_grid,
                                 const Numeric& rtp_pressure,
                                 const Numeric& rtp_temperature,
                                 const Numeric& ConstXsec) {
  const auto n = ConstXsec * number_density(rtp_pressure, rtp_temperature);
  gas_scattering_coef = PropmatVector(f_grid.nelem(), Propmat{n, n, n, n});
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scattering_coefAirSimple(PropmatVector& gas_scattering_coef,
                                  const Vector& f_grid,
                                  const Numeric& rtp_pressure,
                                  const Numeric& rtp_temperature) {
  static constexpr std::array coefficients{
      3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};

  gas_scattering_coef = PropmatVector(f_grid.nelem());

  const auto numden = number_density(rtp_pressure, rtp_temperature);

  for (Index f = 0; f < f_grid.nelem(); f++) {
    const Numeric wavelen = Conversion::freq2wavelen(f_grid[f]) * 1e6;
    Numeric sum = 0;
    Numeric pows = 1;
    for (auto& coef: coefficients) {
      sum += coef * pows;
      pows /= Math::pow2(wavelen);
    }
    gas_scattering_coef[f].A() = numden * 1e-32 * sum / Math::pow4(wavelen);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scattering_matIsotropic(MuelmatVector& gas_scattering_mat,
                                 Vector& gas_scattering_fct_legendre,
                                 const Vector& gas_scattering_los_in,
                                 const Vector& gas_scattering_los_out,
                                 const Index& gas_scattering_output_type) {
  // FIXME: It seems like this method should be split in 2 methods, one for the
  // output type 0 and one for the output type 1. The output type 1 should be
  // called gas_scattering_matIsotropic and the output type 0 should be called
  // gas_scattering_fct_legendreIsotropic. The current implementation is a bit
  // confusing; at the very least there should not be any if
  // (gas_scattering_output_type) statements...

  //out
  if (gas_scattering_output_type) {
    gas_scattering_fct_legendre.resize(1);
    gas_scattering_fct_legendre = 1.;

  } else {
    if (gas_scattering_los_in.nelem() > 0 && gas_scattering_los_out.nelem() > 0) {
      gas_scattering_mat = MuelmatVector(1, Muelmat{1.0});
    } else {
      // set the scattering matrics empty in case the in and out los are empty
      gas_scattering_mat = MuelmatVector();
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void gas_scattering_matRayleigh(MuelmatVector& gas_scattering_mat,
                                Vector& gas_scattering_fct_legendre,
                                const Vector& gas_scattering_los_in,
                                const Vector& gas_scattering_los_out,
                                const Index& gas_scattering_output_type,
                                const Numeric& depolarization_factor) {
  // FIXME: Again, gas_scattering_output_type is not a good way to implement
  // this (see comment above, 2023-05-19, in gas_scattering_matIsotropic).

  ARTS_USER_ERROR_IF(gas_scattering_los_in.nelem() != gas_scattering_los_out.nelem(),
    "The length of the vectors of incoming and outgoing direction must be the same.")



    if (gas_scattering_output_type) {
    gas_scattering_fct_legendre.resize(3);
    gas_scattering_fct_legendre = {1, 0, 0.1};

  } else {

    //if gas_scattering_los_in or gas_scattering_los_out is empty then gas_scattering_mat is empty.
    if (gas_scattering_los_in.nelem()>0 && gas_scattering_los_out.nelem()>0){
      //For the scattering calculation we need the propagation direction of incoming
      //and outgoing radiation. Therefore we have to convert the line of sights to
      //propagation directions.
      Vector in_prop;
      Vector out_prop;

      mirror_los(in_prop, gas_scattering_los_in);
      mirror_los(out_prop, gas_scattering_los_out);

      // calc_scatteringAngle() between gas_scattering_los_in and gas_scattering_los_out
      Numeric za_inc = in_prop[0];
      Numeric za_sca = out_prop[0];

      Numeric aa_inc = 0;
      Numeric aa_sca = 0;
      aa_inc = in_prop[1];
      aa_sca = out_prop[1];

      Numeric theta_rad = scat_angle(za_sca, aa_sca, za_inc, aa_inc);

      // Rayleigh phase matrix in scattering system
      Vector pha_mat_int = calc_rayleighPhaMat(theta_rad);

      // transform the phase matrix
      Matrix pha_mat(4, 4, 0.0);

      // account for depolarization factor
      Numeric delta =
          (1 - depolarization_factor) / (1 + 0.5 * depolarization_factor);
      Numeric delta_prime =
          (1 - 2 * depolarization_factor) / (1 - depolarization_factor);
      Vector depol(6, 0.0);

      depol[5] = 1;
      pha_mat_int[5] *= delta_prime;
      depol[0] = 1;
      depol *= 1 - delta;
      pha_mat_int += depol;

      pha_mat_labCalc(pha_mat,
                      pha_mat_int,
                      gas_scattering_los_out[0],
                      gas_scattering_los_out[1],
                      gas_scattering_los_in[0],
                      gas_scattering_los_in[1],
                      theta_rad);

      MuelmatVector sca_mat_temp(1, Muelmat{pha_mat});

      gas_scattering_mat = sca_mat_temp;
    } else {
      // set the scattering matrics empty in case the in and out los are empty
      gas_scattering_mat = MuelmatVector();
    }
  }
}
