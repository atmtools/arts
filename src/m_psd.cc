/* Copyright (C) 2017
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Jana Mendrok     <jana.mendrok@gmail.com>
   Manfred Brath    <manfred.brath@uni-hamburg.de>
                         
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
  @file   m_psd.cc
  @author Patrick Eriksson, Jana Mendrok, Manfred Brath
  @date   2017-11-05

  \brief  Workspace functions related to particle size distributions.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "lin_alg.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "psd.h"

/*===========================================================================
  === PSDs of Mono type
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdMonoDispersive(Matrix& psd_data,
                       Tensor3& dpsd_data_dx,
                       const Vector& pnd_agenda_input_t,
                       const Matrix& pnd_agenda_input,
                       const ArrayOfString& pnd_agenda_input_names,
                       const ArrayOfString& dpnd_data_dx_names,
                       const ArrayOfArrayOfScatteringMetaData& scat_meta,
                       const Index& species_index,
                       const Numeric& t_min,
                       const Numeric& t_max,
                       const Index& picky) {
  psd_mono_common(psd_data,
                  dpsd_data_dx,
                  "ntot",
                  pnd_agenda_input_t,
                  pnd_agenda_input,
                  pnd_agenda_input_names,
                  dpnd_data_dx_names,
                  scat_meta,
                  species_index,
                  t_min,
                  t_max,
                  picky);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdMonoMass(Matrix& psd_data,
                 Tensor3& dpsd_data_dx,
                 const Vector& pnd_agenda_input_t,
                 const Matrix& pnd_agenda_input,
                 const ArrayOfString& pnd_agenda_input_names,
                 const ArrayOfString& dpnd_data_dx_names,
                 const ArrayOfArrayOfScatteringMetaData& scat_meta,
                 const Index& species_index,
                 const Numeric& t_min,
                 const Numeric& t_max,
                 const Index& picky) {
  psd_mono_common(psd_data,
                  dpsd_data_dx,
                  "mass",
                  pnd_agenda_input_t,
                  pnd_agenda_input,
                  pnd_agenda_input_names,
                  dpnd_data_dx_names,
                  scat_meta,
                  species_index,
                  t_min,
                  t_max,
                  picky);
}

/*===========================================================================
  === PSDs of MGD type
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGamma(Matrix& psd_data,
                      Tensor3& dpsd_data_dx,
                      const Vector& psd_size_grid,
                      const Vector& pnd_agenda_input_t,
                      const Matrix& pnd_agenda_input,
                      const ArrayOfString& pnd_agenda_input_names,
                      const ArrayOfString& dpnd_data_dx_names,
                      const Numeric& n0,
                      const Numeric& mu,
                      const Numeric& la,
                      const Numeric& ga,
                      const Numeric& t_min,
                      const Numeric& t_max,
                      const Index& picky) {
  // Standard checks
  START_OF_PSD_METHODS();

  // Additional (basic) checks
  ARTS_USER_ERROR_IF (nin > 4,
        "The number of columns in *pnd_agenda_input* must "
        "be 0, 1, 2, 3 or 4.");

  // Check fixed parameters
  const Index n0_fixed = (Index) !(std::isnan(n0));
  const Index mu_fixed = (Index) !(std::isnan(mu));
  const Index la_fixed = (Index) !(std::isnan(la));
  const Index ga_fixed = (Index) !(std::isnan(ga));
  //
  ARTS_USER_ERROR_IF (nin + n0_fixed + mu_fixed + la_fixed + ga_fixed != 4,
        "This PSD has four free parameters. This means that "
        "the number\nof columns in *pnd_agenda_input* and the "
        "number of numerics\n(i.e. non-NaN) and among "
        "the GIN arguments n0, mu, la and\nga must add up to "
        "four. And this was found not to be the case.");

  // Create vectors to hold the four MGD and the "extra" parameters
  Vector mgd_pars(4);
  ArrayOfIndex mgd_i_pai = {-1, -1, -1, -1};  // Position in pnd_agenda_input
  {
    Index nhit = 0;
    if (n0_fixed) {
      mgd_pars[0] = n0;
    } else {
      mgd_i_pai[0] = nhit++;
    }
    if (mu_fixed) {
      mgd_pars[1] = mu;
    } else {
      mgd_i_pai[1] = nhit++;
    }
    if (la_fixed) {
      mgd_pars[2] = la;
    } else {
      mgd_i_pai[2] = nhit++;
    }
    if (ga_fixed) {
      mgd_pars[3] = ga;
    } else {
      mgd_i_pai[3] = nhit++;
    }
  }

  // Determine what derivatives to do and their positions
  ArrayOfIndex mgd_do_jac = {
      0,
      0,
      0,
      0,
  };
  ArrayOfIndex mgd_i_jac = {
      -1, -1, -1, -1};  // Position among jacobian quantities
  //
  for (Index i = 0; i < ndx; i++) {
    for (Index j = 0; j < 4; j++) {
      if (dx2in[i] == mgd_i_pai[j]) {
        mgd_do_jac[j] = 1;
        mgd_i_jac[j] = i;
        break;
      }
    }
  }

  // Loop input data and calculate PSDs
  for (Index ip = 0; ip < np; ip++) {
    // Extract MGD parameters
    for (Index i = 0; i < 4; i++) {
      if (mgd_i_pai[i] >= 0) {
        mgd_pars[i] = pnd_agenda_input(ip, mgd_i_pai[i]);
      }
    }
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if n0==0 and no jacobians requested.
    if ((mgd_pars[0] == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
          "This is outside the specified allowed range: [ max(0.,", t_min,
          "), ", t_max, " ]")
      continue;
      // If here, we are ready with this point!
    }

    // Check that la and ga are OK
    ARTS_USER_ERROR_IF (mgd_pars[2] <= 0,
                        "Bad MGD parameter detected: la <= 0");
    ARTS_USER_ERROR_IF (mgd_pars[3] <= 0,
                        "Bad MGD parameter detected: ga <= 0");

    // Calculate PSD and derivatives
    Matrix jac_data(4, nsi);
    //
    mgd_with_derivatives(psd_data(ip, joker),
                         jac_data,
                         psd_size_grid,
                         mgd_pars[0],
                         mgd_pars[1],
                         mgd_pars[2],
                         mgd_pars[3],
                         mgd_do_jac[0],
                         mgd_do_jac[1],
                         mgd_do_jac[2],
                         mgd_do_jac[3]);
    //
    for (Index i = 0; i < 4; i++) {
      if (mgd_do_jac[i]) {
        dpsd_data_dx(mgd_i_jac[i], ip, joker) = jac_data(i, joker);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGammaMass(Matrix& psd_data,
                          Tensor3& dpsd_data_dx,
                          const Vector& psd_size_grid,
                          const Vector& pnd_agenda_input_t,
                          const Matrix& pnd_agenda_input,
                          const ArrayOfString& pnd_agenda_input_names,
                          const ArrayOfString& dpnd_data_dx_names,
                          const Numeric& scat_species_a,
                          const Numeric& scat_species_b,
                          const Numeric& n0,
                          const Numeric& mu,
                          const Numeric& la,
                          const Numeric& ga,
                          const Numeric& t_min,
                          const Numeric& t_max,
                          const Index& picky) {
  // Standard checks
  START_OF_PSD_METHODS();

  // Additional (basic) checks
  ARTS_USER_ERROR_IF (nin < 1 || nin > 4,
        "The number of columns in *pnd_agenda_input* must "
        "be 1, 2, 3 or 4.");
  ARTS_USER_ERROR_IF (scat_species_a <= 0,
                      "*scat_species_a* should be > 0.");
  ARTS_USER_ERROR_IF (scat_species_b <= 0 || scat_species_b >= 5,
                      "*scat_species_b* should be > 0 and < 5.");

  // Check and determine dependent and fixed parameters
  const Index n0_depend = (Index)n0 == -999;
  const Index mu_depend = (Index)mu == -999;
  const Index la_depend = (Index)la == -999;
  const Index ga_depend = (Index)ga == -999;
  //
  ARTS_USER_ERROR_IF (n0_depend + mu_depend + la_depend + ga_depend != 1,
        "One (but only one) of n0, mu, la and ga must be NaN, "
        "to flag that this parameter is the one dependent of "
        "mass content.");
  ARTS_USER_ERROR_IF (mu_depend || ga_depend,
        "Sorry, mu and la are not yet allowed to be the "
        "dependent parameter.");
  //
  const Index n0_fixed = (Index) !(n0_depend || std::isnan(n0));
  const Index mu_fixed = (Index) !(mu_depend || std::isnan(mu));
  const Index la_fixed = (Index) !(la_depend || std::isnan(la));
  const Index ga_fixed = (Index) !(ga_depend || std::isnan(ga));
  //
  ARTS_USER_ERROR_IF (nin + n0_fixed + mu_fixed + la_fixed + ga_fixed != 4,
        "This PSD has four free parameters. This means that "
        "the number\nof columns in *pnd_agenda_input* and the "
        "number of numerics\n(i.e. not -999 or NaN) and among "
        "the GIN arguments n0, mu, la and\nga must add up to "
        "four. And this was found not to be the case.");

  // Create vectors to hold the four MGD and the "extra" parameters
  Vector mgd_pars(4), ext_pars(1);
  ArrayOfIndex mgd_i_pai = {-1, -1, -1, -1};  // Position in pnd_agenda_input
  const ArrayOfIndex ext_i_pai = {0};         // Position in pnd_agenda_input
  {
    Index nhit = 1;  // As mass always occupies first position
    if (n0_fixed) {
      mgd_pars[0] = n0;
    } else if (!n0_depend) {
      mgd_i_pai[0] = nhit++;
    }
    if (mu_fixed) {
      mgd_pars[1] = mu;
    } else if (!mu_depend) {
      mgd_i_pai[1] = nhit++;
    }
    if (la_fixed) {
      mgd_pars[2] = la;
    } else if (!la_depend) {
      mgd_i_pai[2] = nhit++;
    }
    if (ga_fixed) {
      mgd_pars[3] = ga;
    } else if (!ga_depend) {
      mgd_i_pai[3] = nhit++;
    }
  }

  // Determine what derivatives to do and their positions
  ArrayOfIndex mgd_do_jac = {
      0,
      0,
      0,
      0,
  };
  ArrayOfIndex ext_do_jac = {0};
  ArrayOfIndex mgd_i_jac = {
      -1, -1, -1, -1};            // Position among jacobian quantities
  ArrayOfIndex ext_i_jac = {-1};  // Position among jacobian quantities
  //
  for (Index i = 0; i < ndx; i++) {
    if (dx2in[i] == 0)  // That is,  mass is a derivative
    {
      ext_do_jac[0] = 1;
      ext_i_jac[0] = i;
    } else  // Otherwise, either n0, mu, la or ga
    {
      for (Index j = 0; j < 4; j++) {
        if (dx2in[i] == mgd_i_pai[j]) {
          mgd_do_jac[j] = 1;
          mgd_i_jac[j] = i;
          break;
        }
      }
    }
  }

  // Loop input data and calculate PSDs
  for (Index ip = 0; ip < np; ip++) {
    // Extract mass
    ext_pars[0] = pnd_agenda_input(ip, ext_i_pai[0]);
    // Extract core MGD parameters
    for (Index i = 0; i < 4; i++) {
      if (mgd_i_pai[i] >= 0) {
        mgd_pars[i] = pnd_agenda_input(ip, mgd_i_pai[i]);
      }
    }
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if mass==0 and no jacobians requested.
    if ((ext_pars[0] == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
          "This is outside the specified allowed range: [ max(0.,", t_min,
          "), ", t_max, " ]")
      continue;
      // If here, we are ready with this point!
    }

    // Derive the dependent parameter
    //
    Numeric mub1 = 0, eterm = 0, scfac = 0;
    //
    if (n0_depend) {
      mub1 = mgd_pars[1] + scat_species_b + 1;
      eterm = mub1 / mgd_pars[3];
      scfac = (mgd_pars[3] * pow(mgd_pars[2], eterm)) /
              (scat_species_a * tgamma(eterm));
      mgd_pars[0] = scfac * ext_pars[0];
    } else if (la_depend) {
      ARTS_USER_ERROR_IF (ext_pars[0] <= 0,
            "The mass content must be > 0 when la is "
            "the dependent parameter.");
      mub1 = mgd_pars[1] + scat_species_b + 1;
      eterm = mub1 / mgd_pars[3];
      scfac = mgd_pars[3] / (scat_species_a * mgd_pars[0] * tgamma(eterm));
      scfac = pow(scfac, -1 / eterm);
      mgd_pars[2] = scfac * pow(ext_pars[0], -1 / eterm);
    } else {
      ARTS_ASSERT(0);
    }

    // Now when all four MGD parameters are set, check that they were OK from
    // start, or became OK if set
    ARTS_USER_ERROR_IF (mub1 <= 0,
                        "Bad MGD parameter detected: mu + b + 1 <= 0");
    ARTS_USER_ERROR_IF (mgd_pars[2] <= 0,
                        "Bad MGD parameter detected: la <= 0");
    ARTS_USER_ERROR_IF (mgd_pars[3] <= 0,
                        "Bad MGD parameter detected: ga <= 0");

    // Calculate PSS
    Matrix jac_data(4, nsi);
    mgd_with_derivatives(psd_data(ip, joker),
                         jac_data,
                         psd_size_grid,
                         mgd_pars[0],
                         mgd_pars[1],
                         mgd_pars[2],
                         mgd_pars[3],
                         (bool)mgd_do_jac[0] || n0_depend,
                         (bool)mgd_do_jac[1] || mu_depend,
                         (bool)mgd_do_jac[2] || la_depend,
                         (bool)mgd_do_jac[3] || ga_depend);

    // Derivative with respect to mass
    if (ext_do_jac[0]) {
      if (n0_depend) {
        dpsd_data_dx(ext_i_jac[0], ip, joker) = jac_data(0, joker);
        dpsd_data_dx(ext_i_jac[0], ip, joker) *= scfac;
      } else if (la_depend) {
        dpsd_data_dx(ext_i_jac[0], ip, joker) = jac_data(2, joker);
        dpsd_data_dx(ext_i_jac[0], ip, joker) *=
            scfac * (-1 / eterm) * pow(ext_pars[0], -(1 / eterm + 1));
      } else {
        ARTS_ASSERT(0);
      }
    }

    // Derivatives for non-dependent native parameters
    for (Index i = 0; i < 4; i++) {
      if (mgd_do_jac[i]) {
        dpsd_data_dx(mgd_i_jac[i], ip, joker) = jac_data(i, joker);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGammaMassNtot(Matrix& psd_data,
                              Tensor3& dpsd_data_dx,
                              const Vector& psd_size_grid,
                              const Vector& pnd_agenda_input_t,
                              const Matrix& pnd_agenda_input,
                              const ArrayOfString& pnd_agenda_input_names,
                              const ArrayOfString& dpnd_data_dx_names,
                              const Numeric& scat_species_a,
                              const Numeric& scat_species_b,
                              const Numeric& n0,
                              const Numeric& mu,
                              const Numeric& la,
                              const Numeric& ga,
                              const Numeric& t_min,
                              const Numeric& t_max,
                              const Index& picky) {
  psd_mgd_mass_and_something(psd_data,
                             dpsd_data_dx,
                             "Ntot",
                             psd_size_grid,
                             pnd_agenda_input_t,
                             pnd_agenda_input,
                             pnd_agenda_input_names,
                             dpnd_data_dx_names,
                             scat_species_a,
                             scat_species_b,
                             n0,
                             mu,
                             la,
                             ga,
                             t_min,
                             t_max,
                             picky);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGammaMassMeanParticleMass(Matrix& psd_data,
                                          Tensor3& dpsd_data_dx,
                                          const Vector& psd_size_grid,
                                          const Vector& pnd_agenda_input_t,
                                          const Matrix& pnd_agenda_input,
                                          const ArrayOfString& pnd_agenda_input_names,
                                          const ArrayOfString& dpnd_data_dx_names,
                                          const Numeric& scat_species_a,
                                          const Numeric& scat_species_b,
                                          const Numeric& n0,
                                          const Numeric& mu,
                                          const Numeric& la,
                                          const Numeric& ga,
                                          const Numeric& t_min,
                                          const Numeric& t_max,
                                          const Index& picky) {
  psd_mgd_mass_and_something(psd_data,
                             dpsd_data_dx,
                             "mean particle mass",
                             psd_size_grid,
                             pnd_agenda_input_t,
                             pnd_agenda_input,
                             pnd_agenda_input_names,
                             dpnd_data_dx_names,
                             scat_species_a,
                             scat_species_b,
                             n0,
                             mu,
                             la,
                             ga,
                             t_min,
                             t_max,
                             picky);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGammaMassXmean(Matrix& psd_data,
                               Tensor3& dpsd_data_dx,
                               const Vector& psd_size_grid,
                               const Vector& pnd_agenda_input_t,
                               const Matrix& pnd_agenda_input,
                               const ArrayOfString& pnd_agenda_input_names,
                               const ArrayOfString& dpnd_data_dx_names,
                               const Numeric& scat_species_a,
                               const Numeric& scat_species_b,
                               const Numeric& n0,
                               const Numeric& mu,
                               const Numeric& la,
                               const Numeric& ga,
                               const Numeric& t_min,
                               const Numeric& t_max,
                               const Index& picky) {
  psd_mgd_mass_and_something(psd_data,
                             dpsd_data_dx,
                             "mean size",
                             psd_size_grid,
                             pnd_agenda_input_t,
                             pnd_agenda_input,
                             pnd_agenda_input_names,
                             dpnd_data_dx_names,
                             scat_species_a,
                             scat_species_b,
                             n0,
                             mu,
                             la,
                             ga,
                             t_min,
                             t_max,
                             picky);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGammaMassXmedian(Matrix& psd_data,
                                 Tensor3& dpsd_data_dx,
                                 const Vector& psd_size_grid,
                                 const Vector& pnd_agenda_input_t,
                                 const Matrix& pnd_agenda_input,
                                 const ArrayOfString& pnd_agenda_input_names,
                                 const ArrayOfString& dpnd_data_dx_names,
                                 const Numeric& scat_species_a,
                                 const Numeric& scat_species_b,
                                 const Numeric& n0,
                                 const Numeric& mu,
                                 const Numeric& la,
                                 const Numeric& ga,
                                 const Numeric& t_min,
                                 const Numeric& t_max,
                                 const Index& picky) {
  psd_mgd_mass_and_something(psd_data,
                             dpsd_data_dx,
                             "median size",
                             psd_size_grid,
                             pnd_agenda_input_t,
                             pnd_agenda_input,
                             pnd_agenda_input_names,
                             dpnd_data_dx_names,
                             scat_species_a,
                             scat_species_b,
                             n0,
                             mu,
                             la,
                             ga,
                             t_min,
                             t_max,
                             picky);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdModifiedGammaMassSingleMoment(
    Matrix& psd_data,
    Tensor3& dpsd_data_dx,
    const Vector& psd_size_grid,
    const Vector& pnd_agenda_input_t,
    const Matrix& pnd_agenda_input,
    const ArrayOfString& pnd_agenda_input_names,
    const ArrayOfString& dpnd_data_dx_names,
    const Numeric& scat_species_a,
    const Numeric& scat_species_b,
    const Numeric& n_alpha,
    const Numeric& n_b,
    const Numeric& mu,
    const Numeric& gamma,
    const Numeric& t_min,
    const Numeric& t_max,
    const Index& picky) {
  psd_mgd_smm_common(psd_data,
                     dpsd_data_dx,
                     "generic",
                     psd_size_grid,
                     pnd_agenda_input_t,
                     pnd_agenda_input,
                     pnd_agenda_input_names,
                     dpnd_data_dx_names,
                     scat_species_a,
                     scat_species_b,
                     n_alpha,
                     n_b,
                     mu,
                     gamma,
                     t_min,
                     t_max,
                     picky);
}

/*===========================================================================
  === Input: IWC and T
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdDelanoeEtAl14(Matrix& psd_data,
                         Tensor3& dpsd_data_dx,
                         const Vector& psd_size_grid,
                         const Vector& pnd_agenda_input_t,
                         const Matrix& pnd_agenda_input,
                         const ArrayOfString& pnd_agenda_input_names,
                         const ArrayOfString& dpnd_data_dx_names,
                         const Numeric& iwc,
                         const Numeric& n0,
                         const Numeric& dm,
                         const Numeric& rho,
                         const Numeric& alpha,
                         const Numeric& beta,
                         const Numeric& t_min,
                         const Numeric& t_max,
                         const Numeric& dm_min,
                         const Index& picky) {
  // Standard checks
  START_OF_PSD_METHODS();

  // Additional (basic) checks
  ARTS_USER_ERROR_IF (nin > 2,
        "The number of columns in *pnd_agenda_input* must "
        "be 0, 1 or 2");

  // Check and determine dependent and fixed parameters
  const bool n0_depend = (Index)n0 == -999;
  const bool dm_depend = (Index)dm == -999;
  const bool iwc_depend = (Index)iwc == -999;

  // Check fixed parameters
  const bool iwc_fixed = !(std::isnan(iwc)) && !iwc_depend;
  const bool n0_fixed = !(std::isnan(n0)) && !n0_depend;
  const bool dm_fixed = !(std::isnan(dm)) && !dm_depend;

  ARTS_USER_ERROR_IF (!((nin + iwc_fixed + n0_fixed + dm_fixed == 2) ||
        (nin + iwc_fixed + n0_fixed + dm_fixed == 1)),
        "This PSD can have one or two independent parameters, that is \n"
        "the sum of the number of columns in pnd_agenda_input and\n"
        "non-NAN, non-dependent values in iwc, n0, dm must be equal to\n"
        "one or two.");

  ArrayOfIndex i_pai = {-1, -1, -1};  // Position in pnd_agenda_input

  Index nhit = 0;

  if ((n0_depend || dm_depend) && (!iwc_fixed)) {
    i_pai[0] = nhit++;
  }

  if ((!n0_depend) && (!n0_fixed)) {
    i_pai[1] = nhit++;
  }

  if ((!dm_depend) && (!dm_fixed)) {
    i_pai[2] = nhit++;
  }

  // Determine what derivatives to do and their positions
  ArrayOfIndex do_jac = {0, 0, 0};
  ArrayOfIndex i_jac = {-1, -1, -1};  // Position among jacobian quantities

  for (Index i = 0; i < ndx; ++i) {
    for (Index j = 0; j < 3; ++j) {
      if (dx2in[i] == i_pai[j]) {
        do_jac[j] = 1;
        i_jac[j] = i;
        break;
      }
    }
  }

  if (psd_size_grid[0] < std::numeric_limits<Numeric>::epsilon()) {
    ARTS_USER_ERROR_IF (psd_size_grid.nelem() < 2,
          "psd_size_grid has only one element which is 0. This is not allowed.");
  }

  Numeric iwc_p(0.0), n0_p(0.0), dm_p(0.0);
  // Loop input data and calculate PSDs
  for (Index ip = 0; ip < np; ip++) {
    Numeric t = pnd_agenda_input_t[ip];

    // Extract MGD parameters
    if (i_pai[0] >= 0) {
      iwc_p = pnd_agenda_input(ip, i_pai[0]);
    }
    if (i_pai[1] >= 0) {
      n0_p = pnd_agenda_input(ip, i_pai[1]);
    }
    if (i_pai[2] >= 0) {
      dm_p = pnd_agenda_input(ip, i_pai[2]);
    }

    if (n0_depend && dm_depend) {
      n0_p = n0_from_t(t);
      dm_p = dm_from_iwc_n0(iwc_p, n0_p, rho);
    } else if (n0_depend) {
      n0_p = n0_from_iwc_dm(iwc_p, dm_p, rho);
    } else if (dm_depend) {
      dm_p = dm_from_iwc_n0(iwc_p, n0_p, rho);
    }

    // Outside of [t_min,tmax]?
    if ((t < t_min) || (t > t_max)) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
           "This is outside the specified allowed range: [ max(0.,", t_min,
           "), ", t_max, " ]")
      continue;
    }

    if (iwc > 0.0) {
      ARTS_USER_ERROR_IF ((iwc > 0.0) && ((dm_p <= 0.0) || (dm_p < dm_min)),
           "The provided or inferred value of *Dm* (", dm_p, ") is "
           " less than zero or *Dm_min* and this is not allowed. "
           "This means that you have very small or zero values "
           "in *pnd_agenda_input* which is not supported "
           "by this PSD.\n")
    }

    // Calculate PSD and derivatives
    Matrix jac_data(1, nsi);
    Vector x_grid(psd_size_grid);
    x_grid *= 1.0 / dm_p;

    if (x_grid[0] < std::numeric_limits<Numeric>::epsilon()) {
      x_grid[0] = 0.1 * psd_size_grid[1];
    }

    delanoe_shape_with_derivative(
        psd_data(ip, joker), jac_data, x_grid, alpha, beta);
    psd_data(ip, joker) *= n0_p;
    jac_data(0, joker) *= n0_p;

    Vector dndx{jac_data(0, joker)};
    Vector dndn0{psd_data(ip, joker)};
    dndn0 *= (1.0 / n0_p);

    Vector dxddm = x_grid;
    dxddm *= (-1.0 / dm_p);
    Numeric dn0diwc = n0_p / iwc_p;

    if (do_jac[0]) {
      dpsd_data_dx(i_jac[0], ip, joker) = 0.0;

      if (dm_depend) {
        Numeric ddmdiwc = 0.25 * dm_p / iwc_p;
        Vector dndiwc = dxddm;
        dndiwc *= dndx;
        dndiwc *= ddmdiwc;
        dpsd_data_dx(i_jac[0], ip, joker) += dndiwc;
      } else if (n0_depend) {
        Vector dndiwc = dndn0;
        dndiwc *= dn0diwc;
        dpsd_data_dx(i_jac[0], ip, joker) += dndiwc;
      }
    }

    if (do_jac[1]) {
      dpsd_data_dx(i_jac[1], ip, joker) = dndn0;
      if (dm_depend) {
        Vector dndn02 = dndx;
        dndn02 *= dxddm;
        Numeric ddmdn0 = -0.25 / n0_p * dm_p;
        dndn02 *= ddmdn0;
        dpsd_data_dx(i_jac[1], ip, joker) += dndn02;
      }
    }

    if (do_jac[2]) {
      dpsd_data_dx(i_jac[2], ip, joker) = dxddm;
      dpsd_data_dx(i_jac[2], ip, joker) *= dndx;
      if (n0_depend) {
        Vector dndn02 = dndn0;
        Numeric dn0ddm = -4.0 * n0_p / dm_p;
        dndn02 *= dn0ddm;
        dpsd_data_dx(i_jac[2], ip, joker) += dndn02;
      }
    }

    // Ensure that results are zero when IWC is zero.
    if ((!iwc_depend) && (iwc_p == 0.0)) {
      psd_data(0, joker) = 0.0;
      for (size_t i = 0; i < 2; ++i) {
        if (do_jac[i]) {
          dpsd_data_dx(i_jac[i], ip, joker) = 0.0;
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdFieldEtAl07(Matrix& psd_data,
                    Tensor3& dpsd_data_dx,
                    const Vector& psd_size_grid,
                    const Vector& pnd_agenda_input_t,
                    const Matrix& pnd_agenda_input,
                    const ArrayOfString& pnd_agenda_input_names,
                    const ArrayOfString& dpnd_data_dx_names,
                    const Numeric& scat_species_a,
                    const Numeric& scat_species_b,
                    const String& regime,
                    const Numeric& t_min,
                    const Numeric& t_max,
                    const Numeric& t_min_psd,
                    const Numeric& t_max_psd,
                    const Numeric& b_min,
                    const Numeric& b_max,
                    const Index& picky) {
  // Standard checcks
  START_OF_PSD_METHODS();

  // Additional (basic) checks
  ARTS_USER_ERROR_IF (pnd_agenda_input.ncols() != 1,
                      "*pnd_agenda_input* must have one column.");
  ARTS_USER_ERROR_IF (regime != "TR" && regime != "ML",
                      "regime must either be \"TR\" or \"ML\".");
  ARTS_USER_ERROR_IF (scat_species_a <= 0,
                      "*scat_species_a* should be > 0.");
  ARTS_USER_ERROR_IF (scat_species_b < b_min || scat_species_b > b_max,
      "Method called with a mass-dimension-relation exponent b of ",
      scat_species_b, ".\n"
      "This is outside the specified allowed range: [", b_min, ",",
      b_max, "]")

  for (Index ip = 0; ip < np; ip++) {
    // Extract the input variables
    Numeric swc = pnd_agenda_input(ip, 0);
    Numeric t = pnd_agenda_input_t[ip];

    // NaN can be generated for extremly small SWC (such as 1e-78)
    // As a solution set a limit, and consider all below as zero
    if (abs(swc) < 1e-15) {
      swc = 0.0;
    }

    // No calc needed if swc==0 and no jacobians requested.
    if ((swc == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
          "This is outside the specified allowed range: [ max(0.,", t_min,
          "), ", t_max, " ]")
      continue;
      // If here, we are ready with this point!
    }

    // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
    if (t < t_min_psd) {
      t = t_min_psd;
    } else if (t > t_max_psd) {
      t = t_max_psd;
    }

    // Negative swc?
    Numeric psd_weight = 1.0;
    if (swc < 0) {
      psd_weight = -1.0;
      swc *= -1.0;
    }

    // Calculate PSD
    Vector psd_1p(nsi);
    if (swc != 0) {
      psd_snow_F07(psd_1p,
                   psd_size_grid,
                   swc,
                   t,
                   scat_species_a,
                   scat_species_b,
                   regime);
      for (Index i = 0; i < nsi; i++) {
        psd_data(ip, i) = psd_weight * psd_1p[i];
      }
    }

    // Calculate derivative with respect to IWC
    if (ndx) {
      //const Numeric dswc = max( 0.001*swc, 1e-7 );
      const Numeric dswc = 1e-9;
      const Numeric swcp = swc + dswc;
      psd_snow_F07(psd_1p,
                   psd_size_grid,
                   swcp,
                   t,
                   scat_species_a,
                   scat_species_b,
                   regime);
      for (Index i = 0; i < nsi; i++) {
        dpsd_data_dx(0, ip, i) = (psd_1p[i] - psd_data(ip, i)) / dswc;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdMcFarquaharHeymsfield97(Matrix& psd_data,
                                Tensor3& dpsd_data_dx,
                                const Vector& psd_size_grid,
                                const Vector& pnd_agenda_input_t,
                                const Matrix& pnd_agenda_input,
                                const ArrayOfString& pnd_agenda_input_names,
                                const ArrayOfString& dpnd_data_dx_names,
                                const Numeric& scat_species_a,
                                const Numeric& scat_species_b,
                                const Numeric& t_min,
                                const Numeric& t_max,
                                const Numeric& t_min_psd,
                                const Numeric& t_max_psd,
                                const Index& picky,
                                const Index& noisy) {
  // Standard checcks
  START_OF_PSD_METHODS();

  // Extra checks for this PSD
  ARTS_USER_ERROR_IF (pnd_agenda_input.ncols() != 1,
                      "*pnd_agenda_input* must have one column.");
  ARTS_USER_ERROR_IF (noisy && ndx,
        "Jacobian calculations and \"noisy\" can not be "
        "combined.");
  ARTS_USER_ERROR_IF (scat_species_b < 2.9 || scat_species_b > 3.1,
      "This PSD treats pure ice, using Dveq as size grid.\n"
      "This means that *scat_species_b* should be close to 3,\n"
      "but it is outside of the tolerated range of [2.9,3.1].\n"
      "Your value of *scat_species_b* is: ", scat_species_b)
  ARTS_USER_ERROR_IF (scat_species_a < 460 || scat_species_a > 500,
      "This PSD treats pure ice, using Dveq as size grid.\n"
      "This means that *scat_species_a* should be close to 480,\n"
      "but it is outside of the tolerated range of [460,500].\n"
      "Your value of *scat_species_a* is: ", scat_species_a)

  for (Index ip = 0; ip < np; ip++) {
    // Extract the input variables
    Numeric iwc = pnd_agenda_input(ip, 0);
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if iwc==0 and no jacobians requested.
    if ((iwc == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
          "This is outside the specified allowed range: [ max(0.,", t_min,
          "), ", t_max, " ]")
      continue;
      // If here, we are ready with this point!
    }

    // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
    if (t < t_min_psd) {
      t = t_min_psd;
    } else if (t > t_max_psd) {
      t = t_max_psd;
    }

    // Negative iwc?
    Numeric psd_weight = 1.0;
    if (iwc < 0) {
      psd_weight = -1.0;
      iwc *= -1.0;
    }

    // Calculate PSD
    Vector psd_1p(nsi);
    if (iwc != 0) {
      psd_cloudice_MH97(psd_1p, psd_size_grid, iwc, t, noisy);
      for (Index i = 0; i < nsi; i++) {
        psd_data(ip, i) = psd_weight * psd_1p[i];
      }
    }

    // Calculate derivative with respect to IWC
    if (ndx) {
      //const Numeric diwc = max( 0.001*iwc, 1e-9 );
      const Numeric diwc = 1e-9;
      const Numeric iwcp = iwc + diwc;
      psd_cloudice_MH97(psd_1p, psd_size_grid, iwcp, t, noisy);
      for (Index i = 0; i < nsi; i++) {
        dpsd_data_dx(0, ip, i) = (psd_1p[i] - psd_data(ip, i)) / diwc;
      }
    }
  }
}

/*===========================================================================
  === Input: RWC
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdAbelBoutle12(Matrix& psd_data,
                     Tensor3& dpsd_data_dx,
                     const Vector& psd_size_grid,
                     const Vector& pnd_agenda_input_t,
                     const Matrix& pnd_agenda_input,
                     const ArrayOfString& pnd_agenda_input_names,
                     const ArrayOfString& dpnd_data_dx_names,
                     const Numeric& scat_species_a,
                     const Numeric& scat_species_b,
                     const Numeric& t_min,
                     const Numeric& t_max,
                     const Index& picky) {
  psd_mgd_smm_common(psd_data,
                     dpsd_data_dx,
                     "Abel12",
                     psd_size_grid,
                     pnd_agenda_input_t,
                     pnd_agenda_input,
                     pnd_agenda_input_names,
                     dpnd_data_dx_names,
                     scat_species_a,
                     scat_species_b,
                     0,
                     0,
                     0,
                     0,
                     t_min,
                     t_max,
                     picky);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdWangEtAl16(Matrix& psd_data,
                   Tensor3& dpsd_data_dx,
                   const Vector& psd_size_grid,
                   const Vector& pnd_agenda_input_t,
                   const Matrix& pnd_agenda_input,
                   const ArrayOfString& pnd_agenda_input_names,
                   const ArrayOfString& dpnd_data_dx_names,
                   const Numeric& scat_species_a,
                   const Numeric& scat_species_b,
                   const Numeric& t_min,
                   const Numeric& t_max,
                   const Index& picky) {
  psd_mgd_smm_common(psd_data,
                     dpsd_data_dx,
                     "Wang16",
                     psd_size_grid,
                     pnd_agenda_input_t,
                     pnd_agenda_input,
                     pnd_agenda_input_names,
                     dpnd_data_dx_names,
                     scat_species_a,
                     scat_species_b,
                     0,
                     0,
                     0,
                     0,
                     t_min,
                     t_max,
                     picky);
}

/*===========================================================================
  === Input: GWC
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdFieldEtAl19(Matrix& psd_data,
                    Tensor3& dpsd_data_dx,
                    const Vector& psd_size_grid,
                    const Vector& pnd_agenda_input_t,
                    const Matrix& pnd_agenda_input,
                    const ArrayOfString& pnd_agenda_input_names,
                    const ArrayOfString& dpnd_data_dx_names,
                    const Numeric& scat_species_a,
                    const Numeric& scat_species_b,
                    const Numeric& t_min,
                    const Numeric& t_max,
                    const Index& picky) {
  psd_mgd_smm_common(psd_data,
                     dpsd_data_dx,
                     "Field19",
                     psd_size_grid,
                     pnd_agenda_input_t,
                     pnd_agenda_input,
                     pnd_agenda_input_names,
                     dpnd_data_dx_names,
                     scat_species_a,
                     scat_species_b,
                     0,
                     0,
                     0,
                     0,
                     t_min,
                     t_max,
                     picky);
}

/*===========================================================================
  === Input: Atmospheric model PSDs
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdSeifertBeheng06(Matrix& psd_data,
                        Tensor3& dpsd_data_dx,
                        const Vector& psd_size_grid,
                        const Vector& pnd_agenda_input_t,
                        const Matrix& pnd_agenda_input,
                        const ArrayOfString& pnd_agenda_input_names,
                        const ArrayOfString& dpnd_data_dx_names,
                        const String& hydrometeor_type,
                        const Numeric& t_min,
                        const Numeric& t_max,
                        const Index& picky) {
  // Some sizes
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();

  // Checks
  ARTS_USER_ERROR_IF (pnd_agenda_input.ncols() != nin,
        "Length of *pnd_agenda_input_names* and number of "
        "columns in *pnd_agenda_input* must be equal.");
  ARTS_USER_ERROR_IF (pnd_agenda_input.ncols() != 2,
        "*pnd_agenda_input* must have two columns"
        "(mass density and number density).");

  ARTS_USER_ERROR_IF (ndx > 2,
                      "*dpnd_data_dx_names* must have length <=2.");

  // check name tags
  ArrayOfIndex input_idx = {-1, -1};

  for (Index i = 0; i < nin; i++) {
    if ((Index)pnd_agenda_input_names[i].find("mass_density") != String::npos) {
      input_idx[0] = i;  //mass density index
    } else if ((Index)pnd_agenda_input_names[i].find("number_density") !=
               String::npos) {
      input_idx[1] = i;  //number density index
    }
  }

  ARTS_USER_ERROR_IF (input_idx[0] == -1,
                      "mass_density-tag not found ");
  ARTS_USER_ERROR_IF (input_idx[1] == -1,
                      "number_density-tag not found ");

  // look after jacobian tags
  ArrayOfIndex dpnd_data_dx_idx = {-1, -1};

  for (Index i = 0; i < ndx; i++) {
    if ((Index)dpnd_data_dx_names[i].find("mass_density") != String::npos) {
      dpnd_data_dx_idx[0] = i;  //mass density index
    } else if ((Index)dpnd_data_dx_names[i].find("number_density") !=
               String::npos) {
      dpnd_data_dx_idx[1] = i;  //number density index
    }
  }

  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize(np, nsi);
  psd_data = 0.0;
  if (ndx != 0) {
    dpsd_data_dx.resize(ndx, np, nsi);
    dpsd_data_dx = 0.0;
  } else {
    dpsd_data_dx.resize(0, 0, 0);
  }

  for (Index ip = 0; ip < np; ip++) {
    // Extract the input variables
    Numeric WC = pnd_agenda_input(ip, input_idx[0]);
    Numeric N_tot = pnd_agenda_input(ip, input_idx[1]);
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if swc==0 and no jacobians requested.
    if ((WC == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
          "This is outside the specified allowed range: [ max(0.,", t_min,
          "), ", t_max, " ]")
      continue;
      // If here, we are ready with this point!
    }

    // Negative swc?
    Numeric psd_weight = 1.0;
    if (WC < 0) {
      psd_weight = -1.0;
      WC *= -1.0;
    }

    // Calculate PSD and derivatives
    Vector psd_1p(nsi);
    Matrix dpsd_1p(nsi, 2);
    if (WC > 0) {
      psd_SB06(psd_1p, dpsd_1p, psd_size_grid, N_tot, WC, hydrometeor_type);

      for (Index i = 0; i < nsi; i++) {
        psd_data(ip, i) = psd_weight * psd_1p[i];

        for (Index idx = 0; idx < dpnd_data_dx_idx.nelem(); idx++) {
          // with respect to WC

          if (dpnd_data_dx_idx[idx] != -1) {
            dpsd_data_dx(dpnd_data_dx_idx[idx], ip, i) =
                psd_weight * dpsd_1p(i, idx);
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdMilbrandtYau05(Matrix& psd_data,
                       Tensor3& dpsd_data_dx,
                       const Vector& psd_size_grid,
                       const Vector& pnd_agenda_input_t,
                       const Matrix& pnd_agenda_input,
                       const ArrayOfString& pnd_agenda_input_names,
                       const ArrayOfString& dpnd_data_dx_names,
                       const String& hydrometeor_type,
                       const Numeric& t_min,
                       const Numeric& t_max,
                       const Index& picky) {
  // Some sizes
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();

  // Checks
  ARTS_USER_ERROR_IF (pnd_agenda_input.ncols() != nin,
        "Length of *pnd_agenda_input_names* and number of "
        "columns in *pnd_agenda_input* must be equal.");
  ARTS_USER_ERROR_IF (pnd_agenda_input.ncols() != 2,
        "*pnd_agenda_input* must have two columns"
        "(mass density and number density).");

  ARTS_USER_ERROR_IF (ndx > 2,
                      "*dpnd_data_dx_names* must have length <=2.");

  // check name tags
  ArrayOfIndex input_idx = {-1, -1};

  for (Index i = 0; i < nin; i++) {
    if ((Index)pnd_agenda_input_names[i].find("mass_density") != String::npos) {
      input_idx[0] = i;  //mass density index
    } else if ((Index)pnd_agenda_input_names[i].find("number_density") !=
               String::npos) {
      input_idx[1] = i;  //number density index
    }
  }

  ARTS_USER_ERROR_IF (input_idx[0] == -1,
                      "mass_density-tag not found ");
  ARTS_USER_ERROR_IF (input_idx[1] == -1,
                      "number_density-tag not found ");

  // look after jacobian tags
  ArrayOfIndex dpnd_data_dx_idx = {-1, -1};

  for (Index i = 0; i < ndx; i++) {
    if ((Index)dpnd_data_dx_names[i].find("mass_density") != String::npos) {
      dpnd_data_dx_idx[0] = i;  //mass density index
    } else if ((Index)dpnd_data_dx_names[i].find("number_density") !=
               String::npos) {
      dpnd_data_dx_idx[1] = i;  //number density index
    }
  }

  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize(np, nsi);
  psd_data = 0.0;
  if (ndx != 0) {
    dpsd_data_dx.resize(ndx, np, nsi);
    dpsd_data_dx = 0.0;
  } else {
    dpsd_data_dx.resize(0, 0, 0);
  }

  for (Index ip = 0; ip < np; ip++) {
    // Extract the input variables
    Numeric WC = pnd_agenda_input(ip, input_idx[0]);
    Numeric N_tot = pnd_agenda_input(ip, input_idx[1]);
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if wc==0 and no jacobians requested.
    if ((WC == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      ARTS_USER_ERROR_IF (picky,
          "Method called with a temperature of ", t, " K.\n"
          "This is outside the specified allowed range: [ max(0.,", t_min,
          "), ", t_max, " ]")
      continue;
      // If here, we are ready with this point!
    }

    // Negative wc?
    Numeric psd_weight = 1.0;
    if (WC < 0) {
      psd_weight = -1.0;
      WC *= -1.0;
    }

    // Calculate PSD and derivatives
    Vector psd_1p(nsi);
    Matrix dpsd_1p(nsi, 2);
    if (WC > 0) {
      psd_MY05(psd_1p, dpsd_1p, psd_size_grid, N_tot, WC, hydrometeor_type);

      for (Index i = 0; i < nsi; i++) {
        psd_data(ip, i) = psd_weight * psd_1p[i];

        for (Index idx = 0; idx < dpnd_data_dx_idx.nelem(); idx++) {
          // with respect to WC

          if (dpnd_data_dx_idx[idx] != -1) {
            dpsd_data_dx(dpnd_data_dx_idx[idx], ip, i) =
                psd_weight * dpsd_1p(i, idx);
          }
        }
      }
    }
  }
}
