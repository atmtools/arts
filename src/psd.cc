/* Copyright (C) 2017

   Jana Mendrok <jana.mendrok@gmail.com>
   Patrick Eriksson <patrick.eriksson@chalmers.se>

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

/**
  @file   microphysics.cc>
  @author Jana Mendrok, Patrick Eriksson
  @date   2017-11-05

  @brief  Internal functions associated with size distributions
*/

#include "psd.h"

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>
#include <stdexcept>

#include "arts.h"
#include "check_input.h"
#include "cloudbox.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "mc_antenna.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rng.h"
#include "sorting.h"

extern const Numeric PI;
extern const Numeric DENSITY_OF_ICE;
extern const Numeric DENSITY_OF_WATER;

void psd_cloudice_MH97(Vector& psd,
                       const Vector& diameter,
                       const Numeric& iwc,
                       const Numeric& t,
                       const bool noisy) {
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  assert(t > 0.);

  // skip calculation if IWC is 0.0
  if (iwc == 0.0) {
    return;
  }
  assert(iwc > 0.);

  // convert m to microns
  Vector d_um(nD);
  for (Index iD = 0; iD < nD; iD++) d_um[iD] = 1e6 * diameter[iD];
  //convert T from Kelvin to Celsius
  Numeric Tc = t - 273.15;

  //[kg/m3] -> [g/m3] as used by parameterisation
  Numeric ciwc = iwc * 1e3;
  Numeric cdensity = DENSITY_OF_ICE * 1e3;

  Numeric sig_a = 0., sig_b1 = 0.;
  Numeric sig_b2 = 0., sig_m = 0.;
  Numeric sig_aamu = 0., sig_bamu = 0.;
  Numeric sig_abmu = 0., sig_bbmu = 0.;
  Numeric sig_aasigma = 0., sig_basigma = 0;
  Numeric sig_absigma = 0., sig_bbsigma = 0.;

  if (noisy) {
    Rng rng;  //Random Number generator
    Index mc_seed;
    mc_seed = (Index)time(NULL);
    rng.seed(mc_seed, Verbosity());

    sig_a = 0.068, sig_b1 = 0.054;
    sig_b2 = 5.5e-3, sig_m = 0.0029;
    sig_aamu = 0.02, sig_bamu = 0.0005;
    sig_abmu = 0.023, sig_bbmu = 0.5e-3;
    sig_aasigma = 0.02, sig_basigma = 0.5e-3;
    sig_absigma = 0.023, sig_bbsigma = 4.7e-4;

    sig_a = ran_gaussian(rng, sig_a);
    sig_b1 = ran_gaussian(rng, sig_b1);
    sig_b2 = ran_gaussian(rng, sig_b2);
    sig_m = ran_gaussian(rng, sig_m);
    sig_aamu = ran_gaussian(rng, sig_aamu);
    sig_bamu = ran_gaussian(rng, sig_bamu);
    sig_abmu = ran_gaussian(rng, sig_abmu);
    sig_bbmu = ran_gaussian(rng, sig_bbmu);
    sig_aasigma = ran_gaussian(rng, sig_aasigma);
    sig_basigma = ran_gaussian(rng, sig_basigma);
    sig_absigma = ran_gaussian(rng, sig_absigma);
    sig_bbsigma = ran_gaussian(rng, sig_bbsigma);
  }

  // Split IWC in small and large size modes

  // Calculate IWC in each mode
  Numeric a = 0.252 + sig_a;  //g/m^3
  Numeric b1 = 0.837 + sig_b1;
  Numeric IWCs100 = min(ciwc, a * pow(ciwc, b1));
  Numeric IWCl100 = ciwc - IWCs100;

  // Gamma distribution component (small mode)

  Numeric b2 = -4.99e-3 + sig_b2;               //micron^-1
  Numeric m = 0.0494 + sig_m;                   //micron^-1
  Numeric alphas100 = b2 - m * log10(IWCs100);  //micron^-1

  // alpha, and hence dNdD1, becomes NaN if IWC>0.
  // this should be ensured to not happen before.
  //
  // alpha, and hence dNdD1, becomes negative for large IWC.
  // towards this limit, particles anyway get larger than 100um, i.e., outside
  // the size region the small-particle mode gamma distrib is intended for.
  // hence, leave dNdD1 at 0 for those cases.
  Vector dNdD1(nD, 0.);
  if (alphas100 > 0.) {
    Numeric Ns100 = 6 * IWCs100 * pow(alphas100, 5.) /
                    (PI * cdensity * tgamma(5.));  //micron^-5
    for (Index iD = 0; iD < nD; iD++)
      dNdD1[iD] = 1e18 * Ns100 * d_um[iD] *
                  exp(-alphas100 * d_um[iD]);  //micron^-4 -> m^-3 micron^-1
  }

  // Log normal distribution component (large mode)

  // for small IWC, IWCtotal==IWC<100 & IWC>100=0.
  // this will give dNdD2=NaN. avoid that by explicitly setting to 0
  Vector dNdD2(nD, 0.);
  if (IWCl100 > 0.) {
    //FIXME: Do we need to ensure mul100>0 and sigmal100>0?

    Numeric aamu = 5.20 + sig_aamu;
    Numeric bamu = 0.0013 + sig_bamu;
    Numeric abmu = 0.026 + sig_abmu;
    Numeric bbmu = -1.2e-3 + sig_bbmu;
    Numeric amu = aamu + bamu * Tc;
    Numeric bmu = abmu + bbmu * Tc;
    Numeric mul100 = amu + bmu * log10(IWCl100);

    Numeric aasigma = 0.47 + sig_aasigma;
    Numeric basigma = 2.1e-3 + sig_basigma;
    Numeric absigma = 0.018 + sig_absigma;
    Numeric bbsigma = -2.1e-4 + sig_bbsigma;
    Numeric asigma = aasigma + basigma * Tc;
    Numeric bsigma = absigma + bbsigma * Tc;
    Numeric sigmal100 = asigma + bsigma * log10(IWCl100);

    if ((mul100 > 0.) & (sigmal100 > 0.)) {
      Numeric a1 = 6 * IWCl100;  //g/m^3
      Numeric a2_fac = pow(PI, 3. / 2.) * cdensity * sqrt(2) *
                       exp(3 * mul100 + 9. / 2. * pow(sigmal100, 2)) *
                       sigmal100 * pow(1., 3);
      //a2 = a2_fac * d_um; //g/m^3/micron^4
      for (Index iD = 0; iD < nD; iD++)
        dNdD2[iD] = 1e18 * a1 / (a2_fac * d_um[iD]) *
                    exp(-0.5 * pow((log(d_um[iD]) - mul100) / sigmal100, 2));
      //micron^-4 -> m^-3 micron^-1
    }
  }

  for (Index iD = 0; iD < nD; iD++) {
    // FIXME: Do we still need this check here? Non-NaN of each mode should
    // now be ensure by the checks/if-loops for each of the modes separately.
    // I, JM, think (and hope).
    //if ( !std::isnan(dNdD1[iD]) && !std::isnan(dNdD2[iD]) )
    psd[iD] = (dNdD1[iD] + dNdD2[iD]) * 1e6;  // m^-3 m^-1
  }
}

void psd_mgd_mass_and_something(Matrix& psd_data,
                                Tensor3& dpsd_data_dx,
                                const String& something,
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
                                const Index& picky,
                                const Verbosity&) {
  // Standard checks
  START_OF_PSD_METHODS();

  // Additional (basic) checks
  if (nin < 1 || nin > 4)
    throw runtime_error(
        "The number of columns in *pnd_agenda_input* must "
        "be 2, 3 or 4.");
  if (scat_species_a <= 0)
    throw runtime_error("*scat_species_a* should be > 0.");
  if (scat_species_b <= 0 || scat_species_b >= 5)
    throw runtime_error("*scat_species_b* should be > 0 and < 5.");

  // Check and determine dependent and fixed parameters
  const Index n0_depend = (Index)n0 == -999;
  const Index mu_depend = (Index)mu == -999;
  const Index la_depend = (Index)la == -999;
  const Index ga_depend = (Index)ga == -999;
  //
  if (n0_depend + mu_depend + la_depend + ga_depend != 2)
    throw runtime_error(
        "Two (but only two) of n0, mu, la and ga must be NaN, "
        "to flag that these parameters are the ones dependent of "
        "mass content and mean particle size.");
  if (mu_depend || ga_depend)
    throw runtime_error(
        "Sorry, mu and la are not yet allowed to be a "
        "dependent parameter.");
  //
  const Index n0_fixed = (Index) !(n0_depend || std::isnan(n0));
  const Index mu_fixed = (Index) !(mu_depend || std::isnan(mu));
  const Index la_fixed = (Index) !(la_depend || std::isnan(la));
  const Index ga_fixed = (Index) !(ga_depend || std::isnan(ga));
  //
  if (nin + n0_fixed + mu_fixed + la_fixed + ga_fixed != 4)
    throw runtime_error(
        "This PSD has four free parameters. This means that "
        "the number\nof columns in *pnd_agenda_input* and the "
        "number of numerics\n(i.e. not -999 or NaN) and among "
        "the GIN arguments n0, mu, la and\nga must add up to "
        "four. And this was found not to be the case.");

  // Create vectors to hold the four MGD and the "extra" parameters
  Vector mgd_pars(4), ext_pars(2);
  ArrayOfIndex mgd_i_pai = {-1, -1, -1, -1};  // Position in pnd_agenda_input
  const ArrayOfIndex ext_i_pai = {0, 1};      // Position in pnd_agenda_input
  {
    Index nhit = 2;  // As mass and Dm always occupy first position
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
  ArrayOfIndex ext_do_jac = {0, 0};
  ArrayOfIndex mgd_i_jac = {
      -1, -1, -1, -1};                // Position among jacobian quantities
  ArrayOfIndex ext_i_jac = {-1, -1};  // Position among jacobian quantities
  //
  for (Index i = 0; i < ndx; i++) {
    if (dx2in[i] == 0)  // That is, mass is a derivative
    {
      ext_do_jac[0] = 1;
      ext_i_jac[0] = i;
    } else if (dx2in[i] == 1)  // That is, "something" is a derivative
    {
      ext_do_jac[1] = 1;
      ext_i_jac[1] = i;
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
    ext_pars[1] = pnd_agenda_input(ip, ext_i_pai[1]);
    if (ext_pars[1] <= 0) {
      ostringstream os;
      os << "Negative " << something << "found.\nThis is not allowed.";
      throw std::runtime_error(os.str());
    }
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
      if (picky) {
        ostringstream os;
        os << "Method called with a temperature of " << t << " K.\n"
           << "This is outside the specified allowed range: [ max(0.," << t_min
           << "), " << t_max << " ]";
        throw runtime_error(os.str());
      } else {
        continue;
      }  // If here, we are ready with this point!
    }

    // Derive the dependent parameters (see ATD)
    //
    Numeric mu1 = 0, mub1 = 0, eterm = 0, gterm = 0;
    Numeric scfac1 = 0, scfac2 = 0, gab = 0;
    //
    // *** Mean size ***
    if (something == "mean size") {
      if (n0_depend && la_depend) {
        mub1 = mgd_pars[1] + scat_species_b + 1;
        if (mub1 <= 0)
          throw runtime_error("Bad MGD parameter detected: mu + b + 1 <= 0");
        eterm = mub1 / mgd_pars[3];
        // Start by deriving la
        scfac2 = pow(eterm, mgd_pars[3]);
        mgd_pars[2] = scfac2 * pow(ext_pars[1], -mgd_pars[3]);
        // We can now derive n0
        gterm = tgamma(eterm);
        scfac1 =
            (mgd_pars[3] * pow(mgd_pars[2], eterm)) / (scat_species_a * gterm);
        mgd_pars[0] = scfac1 * ext_pars[0];
      } else {
        assert(0);
      }
    }

    // *** Median size ***
    else if (something == "median size") {
      if (n0_depend && la_depend) {
        mub1 = mgd_pars[1] + scat_species_b + 1;
        if (mub1 <= 0)
          throw runtime_error("Bad MGD parameter detected: mu + b + 1 <= 0");
        eterm = mub1 / mgd_pars[3];
        // Start by deriving la
        scfac2 = (mgd_pars[1] + 1 + scat_species_b - 0.327 * mgd_pars[3]) /
                 mgd_pars[3];
        mgd_pars[2] = scfac2 * pow(ext_pars[1], -mgd_pars[3]);
        // We can now derive n0
        gterm = tgamma(eterm);
        scfac1 =
            (mgd_pars[3] * pow(mgd_pars[2], eterm)) / (scat_species_a * gterm);
        mgd_pars[0] = scfac1 * ext_pars[0];
      } else {
        assert(0);
      }
    }

    // *** Mean particle size ***
    else if (something == "mean particle mass") {
      if (n0_depend && la_depend) {
        mu1 = mgd_pars[1] + 1;
        if (mu1 <= 0)
          throw runtime_error("Bad MGD parameter detected: mu + 1 <= 0");
        eterm = (mgd_pars[1] + scat_species_b + 1) / mgd_pars[3];
        gterm = tgamma(eterm);
        // Start by deriving la
        gab = mgd_pars[3] / scat_species_b;
        scfac2 = pow(scat_species_a * gterm / tgamma(mu1 / mgd_pars[3]), gab);
        mgd_pars[2] = scfac2 * pow(ext_pars[1], -gab);
        // We can now derive n0
        scfac1 =
            (mgd_pars[3] * pow(mgd_pars[2], eterm)) / (scat_species_a * gterm);
        mgd_pars[0] = scfac1 * ext_pars[0];
      } else {
        assert(0);
      }
    }

    else if (something == "Ntot") {
      if (n0_depend && la_depend) {
        mu1 = mgd_pars[1] + 1;
        if (mu1 <= 0)
          throw runtime_error("Bad MGD parameter detected: mu + 1 <= 0");
        eterm = (mgd_pars[1] + scat_species_b + 1) / mgd_pars[3];
        gterm = tgamma(eterm);
        // Start by deriving la
        gab = mgd_pars[3] / scat_species_b;
        scfac2 = pow(scat_species_a * gterm / tgamma(mu1 / mgd_pars[3]), gab);
        mgd_pars[2] = scfac2 * pow(ext_pars[1] / ext_pars[0], gab);
        // We can now derive n0
        scfac1 =
            (mgd_pars[3] * pow(mgd_pars[2], eterm)) / (scat_species_a * gterm);
        mgd_pars[0] = scfac1 * ext_pars[0];
      } else {
        assert(0);
      }
    }

    // String something not recognised
    else {
      assert(0);
    }

    // Now when all four MGD parameters are set, check that la and ga are OK
    if (mgd_pars[2] <= 0)
      throw runtime_error("Bad MGD parameter detected: la <= 0");
    if (mgd_pars[3] <= 0)
      throw runtime_error("Bad MGD parameter detected: ga <= 0");

    // Calculate PSD and derivatives
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

    // Derivatives for mass and something
    if (ext_do_jac[0] | ext_do_jac[1]) {
      // *** Mean size ***
      if (something == "mean size") {
        if (n0_depend && la_depend) {
          // Derivative with respect to mass
          if (ext_do_jac[0]) {
            dpsd_data_dx(ext_i_jac[0], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[0], ip, joker) *= scfac1;
          }
          // Derivative with respect to mean size
          if (ext_do_jac[1]) {
            // 1. Term associated with n0
            // Calculated as dpsd/dn0 * dn0/dla * dla/dXm
            dpsd_data_dx(ext_i_jac[1], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                ext_pars[0] * mgd_pars[3] * eterm *
                pow(mgd_pars[2], eterm - 1) / (scat_species_a * gterm);
            // 2. Term associated with la
            // Calculated as dpsd/dla * dla/dXm
            dpsd_data_dx(ext_i_jac[1], ip, joker) += jac_data(2, joker);
            // Apply dla/dXm to sum
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                -mgd_pars[3] * scfac2 * pow(ext_pars[1], -(mgd_pars[3] + 1));
          }
        } else {
          assert(0);
        }
      }

      // *** Median size ***
      else if (something == "median size") {
        if (n0_depend && la_depend) {
          // Derivative with respect to mass
          if (ext_do_jac[0]) {
            dpsd_data_dx(ext_i_jac[0], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[0], ip, joker) *= scfac1;
          }
          // Derivative with respect to median size
          if (ext_do_jac[1]) {
            // 1. Term associated with n0
            // Calculated as dpsd/dn0 * dn0/dla * dla/dXm
            dpsd_data_dx(ext_i_jac[1], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                ext_pars[0] * mgd_pars[3] * eterm *
                pow(mgd_pars[2], eterm - 1) / (scat_species_a * gterm);
            // 2. Term associated with la
            // Calculated as dpsd/dla * dla/dXm
            dpsd_data_dx(ext_i_jac[1], ip, joker) += jac_data(2, joker);
            // Apply dla/dXm to sum
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                -mgd_pars[3] * scfac2 * pow(ext_pars[1], -(mgd_pars[3] + 1));
          }
        } else {
          assert(0);
        }
      }

      // *** Mean particle size ***
      else if (something == "mean particle mass") {
        if (n0_depend && la_depend) {
          // Derivative with respect to mass
          if (ext_do_jac[0]) {
            dpsd_data_dx(ext_i_jac[0], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[0], ip, joker) *= scfac1;
          }
          // Derivative with respect to mean particle size
          if (ext_do_jac[1]) {
            // 1. Term associated with n0
            // Calculated as dpsd/dn0 * dn0/dla * dla/dMm
            dpsd_data_dx(ext_i_jac[1], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                ext_pars[0] * mgd_pars[3] * eterm *
                pow(mgd_pars[2], eterm - 1) / (scat_species_a * gterm);
            // 2. Term associated with la
            // Calculated as dpsd/dla * dla/dMm
            dpsd_data_dx(ext_i_jac[1], ip, joker) += jac_data(2, joker);
            // Apply dla/dMm to sum
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                scfac2 * (-mgd_pars[3] / scat_species_b) *
                pow(ext_pars[1], -(gab + 1));
          } else {
            assert(0);
          }
        }
      }

      else if (something == "Ntot") {
        if (n0_depend && la_depend) {
          // Term part of both derivatives
          const Numeric dn0dla = ext_pars[0] * mgd_pars[3] * eterm *
                                 pow(mgd_pars[2], eterm - 1) /
                                 (scat_species_a * gterm);
          // Derivative with respect to mass
          if (ext_do_jac[0]) {
            // Repeated term
            const Numeric dladw = scfac2 * pow(ext_pars[1], gab) *
                                  (-mgd_pars[3] / scat_species_b) *
                                  pow(ext_pars[0], -(gab + 1));
            // 1. Term associated with n0
            dpsd_data_dx(ext_i_jac[0], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[0], ip, joker) *= scfac1 + dn0dla * dladw;
            // 2. Term associated with la
            Vector term2 = jac_data(2, joker);
            term2 *= dladw;
            // Sum up
            dpsd_data_dx(ext_i_jac[0], ip, joker) += term2;
          }
          // Derivative with respect to Ntot
          if (ext_do_jac[1]) {
            // 1. Term associated with n0
            // Calculated as dpsd/dn0 * dn0/dla * dla/dNtot
            dpsd_data_dx(ext_i_jac[1], ip, joker) = jac_data(0, joker);
            dpsd_data_dx(ext_i_jac[1], ip, joker) *= dn0dla;
            // 2. Term associated with la
            // Calculated as dpsd/dla * dla/dNtot
            dpsd_data_dx(ext_i_jac[1], ip, joker) += jac_data(2, joker);
            // Apply dla/dNtot to sum
            dpsd_data_dx(ext_i_jac[1], ip, joker) *=
                scfac2 * pow(ext_pars[0], -gab) *
                (mgd_pars[3] / scat_species_b) * pow(ext_pars[1], gab - 1);
          }
        } else {
          assert(0);
        }
      }

      // String something not recognised
      else {
        assert(0);
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

void psd_mono_common(Matrix& psd_data,
                     Tensor3& dpsd_data_dx,
                     const String& type,
                     const Vector& pnd_agenda_input_t,
                     const Matrix& pnd_agenda_input,
                     const ArrayOfString& pnd_agenda_input_names,
                     const ArrayOfString& dpnd_data_dx_names,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const Index& species_index,
                     const Numeric& t_min,
                     const Numeric& t_max,
                     const Index& picky,
                     const Verbosity&) {
  // Standard checcks
  const Vector psd_size_grid(1, 0);  // As this WSV is not input for thse WSM
  START_OF_PSD_METHODS();

  // Extra checks for this PSD
  const Index nss = scat_meta.nelem();
  if (nss == 0) throw runtime_error("*scat_meta* is empty!");
  if (nss < species_index + 1) {
    ostringstream os;
    os << "Selected scattering species index is " << species_index
       << " but this "
       << "is not allowed since *scat_meta* has only " << nss << " elements.";
    throw runtime_error(os.str());
  }
  if (scat_meta[species_index].nelem() != 1) {
    ostringstream os;
    os << "This method only works with scattering species consisting of a\n"
       << "single element, but your data do not match this demand.\n"
       << "Selected scattering species index is " << species_index << ".\n"
       << "This species has " << scat_meta[species_index].nelem()
       << " elements.";
    throw runtime_error(os.str());
  }
  //
  if (pnd_agenda_input.ncols() != 1)
    throw runtime_error("*pnd_agenda_input* must have one column.");
  if (nsi != 1)
    throw runtime_error(
        "This method demands that length of "
        "*psd_size_grid* is 1.");

  // Extract particle mass
  Numeric pmass = 0;
  if (type == "mass") {
    pmass = scat_meta[species_index][0].mass;
  }

  for (Index ip = 0; ip < np; ip++) {
    // Extract the input variables
    Numeric x = pnd_agenda_input(ip, 0);
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if n==0 and no jacobians requested.
    if ((x == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      if (picky) {
        ostringstream os;
        os << "Method called with a temperature of " << t << " K.\n"
           << "This is outside the specified allowed range: [ max(0.," << t_min
           << "), " << t_max << " ]";
        throw runtime_error(os.str());
      } else {
        continue;
      }  // If here, we are ready with this point!
    }

    // Set PSD
    //
    if (type == "ntot") {
      psd_data(ip, 0) = x;
      //
      if (ndx) {
        dpsd_data_dx(0, ip, 0) = 1;
      }
    } else if (type == "mass") {
      psd_data(ip, 0) = x / pmass;
      //
      if (ndx) {
        dpsd_data_dx(0, ip, 0) = 1 / pmass;
      }
    } else {
      assert(0);
    }
  }
}

void psd_rain_W16(Vector& psd, const Vector& diameter, const Numeric& rwc) {
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  // skip calculation if RWC is 0.0
  if (rwc == 0.0) {
    return;
  }
  assert(rwc > 0.);

  // a and b relates N0 to lambda N0 = a*lambda^b
  Numeric a = 0.000141;
  Numeric b = 1.49;
  Numeric c1 = DENSITY_OF_WATER * PI / 6;
  Numeric base = c1 / rwc * a * tgamma(4);
  Numeric exponent = 1. / (4 - b);
  Numeric lambda = pow(base, exponent);
  Numeric N0 = a * pow(lambda, b);

  //psd_general_MGD( psd, N0, 0, lambda, 1 );
  N0 *= 1e8;
  lambda *= 100;
  for (Index iD = 0; iD < nD; iD++) {
    psd[iD] = N0 * exp(-lambda * diameter[iD]);
  }
}

void psd_mgd_smm_common(Matrix& psd_data,
                        Tensor3& dpsd_data_dx,
                        const String& psd_name,
                        const Vector& psd_size_grid,
                        const Vector& pnd_agenda_input_t,
                        const Matrix& pnd_agenda_input,
                        const ArrayOfString& pnd_agenda_input_names,
                        const ArrayOfString& dpnd_data_dx_names,
                        const Numeric& scat_species_a,
                        const Numeric& scat_species_b,
                        const Numeric& n_alpha_in,
                        const Numeric& n_b_in,
                        const Numeric& mu_in,
                        const Numeric& gamma_in,
                        const Numeric& t_min,
                        const Numeric& t_max,
                        const Index& picky,
                        const Verbosity&) {
  // Standard checks
  START_OF_PSD_METHODS();

  // All PSDs are defined in terms of hydrometeor mass content
  if (pnd_agenda_input.ncols() != 1)
    throw runtime_error("*pnd_agenda_input* must have one column.");

  // Extra checks for rain PSDs which should be consistent with
  // spherical liquid drops
  if (psd_name == "Abel12" || psd_name == "Wang16"){
    if (scat_species_b < 2.9 || scat_species_b > 3.1) {
      ostringstream os;
      os << "This PSD treats rain, using Dveq as size grid.\n"
	 << "This means that *scat_species_b* should be close to 3,\n"
	 << "but it is outside of the tolerated range of [2.9,3.1].\n"
	 << "Your value of *scat_species_b* is: " << scat_species_b;
      throw runtime_error(os.str());
    }
    if (scat_species_a < 500 || scat_species_a > 540) {
      ostringstream os;
      os << "This PSD treats rain, using Dveq as size grid.\n"
	 << "This means that *scat_species_a* should be close to 520,\n"
	 << "but it is outside of the tolerated range of [500,540].\n"
	 << "Your value of *scat_species_a* is: " << scat_species_a;
      throw runtime_error(os.str());
    }
  }
  // Extra checks for graupel/hail PSDs which assume constant effective density
  //
  if (psd_name == "Field19"){
    if (scat_species_b < 2.8 || scat_species_b > 3.2) {
      ostringstream os;
      os << "This PSD treats graupel, assuming a constant effective density.\n"
	 << "This means that *scat_species_b* should be close to 3,\n"
	 << "but it is outside of the tolerated range of [2.8,3.2].\n"
	 << "Your value of *scat_species_b* is: " << scat_species_b;
      throw runtime_error(os.str());
    }
  }

  for (Index ip = 0; ip < np; ip++) {
    // Extract the input variables
    Numeric water_content = pnd_agenda_input(ip, 0);
    Numeric t = pnd_agenda_input_t[ip];

    // No calc needed if water_content==0 and no jacobians requested.
    if ((water_content == 0.) && (!ndx)) {
      continue;
    }  // If here, we are ready with this point!

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      if (picky) {
        ostringstream os;
        os << "Method called with a temperature of " << t << " K.\n"
           << "This is outside the specified allowed range: [ max(0.," << t_min
           << "), " << t_max << " ]";
        throw runtime_error(os.str());
      } else {
        continue;
      }  // If here, we are ready with this point!
    }

    // Negative wc?
    Numeric psd_weight = 1.0;
    if (water_content < 0) {
      psd_weight = -1.0;
      water_content *= -1.0;
    }

    // PSD settings for different parametrizations
    Numeric gamma = 0.0;
    Numeric n_alpha = 0.0;
    Numeric n_b = 0.0;
    Numeric mu = 0.0;
    if (psd_name == "Abel12"){
      n_alpha = 0.22;
      n_b = 2.2;
      mu = 0.0;
      gamma = 1.0;
    }
    else if (psd_name == "Wang16"){
      // Wang 16 parameters converted to SI units
      n_alpha = 14.764;
      n_b = 1.49;
      mu = 0.0;
      gamma = 1.0;
    }
    else if (psd_name == "Field19"){
      n_alpha = 7.9e9;
      n_b = -2.58;
      mu = 0.0;
      gamma = 1.0;
    }
    else if (psd_name == "generic"){
      n_alpha = n_alpha_in;
      n_b = n_b_in;
      mu = mu_in;
      gamma = gamma_in;
    }
    else {
      assert(0);
    }

    // Calculate PSD
    // Calculate lambda for modified gamma distribution from mass density
    Numeric k = (scat_species_b + mu + 1 - gamma)/gamma;
    Numeric expo = 1.0 / (n_b - k - 1);
    Numeric denom = scat_species_a * n_alpha * tgamma(k + 1);
    Numeric lam = pow(water_content*gamma/denom, expo);
    Numeric n_0 = n_alpha * pow(lam, n_b);
    Vector psd_1p(nsi);
    Matrix jac_data(4, nsi);

    psd_1p = 0.0;
    if (water_content != 0) {
      mgd_with_derivatives(psd_1p, jac_data, psd_size_grid, n_0, mu, lam, gamma,
			   true, // n_0 jacobian
			   false,// mu jacobian
			   true, // lambda jacobian
			   false); // gamma jacobian
    } else {
      assert(0);
    }
    //
    for (Index i = 0; i < nsi; i++) {
      psd_data(ip, i) = psd_weight * psd_1p[i];
    }

    // Calculate derivative with respect to water content
    if (ndx) {
      const Numeric dlam_dwc = pow(gamma/denom, expo) * expo * pow(water_content, expo-1);
      const Numeric dn_0_dwc = n_alpha * n_b * pow(lam, n_b-1) * dlam_dwc;
      for (Index i = 0; i < nsi; i++) {
        dpsd_data_dx(0, ip, i) = psd_weight * (jac_data(0,i)*dn_0_dwc +
					       jac_data(2,i)*dlam_dwc);
      }
    }
  }
}

void psd_snow_F07(Vector& psd,
                  const Vector& diameter,
                  const Numeric& swc,
                  const Numeric& t,
                  const Numeric alpha,
                  const Numeric beta,
                  const String& regime) {
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  assert(t > 0.);

  // skip calculation if SWC is 0.0
  if (swc == 0.0) {
    return;
  }
  assert(swc > 0.);

  Numeric An, Bn, Cn;
  Numeric M2, Mn, M2Mn;
  Numeric base, pp;
  Numeric x, phi23;
  //Numeric dN;

  Vector q(5);

  assert((regime == "TR") || (regime == "ML"));

  //factors of phi23
  if (regime == "TR")
    q = {152., -12.4, 3.28, -0.78, -1.94};
  else  // if( regime=="ML" )
    q = {141., -16.8, 102., 2.07, -4.82};

  //Factors of factors of the moment estimation parametrization
  Vector Aq{13.6, -7.76, 0.479};
  Vector Bq{-0.0361, 0.0151, 0.00149};
  Vector Cq{0.807, 0.00581, 0.0457};

  //convert T from Kelvin to Celsius
  Numeric Tc = t - 273.15;

  // estimate second moment
  M2 = swc / alpha;
  if (beta != 2) {
    // calculate factors of the moment estimation parametrization
    An = exp(Aq[0] + Aq[1] * beta + Aq[2] * pow(beta, 2));
    Bn = Bq[0] + Bq[1] * beta + Bq[2] * pow(beta, 2);
    Cn = Cq[0] + Cq[1] * beta + Cq[2] * pow(beta, 2);

    base = M2 * exp(-Bn * Tc) / An;
    pp = 1. / (Cn);

    M2 = pow(base, pp);
  }

  // order of the moment parametrization
  const Numeric n = 3;

  // calculate factors of the moment estimation parametrization
  An = exp(Aq[0] + Aq[1] * n + Aq[2] * pow(n, 2));
  Bn = Bq[0] + Bq[1] * n + Bq[2] * pow(n, 2);
  Cn = Cq[0] + Cq[1] * n + Cq[2] * pow(n, 2);

  // moment parametrization
  Mn = An * exp(Bn * Tc) * pow(M2, Cn);

  M2Mn = pow(M2, 4.) / pow(Mn, 3.);

  for (Index iD = 0; iD < nD; iD++) {
    // define x
    x = diameter[iD] * M2 / Mn;

    // characteristic function
    phi23 = q[0] * exp(q[1] * x) + q[2] * pow(x, q[3]) * exp(q[4] * x);

    // distribution function
    //dN = phi23*M2Mn;

    //if ( !std::isnan(psd[dN]) )
    //  psd[iD] = dN;

    // set psd directly. Non-NaN should be (and is, hopefully) ensured by
    // checks above (if we encounter further NaN, that should be handled above.
    // which intermediate quantities make problems? at which parametrisation
    // values, ie WC, T, alpha, beta?).
    psd[iD] = phi23 * M2Mn;
  }
}

void psd_SB06(Vector& psd,
              Matrix& dpsd,
              const Vector& mass,
              const Numeric& N_tot,
              const Numeric& WC,
              const String& hydrometeor_type) {
  Numeric N0;
  Numeric Lambda;
  Numeric arg1;
  Numeric arg2;
  Numeric brk;
  Numeric mu;
  Numeric gamma;
  Numeric xmin;
  Numeric xmax;
  Numeric M0min;
  Numeric M0max;
  Numeric M0;
  Numeric M1;
  Numeric c1;
  Numeric c2;
  Numeric L1;
  Numeric mMu;
  Numeric mGamma;
  Numeric brkMu1;

  // Get the coefficients for the right hydrometeor
  if (hydrometeor_type == "cloud_ice")  //Cloud ice water
  {
    mu = 0.;
    gamma = 1. / 3.;
    xmin = 1e-12;
    xmax = 1e-5;
  } else if (hydrometeor_type == "rain")  //Rain
  {
    mu = 0.;
    gamma = 1. / 3.;
    xmin = 2.6e-10;
    xmax = 3e-6;
  } else if (hydrometeor_type == "snow")  //Snow
  {
    mu = 0.;
    gamma = 1. / 2.;
    xmin = 1e-10;
    xmax = 2e-5;
  } else if (hydrometeor_type == "graupel")  //Graupel
  {
    mu = 1.;
    gamma = 1. / 3.;
    xmin = 1e-9;
    xmax = 5e-4;
  } else if (hydrometeor_type == "hail")  //Hail
  {
    mu = 1.;
    gamma = 1. / 3.;
    xmin = 2.6e-10;
    xmax = 5e-4;
  } else if (hydrometeor_type == "cloud_water")  //Cloud liquid water
  {
    mu = 1;
    gamma = 1;
    xmin = 4.2e-15;
    xmax = 2.6e-10;
  } else {
    ostringstream os;
    os << "You use a wrong tag! ";
    throw runtime_error(os.str());
  }

  M0 = N_tot;
  M1 = WC;

  Index nD = mass.nelem();
  psd.resize(nD);
  psd = 0.;

  dpsd.resize(nD, 2);
  dpsd = 0.;

  if (M1 > 0.0) {
    // lower and upper limit check is taken from the ICON code of the two moment
    //scheme

    M0max = M1 / xmax;
    M0min = M1 / xmin;

    //check lower limit of the scheme
    if (M0 > M0min) {
      M0 = M0min;
    }

    //check upper limit of the scheme
    if (M0 < M0max) {
      M0 = M0max;
    }

    //arguments for Gamma function
    arg2 = (mu + 2) / gamma;
    arg1 = (mu + 1) / gamma;

    // results of gamma function
    c1 = tgamma(arg1);
    c2 = tgamma(arg2);

    // variable to shorten the formula
    brk = M0 / M1 * c2 / c1;
    brkMu1 = pow(brk, (mu + 1));

    //Lambda (parameter for modified gamma distribution)
    Lambda = pow(brk, gamma);

    L1 = pow(Lambda, arg1);

    //N0
    N0 = M0 * gamma / tgamma(arg1) * L1;

    // Calculate distribution function
    for (Index iD = 0; iD < nD; iD++) {
      //Distribution function
      psd[iD] = mod_gamma_dist(mass[iD], N0, Lambda, mu, gamma);

      if (std::isnan(psd[iD])) psd[iD] = 0.0;
      if (std::isinf(psd[iD])) psd[iD] = 0.0;

      //Calculate derivatives analytically
      mMu = pow(mass[iD], mu);
      mGamma = pow(mass[iD], gamma);

      // dpsd/dM1
      dpsd(iD, 0) = gamma / c1 * M0 / M1 * mMu * exp(-Lambda * mGamma) *
                    brkMu1 * (-1 - mu + gamma * mGamma * Lambda);

      // dpsd/dM0
      dpsd(iD, 1) = -gamma / c1 * mMu * exp(-Lambda * mGamma) * brkMu1 *
                    (-2 - mu - gamma * mGamma * Lambda);
    }
  } else {
    return;
  }
}

void psd_MY05(Vector& psd,
              Matrix& dpsd,
              const Vector& diameter_max,
              const Numeric N_tot,
              const Numeric WC,
              const String psd_type) {
  Numeric N0;
  Numeric Lambda;
  Numeric arg1;
  Numeric arg2;
  Numeric temp;
  Numeric mu;
  Numeric gamma;
  Numeric alpha;
  Numeric beta;
  Numeric M0;
  Numeric M1;
  Numeric c1;
  Numeric c2;
  Numeric Lmg;
  Numeric DMu;
  Numeric DGamma;

  // Get the coefficients for the right hydrometeor
  if (psd_type == "cloud_ice")  //Cloud ice water
  {
    mu = 0.;
    gamma = 1.;
    alpha = 440.;  //[kg]
    beta = 3;
  } else if (psd_type == "rain")  //Rain
  {
    mu = 0.;
    gamma = 1;
    alpha = 523.5988;  //[kg]
    beta = 3;
  } else if (psd_type == "snow")  //Snow
  {
    mu = 0.;
    gamma = 1;
    alpha = 52.35988;  //[kg]
    beta = 3;
  } else if (psd_type == "graupel")  //Graupel
  {
    mu = 0.;
    gamma = 1;
    alpha = 209.4395;  //[kg]
    beta = 3;
  } else if (psd_type == "hail")  //Hail
  {
    mu = 0.;
    gamma = 1;
    alpha = 471.2389;  //[kg]
    beta = 3;
  } else if (psd_type == "cloud_water")  //Cloud liquid water
  {
    mu = 1;
    gamma = 1;
    alpha = 523.5988;  //[kg]
    beta = 3;
  } else {
    ostringstream os;
    os << "You use a wrong tag! ";
    throw runtime_error(os.str());
  }

  M0 = N_tot;
  M1 = WC;

  Index nD = diameter_max.nelem();
  psd.resize(nD);
  psd = 0.;

  dpsd.resize(nD, 2);
  dpsd = 0.;

  if (M1 > 0.0 && M0 > 0) {
    //arguments for Gamma function
    arg2 = (mu + beta + 1) / gamma;
    arg1 = (mu + 1) / gamma;

    // results of gamma function
    c1 = tgamma(arg1);
    c2 = tgamma(arg2);

    //base of lambda
    temp = alpha * M0 / M1 * c2 / c1;

    //Lambda (parameter for modified gamma distribution)
    Lambda = pow(temp, gamma / beta);

    Lmg = pow(Lambda, arg1);

    //N0
    N0 = M0 * gamma / c1 * Lmg;

    //Distribution function

    // Calculate distribution function
    for (Index iD = 0; iD < nD; iD++) {
      psd[iD] = mod_gamma_dist(diameter_max[iD], N0, Lambda, mu, gamma);

      if (std::isnan(psd[iD])) psd[iD] = 0.0;
      if (std::isinf(psd[iD])) psd[iD] = 0.0;

      //Calculate derivatives analytically
      DMu = pow(diameter_max[iD], mu);
      DGamma = pow(diameter_max[iD], gamma);

      // dpsd/dM1
      dpsd(iD, 0) = (DMu * exp(-DGamma * Lambda) * gamma * M0 * Lmg *
                     (-1 - mu + DGamma * gamma * Lambda) / (M1 * beta * c1));

      // dpsd/dM0
      dpsd(iD, 1) = (DMu * exp(-DGamma * Lambda) * gamma * Lmg *
                     (1 + beta + mu - DGamma * gamma * Lambda) / (beta * c1));
    }

  } else {
    return;
  }
}

Numeric dm_from_iwc_n0(Numeric iwc, Numeric n0, Numeric rho) {
  if (iwc == 0.0) {
    return 1e-9;
  } else {
    return pow(256.0 * iwc / PI / rho / n0, 0.25);
  }
}

Numeric n0_from_iwc_dm(Numeric iwc, Numeric dm, Numeric rho) {
  if (dm > 1e-9) {
    return 256.0 * iwc / PI / rho / pow(dm, 4.0);
  } else {
    return 0.0;
  }
}

Numeric n0_from_t(Numeric t) { return exp(-0.076586 * (t - 273.15) + 17.948); }
