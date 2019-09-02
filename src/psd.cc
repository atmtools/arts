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

/*!
  \file   microphysics.cc
  \author Jana Mendrok, Patrick Eriksson
  \date   2017-11-05
  
  \brief  Internal functions of PSD type
*/

#include "psd.h"

extern const Numeric PI;
extern const Numeric DENSITY_OF_ICE;
extern const Numeric DENSITY_OF_WATER;

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

/*! Calculates particle size distribution of cloud ice using MH97 parametrization.
 *  
 *  Handles a vector of sizes at a time. Implicitly assumes particles of water
 *  ice. Strictly requires IWC and T to be positive, i.e. calling method needs
 *  to ensure this.
 *  
 *  Adapted from the 'old' IWCtopnd_MH97.

    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (supposed to be mass (aka
                       volume) equivalent diameter of pure ice particle) [m]
    \param iwc     atmospheric ice water content [kg/m3]
    \param t       atmospheric temperature [K]
    \param noisy   flag whether to add noise onto PSD parameters according to
                     their reported error statistics
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2017-06-07

*/
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

/*! Calculates particle size distribution of (stratiform) rain using Abel12
 *  parametrization.
 *  
 *  Uses rain water water content. PSD follows an exponential distribution.
 *  Handles a vector of sizes at a time.
 *  
 *  Reference: Abel and Boutle, 2012, "An improved representation of the 
 *  raindrop size distribution for single-moment microphysics schemes".
 *  Ported from CloudArts matlab implementation.

    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (volume equivalent
                       diameter) [m]
    \param rwc     atmospheric rain water content [kg/m3]
  
  \author Bengt Rydberg and Patrick Eriksson
  \date 2017-11-05

*/
void psd_rain_A12(Vector& psd, const Vector& diameter, const Numeric& rwc) {
  Index nD = diameter.nelem();
  psd.resize(nD);
  psd = 0.;

  // skip calculation if RWC is 0.0
  if (rwc == 0.0) {
    return;
  }
  assert(rwc > 0.);

  const Numeric x1 = 0.22;
  const Numeric x2 = 2.20;
  const Numeric c1 = DENSITY_OF_WATER * PI / 6;
  const Numeric lambda = pow(c1 / rwc * x1 * tgamma(4), 1 / (4 - x2));
  const Numeric N0 = x1 * pow(lambda, x2);

  mgd(psd, diameter, N0, 0, lambda, 1);
}

/*! Calculates particle size distribution of (stratiform) rain using Wang16
 *  parametrization.
 *  
 *  Uses rain water water content. PSD follows an exponential distribution.
 *  Handles a vector of sizes at a time.
 *  
 *  Reference: Wang et al., 2016, "Investigation of liquid cloud microphysical
 *  properties of deep convective systems: 1. Parameterization raindrop size
 *  distribution and its application for stratiform rain estimation".
 *  Ported from CloudArts matlab implementation.

    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (volume equivalent
                       diameter) [m]
    \param rwc     atmospheric rain water content [kg/m3]
  
  \author Jana Mendrok, Bengt Rydberg
  \date 2017-06-07

*/
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

/*! Calculates particle size distribution of tropical or midlatitude snow using
 *  F07 parametrization.
 *
 *  Handles a vector of sizes at a time.
 *  Strictly requires SWC and T to be positive and regime to be either "TR" or
 *  "ML", i.e. calling methods need to ensure these.
 *  No further limitations on the allowed temperatures here. Strictly valid it's
 *  only within -60<=t<=0C, the measured t-range the parametrization is based
 *  on. However, this is left to be handled by the calling methods.
 *
 *  Adapted from the 'old' IWCtopnd_F07TR/ML.
 
    \param psd     particle number density per size interval [#/m3*m]
    \param diameter  size of the scattering elements (supposed to be maximum
                       diameter of the ice particles) [m]
    \param swc     atmospheric snow water content [kg/m^3]
    \param t       atmospheric temperature [K]
    \param alpha   mass-dimension relationship scaling factor
                     (m=alpha*(Dmax/D0)^beta) [kg]
    \param beta    mass-dimension relationship exponent [-]
    \param regime  parametrization regime to apply (TR=tropical, ML=midlatitude)
 
 \author Manfred Brath, Jana Mendrok
 \date 2017-06-13
 
 */
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
