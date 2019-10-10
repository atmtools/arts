/* Copyright (C) 2011-2017 Jana Mendrok <jana.mendrok@gmail.com>
                      
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
  \author Jana Mendrok, Daniel Kreyling, Manfred Brath, Patrick Eriksson
  \date   2017-08-01
  
  \brief  Internal functions for microphysics calculations (size distributions etc.)
*/

#include "microphysics.h"

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
#include "psd.h"
#include "rng.h"
#include "sorting.h"

/*! Derives a and b for relationship mass = a * x^b

    The parameters a and b are derived by a fit including all data inside the
    size range [x_fit_start,x_fit_end].

    The vector x must have been checked to have at least 2 elements.

    An error is thrown if less than two data points are found inside
    [x_fit_start,x_fit_end].

    \return a           Derived a parameter.
    \return b           Derived b parameter.
    \param  x           Size grid
    \param  mass        Particle masses
    \param  x_fit_start Start point of x-range to use for fitting
    \param  x_fit_end   Endpoint of x-range to use for fitting
  
  \author Jana Mendrok, Patrick Eriksson
  \date 2017-10-18

*/
void derive_scat_species_a_and_b(Numeric& a,
                                 Numeric& b,
                                 const Vector& x,
                                 const Vector& mass,
                                 const Numeric& x_fit_start,
                                 const Numeric& x_fit_end) {
  const Index nse = x.nelem();
  assert(nse > 1);

  ArrayOfIndex intarr_sort, intarr_unsort(0);
  Vector x_unsorted(nse), m_unsorted(nse);
  Vector q;
  Index nsev = 0;

  for (Index i = 0; i < nse; i++) {
    if (std::isnan(x[i]))
      throw runtime_error("NaN found in selected size grid data.");
    if (std::isnan(mass[i]))
      throw runtime_error("NaN found among particle mass data.");

    if (x[i] >= x_fit_start && x[i] <= x_fit_end) {
      x_unsorted[nsev] = x[i];
      m_unsorted[nsev] = mass[i];
      nsev += 1;
    }
  }

  if (nsev < 2)
    throw runtime_error(
        "Less than two size points found in the range "
        "[x_fit_start,x_fit_end]. It is then not possible "
        "to determine the a and b parameters.");

  get_sorted_indexes(intarr_sort, x_unsorted[Range(0, nsev)]);
  Vector log_x(nsev), log_m(nsev);

  for (Index i = 0; i < nsev; i++) {
    log_x[i] = log(x_unsorted[intarr_sort[i]]);
    log_m[i] = log(m_unsorted[intarr_sort[i]]);
  }

  linreg(q, log_x, log_m);
  a = exp(q[0]);
  b = q[1];
}

/*! Calculates particle size distribution using H11 parametrization.
 *  Each diameter of the scattering elements is a node in the distribution.
 *  One call of this function calculates number density for one scattering
 *  element.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Daniel Kreyling
  \date 2011-10-28

*/
Numeric IWCtopnd_H11(const Numeric diameter_max, const Numeric t) {
  Numeric dNdD;
  Numeric la;
  Numeric mu;

  // convert m to cm
  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t - 273.15;
  //choose parametrization depending on T
  if (T >= -56.) {
    la = 12.13 * exp(-0.055 * T);
  } else {
    la = 0.83 * exp(-0.103 * T);
  }
  if (T >= -68.) {
    mu = -0.57 - 0.028 * T;
  } else {
    mu = -30.93 - 0.472 * T;
  }

  //Distribution function H11

  dNdD = pow(dmax, mu) * exp(-la * dmax);

  if (std::isnan(dNdD)) dNdD = 0.0;
  return dNdD;
}


/*! Calculates particle size and shape distribution using H13 parametrization.
 *  Each diameter of the scattering elements is a node in the distribution.
 *  One call of this function calculates number density for one scattering
 *  element.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Johan Strandgren  
  \date 2013-08-26

*/
Numeric IWCtopnd_H13Shape(const Numeric diameter_max, const Numeric t) {
  Numeric dNdD;
  Numeric la;
  Numeric mu;
  // convert m to cm

  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t - 273.15;
  //choose parametrization depending on T
  if (T >= -58.) {
    la = 9.88 * exp(-0.060 * T);
  } else {
    la = 0.75 * exp(-0.1057 * T);
  }
  if (T >= -61.) {
    mu = -0.59 - 0.030 * T;
  } else {
    mu = -14.09 - 0.248 * T;
  }

  //Distribution function H13Shape

  dNdD = pow(dmax, mu) * exp(-la * dmax);

  if (std::isnan(dNdD)) dNdD = 0.0;
  return dNdD;
}

/*! Calculates area ratio using H13 shape parametrization.
 *  Each scattering element is a node in the aspect ratio distribution.
 *  One call of this function calculates one aspect ratio.  

    \return dNdD particle number density per diameter interval [#/m3/m]
          
    \param diameter_max  maximum diameter of scattering scattering element [m]
    \param t     atmospheric temperature [K]
  
  \author Johan Strandgren  
  \date 2013-08-26

*/
Numeric area_ratioH13(const Numeric diameter_max, const Numeric t) {
  Numeric Ar;
  Numeric alpha;
  Numeric beta;

  // convert m to cm
  Numeric dmax = diameter_max * 1e2;
  //convert T from Kelvin to Celsius
  Numeric T = t - 273.15;
  //Parameterize for all temperatures

  alpha = 0.25 * exp(0.0161 * T);

  beta = -0.25 + 0.0033 * T;

  // Area ratio function depending on temperature

  Ar = alpha * pow(dmax, beta);

  if (std::isnan(Ar)) Ar = 0.0;
  return Ar;
}

/*! Calculates particle size distribution for liquid water clouds using a gamma
 *  parametrization by Hess et al., 1998 (continental stratus).
 *  Each scattering element is a node in the distribution.
 *  One call of this function calculates one particle number density.
 *  Implicitly assumes particles of liquid water.

	\return n  particle number density per radius interval [#/m3*m]
         
	\param lwc atmospheric liquid water content [kg/m3]
	\param r radius of scattering scattering element [m]
  
  \author Daniel Kreyling
  \date 2010-12-16

*/
Numeric LWCtopnd(const Numeric lwc,    //[kg/m^3]
                 const Numeric radius  // [m]
) {
  // skip calculation if LWC is 0.0
  if (lwc == 0.0) {
    return 0.0;
  }
  Numeric rc = 4.7 * 1e-6;  //[um]->[m]
  Numeric alpha = 5.0;
  Numeric gam = 1.05;

  Numeric a4g = (alpha + 4.) / gam;
  Numeric B = (alpha / gam) / pow(rc, gam);
  Numeric A =
      0.75 / PI * lwc / DENSITY_OF_WATER * gam * pow(B, a4g) / tgamma(a4g);
  Numeric dNdr =
      A * (pow(radius, alpha) * exp(-B * pow(radius, gam)));  // [#/m3/m]

  /* alternative implementation
  Numeric rc = 4.7; //micron
	Numeric alpha = 5.0;
	Numeric gam = 1.05;
	
	Numeric B=(alpha/gam)/pow(rc,gam); 
	Numeric A=gam*pow(B,((alpha+1)/gam))/tgamma((alpha+1)/gam);
	Numeric dNdr=A*(pow(radius*1e6,alpha)*exp(-B*pow(radius*1e6,gam)))*1e6; // [#/m3/m]
*/

  if (std::isnan(dNdr)) dNdr = 0.0;
  return dNdr;
}

/*! Calculates the particle number density field for Cloud liquid water according
 *  the modified gamma distribution for cloud water inside Geer and Baordo (2014),
 *  see table A1
 *  One call of this function calculates one particle number density.
 *  Assumptions are: density of particles is constant and particle shape is sphere.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d volume equivalent diameter of scattering particle [m]
 \param lwc liquid water content [kg/m^3]
 
 \author Manfred Brath
 \date 2014-11-28
 
 */
Numeric LWCtopnd_MGD_LWC(const Numeric d,
                         const Numeric rho,
                         const Numeric lwc) {
  Numeric dN;
  Numeric N0;

  // coefficients of modified gamma distri
  const Numeric mu = 2;           //
  const Numeric lambda = 2.13e5;  //
  const Numeric gamma = 1;
  const Numeric fc = 1.4863e30;

  // N0
  N0 = fc * lwc / rho;  //[m^-4]

  //Distribution function
  dN = N0 * pow(d, mu) * exp(-lambda * pow(d, gamma));

  if (std::isnan(dN)) dN = 0.0;
  return dN;
}

/*! Calculates the particle number density field for Cloud ice according
 *  the modified gamma distribution for cloud ice inside Geer and Baordo (2014),
 *  see table A1
 *  One call of this function calculates one particle number density.
 *  Assumptions are: density of particles is constant and particle shape is sphere.
 
 \return dN particle number density per diameter interval [#/m3/m]
 
 \param d volume equivalent diameter of scattering particle [m]
 \param iwc ice water content [kg/m^3]
 
 \author Manfred Brath
 \date 2014-11-28
 
 */
Numeric IWCtopnd_MGD_IWC(const Numeric d,
                         const Numeric rho,
                         const Numeric iwc) {
  Numeric dN;
  Numeric N0;

  // coefficients of modified gamma distri
  const Numeric mu = 2;           //
  const Numeric lambda = 2.05e5;  //
  const Numeric gamma = 1;
  const Numeric fc = 1.1813e30;

  // N0
  N0 = fc * iwc / rho;  //[m^-4]

  //Distribution function
  dN = N0 * pow(d, mu) * exp(-lambda * pow(d, gamma));

  if (std::isnan(dN)) dN = 0.0;
  return dN;
}

/*! Calculates particle size distribution using MP48 parametrization for
 *  precipitating hydrometeors (rain, snow, ...).
 *  Each scattering element is a node in the distribution.
 *  One call of this function calculates one particle number density.  

    \return dNdD particle number density per diameter interval [#/m3*m]
          
    \param PR    precipitation rate [mm/hr]
    \param diameter_melted_equivalent melted equivalent diameter of scattering scattering element [m]
 
  \author Jana Mendrok
  \date 2012-04-04

*/
Numeric PRtopnd_MP48(const Numeric PR,
                     const Numeric diameter_melted_equivalent) {
  // skip calculation if PR is 0.0
  if (PR == 0.0) {
    return 0.0;
  }

  Numeric N0 = 0.08 * 1e8;  // [#/cm3/cm] converted to [#/m3/m]
  Numeric lambda =
      41. * 1e2 * pow(PR, -0.21);  // [1/cm] converted to [1/m] to fit
                                   // diameter_melted_equivalent [m]

  Numeric n = N0 * exp(-lambda * diameter_melted_equivalent);
  return n;
}
