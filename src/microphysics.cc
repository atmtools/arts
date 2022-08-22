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
#include "arts_constants.h"
#include "arts_conversions.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric DENSITY_OF_ICE=Constant::density_of_ice_at_0c;
inline constexpr Numeric DENSITY_OF_WATER=Constant::denity_of_water_at_4c;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);


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


Numeric asymmetry_parameter(ConstVectorView sa_grid,
                            ConstVectorView pfun)
{
  const Index n = sa_grid.nelem();
  
  ARTS_ASSERT(abs(sa_grid[0]-0.0) < 1.0e-3);
  ARTS_ASSERT(abs(sa_grid[n-1]-180.0) < 1.0e-3);
  ARTS_ASSERT(pfun.nelem() == n);

  Vector sa = sa_grid;
  sa *= DEG2RAD;

  // Sine and cosine of scattering angle
  Vector sterm = sa;
  transform(sterm, sin, sterm);
  Vector cterm = sa;
  transform(cterm, cos, cterm);

  // Functions to integrate 
  Vector f1(n), f2(n);
  for (Index i=0; i<n; ++i) {
    f1[i] = sterm[i] * pfun[i];
    f2[i] = cterm[i] * f1[i];
  }

  // We skip some 2, 4 and pi, as they all cancel in the end
  const Numeric normfac = trapz(sa, f1);

  if (normfac < 1e-12) {
    return 0.0;  // If very small scattering cross-section, set to zero
  } else {       // to avoid numerical issues
    return trapz(sa, f2) / normfac;
  }
}


void derive_scat_species_a_and_b(Numeric& a,
                                 Numeric& b,
                                 const Vector& x,
                                 const Vector& mass,
                                 const Numeric& x_fit_start,
                                 const Numeric& x_fit_end) {
  const Index nse = x.nelem();
  ARTS_ASSERT(nse > 1);

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

