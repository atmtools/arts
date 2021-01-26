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

