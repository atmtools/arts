/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   lineshapes.cc
  \brief  Stuff related to lineshape functions.

  This file contains both the lineshape functions themselves and the
  function define_lineshape_data which sets the lineshape lookup
  data. 

  \author Stefan Buehler
  \date   2000-08-21
*/

#include "arts.h"
#include "vecmat.h"
#include "absorption.h"

/*! The Lorentz line shape. This is a quick and dirty implementation.

    \retval ls     The shape function.
    \param  f0     Line center frequency.
    \param  gamma  The pressure broadening parameter.
    \param  sigma  The Doppler broadening parameter. (Not used.)
    \param  f_mono The frequency grid.

    \author Stefan Buehler 16.06.2000 */
void lineshape_lorentz(VECTOR&       ls,
		       Numeric	      f0,
		       Numeric       gamma,
		       Numeric       sigma,
		       const VECTOR& f_mono)
{
  // FIXME: Maybe try if call by reference is faster for f0 and gamma?

  // 1/PI:
  extern const Numeric PI;
  static const Numeric invPI = 1. / PI;

  assert( ls.dim() == f_mono.dim() );

  for ( size_t i=0; i<f_mono.dim(); ++i )
    {
      ls[i] = invPI * gamma / ( pow( f_mono[i]-f0, 2) + pow(gamma,2) );
    }
}


/*! The lookup data for the different lineshapes. */
ARRAY<LineshapeRecord> lineshape_data;


void define_lineshape_data()
{
  // Initialize to empty, just in case.
  lineshape_data.clear();

  lineshape_data.push_back
    (LineshapeRecord
     ("Lorentz",
      "The Lorentz line shape. This is a quick and dirty implementation.",
      -1,
      lineshape_lorentz));

}
