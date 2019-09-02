/* Copyright (C) 2005-2012 Cory Davis <cory@met.ed.ac.uk>
                            
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
  === File description 
  ===========================================================================*/

/*!
  \file   mc_interp.h
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2005-02-28 

  \brief  Interpolation classes and functions created for use within Monte 
  Carlo scattering simulations 

*/
/*===========================================================================
  === External declarations
  ===========================================================================*/

#ifndef mc_interp_h
#define mc_interp_h
#include "array.h"
#include "arts.h"
#include "interpolation.h"
#include "matpackI.h"
#include "optproperties.h"
#include "ppath.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

//! A 2D sequential linear interpolation (SLI) lookup table
/*! This class holds the gridded for 2D SLI as well as the
interpolate member function for retrieving interpolated values. 
 */
class SLIData2 {
 public:
  //grid of x1 values where y is known
  Vector x1a;
  //A vector of x2 values for every x1a
  ArrayOfVector x2a;
  //y values for every x1a, x2a
  ArrayOfVector ya;
  //performs SLI.
  Numeric interpolate(Numeric x1, Numeric x2) const;
  //checks that it is not empty
  //void check() const;
};

ostream& operator<<(ostream& os, const SLIData2& sli);

void interp(MatrixView tia,
            ConstVectorView itw,
            const ArrayOfMatrix& a,
            const GridPos& tc);

void interp(VectorView tia,
            ConstVectorView itw,
            const ArrayOfVector& a,
            const GridPos& tc);

void interp_scat_angle_temperature(  //Output:
    VectorView pha_mat_int,
    Numeric& theta_rad,
    //Input:
    const SingleScatteringData& scat_data_single,
    const Numeric& za_sca,
    const Numeric& aa_sca,
    const Numeric& za_inc,
    const Numeric& aa_inc,
    const Numeric& rtp_temperature);

#endif  // mc_interp_h
