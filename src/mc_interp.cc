/* Copyright (C) 2005 Cory Davis <cory@met.ed.ac.uk>
                            
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
  \file   mc_interp.cc
  \author Cory Davis <cory@met.ed.ac.uk>
  \date   2005-02-28 

  \brief  Interpolation classes and functions created for use within Monte 
  Carlo scattering simulations 

*/
/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "mc_interp.h"
#include "interpolation.h"

//! Perform sequential interpolation
/*!
  \param x1 desired x1 value
  \param x2 desired x2 value

  \return interpolated y value at x1,x2

  \author Cory Davis <cdavis@staffmail.ed.ac.uk>
  \date 2005-02-28
*/

Numeric SLIData2::interpolate(Numeric x1, Numeric x2)
{
  GridPos gp1,gpl,gpr;
  Vector itw1(2),itwl(2),itwr(2);
  Numeric yl,yr;
  //Get interpolation weights for x1 interpolation
  gridpos(gp1,this->x1a,x1);
  interpweights(itw1,gp1);
  //Get y values on bounding x1 grid points for desired x2
  gridpos(gpl,this->x2a[gp1.idx],x2);
  interpweights(itwl,gpl);
  gridpos(gpr,this->x2a[gp1.idx+1],x2);
  interpweights(itwr,gpr);
  yl=interp(itwl,this->ya[gp1.idx],gpl);
  yr=interp(itwr,this->ya[gp1.idx+1],gpr);
  //interpolate these two y values useing x1 interpolation weights
  return itw1[0]*yl+itw1[1]*yr;
}

ostream& operator<< (ostream &os, const SLIData2& /* sli */)
{
  os << "SLIData2    : Output operator not implemented";
  return os;
}
