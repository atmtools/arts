/* Copyright (C) 2003 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                            
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
  ===  File description 
  ===========================================================================*/

/*!
  \file   m_refraction.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2003-01-09

  \brief  Workspace methods releated to refraction.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <iostream>
#include "absorption.h"
#include "arts.h"
#include "auto_md.h"
#include "matpackI.h"



/*===========================================================================
  === WSMs
  ===========================================================================*/

//! refr_indexThayer
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2003-01-09
*/
void refr_indexThayer(
	     Numeric&                    refr_index,
       const Numeric&                    a_pressure,
       const Numeric&                    a_temperature,
       const Vector&                     a_vmr_list,
       const ArrayOfArrayOfSpeciesTag&   gas_species )
{
  refr_index = 1.0 + 77.6e-8 * a_pressure / a_temperature;

  //  for ( Index i=0; i<n; i++ )
  //{
  // e = p_abs[i] * h2o_abs[i];
  //  p = p_abs[i] - e;
  //
  //  refr_index[i] = 1.0 + 77.6e-8 * p / t_abs[i] + 
  //                        72e-8 * e / t_abs[i] +
  //                        3.754e-3 * e / (t_abs[i]*t_abs[i]);
  //}
}
