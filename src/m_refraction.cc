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
#include "refraction.h"



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
  Index   firstH2O = find_first_species_tg( gas_species, 
                                      species_index_from_species_name("H2O") );

  if( firstH2O < 0 )
    throw runtime_error( 
       "Water vapour is a requiered (must be a tag group in *gas_species*)." );

  refr_index_thayer_1974( refr_index, a_pressure, a_temperature, 
                                                        a_vmr_list[firstH2O] );
}



//! refr_indexUnit
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2003-01-18
*/
void refr_indexUnit(
	     Numeric&                    refr_index,
       const Numeric&                    a_pressure,
       const Numeric&                    a_temperature,
       const Vector&                     a_vmr_list,
       const ArrayOfArrayOfSpeciesTag&   gas_species )
{
  refr_index = 1;
}
