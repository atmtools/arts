/* Copyright (C) 2003-2008 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>

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

#include "absorption.h"
#include "arts.h"
#include "check_input.h"
#include "matpackI.h"
#include "messages.h"
#include "refraction.h"
#include "special_interp.h"
#include "abs_species_tags.h"



/*===========================================================================
  === WSMs
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexIR(// WS Output
                  Numeric&       refr_index,
                  // WS Input
                  const Numeric& a_pressure,
                  const Numeric& a_temperature,
                  const Vector&  a_vmr_list,
                  const Verbosity&)
{
  //To suppress warning about unused parameter
  a_vmr_list.nelem();

  refr_index_ir( refr_index, a_pressure, a_temperature );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexThayer(Numeric&        refr_index,
                      const Numeric&  a_pressure,
                      const Numeric&  a_temperature,
                      const Vector&   a_vmr_list,
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const Verbosity&)
{
  if( abs_species.nelem() != a_vmr_list.nelem() )
    throw runtime_error( "The number of tag groups differ between "
                                           "*a_vmr_list* and *abs_species*." );

  Index   firstH2O = find_first_species_tg( abs_species,
                                      species_index_from_species_name("H2O") );

  if( firstH2O < 0 )
    throw runtime_error(
       "Water vapour is a required (must be a tag group in *abs_species*)." );

  refr_index_thayer_1974( refr_index, a_pressure, a_temperature,
                                                        a_vmr_list[firstH2O] );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void refr_indexUnit(Numeric& refr_index,
                    const Verbosity&)
{
  refr_index = 1;
}

