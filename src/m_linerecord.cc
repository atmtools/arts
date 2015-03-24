/* Copyright (C) 2015
 Richard Larsson <ric.larsson@gmail.com>
 
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

//======================================================================
//         Functions for altering the line catalog
//======================================================================

#include "auto_md.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesShiftFrequency(ArrayOfLineRecord& abs_lines, const Numeric& freqeuncy_shift, const Verbosity&)
{
    // Catch use case that is not a use case
    if(abs_lines.nelem()==0)
        throw std::runtime_error("*abs_lines_per_species* is empty.  Is shifting frequency really intended?");
    
    for(Index jj=0;jj<abs_lines.nelem();jj++)
    {
        abs_lines[jj].setF(abs_lines[jj].F()+freqeuncy_shift);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesShiftFrequency(ArrayOfArrayOfLineRecord& abs_lines_per_species, const Numeric& freqeuncy_shift, const Verbosity&)
{
    
    // Catch use case that is not a use case
    if(abs_lines_per_species.nelem()==0)
        throw std::runtime_error("*abs_lines_per_species* is empty.  Is shifting frequency really intended?");
    
    // Simply shift all lines from their original frequency by input *freqeuncy_shift*.
    for(Index ii=0;ii<abs_lines_per_species.nelem();ii++)
    {
        // Get a reference to the current list of lines to save typing:
        ArrayOfLineRecord& ll = abs_lines_per_species[ii];
        for(Index jj=0;jj<ll.nelem();jj++)
        {
            ll[jj].setF(ll[jj].F()+freqeuncy_shift);
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesRelativeLineStrengthShift(ArrayOfArrayOfLineRecord& abs_lines_per_species, const Numeric& relative_line_strength_shift, const Verbosity&)
{
    
    // Catch use case that is not a use case
    if(abs_lines_per_species.nelem()==0)
        throw std::runtime_error("*abs_lines_per_species* is empty.  Is shifting line strength really intended?");
    
    // Simply rescale all lines from their original line strength by input *relative_line_strength_shift*.
    for(Index ii=0;ii<abs_lines_per_species.nelem();ii++)
        for(Index jj=0;jj<abs_lines_per_species[ii].nelem();jj++)
            abs_lines_per_species[ii][jj].setI0(abs_lines_per_species[ii][jj].I0()*(1.0+relative_line_strength_shift));
}