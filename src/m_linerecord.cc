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
        throw std::runtime_error("*abs_lines* is empty.  Is shifting frequency really intended?");
    
    // Shift all line center frequencies
    for(Index jj=0;jj<abs_lines.nelem();jj++)
    {
        abs_lines[jj].setF(abs_lines[jj].F()+freqeuncy_shift);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesRelativeLineStrengthShift(ArrayOfLineRecord& abs_lines, const Numeric& relative_line_strength_shift, const Verbosity&)
{
    
    // Catch use case that is not a use case
    if(abs_lines.nelem()==0)
        throw std::runtime_error("*abs_lines* is empty.  Is shifting line strength really intended?");
    
    // Rescale all line strengths
    for(Index jj=0;jj<abs_lines.nelem();jj++)
    {
        abs_lines[jj].setI0(abs_lines[jj].I0()*(1.0+relative_line_strength_shift));
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
    {
        // Get a reference to the current list of lines to save typing:
        ArrayOfLineRecord& ll = abs_lines_per_species[ii];
        for(Index jj=0;jj<ll.nelem();jj++)
        {
            ll[jj].setI0(ll[jj].I0()*(1.0+relative_line_strength_shift));
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesMatchNLTEQuantumIdentifiers(ArrayOfArrayOfLineRecord& abs_lines_per_species, 
                                                      const ArrayOfQuantumIdentifier& nlte_quantum_identifiers, 
                                                      const ArrayOfArrayOfSpeciesTag& abs_species,
                                                      const Vector& vibrational_energies,
                                                      const Verbosity&)
{
  
  const bool do_ev = vibrational_energies.nelem();
  
  if(do_ev)
  {
    if( vibrational_energies.nelem()!=nlte_quantum_identifiers.nelem() )
    {
      ostringstream os;
      os << "Your vibrational energy levels vector is not the same size as\n"
         << "your *nlte_quantum_identifiers* array.  These must be the same\n"
         << "size and the content should match.\n";
      throw std::runtime_error(os.str());
    }
  }
  
    ArrayOfIndex matches;
    ArrayOfQuantumMatchInfo match_info;

    for (Index qi = 0; qi < nlte_quantum_identifiers.nelem(); qi++)
    {
        for (Index s = 0; s < abs_lines_per_species.nelem(); s++)
        {
          
            // Skip this species if qi is not part of the species represented by this abs_lines
            if(abs_species[s][0].Species() != nlte_quantum_identifiers[qi].Species())
              continue;
            
            ArrayOfLineRecord& species_lines = abs_lines_per_species[s];
            
            // Run internal mathcing routine
            match_lines_by_quantum_identifier(matches, match_info, nlte_quantum_identifiers[qi], species_lines);

            // Use info about mathced lines to tag the relevant parameter
            for (Index i = 0; i < matches.nelem(); i++)
            {
                // For each line record
                LineRecord& lr = species_lines[matches[i]];
                
                // If any of the levels match partially or fully set the right quantum number
                switch (match_info[i].Upper())
                {
                    case QMI_NONE:    break;
                    case QMI_FULL:    
                      lr.SetEvuppIndex(qi); 
                      if(do_ev) 
                        lr.SetEvupp(vibrational_energies[qi]);
                      break;
                    case QMI_PARTIAL: 
                      lr.SetEvuppIndex(qi);
                      if(do_ev) 
                        lr.SetEvupp(vibrational_energies[qi]);
                      break;
                }
                switch (match_info[i].Lower())
                {
                    case QMI_NONE:    break;
                    case QMI_FULL:    
                      lr.SetEvlowIndex(qi);
                      if(do_ev) 
                        lr.SetEvlow(vibrational_energies[qi]);
                      break;;
                    case QMI_PARTIAL: 
                      lr.SetEvlowIndex(qi);
                      if(do_ev) 
                        lr.SetEvlow(vibrational_energies[qi]);
                      break;;
                }
            }
        }
    }
}