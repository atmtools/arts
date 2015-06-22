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
void abs_linesReplaceWithLines(ArrayOfLineRecord& abs_lines, 
                               const ArrayOfLineRecord& replacement_lines, 
                               const Verbosity&)
{
  
  if(replacement_lines.nelem()==0)
    throw std::runtime_error("replacement_lines is empty.\n");
  
  ArrayOfIndex matches;
  ArrayOfQuantumMatchInfo match_info;

  for (Index ri = 0; ri < replacement_lines.nelem(); ri++)
  {
    QuantumIdentifier QI;
    QI.SetSpecies(replacement_lines[ri].Species());
    QI.SetIsotopologue(replacement_lines[ri].Isotopologue());
    QI.SetTransition(replacement_lines[ri].QuantumNumbers().Upper(),replacement_lines[ri].QuantumNumbers().Lower());
    
    
    // Run internal mathcing routine
    match_lines_by_quantum_identifier(matches, match_info, QI, abs_lines);
    
    // We demand that things are formatted the right way and that there are not multiple matches.
    if( matches.nelem()>1 )
    { 
      ostringstream os;
      os << "Multiple matches in comparison.  Something is wrong!\n"
         << "Line is:\n" << replacement_lines[ri]<<std::endl;
      throw std::runtime_error(os.str());
    }
    else if( matches.nelem()==0 )
      { 
        ostringstream os;
        os << "No match found!  Make sure your replacement lines and abs_lines have the same quantum numbers definition.\n"
         << "Line is:\n" << replacement_lines[ri]<<std::endl;
        throw std::runtime_error(os.str());
    }
    
    LineRecord& lr_old = abs_lines[matches[0]];
    
    // If any of the levels match partially or fully set the right quantum number
    switch (match_info[0].Upper())
    {
      case QMI_NONE:
        {
        ostringstream os;
        os << "There are no quantum numbers in your replacement ines so they match to abs_lines.\n"
           << "replacement_line:\n"<<replacement_lines[ri]<<"\nabs_line:\n"<<lr_old<<std::endl;
        throw std::runtime_error(os.str());
        break;
        }
      case QMI_PARTIAL:
        {
        ostringstream os;
        os << "Your replacement lines are only partially defined so they match to the abs_lines.\n"
           << "replacement_line:\n"<<replacement_lines[ri]<<"\nabs_line:\n"<<lr_old<<std::endl;
        throw std::runtime_error(os.str());
        break;
        }
      case QMI_FULL:
        {
        lr_old = replacement_lines[ri];
        break;
        }
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReplaceParameterWithLinesParameter(ArrayOfLineRecord& abs_lines, 
                                                 const ArrayOfLineRecord& replacement_lines, 
                                                 const String& parameter_name,
                                                 const Verbosity&)
{
  
  Index parameter_switch = -1;
  
  if(parameter_name.nelem()==0)
    throw std::runtime_error("parameter_name is empty.\n");
  if(replacement_lines.nelem()==0)
    throw std::runtime_error("replacement_lines is empty.\n");
  else if(parameter_name == "Central Frequency")
    parameter_switch = 0;
  else if(parameter_name == "Line Strength")
    parameter_switch = 1;
  else if(parameter_name == "Pressure Broadening")
    parameter_switch = 2;
  else if(parameter_name == "Line Mixing")
    parameter_switch = 3;
  else if(parameter_name == "Lower State Energy")
    parameter_switch = 4;
  
  
  ArrayOfIndex matches;
  ArrayOfQuantumMatchInfo match_info;

  for (Index ri = 0; ri < replacement_lines.nelem(); ri++)
  {
    const LineRecord& lr = replacement_lines[ri];
    QuantumIdentifier QI;
    QI.SetSpecies(lr.Species());
    QI.SetIsotopologue(lr.Isotopologue());
    QI.SetTransition(lr.QuantumNumbers().Upper(),lr.QuantumNumbers().Lower());
    
    

    // Run internal mathcing routine
    match_lines_by_quantum_identifier(matches, match_info, QI, abs_lines);
    
    // We demand that things are formatted the right way and that there are not multiple matches.
    if( matches.nelem()>1 )
    { 
      ostringstream os;
      os << "Multiple matches in comparison.  Something is wrong!\n"
         << "Line is:\n" << lr<<std::endl;
      throw std::runtime_error(os.str());
    }
    else if( matches.nelem()==0 )
      { 
        ostringstream os;
        os << "No match found!  Make sure your replacement lines and abs_lines have the same quantum numbers definition.\n"
         << "Line is:\n" << lr<<std::endl;
        throw std::runtime_error(os.str());
    }
    
    LineRecord& lr_old = abs_lines[matches[0]];
    
    // If any of the levels match partially or fully set the right quantum number
    switch (match_info[0].Upper())
    {
      case QMI_NONE:
        {
        ostringstream os;
        os << "There are no quantum numbers in your replacement ines so they match to abs_lines.\n"
           << "replacement_line:\n"<<lr<<"\nabs_line:\n"<<lr_old<<std::endl;
        throw std::runtime_error(os.str());
        break;
        }
      case QMI_PARTIAL:
        {
        ostringstream os;
        os << "Your replacement lines are only partially defined so they match to the abs_lines.\n"
           << "replacement_line:\n"<<lr<<"\nabs_line:\n"<<lr_old<<std::endl;
        throw std::runtime_error(os.str());
        break;
        }
      case QMI_FULL:
        switch (parameter_switch)
        {
          case 0: //"Central Frequency":
            lr_old.setF(lr.F());
            break;
          case 1: //"Line Strength":
            lr_old.setI0(lr.I0());
            break;
          case 2: //"Pressure Broadening":
            lr_old.SetPressureBroadeningData(lr.PressureBroadening());
            break;
          case 3: //"Line Mixing":
            lr_old.SetLineMixingData(lr.LineMixing());
            break;
          case 4: //"Lower State Energy":
            lr_old.SetElow(lr.Elow());
            break;
          default:
          {
            ostringstream os;
            os << "Usupported paramter_name\n" << parameter_name
               << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
            break;
          }
            
        }
        break;
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void nlteSetByQuantumIdentifiers(Index& nlte_do,
                                 ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                 const ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
                                 const ArrayOfArrayOfSpeciesTag& abs_species,
                                 const Tensor4&                   t_nlte_field,
                                 const Vector&                    p_grid,
                                 const Vector&                    lat_grid,
                                 const Vector&                    lon_grid,
                                 const Index&                     atmosphere_dim,
                                 const Vector& vibrational_energies,
                                 const Verbosity&)
{
    nlte_do = 1;

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
                        if(lr.EvuppIndex()==-1)
                        {
                            lr.SetEvuppIndex(qi);
                            if(do_ev)
                                lr.SetEvupp(vibrational_energies[qi]);
                        }
                        else
                        {
                            ostringstream os;
                            os << "The linerecord:\n"<<lr<<"\nhad the energy state level of "<<
                            "this quantum identifier:\n"<<nlte_quantum_identifiers[qi]<<
                            "\nset twice by the input quantum identifiers.  All levels must "<<
                            "point at a unique state. " << qi;
                            throw std::runtime_error(os.str());
                        }
                        break;
                    case QMI_PARTIAL:
                        if(lr.EvuppIndex()==-1)
                        {
                            lr.SetEvuppIndex(qi);
                            if(do_ev)
                                lr.SetEvupp(vibrational_energies[qi]);
                        }
                        else
                        {
                            ostringstream os;
                            os << "The linerecord:\n"<<lr<<"\nhad the energy state level of "<<
                            "this quantum identifier:\n"<<nlte_quantum_identifiers[qi]<<
                            "\nset twice by the input quantum identifiers.  All levels must "<<
                            "point at a unique state. " << qi;
                            throw std::runtime_error(os.str());
                        }
                        break;
                }
                switch (match_info[i].Lower())
                {
                    case QMI_NONE:    break;
                    case QMI_FULL:
                        if(lr.EvlowIndex()==-1)
                        {
                            lr.SetEvlowIndex(qi);
                            if(do_ev)
                                lr.SetEvlow(vibrational_energies[qi]);
                        }
                        else
                        {
                            ostringstream os;
                            os << "The linerecord:\n"<<lr<<"\nhad the energy state level of "<<
                            "this quantum identifier:\n"<<nlte_quantum_identifiers[qi]<<
                            "\nset twice by the input quantum identifiers.  All levels must "<<
                            "point at a unique state. " << qi;
                            throw std::runtime_error(os.str());
                        }
                        break;
                    case QMI_PARTIAL:
                        if(lr.EvlowIndex()==-1)
                        {
                            lr.SetEvlowIndex(qi);
                            if(do_ev)
                                lr.SetEvlow(vibrational_energies[qi]);
                        }
                        else
                        {
                            ostringstream os;
                            os << "The linerecord:\n"<<lr<<"\nhad the energy state level of "<<
                            "this quantum identifier:\n"<<nlte_quantum_identifiers[qi]<<
                            "\nset twice by the input quantum identifiers.  All levels must "<<
                            "point at a unique state. " << qi;
                            throw std::runtime_error(os.str());
                        }
                        break;
                }
            }
        }
    }
    
    chk_nlte(t_nlte_field, nlte_quantum_identifiers, abs_lines_per_species,
             p_grid, lat_grid, lon_grid, atmosphere_dim);
}
