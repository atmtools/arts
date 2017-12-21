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
#include "absorption.h"
#include "sorting.h"

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
  
  //if(replacement_lines.nelem()==0)
  //  throw std::runtime_error("replacement_lines is empty.\n");
  
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
void abs_linesChangeParameterForMatchingLines(ArrayOfLineRecord& abs_lines, 
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const QuantumIdentifier& QI, 
                                              const String& parameter_name,
                                              const Numeric& change,
                                              const Index& relative,
                                              const Index& loose_matching,
                                              const Verbosity&)
{
    
    Index parameter_switch = -1;
    
    if(parameter_name.nelem()==0)
        throw std::runtime_error("parameter_name is empty.\n");
    else if(parameter_name == "Central Frequency")
        parameter_switch = 0;
    else if(parameter_name == "Line Strength")
        parameter_switch = 1;
    else if(parameter_name == "Pressure Broadening Self")
        parameter_switch = 2;
    else if(parameter_name == "Pressure Broadening Foreign")
        parameter_switch = 3;
    else if(parameter_name == "Lower State Energy")
        parameter_switch = 4;
    else if(parameter_name == "Pressure Broadening Self Exponent")
        parameter_switch = 5;
    else if(parameter_name == "Pressure Broadening Foreign Exponent")
        parameter_switch = 6;
    
    ArrayOfIndex matches;
    ArrayOfQuantumMatchInfo match_info;
    // Run internal mathcing routine
    match_lines_by_quantum_identifier(matches, match_info, QI, abs_lines);
    // We demand that things are formatted the right way and that there are not multiple matches.
    if( matches.nelem()>1 && loose_matching == 0 )
        throw std::runtime_error("Multiple matches in comparison.  You set loose_matching to not allow this!\n");
    else if( matches.nelem()==0 )
        throw std::runtime_error("No match found!  Make sure your QuantumIdentifier is in abs_lines before using this function.\n(For instance, try to make sure quantum numbers are defined the same way.)\n");
    
    // Broadening species
    const Index this_species = find_first_species_tg( abs_species, QI.Species() );
    ArrayOfIndex broad_spec_locations;
    find_broad_spec_locations(broad_spec_locations, abs_species, this_species );
    // Water index
    const Index h2o_index = find_first_species_tg( abs_species, species_index_from_species_name("H2O") );
    
    bool any=false;
    PressureBroadeningData pb;
    
    for(Index mii =0; mii<matches.nelem(); mii++)
    {
        
        // Skip if there are none in any of the levels
        if((match_info[mii].Upper()==QMI_NONE||match_info[mii].Lower()==QMI_NONE)&&loose_matching==0)
            continue;
        
        // Skip if there are any partials unless we accept loose matching
        if((match_info[mii].Upper()==QMI_PARTIAL||match_info[mii].Lower()==QMI_PARTIAL)&&loose_matching==0)
            continue;
        
        if(!any)
            any = true;
        
        LineRecord& lr = abs_lines[matches[mii]];
        
        switch (parameter_switch)
        {
            case 0: //"Central Frequency":
                if(relative==0)
                    lr.setF(lr.F()+change);
                else 
                    lr.setF(lr.F()*(1.0e0+change));
                break;
            case 1: //"Line Strength":
                if(relative==0)
                    lr.setI0(lr.I0()+change);
                else 
                    lr.setI0(lr.I0()*(1.0e0+change));
                break;
            case 2: //"Pressure Broadening Self":
                pb = lr.PressureBroadening();
                if(relative==0)
                    pb.ChangeSelf(change,this_species,h2o_index,broad_spec_locations);
                else
                    pb.ChangeSelfRelative(change,this_species,h2o_index,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
                break;
            case 3: //"Pressure Broadening Foreign":
                pb = lr.PressureBroadening();
                if(relative==0)
                    pb.ChangeForeign(change,broad_spec_locations);
                else
                    pb.ChangeForeignRelative(change,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
                break;
            case 4: //"Lower State Energy":
                if(relative==0)
                    lr.SetElow(lr.Elow()+change);
                else 
                    lr.SetElow(lr.Elow()*(1.0e0+change));
                break;
            case 5: //"Pressure Broadening Self":
                pb = lr.PressureBroadening();
                if(relative==0)
                    pb.ChangeSelfExponent(change,this_species,h2o_index,broad_spec_locations);
                else
                    pb.ChangeSelfExponentRelative(change,this_species,h2o_index,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
                break;
            case 6: //"Pressure Broadening Foreign":
                pb = lr.PressureBroadening();
                if(relative==0)
                    pb.ChangeForeignExponent(change,broad_spec_locations);
                else
                    pb.ChangeForeignExponentRelative(change,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
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
    }
    
    if(!any)
        throw std::runtime_error("You have no matches.  This is not accepted as a valid use case.  (Is your matching information correct?)\n");
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetParameterForMatchingLines(ArrayOfLineRecord& abs_lines, 
                                           const ArrayOfArrayOfSpeciesTag& abs_species,
                                           const QuantumIdentifier& QI, 
                                           const String& parameter_name,
                                           const Numeric& new_value,
                                           const Index& loose_matching,
                                           const Verbosity&)
{
    
    Index parameter_switch = -1;
    
    if(parameter_name.nelem()==0)
        throw std::runtime_error("parameter_name is empty.\n");
    else if(parameter_name == "Central Frequency")
        parameter_switch = 0;
    else if(parameter_name == "Line Strength")
        parameter_switch = 1;
    else if(parameter_name == "Pressure Broadening Self")
        parameter_switch = 2;
    else if(parameter_name == "Pressure Broadening Foreign")
        parameter_switch = 3;
    else if(parameter_name == "Lower State Energy")
        parameter_switch = 4;
    else if(parameter_name == "Pressure Broadening Self Exponent")
        parameter_switch = 5;
    else if(parameter_name == "Pressure Broadening Foreign Exponent")
        parameter_switch = 6;
    
    ArrayOfIndex matches;
    ArrayOfQuantumMatchInfo match_info;
    // Run internal mathcing routine
    match_lines_by_quantum_identifier(matches, match_info, QI, abs_lines);
    // We demand that things are formatted the right way and that there are not multiple matches.
    if( matches.nelem()>1 && loose_matching == 0 )
        throw std::runtime_error("Multiple matches in comparison.  You set loose_matching to not allow this!\n");
    else if( matches.nelem()==0 )
        throw std::runtime_error("No match found!  Make sure your QuantumIdentifier is in abs_lines before using this function.\n(For instance, try to make sure quantum numbers are defined the same way.)\n");
    
    // Broadening species
    const Index this_species = find_first_species_tg( abs_species, QI.Species() );
    ArrayOfIndex broad_spec_locations;
    find_broad_spec_locations(broad_spec_locations, abs_species, this_species );
    // Water index
    const Index h2o_index = find_first_species_tg( abs_species, species_index_from_species_name("H2O") );
    
    bool any=false;
    PressureBroadeningData pb;
    
    for(Index mii =0; mii<matches.nelem(); mii++)
    {
        // Skip if there are none in any of the levels
        if((match_info[mii].Upper()==QMI_NONE||match_info[mii].Lower()==QMI_NONE)&&loose_matching==0)
            continue;
        
        // Skip if there are any partials unless we accept loose matching
        if((match_info[mii].Upper()==QMI_PARTIAL||match_info[mii].Lower()==QMI_PARTIAL)&&loose_matching==0)
            continue;
        
        if(!any)
            any = true;
        
        LineRecord& lr = abs_lines[matches[mii]];
        
        switch (parameter_switch)
        {
            case 0: //"Central Frequency":
                lr.setF(new_value);
                break;
            case 1: //"Line Strength":
                lr.setI0(new_value);
                break;
            case 2: //"Pressure Broadening Self":
                pb = lr.PressureBroadening();
                pb.SetSelf(new_value,this_species,h2o_index,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
                break;
            case 3: //"Pressure Broadening Foreign":
                pb = lr.PressureBroadening();
                pb.SetForeign(new_value,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
                break;
            case 4: //"Lower State Energy":
                lr.SetElow(new_value);
                break;
            case 5: //"Pressure Broadening Self Exponent":
                pb = lr.PressureBroadening();
                pb.SetSelfExponent(new_value,this_species,h2o_index,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
                break;
            case 6: //"Pressure Broadening Foreign Exponent":
                pb = lr.PressureBroadening();
                pb.SetForeignExponent(new_value,broad_spec_locations);
                lr.SetPressureBroadeningData(pb);
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
    }
    
    if(!any)
        throw std::runtime_error("You have no matches.  This is not accepted as a valid use case.  (Is your matching information correct?)\n");
}


/* Workspace method: Doxygen documentation will be auto-generated */
void nlteSetByQuantumIdentifiers(Index&                           nlte_do,
                                 ArrayOfArrayOfLineRecord&        abs_lines_per_species,
                                 const ArrayOfQuantumIdentifier&  nlte_quantum_identifiers,
                                 const ArrayOfArrayOfSpeciesTag&  abs_species,
                                 const Tensor4&                   t_nlte_field,
                                 const Vector&                    p_grid,
                                 const Vector&                    lat_grid,
                                 const Vector&                    lon_grid,
                                 const Index&                     atmosphere_dim,
                                 const Vector&                    vibrational_energies,
                                 const Index&                     population_type,
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
    
    if(not (population_type < (Index) LinePopulationType::End) and not (population_type > 0))
    {
      throw std::runtime_error("Cannot understand given population type");
    }
    
    // All energies must be positive
    for(Index ii=0; ii< vibrational_energies.nelem();ii++)
      if(vibrational_energies[ii]<0)
        {
          ostringstream os;
          os << "Some of your vibrational energy levels are negative.  They should be positive.\n"
             << "Your vibrational levels are:\n" <<vibrational_energies;
          throw std::runtime_error(os.str());
        }
        
    
    #pragma omp parallel for        \
    if (!arts_omp_in_parallel()) 
    for (Index qi = 0; qi < nlte_quantum_identifiers.nelem(); qi++)
    {
        ArrayOfIndex matches;
        ArrayOfQuantumMatchInfo match_info;
        
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
                        if(lr.NLTEUpperIndex()==-1)
                        {
                            lr.SetNLTEUpperIndex(qi);
                            if(do_ev)
                              lr.SetEvupp(vibrational_energies[qi]);
                            lr.SetLinePopulationTypeFromIndex(population_type);
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
                        if(lr.NLTEUpperIndex()==-1)
                        {
                            lr.SetNLTEUpperIndex(qi);
                            if(do_ev)
                              lr.SetEvupp(vibrational_energies[qi]);
                            lr.SetLinePopulationTypeFromIndex(population_type);
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
                        if(lr.NLTELowerIndex()==-1)
                        {
                            lr.SetNLTELowerIndex(qi);
                            if(do_ev)
                              lr.SetEvlow(vibrational_energies[qi]);
                            lr.SetLinePopulationTypeFromIndex(population_type);
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
                        if(lr.NLTELowerIndex()==-1)
                        {
                            lr.SetNLTELowerIndex(qi);
                            if(do_ev)
                              lr.SetEvlow(vibrational_energies[qi]);
                            lr.SetLinePopulationTypeFromIndex(population_type);
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


/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromabs_linesSet(Vector& f_grid,
                            const ArrayOfLineRecord& abs_lines,
                            const Numeric&           half_width,
                            const Index&             nr_f_per_line,
                            const Index&             line_nr,
                            const Verbosity&         verbosity)
{
    CREATE_OUT2;
    
    const Index lines = abs_lines.nelem();
    
    if(line_nr<0)
        f_grid.resize(lines*nr_f_per_line);
    else
        f_grid.resize(nr_f_per_line);
    
    out2 << "  Creating f_grid vector of length "<<f_grid.nelem();
    
    if(lines == 0)
        throw std::runtime_error("You need at least one line to run this code.\n");
    if(line_nr>=lines)
        throw std::runtime_error("You specified a line number that is outside the range of abs_lines.\n");
    if(nr_f_per_line<1)
        throw std::runtime_error("You need more than 0 frequencies per line to execute this function.\n");
    
    // Helper variable to ensure that there are no overlaps
    Numeric f_max = 0.0;
    
    if(line_nr>=0) // there is a line, then set frequency for a single line
    {
        if((abs_lines[line_nr].F()-half_width)<f_max)
            throw std::runtime_error("Frequencies below 0 Hz are not supported by this function.\n");
        VectorNLinSpace(f_grid,nr_f_per_line,abs_lines[line_nr].F()-half_width,abs_lines[line_nr].F()+half_width,verbosity);
    }
    else // if there are many lines, then set frequency from the many lines
    {
        Vector tmp;
        for(Index ii=0;ii<lines;ii++)
        {
            if((abs_lines[ii].F()-half_width)<f_max)
                throw std::runtime_error("Frequency overlaps are not supported by this function.\n");
            else
                f_max = abs_lines[ii].F()+half_width;
            VectorNLinSpace(tmp,nr_f_per_line,abs_lines[ii].F()-half_width,abs_lines[ii].F()+half_width,verbosity);
            f_grid[Range(ii*nr_f_per_line,nr_f_per_line)]=tmp;
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromabs_lines_per_speciesSetFromSpeciesTag(Vector& f_grid,
                                                      const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                                      const ArrayOfArrayOfSpeciesTag& abs_species,
                                                      const Numeric&                  half_width,
                                                      const Index&                    nr_f_per_line,
                                                      const String&                   species_tag,
                                                      const Verbosity&                verbosity)
{
    CREATE_OUT2;
    
    if(species_tag == "")
        throw std::runtime_error("You need at least one tag in this code.\n");
    
    ArrayOfSpeciesTag st;
    array_species_tag_from_string(st,species_tag);
    
    for(Index ii=0; ii<abs_species.nelem(); ii++)
    {
        bool test = false;
        if(abs_species[ii].nelem()==st.nelem())
        {
            for(Index jj=0;  jj<st.nelem(); jj++)
            {
                if(st[jj]==abs_species[ii][jj])
                    test=true;
                else
                {
                    test=false;
                    break;
                }
            }
            if(test)
            {
                f_gridFromabs_linesSet(f_grid,abs_lines_per_species[ii],half_width,nr_f_per_line,-1,verbosity);
                return;
            }
        }
    }
    
    throw std::runtime_error("No frequency set for the given species_tag.\n");
}
