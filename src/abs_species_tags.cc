/* Copyright (C) 2002-2012 Stefan Buehler  <sbuehler@ltu.se>

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
  \file   abs_species_tags.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue May 31 17:18:22 2005
  
  \brief  Stuff related to absorption species tags.
  
  This file contains functions related to SpeciesTags. It is better to
  separate this from the other absorption functions, since this part
  was actually improved in arts-1-1 and should be kept, whereas most
  other absorption stuff is back-ported from arts-1-0.
*/

#include "arts.h"
#include <cfloat>
#include <map>
#include "absorption.h"
#include "auto_md.h"
#include "abs_species_tags.h"
#include "global_data.h"

//! Constructor from a tag definition String. 
/*! 
  For examples see documentation of member function Name(). 

  \param def String containing tag definition.
  
  \exception runtime_error The given String could not be mapped to
  a sensible tag description.
*/
SpeciesTag::SpeciesTag(String def) 
{
  // Save input string for error messages:
  String def_original = def;
    
  // Species lookup data:
  using global_data::species_data;
  // Name of species and isotopologue (aux variables):
  String name, isoname;
  // Aux index:
  Index n;

  // Set default values for isotopologue
  misotopologue = -1;

  // Set frequency limits to default values (no limits):
  mlf = -1;
  muf = -1;

  // Set CIA species to -1 by default
  mcia_second = -1;

  // Set type to normal LBL species by default
  mtype = TYPE_PLAIN;
  
  // Set line mixing to off by default
  mline_mixing = LINE_MIXING_OFF;

  // We cannot set a default value for the isotopologue, because the
  // default should be `ALL' and the value for `ALL' depends on the
  // species. 
    

  // Extract the species name:
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      name = def.substr(0,n);      // Extract before '-'
      def.erase(0,n+1);             // Remove from def
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In that
      // case we assume that it contains just a species name and
      // nothing else
      name = def;
      def  = "";
    }

  // Obtain species index from species name.
  // (This will also remove possible whitespace.)
  mspecies = species_index_from_species_name( name );

  // Remove whitespace
  name.trim();

  // Check if species name contains the special tag for
  // Faraday Rotation
  if (name == "free_electrons")
    {
      mtype = TYPE_FREE_ELECTRONS;
      return;
    }

  // Check if species name contains the special tag for
  // Particles
  if (name == "particles")
    {
      mtype = TYPE_PARTICLES;
      return;
    }

  if ( 0 > mspecies )
    {
      ostringstream os;
      os << "Species \"" << name << "\" is not a valid species.";
      throw runtime_error(os.str());
    }

  // Get a reference to the relevant Species Record:
  const SpeciesRecord& spr = species_data[mspecies];

  if ( 0 == def.nelem() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopologues and no frequency limits.
      // Frequency defaults are already set. Set isotopologue defaults:
      misotopologue = spr.Isotopologue().nelem();
      // This means all isotopologues.

      return;
    }
    
  // Extract the isotopologue name/Zeeman flag:
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      isoname = def.substr(0,n);    // Extract before '-'
      def.erase(0,n+1);             // Remove from def

      if ("Z" == isoname)
        {
          mtype = TYPE_ZEEMAN;
          // Zeeman flag was present, now extract the isotopologue name:
          n    = def.find('-');    // find the '-'
          if (n != def.npos )
            {
              isoname = def.substr(0,n);    // Extract before '-'
              def.erase(0,n+1);             // Remove from def
            }
          else
            {
              // n==def.npos means that def does not contain a '-'. In that
              // case we assume that it contains just the isotopologue name and
              // nothing else.
              isoname = def;
              def  = "";
            }
        }

      if ("HXSEC" == isoname)
        {
          mtype = TYPE_HITRAN_XSEC;
          // Hitran Xsec flag was present, now extract the isotopologue name:
          n    = def.find('-');    // find the '-'
          if (n != def.npos )
            {
              isoname = def.substr(0,n);    // Extract before '-'
              def.erase(0,n+1);             // Remove from def
            }
          else
            {
              // n==def.npos means that def does not contain a '-'. In that
              // case we assume that it contains just the isotopologue name and
              // nothing else.
              isoname = def;
              def  = "";
            }
        }

        if ("LM" == isoname)
        {
            mline_mixing = LINE_MIXING_ON;

            // Line mixing flag was present, now extract the isotopologue name:
            n    = def.find('-');    // find the '-'
            if (n != def.npos )
            {
                isoname = def.substr(0,n);    // Extract before '-'
                def.erase(0,n+1);             // Remove from def
            }
            else
            {
                // n==def.npos means that def does not contain a '-'. In that
                // case we assume that it contains just the isotopologue name and
                // nothing else.
                isoname = def;
                def  = "";
            }
        }
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In that
      // case we assume that it contains just the isotopologue name or
      // Zeeman flag and nothing else.
      isoname = def;
      def  = "";
      if ("Z" == isoname)
        {
          mtype = TYPE_ZEEMAN;
          // This means that there is nothing else to parse. Apparently
          // the user wants all isotopologues and no frequency limits.
          misotopologue = spr.Isotopologue().nelem();
          return;
        }
      if ("LM" == isoname)
        {
          mline_mixing = LINE_MIXING_ON;
          // This means that there is nothing else to parse. Apparently
          // the user wants all isotopologues and no frequency limits.
          misotopologue = spr.Isotopologue().nelem();
          return;
        }
      if ("HXSEC" == isoname)
        {
          mtype = TYPE_HITRAN_XSEC;
          // This means that there is nothing else to parse. Apparently
          // the user wants all isotopologues and no frequency limits.
          misotopologue = spr.Isotopologue().nelem();
          return;
        }
    }

  // Check for joker:
  if ( "*" == isoname )
    {
      // The user wants all isotopologues. Set this accordingly:
      misotopologue = spr.Isotopologue().nelem();
    }
//  else if ( "nl" == isoname )     // Check for "nl":
//    {
//      // The user wants no lines at all. Set this accordingly:
//      misotopologue = -1;
//    }
  else if ( "CIA" == isoname )     // Check for "cia":
    {
      // The user wants this to use the CIA catalog:
      mtype = TYPE_CIA;
      misotopologue = -1;

      // We have to read in the second species, and the dataset index
      n    = def.find('-');    // find the '-'

      if (n == def.npos )
        {
          ostringstream os;
          os << "Invalid species tag " << def_original << ".\n"
             << "I am missing a minus sign (and a dataset index after that.)";
          throw runtime_error(os.str());
        }

      String otherspec = def.substr(0,n);    // Extract before '-'
      def.erase(0,n+1);                      // Remove from def

    
      mcia_second = species_index_from_species_name(otherspec);
      
      if ( 0 > mcia_second )
        {
          ostringstream os;
          os << "CIA species \"" << otherspec << "\" is not a valid species.";
          throw runtime_error(os.str());
        }

      // Convert remaining def to dataset index.
      
      // Check that everything remaining is just numbers.
      for (Index i=0; i<def.nelem(); ++i)
          if (!isdigit(def[i])) {
              ostringstream os;
              os << "Invalid species tag " << def_original << ".\n"
                 << "The tag should end with a dataset index";
              throw runtime_error(os.str());
          }
      
      // Do the conversion from string to index:
      istringstream is(def);
      is >> mcia_dataset;

      def = "";
    }
  else
    {
      // Make an array containing the isotopologue names:
      ArrayOfString ins;
      for ( Index i=0; i<spr.Isotopologue().nelem(); ++i )
        ins.push_back( spr.Isotopologue()[i].Name() );

      misotopologue = find_first (ins, isoname);

      // Check if we found a matching isotopologue:
      if ( misotopologue < 0 ) 
        {
          ostringstream os;
          os << "Isotopologue " << isoname << " is not a valid isotopologue or "
             << "absorption model for species " << name << ".\n"
             << "Valid options are:\n";
          for ( Index i=0; i<ins.nelem(); ++i )
            os << name << "-" << ins[i] << "\n";
          throw runtime_error(os.str());
        }
      
      // Check if the found isotopologue represents a predefined model
      // (continuum or full absorption model) and set the type accordingly:
      if ( !isdigit(isoname[0]) )
          mtype = TYPE_PREDEF;
    }

  if ( 0 == def.nelem() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants no frequency limits.  Frequency defaults are
      // already set, so we can return directly.

      return;
    }

  if (def[0] != '*' && !isdigit(def[0]))
    {
      ostringstream os;
      os << "Expected frequency limits, but got \"" << def << "\"";
      throw runtime_error(os.str());
    }

  // Look for the two frequency limits:

  // Extract first frequency
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      // Frequency as a String:
      String fname;
      fname = def.substr(0,n);              // Extract before '-'
      def.erase(0,n+1);               // Remove from def

      // Check for joker:
      if ( "*" == fname )
        {
          // The default for mlf is already -1, meaning `ALL'.
          // So there is nothing to do here.
        }
      else if (!isdigit(fname[0]))
        {
          ostringstream os;
          os << "Expected frequency limit, but got \"" << fname << "\"";
          throw runtime_error(os.str());
        }
      else
        {
          // Convert to Numeric:
          char *endptr;
          mlf = strtod(fname.c_str(), &endptr);
          if (endptr != fname.c_str() + fname.nelem())
            {
              ostringstream os;
              os << "Error parsing frequency limit \"" << fname << "\"";
              throw runtime_error(os.str());
            }
        }
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In this
      // case that is not allowed!
      throw runtime_error("You must either specify both frequency limits\n"
                          "(at least with jokers), or none.");
    }


  // Now there should only be the upper frequency left in def.
  // Check for joker:
  if ( "*" == def )
    {
      // The default for muf is already -1, meaning `ALL'.
      // So there is nothing to do here.
    }
  else if (!isdigit(def[0]))
    {
      ostringstream os;
      os << "Expected frequency limit, but got \"" << def << "\"";
      throw runtime_error(os.str());
    }
  else
    {
      // Convert to Numeric:
      char *endptr;
      muf = strtod(def.c_str(), &endptr);
      if (endptr != def.c_str() + def.nelem())
        {
          ostringstream os;
          os << "Error parsing frequency limit \"" << def << "\"";
          throw runtime_error(os.str());
        }
    }
}


//! Return the full name of the tag.
/*! 
  Examples:

  \verbatim
  O3-*-*-*         : All O3 lines
  O3-nl            : O3, but without any lines
  O3-666-*-*       : All O3-666 lines
  O3-*-500e9-501e9 : All O3 lines between 500 and 501 GHz.
  \endverbatim

  \return The tag name as a string.
*/
String SpeciesTag::Name() const 
{
    // Species lookup data:
    using global_data::species_data;
    // A reference to the relevant record of the species data:
    const  SpeciesRecord& spr = species_data[mspecies];
    // For return value:
    ostringstream os;
    
    // First the species name:
    os << spr.Name() << "-";
    
    // Is this a CIA tag?
    if (mtype==TYPE_CIA)
      {
        os << "CIA-"
           << species_name_from_species_index(mcia_second) << "-"
           << mcia_dataset;
        
      }
    else if (mtype == TYPE_FREE_ELECTRONS || mtype == TYPE_PARTICLES)
      {
        os << spr.Name();
      }
        // Hitran Xsec flag.
    else if (mtype == TYPE_HITRAN_XSEC)
      {
        os << "HXSEC";
      }
    else
      {
        // Zeeman flag.
        if (mtype == TYPE_ZEEMAN) os << "Z-";


        // Line Mixing Type
        if (mline_mixing != LINE_MIXING_OFF)
        {
            os << "LM-";
        }

        // Now the isotopologue. Can be a single isotopologue or ALL.
        if ( misotopologue == spr.Isotopologue().nelem() )
          {
            // One larger than allowed means all isotopologues!
            os << "*-";
          }
        else if ( misotopologue == -1 )
          {
            // -1 means no lines!
            os << "nl-";
          }
        else
          {
            os << spr.Isotopologue()[misotopologue].Name() << "-";
          }

        // Now the frequency limits, if there are any. For this we first
        // need to determine the floating point precision.
        
        // Determine the precision, depending on whether Numeric is double
        // or float:
        int precision;
#ifdef USE_FLOAT
        precision = FLT_DIG;
#else
#ifdef USE_DOUBLE
        precision = DBL_DIG;
#else
#error Numeric must be double or float
#endif
#endif
        
        if ( 0 > mlf )
          {
            // mlf < 0 means no lower limit.
            os << "*-";
          }
        else
          {
            os << setprecision(precision);
            os << mlf << "-";
          }
        
        if ( 0 > muf )
          {
            // muf < 0 means no upper limit.
            os << "*";
          }
        else
          {
            os << setprecision(precision);
            os << muf;
          }
      }
    return os.str();
}

ostream& operator << (ostream& os, const SpeciesTag& ot)
{
  return os << ot.Name();
}

//! Return the name of a tag group as a string. 
/*!
   A tag group consists of several elementary SpeciesTags. This function
   returns a String with the name of the entire tag group. This is
   nice for informational output messages, for example in the
   absorption routines.

   \param  tg  The tag group in question.
   \return The full name of the tag group, as it could occur in the controlfile.

   \author Stefan Buehler
   \date   2001-03-13
*/
String get_tag_group_name( const ArrayOfSpeciesTag& tg )
{
  String name = "";
  Index i;
  
  for ( i=0; i<tg.nelem()-1; ++i )
    {
      name += tg[i].Name() + ", ";
    }
  name += tg[i].Name();

  return name;
}

//! Return the species of a tag group as a string
/*!
   A tag group consists of several elementary SpeciesTags, which must
   all belong to the same molecular species. This function
   returns a string with the name of the species. This is
   nice for informational output messages, for example in the
   absorption routines.

   E.g., if the tag group is: "H2O-161, H2O-181", then the function
   will return "H2O".

   It also does a safety check that really all tags belong to the same
   species.

   \param  tg  The tag group in question.
   \return The name of the species, as it should be used in gridded
           atmospheric fields, for example.

   \author Stefan Buehler
   \date   2009-06-11
*/
String get_species_name( const ArrayOfSpeciesTag& tg )
{
  // Get species index of first tag:
  Index spec_ind = tg[0].Species();

  // As a safety check, make sure that all other species indices are
  // the same:
  for (Index i=1; i<tg.nelem(); ++i)
    {
      //      out1 << spec_ind << " " << tg[i].Species() << "\n";
      if (tg[i].Species() != spec_ind)
        {
          ostringstream os;
          os << "All tags in a tag group must belong to the same species!\n"
             << "The offending tag group is: " << get_tag_group_name(tg);
          throw runtime_error( os.str() );
        }
    }

  return species_name_from_species_index( spec_ind );
}

//! Find first occurrence of species in tag groups.
/*! 
  The species to look for must be specified by its species index, not
  by the name. Use the helper function to get the species index from
  the species name if necessary.

  \see species_index_from_species_name.

  \param tgs The species tags to search in.
  \param spec The species index of the species to look for.
  
  \return The index of spec in tgs, -1 if not found.

  \author Stefan Buehler
  \date   2003-01-13
*/
Index find_first_species_tg( const ArrayOfArrayOfSpeciesTag& tgs,
                             const Index& spec )
{
  return find_next_species_tg(tgs, spec, 0);
}

//! Find next occurrence of species in tag groups.
/*! 
  The species to look for must be specified by its species index, not
  by the name. Use the helper function to get the species index from
  the species name if necessary.

  \see species_index_from_species_name.

  \param tgs The species tags to search in.
  \param spec The species index of the species to look for.
  \param start The index position at which to start the 
               search (0 would be the very beginning).

  \return The index of spec in tgs, -1 if not found.

  \author Stefan Buehler
  \date   2007-11-16
*/
Index find_next_species_tg( const ArrayOfArrayOfSpeciesTag& tgs,
                            const Index& spec,
                            const Index& start )
{
  for ( Index i=start; i<tgs.nelem(); ++i )
    {
      // We compare the given species index spec to the index of the
      // first element in each tag group. (All elements of a group
      // must belong to the same species.)
      //
      // If they match, then this i is the index of the tag group we
      // want. 
      if ( spec == tgs[i][0].Species() )
        return i;
    }

  // If we get here, then spec did not match the species of any of the
  // tag groups.
  return -1;
}


//! Converts a String to ArrayOfSpeciesTag
/*!
   This function is used when preparing strings read from e.g. control
   files to be stored as SpeciesTag in abs_species.
   
   Note: This is originally a part of abs_speciesSet.
   
   \param tags  Array of SpeciesTag.
   \param names String with species.
   
   \author Mattias Ekstroem
   \date   2004-09-30
*/   
void array_species_tag_from_string( ArrayOfSpeciesTag& tags, 
                                    const String& names )
{
  // There can be a comma separated list of tag definitions, so we
  // need to break the String apart at the commas.
  ArrayOfString tag_def;

  bool go_on = true;
  String these_names = names;
  while (go_on)
  {
    //          Index n = find_first( these_names, ',' );
    Index n = these_names.find(',');
    if ( n == these_names.npos ) // Value npos indicates not found.
    {
      // There are no more commas.
      //              cout << "these_names: (" << these_names << ")\n";
      tag_def.push_back(these_names);
      go_on = false;
    }
    else
    {
      tag_def.push_back(these_names.substr(0,n));
      these_names.erase(0,n+1);
    }
  }
  // tag_def now holds the different tag Strings for this group.
  
  // Set size to zero, in case the method has been called before.
  tags.resize(0);
    
  for ( Index s=0; s<tag_def.nelem(); ++s )
  {
    SpeciesTag this_tag(tag_def[s]);

    // Safety checks:
    if (s>0)
      {
        // Tags inside a group must belong to the same species.
        if ( tags[0].Species() != this_tag.Species() )
            throw runtime_error("Tags in a tag group must belong to the same species.");

        // Zeeman tags and plain line by line tags must not be mixed. (Because
        // there can be only one line list per tag group.)
        if (
            ((tags[0].Type()==SpeciesTag::TYPE_ZEEMAN) &&
             (this_tag.Type()==SpeciesTag::TYPE_PLAIN))
            ||
            ((tags[0].Type()==SpeciesTag::TYPE_PLAIN) &&
             (this_tag.Type()==SpeciesTag::TYPE_ZEEMAN))
            )
            throw runtime_error("Zeeman tags and plain line-by-line tags must "
                                "not be mixed in the same tag group.");
      }
    
    tags.push_back(this_tag);
  }
}


//! Check the correctness of abs_species.
/*!
   Checks on the correctness of the tags will be performed,
   e.g. free_electrons and particles species are only allowed once in
   abs_species.

   \param tags  Array of Array of SpeciesTag.

   \author Oliver Lemke
   \date   2013-04-23
*/
void check_abs_species( const ArrayOfArrayOfSpeciesTag& abs_species )
{
    Index num_free_electrons = 0;
    for ( Index i=0; i<abs_species.nelem(); ++i )
    {
        bool has_free_electrons = false;
        bool has_particles = false;
        bool has_hitran_xsec = false;
        for ( Index s=0; s<abs_species[i].nelem(); ++s )
        {
            if (abs_species[i][s].Type() == SpeciesTag::TYPE_FREE_ELECTRONS)
            {
                num_free_electrons++;
                has_free_electrons = true;
            }

            if (abs_species[i][s].Type() == SpeciesTag::TYPE_PARTICLES)
            {
              has_particles = true;
            }

            if (abs_species[i][s].Type() == SpeciesTag::TYPE_HITRAN_XSEC)
            {
                has_hitran_xsec = true;
            }
        }

        if (abs_species[i].nelem() > 1 && has_free_electrons)
            throw std::runtime_error("'free_electrons' must not be combined "
                                "with other tags in the same group.");
        if (num_free_electrons > 1)
            throw std::runtime_error("'free_electrons' must not be defined "
                                "more than once.");

        if (abs_species[i].nelem() > 1 && has_particles)
            throw std::runtime_error("'particles' must not be combined "
                                     "with other tags in the same group.");

        if (abs_species[i].nelem() > 1 && has_hitran_xsec)
          throw std::runtime_error("'hitran_xsec' must not be combined "
                                   "with other tags in the same group.");
    }
}


/**
 Is this a Zeeman tag group?

 It is not enough to just look at the first tag to find out if this is a Zeeman
 tag group, since there could be some predefined absorption or CIA tag first.
 
 \return true if this is a Zeeman tag group, otherwise false.
 \param  tg  The tag group to check
 
 \author Stefan Buehler
 \date   2013-01-15 */
bool is_zeeman(const ArrayOfSpeciesTag& tg)
{
    for ( Index s=0; s<tg.nelem(); ++s )
        if (tg[s].Type()==SpeciesTag::TYPE_ZEEMAN)
            return true;
        
    return false;
}


/** Returns the index of the tag group tg2 within the array of tag
    groups tgs1. Slightly modified copy of get_tagindex_for_Strings.
   
    @exception runtime_error  Could not find tg2 in tgs1.
   
    \retval tgs1_index     Index in tgs1 for tg2
    \param  tgs1           The tags groups to search in.
    \param  tg2            The tag group for which the index shall be found.
    
    \author Patrick Eriksson, Axel von Engeln, and Stefan Buehler
    \date 2001-01-31 */
void get_tag_group_index_for_tag_group( Index&               tgs1_index, 
                                        const ArrayOfArrayOfSpeciesTag&      tgs1, 
                                        const ArrayOfSpeciesTag&  tg2 )
{
  bool found = false;

  for ( Index i=0;
        i<tgs1.nelem() && !found;
        ++i )
    {
      // Is at least the size correct?
      if ( tg2.nelem() == tgs1[i].nelem() )
        {
          bool ok = true;

          for ( Index j=0; j<tg2.nelem(); ++j )
            {
              if ( tg2[j].Name() != tgs1[i][j].Name() )
                ok = false;
            }

          if ( ok )
            {
              found = true;
              tgs1_index = i;
            }
        }
    }

  if ( !found )
    {
      ostringstream os;
      os << "The tag String \"" << tg2 << 
        "\" does not match any of the given tags.\n";
      throw runtime_error(os.str());
    }
}

