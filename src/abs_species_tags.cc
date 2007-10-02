/* Copyright (C) 2002-2007 Stefan Buehler  <sbuehler@ltu.se>

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

//! Constructor from a tag definition String. 
/*! 
  For examples see documentation of member function Name(). 

  \param def String containing tag definition.
  
  \exception runtime_error The given String could not be mapped to
  a sensible tag description.
*/
SpeciesTag::SpeciesTag(String def) 
{
  // Species lookup data:
  extern const Array<SpeciesRecord> species_data;
  // Name of species and isotope (aux variables):
  String name, isoname;
  // Aux index:
  Index n;

  // Set frequency limits to default values (no limits):
  mlf = -1;
  muf = -1;

  // We cannot set a default value for the isotope, because the
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
      // the user wants all isotopes and no frequency limits.
      // Frequency defaults are already set. Set isotope defaults:
      misotope = spr.Isotope().nelem();
      // This means all isotopes.

      return;
    }
    
  // Extract the isotope name:
  n    = def.find('-');    // find the '-'
  if (n != def.npos )
    {
      isoname = def.substr(0,n);          // Extract before '-'
      def.erase(0,n+1);             // Remove from def
    }
  else
    {
      // n==def.npos means that def does not contain a '-'. In that
      // case we assume that it contains just the isotope name and
      // nothing else.
      isoname = def;
      def  = "";
    }
    
  // Check for joker:
  if ( "*" == isoname )
    {
      // The user wants all isotopes. Set this accordingly:
      misotope = spr.Isotope().nelem();
    }
  else if ( "nl" == isoname )     // Check for "nl":
    {
      // The user wants no lines at all. Set this accordingly:
      misotope = -1;
    }
  else
    {
      // Make an array containing the isotope names:
      ArrayOfString ins;
      for ( Index i=0; i<spr.Isotope().nelem(); ++i )
        ins.push_back( spr.Isotope()[i].Name() );

      misotope = find_first (ins, isoname);

      // Check if we found a matching isotope:
      if ( misotope >= ins.nelem() ) 
        {
          ostringstream os;
          os << "Isotope " << isoname << " is not a valid isotope for "
             << "species " << name << ".\n"
             << "Valid isotopes are:";
          for ( Index i=0; i<ins.nelem(); ++i )
            os << " " << ins[i];
          throw runtime_error(os.str());
        }
    }

  if ( 0 == def.nelem() )
    {
      // This means that there is nothing else to parse. Apparently
      // the user wants no frequency limits.  Frequency defaults are
      // already set, so we can return directly.

      return;
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
      else
        {
          // Convert to Numeric:
          istringstream is(fname);
          is >> mlf;
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
  else
    {
      // Convert to Numeric:
      istringstream is(def);
      is >> muf;
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
  extern const Array<SpeciesRecord> species_data;
  // A reference to the relevant record of the species data:
  const  SpeciesRecord& spr = species_data[mspecies];
  // For return value:
  ostringstream os;

  // First the species name:
  os << spr.Name() << "-";

  // Now the isotope. Can be "nl", a single isotope or ALL.
  if ( misotope == spr.Isotope().nelem() )
    {
      // One larger than allowed means all isotopes!
      os << "*-";
    }
  else if ( misotope == -1 )
    {
      // -1 means no lines!
      os << "nl-";
    }
  else
    {
      os << spr.Isotope()[misotope].Name() << "-";
    }

  // Now the frequency limits, if there are any. For this we first
  // need to determine the floating point precision.

  // Determine the precision, depending on whether Numeric is double
  // or float:  
  Index precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: assert(false); arts_exit (); // It must be either float or double.
  }

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
  for ( Index i=0; i<tgs.nelem(); ++i )
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

//! Return species index for given species name.
/*! 
  This is useful in connection with other functions that need a species
  index.

  \see find_first_species_tg.

  \param name Species name.

  \return Species index, -1 means not found.

  \author Stefan Buehler
  \date   2003-01-13
*/
Index species_index_from_species_name( String name )
{
  // The species map. This is used to find the species id.
  extern map<String, Index> SpeciesMap;

  // For the return value:
  Index mspecies;
  
  // Remove leading whitespace, if there is any:
  while ( 0 != name.nelem() && (
                                ' '  == name[0] ||
                                '\t' == name[0] ||
                                '\n' == name[0] ||
                                '\r' == name[0]
                                )
          )    name.erase(0,1);

  // Remove trailing whitespace, if there is any:
  while ( 0 != name.nelem() && (
                                ' '  == name[name.nelem()-1] ||
                                '\t' == name[name.nelem()-1] ||
                                '\n' == name[name.nelem()-1] ||
                                '\r' == name[name.nelem()-1]
                                )
          )    name.erase(name.nelem()-1);

  //  cout << "name / def = " << name << " / " << def << endl;

  // Look for species name in species map:
  map<String, Index>::const_iterator mi = SpeciesMap.find(name);
  if ( mi != SpeciesMap.end() )
    {
      // Ok, we've found the species. Set mspecies.
      mspecies = mi->second;
    }
  else
    {
      // The species does not exist!
      mspecies = -1;
    }

  return mspecies;
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

    // Safety check: For s>0 check that the tags belong to the same species.
    if (s>0)
    if ( tags[0].Species() != this_tag.Species() )
      throw runtime_error("Tags in a tag group must belong to the same species.");

    tags.push_back(this_tag);
  }
}
   
/** 
   Returns the index among some tag groups for an array of tag Strings. 
   
   \verbatim
   For example, if tags1 correspond to the definition
     ["O3","H2O-161,H2O-162"]
   and the tag Strings are
     ["H2O-161,H2O-162","O3"]
   the tags1_index becomes
     [2,1]
   \endverbatim

   @exception runtime_error  Some String is not a valid tag item.
   @exception runtime_error  Not all Strings are not found among the tags.

   \retval tags1_index     Index in tags1 for tags2_Strings
   \param  tags1           The tags to search in.
   \param  tags2_Strings   The tag Strings for which indeces shall be found.

   \author Patrick Eriksson 
   \date 2000-12-06
*/
void get_tagindex_for_Strings( 
                              ArrayOfIndex&   tags1_index, 
                              const ArrayOfArrayOfSpeciesTag&      tags1, 
                              const ArrayOfString&  tags2_Strings )
{
  const Index   n1 = tags1.nelem();
  const Index   n2 = tags2_Strings.nelem();
  ArrayOfArrayOfSpeciesTag   tags2;                // Internal tag names for tag_Strings
  Index   i1, i2, nj, j, found, ok;

  tags1_index.resize(n2);
  //  cout << "tags2_Strings: " << tags2_Strings << "\n";
  SpeciesSet( tags2, "tags2", tags2_Strings );

  for ( i2=0; i2<n2; i2++ )
  {
    found = 0;
    for ( i1=0; (i1<n1) && !found; i1++ )
    {
      nj = tags2[i2].nelem(); 
      if ( nj  == tags1[i1].nelem() )
      {
        ok = 1;
        for ( j=0; j<nj; j++ )
        {
          if ( tags2[i2][j].Name() != tags1[i1][j].Name() )
            ok = 0;
        }
        if ( ok )
        {
          found = 1;
          tags1_index[i2] = i1;
        }
      }
    } 
    if ( !found )
    {
      ostringstream os;
      os << "The tag String \"" << tags2_Strings[i2] << 
            "\" does not match any of the given tags.\n";
      throw runtime_error(os.str());
    }
  }
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

