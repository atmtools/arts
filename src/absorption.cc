/* Copyright (C) 2002 Stefan Buehler  <sbuehler@uni-bremen.de>

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
  \file   absorption.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Thu May 30 15:57:39 2002
  
  \brief  Physical absorption routines. 
*/

#include <cfloat>
#include <map>
#include "absorption.h"


/** The map associated with species_data. */
std::map<String, Index> SpeciesMap;


//----------------------------------------------------------------------
// Functions should be moved over here one by one from the old file
// old_absorption.cc.
//----------------------------------------------------------------------

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
  // The species map. This is used to find the species id.
  extern std::map<String, Index> SpeciesMap;
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
      ostringstream os;
      os << "Species " << name << " is not a valid species.";
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

      // Use the find algorithm from the STL to find the isotope ID. It
      // returns an iterator, so to get the index we take the
      // difference to the begin() iterator.
      misotope = find( ins.begin(),
                       ins.end(),
                       isoname ) - ins.begin();

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
      throw runtime_error("You must either speciefy both frequency limits\n"
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
  default: assert(false); exit(1); // It must be either float or double.
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

//! Print the name of a tag group. 
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


