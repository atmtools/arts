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
  \file   abs_species_tags.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue May 31 17:18:22 2005
  
  \brief  Header file for stuff related to absorption species tags.
  
  This file contains class definition and function headers related to
  SpeciesTags. It is better to separate this from the other absorption
  functions, since this part was actually improved in arts-1-1 and
  should be kept, whereas most other absorption stuff is back-ported
  from arts-1-0.
*/

#ifndef abs_species_h
#define abs_species_h

#include <stdexcept>
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "make_array.h"
#include "bifstream.h"

/** A tag group can consist of the sum of several of these.

    \author Stefan Buehler */
class SpeciesTag {
public:
  /** Default constructor. */
  SpeciesTag() { /* Nothing to be done here. */ }

  // Documentation is with implementation.
  SpeciesTag(String def); 

  // Documentation is with implementation.
  String Name() const;
    
  /** Molecular species index. */
  Index Species() const { return mspecies; }

  /** Isotopic species index.
      If this is equal to the number of isotopes (one more than
      allowed) it means all isotopes of this species. */ 
  Index Isotope() const { return misotope; }

  /** The lower line center frequency in Hz.
      If this is <0 it means no lower limit. */
  Numeric Lf() const { return mlf; }

  /** The upper line center frequency in Hz:
      If this is <0 it means no upper limit. */
  Numeric Uf() const { return muf; }

  //! Comparison operator for species tags.
  /*!
    This returns false as soon as a singe discrepancy is
    detected. Otherwise it returns true at the end.
  
    \param other The other tag to compare to.
  
    \return true if the two tags are equal.
    
    \author Stefan Buehler
    \date   2002-11-29
  */
  bool operator==(const SpeciesTag& other) const
  {
    if ( other.mspecies != mspecies ) return false;
    if ( other.misotope != misotope ) return false;
    if ( other.mlf      != mlf      ) return false;
    if ( other.muf      != muf      ) return false;
    return true;
  }

private:

  //! Molecular species index.
  Index mspecies;

  //! Isotopic species index.
  /*!
    If this is equal to the number of isotopes (one more than
    allowed) it means all isotopes of this species. If it is <0 it
    means no isotope (no lines), corresponding to "H2O-nl" */
  Index misotope;

  //! The lower limit line center frequency in Hz.
  /*! If this is <0 it means no lower limit. */
  Numeric mlf;

  //! The upper line center frequency in Hz.
  /*! If this is <0 it means no upper limit. */
  Numeric muf;
};


/** Output operator for SpeciesTag. 

    \author Stefan Buehler */
ostream& operator << (ostream& os, const SpeciesTag& ot);


/** A tag group is an array of SpeciesTags. This corresponds to one
    "species" in the controlfile. Example: "O3-666, O3-668"

    \author Stefan Buehler */
typedef  Array<SpeciesTag> ArrayOfSpeciesTag;

/** Contains the available tag groups. Contrary to the Bredbeck
    definition, tag groups may only consist of tags belonging to the
    same species. The reason for this is that there is one VMR profile
    associated with each tag group.

    \author Stefan Buehler */
typedef  Array<ArrayOfSpeciesTag> ArrayOfArrayOfSpeciesTag;


//======================================================================
//             Functions related to species and tags
//======================================================================

String get_tag_group_name( const Array<SpeciesTag>& tg );

Index species_index_from_species_name( String name );

Index find_first_species_tg( const ArrayOfArrayOfSpeciesTag& tgs,
                             const Index& spec );

void array_species_tag_from_string( ArrayOfSpeciesTag& tags,
                                    const String& names );


//=====================================================================
//           Definitions for species_index
//=====================================================================

#define SPECIES_INDEX_N2 0
#define SPECIES_INDEX_O2 1
#define SPECIES_INDEX_H2O 2
#define SPECIES_INDEX_O3 3
#define SPECIES_INDEX_CO2 4
#define SPECIES_INDEX_COUNT 5


//--------------------------------------------------------------------------------
// Functions from ARTS-1-0. Are they still needed?
//--------------------------------------------------------------------------------


void get_tagindex_for_Strings( 
              ArrayOfIndex&   tags1_index, 
        const ArrayOfArrayOfSpeciesTag&      tags1, 
        const ArrayOfString&  tags2_Strings );

void get_tag_group_index_for_tag_group( 
              Index&         tags1_index, 
        const ArrayOfArrayOfSpeciesTag&      tags1, 
        const Array<SpeciesTag>&  tags2 );


#endif // abs_species_h
