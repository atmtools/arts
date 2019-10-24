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
#include "array.h"
#include "bifstream.h"
#include "matpackI.h"
#include "mystring.h"

/** A tag group can consist of the sum of several of these.

    \author Stefan Buehler */
class SpeciesTag {
 public:
  /** Default constructor. */
  SpeciesTag()
      : mspecies(-1),
        misotopologue(-1),
        mlf(0.),
        muf(0.),
        mtype(-1),
        mcia_second(-1),
        mcia_dataset(-1) { /* Nothing to be done here. */
  }

  // Documentation is with implementation.
  SpeciesTag(String def);

  // Documentation is with implementation.
  String Name() const;

  /** Molecular species index. */
  Index Species() const { return mspecies; }

  /** Molecular species index. */
  Index BathSpecies() const { return mcia_second; }

  /** Name of main species */
  String SpeciesNameMain() const;

  /** Mass of main species */
  Numeric SpeciesMass() const;

  /** Check if the species is same as SpeciesTag(s).Species() */
  bool IsSpecies(const String& s) const;

  /** Check if the isotopologue is same as SpeciesTag(s).Isotopologue() */
  bool IsIsotopologue(const String& s) const;

  /** Isotopologue species index.
      If this is equal to the number of isotopologues (one more than
      allowed) it means all isotopologues of this species. */
  Index Isotopologue() const { return misotopologue; }

  /** The lower line center frequency in Hz.
      If this is <0 it means no lower limit. */
  Numeric Lf() const { return mlf; }

  /** The upper line center frequency in Hz:
      If this is <0 it means no upper limit. */
  Numeric Uf() const { return muf; }

  /** Species index of the 2nd CIA species */
  Index CIASecond() const { return mcia_second; }

  /** CIA dataset index inside this CIA file. */
  Index CIADataset() const { return mcia_dataset; }

  //! Comparison operator for species tags.
  /*!
    This returns false as soon as a singe discrepancy is
    detected. Otherwise it returns true at the end.
  
    \param other The other tag to compare to.
  
    \return true if the two tags are equal.
    
    \author Stefan Buehler
    \date   2002-11-29
  */
  bool operator==(const SpeciesTag& other) const {
    if (other.mspecies != mspecies) return false;
    if (other.misotopologue != misotopologue) return false;
    if (other.mlf != mlf) return false;
    if (other.muf != muf) return false;
    if (other.mtype != mtype) return false;
    if (mtype == TYPE_CIA && (other.mcia_second != mcia_second ||
                              other.mcia_dataset != mcia_dataset))
      return false;
    return true;
  }

  /** Enum for type of this tag.

  See private member mtype for more explanations.   */
  enum {
    TYPE_PLAIN,
    TYPE_ZEEMAN,
    TYPE_PREDEF,
    TYPE_CIA,
    TYPE_FREE_ELECTRONS,
    TYPE_PARTICLES,
    TYPE_HITRAN_XSEC
  };

  /** Return the type of this tag.
   
   See private member mtype for more explanations.   */
  Index Type() const { return mtype; }

 private:
  //! Molecular species index.
  Index mspecies;

  //! Isotopologue species index.
  /*!
    If this is equal to the number of isotopologues (one more than
    allowed) it means all isotopologues of this species. If it is <0 it
    means no isotopologue (no lines), corresponding to "H2O-nl" */
  Index misotopologue;

  //! The lower limit line center frequency in Hz.
  /*! If this is <0 it means no lower limit. */
  Numeric mlf;

  //! The upper line center frequency in Hz.
  /*! If this is <0 it means no upper limit. */
  Numeric muf;

  /** Type of this tag.
   
   The type can be:
   <PRE>
   TYPE_PLAIN:          A normal line-by-line tag
   TYPE_ZEEMAN:         A line-by-line tag with Zeeman calculation
   TYPE_PREDEF:         A tag for a predefined absorption model (continuum or
                         full absorption model)
   TYPE_CIA:            A HITRAN collision induces absorption (CIA) tag
   TYPE_FREE_ELECTRONS: A free electrons tag
   TYPE_PARTICLES:      A particle tag
   TYPE_HITRAN_XSEC:    A HITRAN absorption cross section tag
   </PRE>
   */
  Index mtype;

  //! 2nd CIA species index.
  /*! Contains the species index of the second CIA species that should be used for this tag. */
  Index mcia_second;

  //! CIA dataset index.
  /*! A CIA file contains several datasets. This index specifies which one we want. */
  Index mcia_dataset;
};

/** Output operator for SpeciesTag. 

    \author Stefan Buehler */
ostream& operator<<(ostream& os, const SpeciesTag& ot);

/** A tag group is an array of SpeciesTags. This corresponds to one
    "species" in the controlfile. Example: "O3-666, O3-668"

    \author Stefan Buehler */
typedef Array<SpeciesTag> ArrayOfSpeciesTag;

/** Contains the available tag groups. Contrary to the Bredbeck
    definition, tag groups may only consist of tags belonging to the
    same species. The reason for this is that there is one VMR profile
    associated with each tag group.

    \author Stefan Buehler */
typedef Array<ArrayOfSpeciesTag> ArrayOfArrayOfSpeciesTag;

//======================================================================
//             Functions related to species and tags
//======================================================================

String get_tag_group_name(const ArrayOfSpeciesTag& tg);

String get_species_name(const ArrayOfSpeciesTag& tg);

Index find_first_species_tg(const ArrayOfArrayOfSpeciesTag& tgs,
                            const Index& spec);

Index find_next_species_tg(const ArrayOfArrayOfSpeciesTag& tgs,
                           const Index& spec,
                           const Index& start);

void array_species_tag_from_string(ArrayOfSpeciesTag& tags,
                                   const String& names);

void check_abs_species(const ArrayOfArrayOfSpeciesTag& tags);

bool is_zeeman(const ArrayOfSpeciesTag& tg);

//--------------------------------------------------------------------------------
// Functions from ARTS-1-0. Are they still needed?
//--------------------------------------------------------------------------------

void get_tag_group_index_for_tag_group(Index& tags1_index,
                                       const ArrayOfArrayOfSpeciesTag& tags1,
                                       const Array<SpeciesTag>& tags2);

#endif  // abs_species_h
