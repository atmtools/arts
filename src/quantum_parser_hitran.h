/* Copyright (C) 2015, The ARTS Developers.

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

/** \file

  Parser for quantum numbers from HITRAN 2004 and later.

  \author Oliver Lemke
*/

#ifndef quantum_parser_hitran_h
#define quantum_parser_hitran_h

#include "quantum_parser.h"

//! Parser for quantum number strings from HITRAN 2004 and later
class QuantumParserHITRAN2004 {
 public:
  //! Constructor initializing the parser
  QuantumParserHITRAN2004();

  //! Parse quantum numbers from string
  /**
     \param[in,out] qid (out) Quantum numbers extracted from the
                 string (in) must have defined Species()
     \param[in]  quantum_string Quantum number string as found
                 in the HITRAN catalog (Length: 15*4 characters).
    */
  void Parse(QuantumIdentifier& qid, const String& quantum_string) const;

 private:
  // Internal type definitions

  class QuantumClassGroup {
   public:
    QuantumClassGroup() : iclass(-1), igroup(-1){};

    Index iclass;
    Index igroup;
  };

  typedef enum {
    CI_CLASS1 = 0,
    CI_CLASS2,
    CI_CLASS3,
    CI_CLASS4,
    CI_CLASS5,
    CI_CLASS6,
    CI_CLASS7,
    CI_CLASS8,
    CI_CLASS9,
    CI_CLASS10,
    CI_FINAL
  } ClassIds;

  typedef enum {
    GI_GROUP1 = 0,
    GI_GROUP2,
    GI_GROUP3,
    GI_GROUP4,
    GI_GROUP5,
    GI_GROUP6,
    GI_GROUP6OH,
    GI_FINAL,
    GI_UNDEFINED
  } GroupIds;

  typedef struct {
    Array<QuantumFieldDescription> upper;
    Array<QuantumFieldDescription> lower;
  } QuantumGroup;

  typedef Array<QuantumFieldDescription> QuantumClass;

  // Internal functions

  void SetClassGroup(const String& species_name,
                     const ClassIds iclass,
                     const GroupIds igroup);

  // Member variables

  Array<QuantumClass> mclass;
  Array<QuantumGroup> mgroup;
  Array<QuantumClassGroup> mspecies;
};

enum class
    QuantumNumberTypeLabelsHitran : Index {  // Note these are from comparing HITRAN par-format with new labeled format...
      O2_X_is_X,  // Lambda in state 0, spin is 1
      O2_X_is_a,  // Lambda in state 2, spin is 0
      O2_X_is_b,  // Lambda in state 0, spin is 0
      NO_X_is_X,  // Lambda in state 1, spin is 0.5, Omega is Omega
      OH_X_is_X,  // Lambda in state 1, spin is 0.5, Omega is Omega
      OH_X_is_A,  // Lambda in state 0, spin is 0.5, if value at Omega is 2 then N = J + S, if value at Omega is 1 then N = J - S,  Omega is itself undefined.  This is a poorly documented HITRAN feature...
      ClO_X_is_X,  // Lambda in state 1, spin is 0.5, Omega is Omega
    };

#endif /* quantum_parser_hitran_h */
