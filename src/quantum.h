/* Copyright (C) 2013
   Oliver Lemke <olemke@core-dump.info>

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

/**
 * @file quantum.cc
 * @author Oliver Lemke
 * @date 2013-04-12
 * 
 * @brief Classes to handle quantum numbers
 */

#ifndef quantum_h
#define quantum_h

#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include "array.h"
#include "enums.h"
#include "interpolation.h"
#include "matpack.h"
#include "mystring.h"
#include "rational.h"
#include "isotopologues.h"

/** Enum for Quantum Numbers used for indexing */
ENUMCLASS(QuantumNumberType, char,
  J,       // Total angular momentum
  dJ,      // Delta total angular momentum
  M,       // Projection of J along magnetic field
  N,       // J minus spin
  dN,      // Delta J minus spin
  S,       // Spin angular momentum (from electrons) NOTE: S_global for HITRAN S
  tau,
  n,
  F,       // J + nuclear spin
  Ka,      //(This is a projection of J along one axis)
  Kc,      //(This is a projection of J along another axis)
  Omega,   // This is an absolute projection of J and S
  i,       //(Is related to Omega)
  Lambda,  // This is Sigma or Pi or Lambda states (as seen in literature)
  alpha,   // Alpha from HITRAN
  Sym,     // Symmetry expression
  parity,  // parity value (+/-)
  kronigParity,  // ???
  v1,      // Vibrational mode 1
  v2,      // Vibrational mode 2
  v3,      // Vibrational mode 3
  v4,      // Vibrational mode 4
  v5,      // Vibrational mode 5
  v6,      // Vibrational mode 6
  v7,
  v8,
  v9,
  v10,
  v11,
  v12,
  l1,      // The absolute sum of l_j for v_j
  l2,      // Vibrational angular momentum associated with v2
  l3,
  l4,
  l5,
  l6,
  l7,
  l8,
  l9,
  l10,
  l11,
  l12,
  pm,      // Symmetry type for l=0
  r,       // Rank of the level within a set of the same vibrational symmetry
  S_global,  // Symmetry of the level
  ElectronState,  // Electronic state
  n_global,  // Torosional quanta
  C,         // Another symmetry expression
  Hund  // Flag for Hund case type.  This flag lets Zeeman know what to expect
)  // QuantumNumberType

constexpr QuantumNumberType string2quantumnumbertype(const std::string_view s) {
  QuantumNumberType out = toQuantumNumberType(s);
  if (out == QuantumNumberType::FINAL) {
    if (s.find("F#") < s.length()) out = QuantumNumberType::F;  // HITRAN has many names for F
    else if (s == "K") out = QuantumNumberType::Ka;  // HITRAN name
    else if (s == "v") out = QuantumNumberType::v1;  // HITRAN name
    else if (s == "l") out = QuantumNumberType::l1;  // HITRAN name
    else if (s == "ElecStateLabel") out = QuantumNumberType::ElectronState;  // HITRAN name
  }
  return out;
}

/** Enum for Hund cases */
enum class Hund : Index { CaseA = int('a'), CaseB = int('b') };

/** Container class for Quantum Numbers */
class QuantumNumbers {
 public:
  typedef std::array<Rational, Index(QuantumNumberType::FINAL)>
      QuantumContainer;

  /** Initializer with undefined values */
  constexpr QuantumNumbers() noexcept
      : mqnumbers({RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 3
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 6
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 9
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 12
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 15
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 18
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 21
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 24
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 27
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 30
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 33
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 36
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 39
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 42
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 45
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 48
                   RATIONAL_UNDEFINED}) {}
  
  // Temporary initialization until there is a more reliable way to setup 
  // a full initialization for all quantum numbers but only choose a few
  // select ones based on the problem at hand
  constexpr QuantumNumbers(Rational J, Rational N, Rational v) noexcept
      : mqnumbers({J, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 3
                   N, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 6
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 9
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 12
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 15
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 18
                   v, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 21
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 24
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 27
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 30
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 33
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 36
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 39
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 42
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 45
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,  // 48
                   RATIONAL_UNDEFINED}) {}

  /** Access operator
   * 
   * @param[in] qn Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  constexpr const Rational& operator[](const Index qn) const noexcept {
    return mqnumbers[qn];
  }

  /** Access operator
   * 
   * @param qn[in] Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  constexpr const Rational& operator[](const QuantumNumberType qn) const noexcept {
    return mqnumbers[Index(qn)];
  }
  
  /** Access operator
   * 
   * @param[in] qn Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  constexpr Rational& operator[](const Index qn) noexcept {
    return mqnumbers[qn];
  }
  
  /** Access operator
   * 
   * @param qn[in] Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  constexpr Rational& operator[](const QuantumNumberType qn) noexcept {
    return mqnumbers[Index(qn)];
  }

  /** Set quantum number at position
   * 
   * @param[in] qn Index Pos to set at
   * @param[in] r Rational to set
   */
  constexpr void Set(Index qn, Rational r) {
    mqnumbers[qn] = r;
  }

  /** Set quantum number at position
   * 
   * @param[in] qn Index Pos to set at
   * @param[in] r Rational to set
   */
  constexpr void Set(QuantumNumberType qn, Rational r) {
    mqnumbers[Index(qn)] = r;
  }

  /** Set quantum number at position
   * 
   * @param[in] qn String Pos to set at by name
   * @param[in] r Rational to set
   */
  constexpr void Set(const std::string_view qn, Rational r)  {
    mqnumbers[Index(string2quantumnumbertype(qn))] = r;
  }

  /** Get the numbers
   * 
   * @return const QuantumContainer& All the numbers
   */
  constexpr const QuantumContainer& GetNumbers() const { return mqnumbers; }

  /** The number of defined quantum numbers
   * 
   * @return Index Count of defined numbers
   */
  constexpr Index nNumbers() const {
    Index out=0;
    for (auto& qn: mqnumbers) if (qn.isDefined()) out++;
    return out;
  }
  
  /** Check if there's any quantum numbers defined */
  constexpr bool Any() const {
    for (auto& qn: mqnumbers) if (qn.isDefined()) return true;
    return false;
  }
  
  /** Returns this as a string */
  String toString() const;

 private:
  QuantumContainer mqnumbers;
};

/** Check if all quantum numbers are the same between a and b
 * 
 * @param[in] a One set of quantum numbers
 * @param[in] b Another set of quantum numbers
 * @return true If all quantum numbers match
 * @return false Otherwise
 */
constexpr bool operator==(const QuantumNumbers& a, const QuantumNumbers& b) {
  for (Index i=0; i<Index(QuantumNumberType::FINAL); i++) {
    Rational ra = a.GetNumbers()[i];
    Rational rb = b.GetNumbers()[i];
    if (ra.isDefined() or rb.isDefined()) {
      if (ra not_eq rb) {
        return false;
      }
    }
  }
  return true;
}

/** Opposite of operator==()
 * 
 * @param[in] a One set of quantum numbers
 * @param[in] b Another set of quantum numbers
 * @return false If all quantum numbers match
 * @return true Otherwise
 */
constexpr bool operator!=(const QuantumNumbers& a, const QuantumNumbers& b) {
  return not (a == b);
}

/** Input operator */
std::istream& operator>>(std::istream& is, QuantumNumbers& qn);

/** Output operator */
std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn);

typedef Array<QuantumNumbers> ArrayOfQuantumNumbers;

namespace Quantum {
ENUMCLASS (IdentifierType, unsigned char,
           None,
           All,
           Transition,
           EnergyLevel
          )

/** Class to identify and match lines by their quantum numbers
 * 
 * Describes a transition, an energy level, all numbers or none
 * as matchable.  Useful to match to line energy levels and
 * to identify lines themselves.
 * 
 * For transitions, the QI contains upper and lower quantum numbers.
 * For energy levels, it only holds one set of quantum numbers which
 * are then matched against the upper and lower qns of the lines. For
 * all and none no quantum numbers are considered
 * 
 * File format:
 *   Transition:   SPECIES_NAME-ISOTOPE TR UP QUANTUMNUMBERS LO QUANTUMNUMBERS
 *   Energy level: SPECIES_NAME-ISOTOPE EN QUANTUMNUMBERS
 *   All lines:    SPECIES_NAME-ISOTOPE ALL
 *   No lines:    SPECIES_NAME-ISOTOPE NONE
 * 
 * Example written out:
 *   H2O-161 TR UP J 0/1 v1 2/3 LO J 1/1 v2 1/2
 *   H2O-161 EN J 0/1 v1 2/3
 *   H2O-161 ALL
 *   H2O-161 NONE
 */
struct Identifier {
  IdentifierType type;  // Type of ID
  Index spec_ind;       // Index to valid IsotopeRecord in Isotopologues, or -1
  QuantumNumbers upp;   // Upper quantum numbers, or energy level quantum numbers
  QuantumNumbers low;   // Lower quantum numbers
  
  constexpr Index SpecInd(const Species::IsotopeRecord& ir) const ARTS_NOEXCEPT {
    Index ind = Species::find_species_index(ir);
    ARTS_ASSERT(ind >= 0 and not Species::is_predefined_model(ir),
                "Must be valid, non-joker, isotopologue.  Is: ", ir)
    return ind;
  }
  
  //! Default to nothing
  constexpr Identifier() noexcept :
    type(IdentifierType::None), spec_ind(-1), upp(), low() {}
  
  constexpr Identifier(const Species::IsotopeRecord& ir, IdentifierType t=IdentifierType::All) ARTS_NOEXCEPT :
    type(t), spec_ind(SpecInd(ir)), upp(), low() {}
  
  constexpr Identifier(const Species::IsotopeRecord& ir,
                       QuantumNumbers el) ARTS_NOEXCEPT :
    type(IdentifierType::EnergyLevel), spec_ind(SpecInd(ir)),
    upp(std::move(el)), low() {}
  
  constexpr Identifier(const Species::IsotopeRecord& ir,
                       QuantumNumbers upper, QuantumNumbers lower) ARTS_NOEXCEPT :
    type(IdentifierType::Transition), spec_ind(SpecInd(ir)),
    upp(std::move(upper)), low(std::move(lower)) {}
  
  Identifier(const Species::IsotopeRecord& ir, 
             const std::vector<QuantumNumberType>& keys,
             const std::vector<Rational>& upper,
             const std::vector<Rational>& lower) ARTS_NOEXCEPT :
    type(IdentifierType::Transition), spec_ind(SpecInd(ir)), upp(), low() {
    const std::size_t n = keys.size();
    ARTS_ASSERT(n == upper.size() and n == lower.size())
    for (std::size_t i=0; i<n; i++) {
      upp[keys[i]] = upper[i];
      low[keys[i]] = lower[i];
    }
  }
  
  void SetFromString(String str);
  
  explicit Identifier(String x) { SetFromString(x); }
  
  constexpr const Species::IsotopeRecord& Isotopologue() const noexcept {
    return Species::Isotopologues[spec_ind];
  }
  
  constexpr void Isotopologue(const Species::IsotopeRecord& ir) ARTS_NOEXCEPT {
    spec_ind = SpecInd(ir);
  }
  
  constexpr Species::Species Species() const noexcept {
    return Isotopologue().spec;
  }
  
  constexpr const QuantumNumbers& Upper() const ARTS_NOEXCEPT {
    ARTS_ASSERT(IdentifierType::Transition == type)
    return upp;
  }
  
  constexpr QuantumNumbers& Upper() ARTS_NOEXCEPT {
    ARTS_ASSERT(IdentifierType::Transition == type)
    return upp;
  }
  
  constexpr const QuantumNumbers& Lower() const ARTS_NOEXCEPT {
    ARTS_ASSERT(IdentifierType::Transition == type)
    return low;
  }
  
  constexpr QuantumNumbers& Lower() ARTS_NOEXCEPT {
    ARTS_ASSERT(IdentifierType::Transition == type)
    return low;
  }
  
  constexpr const QuantumNumbers& Level() const ARTS_NOEXCEPT {
    ARTS_ASSERT(IdentifierType::EnergyLevel == type)
    return upp;
  }
  
  constexpr QuantumNumbers& Level() ARTS_NOEXCEPT {
    ARTS_ASSERT(IdentifierType::EnergyLevel == type)
    return upp;
  }
  
  constexpr bool operator==(const Identifier& other) const noexcept {
    if (other.type not_eq type) return false;
    if (not Species::same_or_joker(other.Isotopologue(), Isotopologue())) return false;
    if (other.upp not_eq upp) return false;
    if (other.low not_eq low) return false;
    return true;
  }
  
  constexpr bool operator!=(const Identifier& other) const noexcept {
    return not operator==(other);
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Identifier& id) {
    switch (id.type) {
      case IdentifierType::None: return os << "NONE"; break;
      case IdentifierType::All:
        return os << Species::Isotopologues[id.spec_ind].FullName() << " ALL";
      case IdentifierType::EnergyLevel:
        return os << Species::Isotopologues[id.spec_ind].FullName() << " EN " << id.upp;
      case IdentifierType::Transition:
        return os << Species::Isotopologues[id.spec_ind].FullName() << " TR"
                     " UP " << id.upp <<
                     " LO " << id.low;
      case IdentifierType::FINAL: { /* Leave last */
      }
    }
    
    return os;
  }
  
  /** Return a quantum identifer as if it wants to match to lower energy level */
  constexpr Identifier LowerId() const ARTS_NOEXCEPT {
    return Identifier(Species::Isotopologues[spec_ind], low);
  };
  
  /** Return a quantum identifer as if it wants to match to lower energy level */
  constexpr Identifier UpperId() const ARTS_NOEXCEPT {
    return Identifier(Species::Isotopologues[spec_ind], upp);
  };
};
}

using QuantumIdentifierType = Quantum::IdentifierType;

/*! Identifier for species, energy levels, or transitions */
using QuantumIdentifier = Quantum::Identifier;

/** List of QuantumIdentifier */
typedef Array<QuantumIdentifier> ArrayOfQuantumIdentifier;

/** Check for valid quantum number name
 * 
 * @param[in] name Parameter
 * @return true If the parameter exist
 * @return false Otherwise
 */
constexpr bool IsValidQuantumNumberName(const std::string_view name) {
  return good_enum(string2quantumnumbertype(name));
}

/** Check for valid quantum number name and throws if it is invalid
 * 
 * @param[in] name Parameter
 */
void ThrowIfQuantumNumberNameInvalid(String name);

/** Updates the quantum identifier based on a lists of strings
 * 
 * The input lists of strings should be paired as {key, value}
 * 
 * \param[in,out] qid Identifier to update
 * \param[in] upper_list List of strings to update upper state
 * \param[in] lower_list List of strings to update lower state
 */
void update_id(QuantumIdentifier& qid, const std::vector<std::array<String, 2> >& upper_list, const std::vector<std::array<String, 2> >& lower_list);

#endif
