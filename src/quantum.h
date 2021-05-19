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
#include "species.h"

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
  
  /** Returns this as a string */
  String toString() const;

 private:
  QuantumContainer mqnumbers;
};

typedef Array<QuantumNumbers> ArrayOfQuantumNumbers;

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
class QuantumIdentifier {
 public:
  /** Ways to identify quantum numbers */
  typedef enum : Index { TRANSITION, ENERGY_LEVEL, ALL, NONE } QType;

  /** Initialize with no matches */
  constexpr QuantumIdentifier() noexcept
      : mqtype(QuantumIdentifier::NONE), mspecies(-1), miso(-1) {}

  /** Initialize with no quantum numbers defined but known species and matching type
   * 
   * @param[in] qt Way to identify quantum numbers
   * @param[in] species Species index-mapped
   * @param[in] iso Isotopologue index-mapped
   */
  constexpr QuantumIdentifier(const QuantumIdentifier::QType& qt,
                              const Index species,
                              const Index iso) noexcept
      : mqtype(qt), mspecies(species), miso(iso) {}

  /** Initialize with transition identifier type
   * 
   * @param[in] species Species index-mapped
   * @param[in] iso Isotopologue index-mapped
   * @param[in] upper Upper state quantum numbers
   * @param[in] lower Lower state quantum numbers
   */
  constexpr QuantumIdentifier(const Index spec,
                              const Index isot,
                              const QuantumNumbers& upper,
                              const QuantumNumbers& lower) noexcept
      : mqtype(QuantumIdentifier::TRANSITION),
        mspecies(spec),
        miso(isot),
        mqm({upper, lower}) {}

  /** Initialize with transition identifier type
   * 
   * @param[in] species Species index-mapped
   * @param[in] iso Isotopologue index-mapped
   * @param[in] ids List of quantum number keys
   * @param[in] upper Upper state quantum numbers
   * @param[in] lower Lower state quantum numbers
   */
  QuantumIdentifier(const Index spec,
                    const Index isot,
                    const std::vector<QuantumNumberType>& keys,
                    const std::vector<Rational>& upper,
                    const std::vector<Rational>& lower);

  /** Initialize with energy level identifier type
   * 
   * @param[in] species Species index-mapped
   * @param[in] iso Isotopologue index-mapped
   * @param[in] qnr Quantum numbers
   */
  constexpr QuantumIdentifier(const Index spec,
                              const Index isot,
                              const QuantumNumbers& qnr) noexcept
      : mqtype(QuantumIdentifier::ENERGY_LEVEL),
        mspecies(spec),
        miso(isot),
        mqm({qnr}) {}

  /** Construct a new Quantum Identifier object from text
   * 
   * @param[in] x Text
   */
  explicit QuantumIdentifier(String x) { SetFromString(x); }

  /** Upper level index */
  static constexpr Index TRANSITION_UPPER_INDEX = 0;

  /** Lower level index */
  static constexpr Index TRANSITION_LOWER_INDEX = 1;

  /** Energy level index */
  static constexpr Index ENERGY_LEVEL_INDEX = 0;

  /** Set the Species
   * 
   * @param[in] sp Species index-mapped
   */
  constexpr void Species(Index sp) { mspecies = sp; }

  /** Set the Isotopologue
   * 
   * @param[in] is Isotopologue index-mapped
   */
  constexpr void Isotopologue(Index iso) { miso = iso; }

  /** Set to transition type identifier
   * 
   * @param[in] upper Upper state quantum numbers
   * @param[in] lower Lower state quantum numbers
   */
  constexpr void SetTransition(const QuantumNumbers& upper, const QuantumNumbers& lower) {
    mqtype = QuantumIdentifier::TRANSITION;
    mqm[TRANSITION_UPPER_INDEX] = upper;
    mqm[TRANSITION_LOWER_INDEX] = lower;
  }

  /** Set to energy level identifier
   * 
   * @param[in] q Quantum numbers
   */
  constexpr void SetEnergyLevel(const QuantumNumbers& q) {
    mqtype = QuantumIdentifier::ENERGY_LEVEL;
    mqm[ENERGY_LEVEL_INDEX] = q;
  }

  /** Set to All identifier */
  constexpr void SetAll() { mqtype = QuantumIdentifier::ALL; };
  
  /** Set to NONE identifier */
  constexpr void SetNone() { mqtype = QuantumIdentifier::NONE; };

  /** Set key to transition type */
  constexpr void SetTransition() { mqtype = QuantumIdentifier::TRANSITION; };

  /** Set from a String object
   * 
   * @param[in] str The string to set this from
   */
  void SetFromString(String str);

  /** Set CO2 transition from String objects
   * 
   * @param[in] upper Upper state quantum numbers
   * @param[in] lower Lower state quantum numbers
   * @param[in] iso Isotopologue by string
   */
  void SetFromStringForCO2Band(String upper, String lower, String iso);

  /** @return constexpr QType Type of identifier */
  constexpr QType Type() const { return mqtype; }
  
  /** Set Type */
  constexpr void Type(QType x) { mqtype = x; }
  
  /** Checks if input is a valid Type */
  static constexpr bool validIndexForType(Index x) noexcept {
    return x == Index(TRANSITION) or x == Index(ENERGY_LEVEL) or x == Index(ALL) or x == Index(NONE);
  }
  
  /** @return QType if string is a Type or -1 if not */
  static constexpr QType string2Type(const std::string_view str) noexcept  {
    if ("ENERGY_LEVEL" == str) {
      return QuantumIdentifier::ENERGY_LEVEL;
    } else if ("TRANSITION" == str) {
      return QuantumIdentifier::TRANSITION;
    } else if ("ALL" == str) {
      return QuantumIdentifier::ALL;
    } else if ("NONE" == str) {
      return QuantumIdentifier::NONE;
    } else {
      return QType(-1);
    }
  }

  /** @return QType as String */
  constexpr std::string_view TypeStr() const noexcept {
    switch (mqtype) {
      case QuantumIdentifier::TRANSITION:   return "TR";
      case QuantumIdentifier::ENERGY_LEVEL: return "EN";
      case QuantumIdentifier::ALL:          return "ALL";
      case QuantumIdentifier::NONE:         return "NONE";
      default: return "There's an error";
    }
  }

  /** Return the Species by name */
  String SpeciesName() const;

  /** Return the Species by index */
  constexpr Index Species() const noexcept { return mspecies; }
  
  /** Return the Species mass */
  Numeric SpeciesMass() const;

  /** Return the Species by index reference */
  constexpr Index& Species() noexcept { return mspecies; }

  /** Return the Isotopologue by index */
  constexpr Index Isotopologue() const noexcept { return miso; }

  /** Return the Isotopologue by index reference */
  constexpr Index& Isotopologue() noexcept { return miso; }

  /** Return the quantum numbers array const reference */
  constexpr const std::array<QuantumNumbers, 2>& QuantumMatch() const noexcept { return mqm; }

  /** Return the quantum numbers array reference */
  constexpr std::array<QuantumNumbers, 2>& QuantumMatch() noexcept { return mqm; }

  /** Return a quantum identifer as if it wants to match to upper energy level */
  constexpr QuantumIdentifier UpperQuantumId() const ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == TRANSITION);
    return QuantumIdentifier(mspecies, miso, mqm[TRANSITION_UPPER_INDEX]);
  };

  /** Return a quantum identifer as if it wants to match to lower energy level */
  constexpr QuantumIdentifier LowerQuantumId() const ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == TRANSITION);
    return QuantumIdentifier(mspecies, miso, mqm[TRANSITION_LOWER_INDEX]);
  };

  /** Return the upper quantum numbers by const reference */
  constexpr const QuantumNumbers& UpperQuantumNumbers() const noexcept {
    return mqm[TRANSITION_UPPER_INDEX];
  };

  /** Return the lower quantum numbers by const reference */
  constexpr const QuantumNumbers& LowerQuantumNumbers() const noexcept {
    return mqm[TRANSITION_LOWER_INDEX];
  };

  /** Return a upper quantum number by copy */
  constexpr Rational UpperQuantumNumber(QuantumNumberType X) const ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == TRANSITION);
    return mqm[TRANSITION_UPPER_INDEX][X];
  };

  /** Return a lower quantum number by copy */
  constexpr Rational LowerQuantumNumber(QuantumNumberType X) const ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == TRANSITION);
    return mqm[TRANSITION_LOWER_INDEX][X];
  };
  
  /** Return a upper quantum number by copy */
  constexpr Rational& UpperQuantumNumber(QuantumNumberType X) ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == TRANSITION);
    return mqm[TRANSITION_UPPER_INDEX][X];
  };
  
  /** Return a lower quantum number by copy */
  constexpr Rational& LowerQuantumNumber(QuantumNumberType X) ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == TRANSITION);
    return mqm[TRANSITION_LOWER_INDEX][X];
  };

  /** Return the energy level quantum numbers by const reference */
  constexpr const QuantumNumbers& EnergyLevelQuantumNumbers() const noexcept {
    return mqm[ENERGY_LEVEL_INDEX];
  }
  
  /** Return a energy level quantum number by copy */
  constexpr Rational EnergyLevelQuantumNumber(QuantumNumberType X) const ARTS_NOEXCEPT {
    ARTS_ASSERT(mqtype == ENERGY_LEVEL);
    return mqm[ENERGY_LEVEL_INDEX][X];
  };

  /** Return the upper quantum numbers by reference */
  constexpr QuantumNumbers& UpperQuantumNumbers() noexcept {
    return mqm[TRANSITION_UPPER_INDEX];
  };

  /** Return the lower quantum numbers by reference */
  constexpr QuantumNumbers& LowerQuantumNumbers() noexcept {
    return mqm[TRANSITION_LOWER_INDEX];
  };

  /** Return the energy level quantum numbers by reference */
  constexpr QuantumNumbers& EnergyLevelQuantumNumbers() noexcept {
    return mqm[ENERGY_LEVEL_INDEX];
  }

  /** Check if there are any quantum numbers defined */
  constexpr bool any_quantumnumbers() const noexcept {
    Index qni = 0;
    switch (mqtype) {
      case QuantumIdentifier::ENERGY_LEVEL:
      case QuantumIdentifier::TRANSITION:
        for (const auto& qns : mqm) do {
          if (not qns[qni].isUndefined()) return true;
          qni++;
        } while (qni not_eq Index(QuantumNumberType::FINAL));
        break;
      case QuantumIdentifier::ALL:
      case QuantumIdentifier::NONE:
        break;
    }
    return false;
  }

  /** Check if *this is a energy level type of identifier */
  constexpr bool IsEnergyLevelType() const noexcept { return mqtype == ENERGY_LEVEL; }

 private:
  QType mqtype;
  Index mspecies;
  Index miso;
//   Index misotopologue_index;
  std::array<QuantumNumbers, 2> mqm;
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

/** Is everything the same between the identifiers
 * 
 * May throw if different Qtypes are compared.
 * 
 * @param[in] a One identifier
 * @param[in] b Another identifier
 * @return true If all doubly defined quantum numbers match
 * @return false Otherwise
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr bool operator==(const QuantumIdentifier& a, const QuantumIdentifier& b) {
  if (!(a.Isotopologue() == b.Isotopologue() && a.Species() == b.Species() &&
        a.Type() == b.Type()))
    return false;
  
  if (a.Type() == QuantumIdentifier::ENERGY_LEVEL)
    return a.QuantumMatch()[a.ENERGY_LEVEL_INDEX] == b.QuantumMatch()[b.ENERGY_LEVEL_INDEX];
  else if (a.Type() == QuantumIdentifier::TRANSITION)
    return a.QuantumMatch()[a.TRANSITION_LOWER_INDEX] == b.QuantumMatch()[b.TRANSITION_LOWER_INDEX] and
           a.QuantumMatch()[a.TRANSITION_UPPER_INDEX] == b.QuantumMatch()[b.TRANSITION_UPPER_INDEX];
  else if (a.Type() == QuantumIdentifier::ALL)
    return true;
  else if (a.Type() == QuantumIdentifier::NONE)
    return false;
  else {
    ARTS_ASSERT(false, "Programmer error --- added type is missing");
  }
}
#pragma GCC diagnostic pop

/** Is anything different between the identifiers
 * 
 * May throw if different Qtypes are compared.
 * 
 * @param[in] a One identifier
 * @param[in] b Another identifier
 * @return true If some quantum numbers mismatch
 * @return false Otherwise
 */
constexpr bool operator!=(const QuantumIdentifier& a, const QuantumIdentifier& b) {
  return not operator==(a, b);
}

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

/** Input operator */
std::istream& operator>>(std::istream& is, QuantumNumbers& qn);

/** Output operator */
std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn);

/** Output operator */
std::ostream& operator<<(std::ostream& os, const QuantumIdentifier& qi);

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
