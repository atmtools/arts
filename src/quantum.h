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
#include "interpolation.h"
#include "matpack.h"
#include "mystring.h"
#include "rational.h"

/** Enum for Quantum Numbers used for indexing
 *
 * If you add anything here, remember to also adapt
 * operator<<(ostream&, const QuantumNumbers&) and
 * operator>>(istream&, QuantumNumbers&)
 * to handle the added numbers.
 */
enum class QuantumNumberType : Index {
  J = 0,   // Total angular momentum
  dJ,      // Delta total angular momentum
  M,       // Projection of J along magnetic field
  N,       // J minus spin
  dN,      // Delta J minus spin
  S,       // Spin angular momentum (from electrons) NOTE: S_global for HITRAN S
  F,       // J + nuclear spin
  K,       //(This is a projection of J along one axis)
  Ka,      //(This is a projection of J along one axis)
  Kc,      //(This is a projection of J along another axis)
  Omega,   // This is an absolute projection of J and S
  i,       //(Is related to Omega)
  Lambda,  // This is Sigma or Pi or Lambda states (as seen in literature)
  alpha,   // Alpha from HITRAN
  Sym,     // Symmetry expression
  parity,  // parity value (+/-)
  v1,      // Vibrational mode 1
  v2,      // Vibrational mode 2
  l2,      // Vibrational angular momentum associated with v2
  v3,      // Vibrational mode 3
  v4,      // Vibrational mode 4
  v5,      // Vibrational mode 5
  v6,      // Vibrational mode 6
  l,       // The absolute sum of l_j for v_j
  pm,      // Symmetry type for l=0
  r,       // Rank of the level within a set of the same vibrational symmetry
  S_global,  // Symmetry of the level
  X,         // Electronic state
  n_global,  // Torosional quanta
  C,         // Another symmetry expression
  Hund,  // Flag for Hund case type.  This flag lets Zeeman know what to expect
  FINAL_ENTRY  // We need this to determine the number of elements in this enum
};

inline QuantumNumberType string2quantumnumbertype(const String& s) {
  #define INPUT_QUANTUM(ID) \
  if (s == #ID) return QuantumNumberType::ID
  INPUT_QUANTUM(J);
  else INPUT_QUANTUM(dJ);
  else INPUT_QUANTUM(M);
  else INPUT_QUANTUM(N);
  else INPUT_QUANTUM(dN);
  else INPUT_QUANTUM(S);
  else INPUT_QUANTUM(F);
  else INPUT_QUANTUM(K);
  else INPUT_QUANTUM(Ka);
  else INPUT_QUANTUM(Kc);
  else INPUT_QUANTUM(Omega);
  else INPUT_QUANTUM(i);
  else INPUT_QUANTUM(Lambda);
  else INPUT_QUANTUM(alpha);
  else INPUT_QUANTUM(Sym);
  else INPUT_QUANTUM(parity);
  else INPUT_QUANTUM(v1);
  else INPUT_QUANTUM(v2);
  else INPUT_QUANTUM(l2);
  else INPUT_QUANTUM(v3);
  else INPUT_QUANTUM(v4);
  else INPUT_QUANTUM(v5);
  else INPUT_QUANTUM(v6);
  else INPUT_QUANTUM(l);
  else INPUT_QUANTUM(pm);
  else INPUT_QUANTUM(r);
  else INPUT_QUANTUM(S_global);
  else INPUT_QUANTUM(X);
  else INPUT_QUANTUM(n_global);
  else INPUT_QUANTUM(C);
  else INPUT_QUANTUM(Hund);
  else return QuantumNumberType::FINAL_ENTRY;
  #undef INPUT_QUANTUM
}

inline String quantumnumbertype2string(QuantumNumberType s) {
  #define INPUT_QUANTUM(ID) \
  if (s == QuantumNumberType::ID) return #ID
  INPUT_QUANTUM(J);
  else INPUT_QUANTUM(dJ);
  else INPUT_QUANTUM(M);
  else INPUT_QUANTUM(N);
  else INPUT_QUANTUM(dN);
  else INPUT_QUANTUM(S);
  else INPUT_QUANTUM(F);
  else INPUT_QUANTUM(K);
  else INPUT_QUANTUM(Ka);
  else INPUT_QUANTUM(Kc);
  else INPUT_QUANTUM(Omega);
  else INPUT_QUANTUM(i);
  else INPUT_QUANTUM(Lambda);
  else INPUT_QUANTUM(alpha);
  else INPUT_QUANTUM(Sym);
  else INPUT_QUANTUM(parity);
  else INPUT_QUANTUM(v1);
  else INPUT_QUANTUM(v2);
  else INPUT_QUANTUM(l2);
  else INPUT_QUANTUM(v3);
  else INPUT_QUANTUM(v4);
  else INPUT_QUANTUM(v5);
  else INPUT_QUANTUM(v6);
  else INPUT_QUANTUM(l);
  else INPUT_QUANTUM(pm);
  else INPUT_QUANTUM(r);
  else INPUT_QUANTUM(S_global);
  else INPUT_QUANTUM(X);
  else INPUT_QUANTUM(n_global);
  else INPUT_QUANTUM(C);
  else INPUT_QUANTUM(Hund);
  throw std::runtime_error("Bad quantum number type");
  #undef INPUT_QUANTUM
}

/** Enum for Hund cases */
enum class Hund : Index { CaseA = 0, CaseB = 1 };

/** Container class for Quantum Numbers */
class QuantumNumbers {
 public:
  typedef std::array<Rational, Index(QuantumNumberType::FINAL_ENTRY)>
      QuantumContainer;

  /** Initializer with undefined values */
  constexpr QuantumNumbers() noexcept
      : mqnumbers({RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED, RATIONAL_UNDEFINED, RATIONAL_UNDEFINED,
                   RATIONAL_UNDEFINED}) {}

  /** Access operator
   * 
   * @param[in] qn Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  constexpr Rational operator[](const Index qn) const noexcept {
    return mqnumbers[qn];
  }

  /** Access operator
   * 
   * @param qn[in] Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  constexpr Rational operator[](const QuantumNumberType qn) const noexcept {
    return mqnumbers[Index(qn)];
  }
  
  /** Access operator
   * 
   * @param[in] qn Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  Rational& operator[](const Index qn) noexcept {
    return mqnumbers[qn];
  }
  
  /** Access operator
   * 
   * @param qn[in] Index Pos to access
   * @return constexpr Rational Copy of value at pos
   */
  Rational& operator[](const QuantumNumberType qn) noexcept {
    return mqnumbers[Index(qn)];
  }

  /** Set quantum number at position
   * 
   * @param[in] qn Index Pos to set at
   * @param[in] r Rational to set
   */
  void Set(Index qn, Rational r) {
    assert(qn < Index(QuantumNumberType::FINAL_ENTRY));
    mqnumbers[qn] = r;
  }

  /** Set quantum number at position
   * 
   * @param[in] qn Index Pos to set at
   * @param[in] r Rational to set
   */
  void Set(QuantumNumberType qn, Rational r) {
    assert(qn != QuantumNumberType::FINAL_ENTRY);
    mqnumbers[Index(qn)] = r;
  }

  /** Set quantum number at position
   * 
   * @param[in] qn String Pos to set at by name
   * @param[in] r Rational to set
   */
  void Set(String qn, Rational r) {
    mqnumbers[Index(string2quantumnumbertype(qn))] = r;
  }

  /** Get the numbers
   * 
   * @return const QuantumContainer& All the numbers
   */
  const QuantumContainer& GetNumbers() const { return mqnumbers; }

  /** The number of defined quantum numbers
   * 
   * @return Index Count of defined numbers
   */
  Index nNumbers() const {
    return std::accumulate(
        mqnumbers.cbegin(), mqnumbers.cend(), 0, [](Index i, Rational r) {
          return r.isUndefined() ? i : i + 1;
        });
  }

  /** Compare Quantum Numbers
   * Ignores any undefined numbers in the comparison
   *
   * @param[in] qn  Quantum Numbers to compare to
   *
   * @return true For a match
   * @return false Otherwise
   */
  bool Compare(const QuantumNumbers& qn) const;

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
  typedef enum { TRANSITION, ENERGY_LEVEL, ALL, NONE } QType;

  /** Initialize with no matches */
  constexpr QuantumIdentifier() noexcept
      : mqtype(QuantumIdentifier::NONE), mspecies(-1), miso(-1) {}

  /** Initialize with no quantum numbers defined but known species and matching type
   * 
   * @param[in] qt Way to identify quantum numbers
   * @param[in] species Species index-mapped
   * @param[in] iso Isotopologue index-mapped
   */
  constexpr QuantumIdentifier(const QuantumIdentifier::QType qt,
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
                    const std::vector<Rational>& lower)
      : mqtype(QuantumIdentifier::TRANSITION),
        mspecies(spec),
        miso(isot) {
          for(size_t i=0; i<keys.size(); i++) {
            mqm[TRANSITION_UPPER_INDEX][keys[i]] = upper[i];
            mqm[TRANSITION_LOWER_INDEX][keys[i]] = lower[i];
          }
        }

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
  QuantumIdentifier(String x) { SetFromString(x); }

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
  void SetSpecies(const Index& sp) { mspecies = sp; }

  /** Set the Isotopologue
   * 
   * @param[in] is Isotopologue index-mapped
   */
  void SetIsotopologue(const Index& iso) { miso = iso; }

  /** Set to transition type identifier
   * 
   * @param[in] upper Upper state quantum numbers
   * @param[in] lower Lower state quantum numbers
   */
  void SetTransition(const QuantumNumbers& upper, const QuantumNumbers& lower);

  /** Set tp energy level identifier
   * 
   * @param[in] q Quantum numbers
   */
  void SetEnergyLevel(const QuantumNumbers& q);

  /** Set to All identifier */
  void SetAll();

  /** Set key to transition type */
  void SetTransition() { mqtype = QuantumIdentifier::TRANSITION; };

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

  /** @return QType as String */
  String TypeStr() const;

  /** Return the Species by name */
  String SpeciesName() const;

  /** Return the Species by index */
  constexpr Index Species() const { return mspecies; }

  /** Return the Species by index reference */
  Index& Species() { return mspecies; }

  /** Return the Isotopologue by index */
  constexpr Index Isotopologue() const { return miso; }

  /** Return the Isotopologue by index reference */
  Index& Isotopologue() { return miso; }

  /** Return the quantum numbers array const reference */
  const std::array<QuantumNumbers, 2>& QuantumMatch() const { return mqm; }

  /** Return the quantum numbers array reference */
  std::array<QuantumNumbers, 2>& QuantumMatch() { return mqm; }

  /** Return a quantum identifer as if it wants to match to upper energy level */
  constexpr QuantumIdentifier UpperQuantumId() const noexcept {
    return QuantumIdentifier(mspecies, miso, mqm[TRANSITION_UPPER_INDEX]);
  };

  /** Return a quantum identifer as if it wants to match to lower energy level */
  constexpr QuantumIdentifier LowerQuantumId() const noexcept {
    return QuantumIdentifier(mspecies, miso, mqm[TRANSITION_LOWER_INDEX]);
  };

  /** Return the upper quantum numbers by const reference */
  const QuantumNumbers& UpperQuantumNumbers() const noexcept {
    return mqm[TRANSITION_UPPER_INDEX];
  };

  /** Return the lower quantum numbers by const reference */
  const QuantumNumbers& LowerQuantumNumbers() const noexcept {
    return mqm[TRANSITION_LOWER_INDEX];
  };

  /** Return a upper quantum number by copy */
  constexpr Rational UpperQuantumNumber(QuantumNumberType X) const noexcept {
    return mqm[TRANSITION_UPPER_INDEX][X];
  };

  /** Return a lower quantum number by copy */
  constexpr Rational LowerQuantumNumber(QuantumNumberType X) const noexcept {
    return mqm[TRANSITION_LOWER_INDEX][X];
  };
  
  /** Return a upper quantum number by copy */
  Rational& UpperQuantumNumber(QuantumNumberType X) noexcept {
    return mqm[TRANSITION_UPPER_INDEX][X];
  };
  
  /** Return a lower quantum number by copy */
  Rational& LowerQuantumNumber(QuantumNumberType X) noexcept {
    return mqm[TRANSITION_LOWER_INDEX][X];
  };

  /** Return the energy level quantum numbers by const reference */
  const QuantumNumbers& EnergyLevelQuantumNumbers() const noexcept {
    return mqm[ENERGY_LEVEL_INDEX];
  }
  
  /** Return a energy level quantum number by copy */
  constexpr Rational EnergyLevelQuantumNumber(QuantumNumberType X) const noexcept {
    return mqm[ENERGY_LEVEL_INDEX][X];
  };

  /** Return the upper quantum numbers by reference */
  QuantumNumbers& UpperQuantumNumbers() {
    return mqm[TRANSITION_UPPER_INDEX];
  };

  /** Return the lower quantum numbers by reference */
  QuantumNumbers& LowerQuantumNumbers() {
    return mqm[TRANSITION_LOWER_INDEX];
  };

  /** Return the energy level quantum numbers by reference */
  QuantumNumbers& EnergyLevelQuantumNumbers() {
    return mqm[ENERGY_LEVEL_INDEX];
  }

  /** Return if this is in other
   * 
   * All quantum numbers defined in *this must be the
   * same as the quantum numbers in other for a call
   * to In() to return true.  The numbers in other
   * must not all be defined in *this.
   * 
   * @param[in] other Another quantum identifier
   * @return true If the above description holds
   * @return false Otherwise
   */
  bool In(const QuantumIdentifier& other) const;

  /** Return if this is in other's lower energy state
   * 
   * All quantum numbers defined in *this must be the
   * same as the quantum numbers in other for a call
   * to In() to return true.  The numbers in other
   * must not all be defined in *this.
   * 
   * @param[in] other Another quantum identifier
   * @return true If the above description holds
   * @return false Otherwise
   */
  bool InLower(const QuantumIdentifier& other) const;

  /** Return if this is in other's upper energy state
   * 
   * All quantum numbers defined in *this must be the
   * same as the quantum numbers in other for a call
   * to In() to return true.  The numbers in other
   * must not all be defined in *this.
   * 
   * @param[in] other Another quantum identifier
   * @return true If the above description holds
   * @return false Otherwise
   */
  bool InUpper(const QuantumIdentifier& other) const;

  /** Check if there are any quantum numbers defined */
  bool any_quantumnumbers() const;

  /** Check if *this is a energy level type of identifier */
  bool IsEnergyLevelType() const { return mqtype == ENERGY_LEVEL; }

 private:
  QType mqtype;
  Index mspecies;
  Index miso;
  std::array<QuantumNumbers, 2> mqm;
};

/** Is everything the same between the identifiers
 * 
 * May throw if different Qtypes are compared.
 * 
 * @param[in] a One identifier
 * @param[in] b Another identifier
 * @return true If all quantum numbers match
 * @return false Otherwise
 */
inline bool operator==(const QuantumIdentifier& a, const QuantumIdentifier& b) {
  if (!(a.Isotopologue() == b.Isotopologue() && a.Species() == b.Species() &&
        a.Type() == b.Type()))
    return false;

  if (a.Type() == QuantumIdentifier::ENERGY_LEVEL)
    return a.QuantumMatch()[a.ENERGY_LEVEL_INDEX].Compare(
        b.QuantumMatch()[b.ENERGY_LEVEL_INDEX]);
  else if (a.Type() == QuantumIdentifier::TRANSITION)
    return a.QuantumMatch()[a.TRANSITION_LOWER_INDEX].Compare(
               b.QuantumMatch()[b.TRANSITION_LOWER_INDEX]) &&
           a.QuantumMatch()[a.TRANSITION_UPPER_INDEX].Compare(
               b.QuantumMatch()[b.TRANSITION_UPPER_INDEX]);
  else if (a.Type() == QuantumIdentifier::ALL)
    return true;
  else if (a.Type() == QuantumIdentifier::NONE)
    return false;
  else
    throw std::runtime_error("Programmer error --- added type is missing");
}

/** Is anything different between the identifiers
 * 
 * May throw if different Qtypes are compared.
 * 
 * @param[in] a One identifier
 * @param[in] b Another identifier
 * @return true If some quantum numbers mismatch
 * @return false Otherwise
 */
inline bool operator!=(const QuantumIdentifier& a, const QuantumIdentifier& b) {
  return not operator==(a, b);
}

/** Check if all quantum numbers are the same between a and b
 * 
 * @param[in] a One set of quantum numbers
 * @param[in] b Another set of quantum numbers
 * @return true If all quantum numbers match
 * @return false Otherwise
 */
inline bool operator==(const QuantumNumbers& a, const QuantumNumbers& b) {
  return a.Compare(b) and b.Compare(a);
}

typedef Array<QuantumIdentifier> ArrayOfQuantumIdentifier;

/** Check for valid quantum number name
 * 
 * @param[in] name Parameter
 * @return true If the parameter exist
 * @return false Otherwise
 */
bool IsValidQuantumNumberName(String name);

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

// Output from EnergyLevelMap
struct Output2{
  Numeric r_low;
  Numeric r_upp;
};

// Output from EnergyLevelMap
struct Output4{
  Numeric E_low;
  Numeric E_upp;
  Numeric T_low;
  Numeric T_upp;
};

enum class EnergyLevelMapType {
  Tensor3_t,
  Vector_t,
  Numeric_t,
  None_t,
};

EnergyLevelMapType string2energylevelmaptype(const String& s);

String energylevelmaptype2string(EnergyLevelMapType type);

class EnergyLevelMap {
private:
  EnergyLevelMapType mtype;
  ArrayOfQuantumIdentifier mlevels;
  Vector mvib_energy;
  Tensor4 mvalue;
  
public:
  void ThrowIfNotOK() const {
    if (not (mvalue.nbooks() == mlevels.nelem() and 
      (mvib_energy.nelem() == mlevels.nelem() or mvib_energy.nelem() == 0)))
      throw std::runtime_error("Bad dimensions");
    if (mtype == EnergyLevelMapType::Tensor3_t) {}
    else if (mtype == EnergyLevelMapType::Vector_t) {
      if (mvalue.npages() not_eq 1 or mvalue.nrows() not_eq 1)
        throw std::runtime_error("Bad dimensions for vector type");
    }
    else if (mtype == EnergyLevelMapType::Numeric_t) {
      if (mvalue.npages() not_eq 1 or mvalue.nrows() not_eq 1 or mvalue.ncols() not_eq 1)
        throw std::runtime_error("Bad dimensions for numeric type");
    }
    else if (mtype == EnergyLevelMapType::None_t) {
      if (mvalue.npages() not_eq 0 or mvalue.nrows() not_eq 0 or mvalue.ncols() not_eq 0)
        throw std::runtime_error("Bad dimensions for none type");
    }
    
    for (auto& e: mvib_energy) if (e < 0) throw std::runtime_error("Bad energies");
  }
  
  EnergyLevelMap() : mtype(EnergyLevelMapType::None_t), mlevels(0),
  mvib_energy(0), mvalue(0, 0, 0, 0) {ThrowIfNotOK();}
  
  EnergyLevelMap(EnergyLevelMapType new_type, Index pages, Index rows,
                 Index cols, const EnergyLevelMap& old) : 
  mtype(new_type), mlevels(old.mlevels), mvib_energy(old.mvib_energy),
  mvalue(old.mlevels.nelem(), pages, rows, cols) {ThrowIfNotOK();};
  
  // Create Tensor3_t from the raw inputs
  EnergyLevelMap(const Tensor4& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies=Vector(0));
  
  // Create Vector_t from Tensor3_t
  EnergyLevelMap InterpToGridPos(Index atmosphere_dim, const ArrayOfGridPos& p, const ArrayOfGridPos& lat, const ArrayOfGridPos& lon) const;
  
  // Create Numeric_t from Vector_t
  EnergyLevelMap operator[](Index ip) const;
  
  // Create Numeric_t from Tensor3_t
  EnergyLevelMap operator()(Index ip, Index ilat, Index ilon) const;
  
  /** Energy level type */
  EnergyLevelMapType Type() const noexcept {return mtype;}
  
  /** Energy level type */
  const ArrayOfQuantumIdentifier& Levels() const noexcept {return mlevels;}
  
  /** Energy level type */
  const Vector& Energies() const noexcept {return mvib_energy;}
  
  /** Energy level type */
  const Tensor4& Data() const noexcept {return mvalue;}
  
  /** Energy level type */
  EnergyLevelMapType& Type() noexcept {return mtype;}
  
  /** Energy level type */
  ArrayOfQuantumIdentifier& Levels() noexcept {return mlevels;}
  
  /** Energy level type */
  Vector& Energies() noexcept {return mvib_energy;}
  
  /** Energy level type */
  Tensor4& Data() noexcept {return mvalue;}
  
  //////////////////////
  // Numeric_t access //
  //////////////////////
  
  /** Get the output required for Population::NLTE
   * 
   * @param[in] transition A line-by-line transition
   * @return Upper and lower level distributions
   */
  Output2 get_ratio_params(const QuantumIdentifier& transition) const;
  
  /** Get the output required for Population::NLTE-VibrationalTemperatures
   * 
   * @param[in] transition A line-by-line transition
   * @return Upper and lower level distributions and energies
   */
  Output4 get_vibtemp_params(const QuantumIdentifier& transition, const Numeric T) const;
};

std::ostream& operator<<(std::ostream& os, const EnergyLevelMap& elm);

#endif
