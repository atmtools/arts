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

/** \file
    Classes to handle quantum numbers.

    \author Oliver Lemke
*/

#ifndef quantum_h
#define quantum_h

#include <map>
#include <iostream>
#include <stdexcept>
#include "matpack.h"
#include "rational.h"
#include "mystring.h"
#include "array.h"


//! Enum for Quantum Numbers used for indexing
/*
 If you add anything here, remember to also adapt
 operator<<(ostream&, const QuantumNumbers&) and
 operator>>(istream&, QuantumNumbers&)
 to handle the added numbers.
 */
enum class QuantumNumberType : Index {
    J=0,         // Total angular momentum
    dJ,          // Delta total angular momentum
    M,           // Projection of J along magnetic field
    N,           // J minus spin
    dN,          // Delta J minus spin
    S,           // Spin angular momentum (from electrons) NOTE: S_global for HITRAN S
    F,           // J + nuclear spin
    K,           //(This is a projection of J along one axis)
    Ka,          //(This is a projection of J along one axis)
    Kc,          //(This is a projection of J along another axis)
    Omega,       // This is an absolute projection of J and S
    i,           //(Is related to Omega)
    Lambda,      // This is Sigma or Pi or Lambda states (as seen in literature)
    alpha,       // Alpha from HITRAN // FIXME richard
    Sym,         // Symmetry expression
    parity,      // parity value (+/-)
    v1,          // Vibrational mode 1
    v2,          // Vibrational mode 2
    l2,          // Vibrational angular momentum associated with v2
    v3,          // Vibrational mode 3
    v4,          // Vibrational mode 4
    v5,          // Vibrational mode 5
    v6,          // Vibrational mode 6
    l,           // The absolute sum of l_j for v_j
    pm,          // Symmetry type for l=0
    r,           // Rank of the level within a set of the same vibrational symmetry
    S_global,    // Symmetry of the level
    X,           // Electronic state
    n_global,    // Torosional quanta
    C,           // Another symmetry expression
    Hund,        // Flag for Hund case type.  This flag lets Zeeman know what to expect
    FINAL_ENTRY  // We need this to determine the number of elements in this enum
};


//! Enum for details about matched quantum numbers.
typedef enum {
    QMI_NONE = 0,
    QMI_FULL = 1,
    QMI_PARTIAL = 2,
} QuantumMatchInfoEnum;

enum class Hund : Index {
    CaseA=0,
    CaseB=1
};

//! Class that holds details for matching info on upper and lower quantum numbers.
class QuantumMatchInfo
{
public:
    void SetUpper(const QuantumMatchInfoEnum qmie) { mupper = qmie; }
    void SetLower(const QuantumMatchInfoEnum qmie) { mlower = qmie; }

    const QuantumMatchInfoEnum& Upper() const { return mupper; }
    const QuantumMatchInfoEnum& Lower() const { return mlower; }

    QuantumMatchInfoEnum& Upper() { return mupper; }
    QuantumMatchInfoEnum& Lower() { return mlower; }

private:
    QuantumMatchInfoEnum mupper;
    QuantumMatchInfoEnum mlower;
};

typedef Array<QuantumMatchInfo> ArrayOfQuantumMatchInfo;


//! Container class for Quantum Numbers
class QuantumNumbers
{
public:
    typedef Array<Rational> QuantumContainer;

    QuantumNumbers()
    {
        mqnumbers.resize(Index(QuantumNumberType::FINAL_ENTRY));
        for (Index i = 0; i < mqnumbers.nelem(); i++)
            mqnumbers[i] = RATIONAL_UNDEFINED;
    }

    //! Return copy of quantum number
    const Rational operator[](const Index qn) const
    {
        assert(qn < Index(QuantumNumberType::FINAL_ENTRY));
        return mqnumbers[qn];
    }
    
    //! Return copy of quantum number
    const Rational operator[](const QuantumNumberType qn) const
    {
      assert(qn != QuantumNumberType::FINAL_ENTRY);
      return mqnumbers[Index(qn)];
    }

    //! Set quantum number
    void Set(Index qn, Rational r)
    {
        assert(qn < Index(QuantumNumberType::FINAL_ENTRY));
        mqnumbers[qn] = r;
    }
    
    //! Set quantum number
    void Set(QuantumNumberType qn, Rational r)
    {
      assert(qn != QuantumNumberType::FINAL_ENTRY);
      mqnumbers[Index(qn)] = r;
    }

    //! Set quantum number
    void Set(String name, Rational r)
    {
        // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
if (name == #ID) this->Set(QuantumNumberType::ID, r)

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
        else
        {
            std::ostringstream os;
            os << "Unknown quantum number: " << name << " (" << r << ").";
            throw std::runtime_error(os.str());
        }

#undef INPUT_QUANTUM
    }
    
    const QuantumContainer& GetNumbers() const { return mqnumbers; }

    Index nNumbers() const
    {
        Index n = 0;
        for (QuantumContainer::const_iterator it = mqnumbers.begin();
             it != mqnumbers.end(); it++)
            if (!(*it).isUndefined())
                n++;
        return n;
    }

    //! Compare Quantum Numbers
    /**
     Ignores any undefined numbers in the comparison

     \param[in] qn  Quantum Numbers to compare to

     \returns True for match
     */
    bool Compare(const QuantumNumbers& qn) const;

    bool CompareDetailed(QuantumMatchInfoEnum& imatch, const QuantumNumbers& qn) const;
    
    bool operator<(const QuantumNumbers& other) const;
    
    bool operator>(const QuantumNumbers& other) const;

private:
    QuantumContainer mqnumbers;
};


typedef Array<QuantumNumbers> ArrayOfQuantumNumbers;


//! Record containing upper and lower quantum numbers
class QuantumNumberRecord
{
public:
    QuantumNumberRecord() : mqns(2) {}
    QuantumNumberRecord(const ArrayOfQuantumNumbers& qns) : mqns(qns) {assert(mqns.nelem() == 2);}
    
    enum class Level : Index { Upper=0, Lower=1 };
  
    //! Set lower quantum number
    void SetLower(const Index i, const Rational r) {  mqns[Index(Level::Lower)].Set(i, r); }
    void SetLower(const QuantumNumberType i, const Rational r) { mqns[Index(Level::Lower)].Set(i, r); }
    void SetLower(const String i, const Rational r) { mqns[Index(Level::Lower)].Set(i, r); }

    //! Set upper quantum number
    void SetUpper(const Index i, const Rational r) { mqns[Index(Level::Upper)].Set(i, r); }
    void SetUpper(const QuantumNumberType i, const Rational r) { mqns[Index(Level::Upper)].Set(i, r); }
    void SetUpper(const String i, const Rational r) { mqns[Index(Level::Upper)].Set(i, r); }

    //! Get lower quantum number
    Rational Lower(Index i) const { return mqns[Index(Level::Lower)][i]; }
    Rational Lower(QuantumNumberType i) const { return mqns[Index(Level::Lower)][Index(i)]; }

    //! Get upper quantum number
    Rational Upper(Index i) const { return mqns[Index(Level::Upper)][i]; }
    Rational Upper(QuantumNumberType i) const { return mqns[Index(Level::Upper)][Index(i)]; }

    //! Get lower quantum numbers
    QuantumNumbers& Lower() { return mqns[Index(Level::Lower)]; }

    //! Get lower quantum numbers
    const QuantumNumbers& Lower() const { return mqns[Index(Level::Lower)]; }

    //! Get upper quantum numbers
    QuantumNumbers& Upper() { return mqns[Index(Level::Upper)]; }

    //! Get upper quantum numbers
    const QuantumNumbers& Upper() const { return mqns[Index(Level::Upper)]; }
    
    const ArrayOfQuantumNumbers& Data() const {return mqns;}
    ArrayOfQuantumNumbers& Data() {return mqns;}

private:
    ArrayOfQuantumNumbers mqns;
};


//! Class to identify and match lines by their quantum numbers
/*!
 Describes either a transition or an energy level and can be used
 to find matching lines.
 
 For transitions, the QI contains upper and lower quantum numbers.
 For energy levels, it only holds one set of quantum numbers which
 are then matched against the upper and lower qns of the lines.
 
 File format:
 
 Transition:   SPECIES_NAME-ISOTOPE TR UP QUANTUMNUMBERS LO QUANTUMNUMBERS
 Energy level: SPECIES_NAME-ISOTOPE EN QUANTUMNUMBERS
 All lines:    SPECIES_NAME-ISOTOPE ALL
 
 H2O-161 TR UP J 0/1 v1 2/3 LO J 1/1 v2 1/2
 H2O-161 EN J 0/1 v1 2/3
 H2O-161 ALL

*/

class QuantumIdentifier
{
public:
    //! Identify lines by transition or energy level
    typedef enum {
        TRANSITION,
        ENERGY_LEVEL,
        ALL,
        NONE 
    } QType;
    
    QuantumIdentifier() : mqtype(NONE), mspecies(-1), miso(-1), mqm(0) {}
    
    QuantumIdentifier(const QuantumIdentifier::QType qt, const Index species, const Index iso) : mspecies(species), miso(iso) { SetType(qt); }

    QuantumIdentifier(const Index spec, const Index isot, const QuantumNumbers& upper, const QuantumNumbers& lower)
    : mqtype(QuantumIdentifier::TRANSITION), mspecies(spec), miso(isot), mqm(2) { mqm[TRANSITION_LOWER_INDEX] = lower; mqm[TRANSITION_UPPER_INDEX] = upper; }
    
    QuantumIdentifier(const Index spec, const Index isot, const QuantumNumbers& qnr)
    : mspecies(spec), miso(isot)
    {
      SetEnergyLevel(qnr);
    }
    
    QuantumIdentifier(String x) {SetFromString(x);}

    static const Index TRANSITION_UPPER_INDEX = 0;
    static const Index TRANSITION_LOWER_INDEX = 1;
    static const Index ENERGY_LEVEL_INDEX = 0;

    void SetType(const QuantumIdentifier::QType qt)
    {
        mqtype = qt;
        switch (qt) {
            case QuantumIdentifier::TRANSITION:
                mqm.resize(2);
                break;
            case QuantumIdentifier::ENERGY_LEVEL:
                mqm.resize(1);
                break;
            case QuantumIdentifier::ALL:
            case QuantumIdentifier::NONE:
              mqm.resize(0);
              break;
            default:
                break;
        }
    }

    void SetSpecies(const Index &sp) { mspecies = sp; }
    void SetIsotopologue(const Index &iso) { miso = iso; }
    void SetTransition(const QuantumNumbers& upper, const QuantumNumbers& lower);
    void SetEnergyLevel(const QuantumNumbers& q);
    void SetAll();
    void SetFromString(String str);
    void SetFromStringForCO2Band(String upper, String lower, String iso);

    QType Type() const { return mqtype; }
    String TypeStr() const;
    String  SpeciesName() const;
    Index  Species() const { return mspecies; }
    Index& Species() { return mspecies; }
    Index  Isotopologue() const { return miso; }
    Index& Isotopologue() { return miso; }
    const ArrayOfQuantumNumbers& QuantumMatch() const { return mqm; }
    ArrayOfQuantumNumbers& QuantumMatch() { return mqm; }
    
    QuantumIdentifier UpperQuantumId() const {assert(mqtype==TRANSITION); QuantumIdentifier qi(mspecies, miso, mqm[TRANSITION_UPPER_INDEX]); return qi;};
    QuantumIdentifier LowerQuantumId() const {assert(mqtype==TRANSITION); QuantumIdentifier qi(mspecies, miso, mqm[TRANSITION_LOWER_INDEX]); return qi;};
    
    const QuantumNumbers& UpperQuantumNumbers() const {assert(mqtype==TRANSITION); return mqm[TRANSITION_UPPER_INDEX];};
    const QuantumNumbers& LowerQuantumNumbers() const {assert(mqtype==TRANSITION); return mqm[TRANSITION_LOWER_INDEX];};
    QuantumNumbers& UpperQuantumNumbers() {assert(mqtype==TRANSITION); return mqm[TRANSITION_UPPER_INDEX];};
    QuantumNumbers& LowerQuantumNumbers() {assert(mqtype==TRANSITION); return mqm[TRANSITION_LOWER_INDEX];};
    
    //! Tests if RHS contains LHS some how
    bool In(const QuantumIdentifier& other) const;
    bool InLower(const QuantumIdentifier& other) const;
    bool InUpper(const QuantumIdentifier& other) const;
    
    //! Tests if there are any defined quantum numbers
    bool any_quantumnumbers() const;

private:
    QType mqtype;
    Index mspecies;
    Index miso;
    ArrayOfQuantumNumbers mqm;
};

inline bool operator==(const QuantumIdentifier& a,const QuantumIdentifier& b){
    if(!(a.Isotopologue()==b.Isotopologue() && 
        a.Species()==b.Species() &&
        a.Type()==b.Type()))
        return false;
    
    if(a.Type()==QuantumIdentifier::ENERGY_LEVEL)
        return a.QuantumMatch()[a.ENERGY_LEVEL_INDEX].Compare(b.QuantumMatch()[b.ENERGY_LEVEL_INDEX]);
    else if(a.Type()==QuantumIdentifier::TRANSITION)
        return a.QuantumMatch()[a.TRANSITION_LOWER_INDEX].Compare(b.QuantumMatch()[ b.TRANSITION_LOWER_INDEX ]) && 
        a.QuantumMatch()[a.TRANSITION_UPPER_INDEX].Compare(b.QuantumMatch()[ b.TRANSITION_UPPER_INDEX ]) ;
    else if(a.Type()==QuantumIdentifier::ALL)
      return true;
    else if(a.Type()==QuantumIdentifier::NONE)
      return false;
    else
      throw std::runtime_error("Programmer error --- added type is missing");
}

typedef Array<QuantumIdentifier> ArrayOfQuantumIdentifier;


//! Check for valid quantum number name
bool IsValidQuantumNumberName(String name);

//! Throws runtime error if quantum number name is invalid
void ThrowIfQuantumNumberNameInvalid(String name);

std::istream& operator>>(std::istream& is, QuantumNumbers& qn);
std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn);

std::ostream& operator<<(std::ostream& os, const QuantumNumberRecord& qr);

std::istream& operator>>(std::istream& is, QuantumIdentifier& qi);
std::ostream& operator<<(std::ostream& os, const QuantumIdentifier& qi);

#endif

