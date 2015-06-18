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
typedef enum {
    QN_J=0,         // Total angular momentum
    QN_dJ,          // Delta total angular momentum
    QN_M,           // Projection of J along magnetic field
    QN_N,           // J minus spin
    QN_dN,          // Delta J minus spin
    QN_S,           // Spin angular momentum (from electrons) NOTE: S_global for HITRAN S
    QN_F,           // J + nuclear spin
    QN_K,           //(This is a projection of J along one axis)
    QN_Ka,          //(This is a projection of J along one axis)
    QN_Kc,          //(This is a projection of J along another axis)
    QN_Omega,       // This is an absolute projection of J and S
    QN_i,           //(Is related to Omega)
    QN_alpha,       // Alpha from HITRAN // FIXME richard
    QN_Sym,         // Symmetry expression
    QN_v1,          // Vibrational mode 1
    QN_v2,          // Vibrational mode 2
    QN_l2,          // Vibrational angular momentum associated with v2
    QN_v3,          // Vibrational mode 3
    QN_v4,          // Vibrational mode 4
    QN_v5,          // Vibrational mode 5
    QN_v6,          // Vibrational mode 6
    QN_l,           // The absolute sum of l_j for v_j
    QN_pm,          // Symmetry type for l=0
    QN_r,           // Rank of the level within a set of the same vibrational symmetry
    QN_S_global,    // Symmetry of the level
    QN_X,           // Electronic state
    QN_n_global,    // Torosional quanta
    QN_C,           // Another symmetry expression
    QN_FINAL_ENTRY  // We need this to determine the number of elements in this enum
} QuantumIds;


//! Enum for details about matched quantum numbers.
typedef enum {
    QMI_NONE = 0,
    QMI_FULL = 1,
    QMI_PARTIAL = 2,
} QuantumMatchInfoEnum;


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
        mqnumbers.resize(QN_FINAL_ENTRY);
        for (Index i = 0; i < mqnumbers.nelem(); i++)
            mqnumbers[i] = RATIONAL_UNDEFINED;
    }

    //! Return copy of quantum number
    const Rational operator[](const Index qn) const
    {
        assert(qn < QN_FINAL_ENTRY);
        return mqnumbers[qn];
    }

    //! Set quantum number
    void Set(Index qn, Rational r)
    {
        assert(qn < QN_FINAL_ENTRY);
        mqnumbers[qn] = r;
    }

    //! Set quantum number
    void Set(String name, Rational r)
    {
        // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
if (name == #ID) this->Set(QN_ ## ID, r)

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
        else INPUT_QUANTUM(alpha);
        else INPUT_QUANTUM(Sym);
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

private:
    QuantumContainer mqnumbers;
};


//! Record containing upper and lower quantum numbers
class QuantumNumberRecord
{
public:
    //! Set lower quantum number
    void SetLower(const Index i, const Rational r) { mqn_lower.Set(i, r); }

    //! Set upper quantum number
    void SetUpper(const Index i, const Rational r) { mqn_upper.Set(i, r); }

    //! Get lower quantum number
    Rational Lower(Index i) const { return mqn_lower[i]; }

    //! Get upper quantum number
    Rational Upper(Index i) const { return mqn_upper[i]; }

    //! Get lower quantum numbers
    QuantumNumbers& Lower() { return mqn_lower; }

    //! Get lower quantum numbers
    const QuantumNumbers& Lower() const { return mqn_lower; }

    //! Get upper quantum numbers
    QuantumNumbers& Upper() { return mqn_upper; }

    //! Get upper quantum numbers
    const QuantumNumbers& Upper() const { return mqn_upper; }

private:
    //! Upper state quantum numbers
    QuantumNumbers mqn_upper;
    //! Lower state quantum numbers
    QuantumNumbers mqn_lower;
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
 
 H2O-161 TR UP J 0/1 v1 2/3 LO J 1/1 v2 1/2
 H2O-161 EN J 0/1 v1 2/3

*/
class QuantumIdentifier
{
public:
    //! Identify lines by transition or energy level
    typedef enum {
        TRANSITION,
        ENERGY_LEVEL
    } QType;

    typedef Array<QuantumNumbers> QuantumMatchCriteria;

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
            default:
                break;
        }
    }

    void SetSpecies(const Index &sp) { mspecies = sp; }
    void SetIsotopologue(const Index &iso) { miso = iso; }
    void SetTransition(const QuantumNumbers q1, const QuantumNumbers q2);
    void SetEnergyLevel(const QuantumNumbers q);
    void SetFromString(String str);

    QType Type() const { return mqtype; }
    String TypeStr() const;
    Index Species() const { return mspecies; }
    Index Isotopologue() const { return miso; }
    const QuantumMatchCriteria& QuantumMatch() const { return mqm; }
    QuantumMatchCriteria& QuantumMatch() { return mqm; }

private:
    QType mqtype;
    Index mspecies;
    Index miso;
    QuantumMatchCriteria mqm;
};


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

