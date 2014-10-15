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
#include "matpack.h"
#include "rational.h"


//! Enum for Quantum Numbers used for indexing
/*
 If you add anything here, remember to also adapt
 operator<<(ostream&, const QuantumNumbers&) and
 operator>>(istream&, QuantumNumbers&)
 to handle the added numbers.
 */
typedef enum {
    QN_J=0,         // Total
    QN_M,           // Projection of J
    QN_N,           // Total-Spin
    QN_S,           // Electronic spin
    QN_F,           // Total + nuclear spin
    QN_Omega,       // Absolute of projection of total + projection of spin
    QN_K1,          // Either K or Ka in HITRAN (This is a projection of J along one axis)
    QN_K2,          // Either Kb in HITRAN (This is a projection of J along another axis)
    QN_v1,          // Vibrational mode 1
    QN_v2,          // Vibrational mode 2
    QN_v3,          // Vibrational mode 3
    QN_FINAL_ENTRY  // We need this to determine the number of elements in this enum
} QuantumIds;


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

    const QuantumContainer& GetNumbers() const { return mqnumbers; }

    //! Compare Quantum Numbers
    /**
     Ignores any undefined numbers in the comparison

     \param[in] qn  Quantum Numbers to compare to

     \returns True for match
     */
    bool Compare(const QuantumNumbers& qn) const;

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


std::istream& operator>>(std::istream& is, QuantumNumbers& qn);
std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn);

std::ostream& operator<<(std::ostream& os, const QuantumNumberRecord& qr);

#endif

