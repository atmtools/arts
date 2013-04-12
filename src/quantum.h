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

#include "arts.h"
#include "rational.h"

#define QUANTUM_UNDEFINED Rational(0, 0)

//! Enum for Quantum Numbers used for indexing
typedef enum {
    QN_J,           // Total
    QN_N,           // Total-Spin
    QN_S,           // Electronic spin
    QN_v1,          // Vibrational mode 1
    QN_v2,          // Vibrational mode 2
    QN_v3,          // Vibrational mode 3
    QN_FINAL_ENTRY  // We need this to determine the number of elements in this enum
} QuantumNames;


//! Container class for Quantum Numbers
class QuantumNumbers
{
public:
    QuantumNumbers()
    {
    }

    QuantumNumbers(Index nquantum) : mqnumbers(nquantum, QUANTUM_UNDEFINED)
    {
    }

    //! Return copy of quantum number
    Rational operator[](const Index qn) const { return mqnumbers[qn]; }

    //! Return reference to quantum number
    Rational& operator[](const Index qn) { return mqnumbers[qn]; }

private:
    ArrayOfRational mqnumbers;
};


//! Record containing upper and lower quantum numbers
class QuantumNumberRecord
{
public:
    QuantumNumberRecord() : mqn_upper(QN_FINAL_ENTRY),
                            mqn_lower(QN_FINAL_ENTRY)
    {
    }

    //! Set lower quantum number
    void SetLower(const Index i, const Rational r) { mqn_lower[i] = r; }

    //! Set upper quantum number
    void SetUpper(const Index i, const Rational r) { mqn_upper[i] = r; }

    //! Get lower quantum numbers
    QuantumNumbers& Lower() { return mqn_lower; }

    //! Get lower quantum numbers
    const QuantumNumbers& Lower() const { return mqn_lower; }

    //! Get upper quantum numbers
    QuantumNumbers& Upper() { return mqn_upper; }

    //! Get upper quantum numbers
    const QuantumNumbers& Upper() const { return mqn_upper; }

private:
    QuantumNumbers mqn_upper;
    QuantumNumbers mqn_lower;
};

#endif

