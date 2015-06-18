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

#include <stdexcept>
#include "quantum.h"
#include "global_data.h"
#include "absorption.h"


bool QuantumNumbers::Compare(const QuantumNumbers& qn) const
{
    const QuantumContainer& qnumbers2 = qn.GetNumbers();

    bool match = true;

    Index qnri = 0;

    // Compare all quantum numbers in mqnumbers and qnumbers2
    while (match && qnri != QN_FINAL_ENTRY)
    {
        // If one of the two numbers is undefined, it is considered as
        // a match.
        if (!mqnumbers[qnri].isUndefined()
            && !qnumbers2[qnri].isUndefined()
            && mqnumbers[qnri] != qnumbers2[qnri])
            match = false;

        qnri++;
    }

    return match;
}


bool QuantumNumbers::CompareDetailed(QuantumMatchInfoEnum& imatch, const QuantumNumbers& qn) const
{
    const QuantumContainer& qnumbers2 = qn.GetNumbers();

    bool match = true;

    Index qnri = 0;

    imatch = QMI_FULL;

    // Compare all quantum numbers in mqnumbers and qnumbers2
    while (match && qnri != QN_FINAL_ENTRY)
    {
        // If one of the two numbers is undefined, it is considered as
        // a match.
        if (   (!mqnumbers[qnri].isUndefined() && qnumbers2[qnri].isUndefined())
            || (mqnumbers[qnri].isUndefined() && !qnumbers2[qnri].isUndefined()))
        {
            imatch = QMI_PARTIAL;
        }
        else  if (!mqnumbers[qnri].isUndefined()
                  && !qnumbers2[qnri].isUndefined()
                  && mqnumbers[qnri] != qnumbers2[qnri])
        {
            match = false;
            imatch = QMI_PARTIAL;
        }

        qnri++;
    }

    if (!match) imatch = QMI_NONE;

    return match;
}


bool IsValidQuantumNumberName(String name)
{
    bool valid = false;
    // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
if (name == #ID) valid = true

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
#undef INPUT_QUANTUM
    return valid;
}

void ThrowIfQuantumNumberNameInvalid(String name)
{
    if (!IsValidQuantumNumberName(name))
    {
        ostringstream os;
        os << "Invalid quantum number: " << name;
        throw std::runtime_error(os.str());
    }
}

std::istream& operator>>(std::istream& is, QuantumNumbers& qn)
{
    String name;
    Rational r;

    is >> name >> r;

    // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
    if (name == #ID) qn.Set(QN_ ## ID, r)

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

    return is;
}


std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn)
{
    bool first = true;
    // Define a helper macro to save some typing.
#define OUTPUT_QUANTUM(ID) \
    if (!qn[QN_ ## ID].isUndefined()) \
      { if (!first) os << " "; first = false; os << #ID << " " << qn[QN_ ## ID]; }

        OUTPUT_QUANTUM(J);
        OUTPUT_QUANTUM(dJ);
        OUTPUT_QUANTUM(M);
        OUTPUT_QUANTUM(N);
        OUTPUT_QUANTUM(dN);
        OUTPUT_QUANTUM(S);
        OUTPUT_QUANTUM(F);
        OUTPUT_QUANTUM(K);
        OUTPUT_QUANTUM(Ka);
        OUTPUT_QUANTUM(Kc);
        OUTPUT_QUANTUM(Omega);
        OUTPUT_QUANTUM(i);
        OUTPUT_QUANTUM(alpha);
        OUTPUT_QUANTUM(Sym);
        OUTPUT_QUANTUM(v1);
        OUTPUT_QUANTUM(v2);
        OUTPUT_QUANTUM(l2);
        OUTPUT_QUANTUM(v3);
        OUTPUT_QUANTUM(v4);
        OUTPUT_QUANTUM(v5);
        OUTPUT_QUANTUM(v6);
        OUTPUT_QUANTUM(l);
        OUTPUT_QUANTUM(pm);
        OUTPUT_QUANTUM(r);
        OUTPUT_QUANTUM(S_global);
        OUTPUT_QUANTUM(X);
        OUTPUT_QUANTUM(n_global);
        OUTPUT_QUANTUM(C);

#undef OUTPUT_QUANTUM

    return os;
}


String QuantumIdentifier::TypeStr() const {
    String t;
    switch (mqtype) {
        case QuantumIdentifier::TRANSITION:
            t = "TR";
            break;
        case QuantumIdentifier::ENERGY_LEVEL:
            t = "EN";
            break;
        default:
            assert(0);
            break;
    }
    return t;
}


void QuantumIdentifier::SetTransition(const QuantumNumbers q1, const QuantumNumbers q2)
{
    mqtype = QuantumIdentifier::TRANSITION;
    mqm.resize(2);
    mqm[TRANSITION_UPPER_INDEX] = q1;
    mqm[TRANSITION_LOWER_INDEX] = q2;
}


void QuantumIdentifier::SetEnergyLevel(const QuantumNumbers q)
{
    mqtype = QuantumIdentifier::ENERGY_LEVEL;
    mqm.resize(1);
    mqm[ENERGY_LEVEL_INDEX] = q;
}


void QuantumIdentifier::SetFromString(String str)
{
    // Global species lookup data:
    using global_data::species_data;

    // We need a species index sorted by Arts identifier. Keep this in a
    // static variable, so that we have to do this only once.  The ARTS
    // species index is ArtsMap[<Arts String>].
    static map<String, SpecIsoMap> ArtsMap;

    // Remember if this stuff has already been initialized:
    static bool hinit = false;

    if ( !hinit )
    {
        for ( Index i=0; i<species_data.nelem(); ++i )
        {
            const SpeciesRecord& sr = species_data[i];

            for ( Index j=0; j<sr.Isotopologue().nelem(); ++j)
            {
                SpecIsoMap indicies(i,j);
                String buf = sr.Name()+"-"+sr.Isotopologue()[j].Name();

                ArtsMap[buf] = indicies;
            }
        }
        hinit = true;
    }

    std::istringstream is(str);
    String token;

    is >> token;

    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(token);
    if ( i == ArtsMap.end() )
    {
        ostringstream os;
        os << "ARTS Tag: " << token << " is unknown.";
        throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;
    SetSpecies(id.Speciesindex());
    SetIsotopologue(id.Isotopologueindex());

    is >> token;
    if (token == "TR")
    {
        SetType(QuantumIdentifier::TRANSITION);
        is >> token;
        if (token != "UP")
        {
            std::ostringstream os;
            os << "Expected 'UP', but got: " << token;
            throw std::runtime_error(os.str());
        }

        is >> token;
        Rational r;
        while (is)
        {
            ThrowIfQuantumNumberNameInvalid(token);
            is >> r;
            QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX].Set(token, r);
            is >> token;
            if (token == "LO") break;
        }

        if (!is)
        {
            std::ostringstream os;
            os << "Premature end of data, expected 'LO'.";
            throw std::runtime_error(os.str());
        }
        is >> token;
        while (is)
        {
            ThrowIfQuantumNumberNameInvalid(token);
            is >> r;
            QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX].Set(token, r);
            is >> token;
        }
    }
    else if (token == "EN")
    {
        SetType(QuantumIdentifier::ENERGY_LEVEL);

        is >> token;
        Rational r;
        while (is)
        {
            ThrowIfQuantumNumberNameInvalid(token);
            is >> r;
            QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX].Set(token, r);
            is >> token;
        }
    }
    else
    {
        std::ostringstream os;
        os << "Error parsing QuantumIdentifier. Expected TR or EN, but got: " << token << "\n"
        << "QI: " << str;
        throw std::runtime_error(os.str());
    }
}


std::ostream& operator<<(std::ostream& os, const QuantumNumberRecord& qr)
{
    os << "Upper: " << qr.Upper() << " ";
    os << "Lower: " << qr.Lower();

    return os;
}


std::ostream& operator<<(std::ostream& os, const QuantumIdentifier& qi)
{
    using global_data::species_data;

    os << species_data[qi.Species()].Name() << "-"
    << species_data[qi.Species()].Isotopologue()[qi.Isotopologue()].Name()
    << " ";

    if (qi.Type() == QuantumIdentifier::TRANSITION)
    {
        os << "TR UP " << qi.QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX];
        os << " LO " << qi.QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX];
    }
    else if (qi.Type() == QuantumIdentifier::ENERGY_LEVEL)
    {
        os << "EN " << qi.QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL_INDEX];
    }
    else
    {
        assert(0);
    }

    return os;
}

