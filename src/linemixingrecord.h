/* Copyright 2013, The ARTS Developers.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/** \file
 LineMixingRecord class for storing line mixing data.

 \author Oliver Lemke
 */

#ifndef linemixingrecord_h
#define linemixingrecord_h

#include <iostream>
#include "matpackI.h"
#include "quantum.h"
#include "abs_species_tags.h"


class LineMixingRecord
{
public:
    LineMixingRecord() : mspecies(-1), misotopologue(-1) {};
    LineMixingRecord(Index species, Index iso) : mspecies(species), misotopologue(iso) {};

    Index Species() const { return mspecies; }
    Index Isotopologue() const { return misotopologue; }

    QuantumNumberRecord& Quantum() { return mquantum; }
    const QuantumNumberRecord& Quantum() const { return mquantum; }

    Vector& Data() { return mdata; }
    const Vector& Data() const { return mdata; }

private:
    /** Species index */
    Index mspecies;
    /** Isotopologue index */
    Index misotopologue;
    /** Quantum Numbers */
    QuantumNumberRecord mquantum;
    /** Line mixing data */
    Vector mdata;
};

std::ostream& operator<<(std::ostream& os, const LineMixingRecord& lmr);

typedef Array<LineMixingRecord> ArrayOfLineMixingRecord;
typedef Array<ArrayOfLineMixingRecord> ArrayOfArrayOfLineMixingRecord;

#endif /* linemixingrecord_h */
