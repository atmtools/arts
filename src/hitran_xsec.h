/* Copyright (C) 2018 Oliver Lemke <oliver.lemke@uni-hamburg.de>

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

/*!
  \file   hitran_xsec.h
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/

#ifndef HITRAN_XSEC_H
#define HITRAN_XSEC_H

#include "arts.h"

#include "array.h"
#include "mystring.h"
#include "matpackI.h"
#include "bifstream.h"
#include "messages.h"

class XsecRecord
{
public:
    /** Return species index */
    Index Species() const { return mspecies; };

    /** Return species name */
    String SpeciesName() const;

    /** Set species name */
    void SetSpecies(const Index species) { mspecies = species; };

    /** Get coefficients */
    ConstVectorView Coeffs() const { return mcoeffs; };

    /** Get reference pressures */
    ConstVectorView RefPressure() const { return mrefpressure; };

    /** Get reference temperatures */
    ConstVectorView RefTemperature() const { return mreftemperature; };

    /** Get frequency grids of cross sections */
    const ArrayOfVector& Fgrids() const { return mfgrids; };

    /** Get cross sections */
    const ArrayOfVector& Xsecs() const { return mxsecs; };

    /** Get slope of temperature fit */
    const ArrayOfVector& TemperatureSlope() const { return mtslope; };

    /** Get intersect of temperature fit */
    const ArrayOfVector& TemperatureIntersect() const { return mtintersect; };


    /** Interpolate cross section data.

     Interpolate Xsec data to given frequency vector and given scalar pressure.
     Uses third order interpolation in both coordinates, if grid length allows,
     otherwise lower order or no interpolation.

     \param[out] result     Xsec value for given frequency grid and temperature.
     \param[in] f_grid      Frequency grid.
     \param[in] pressure    Scalar pressure.
     \param[in] temperature Scalar temperature.
     \param[in] apply_tfit  Set to 0 to not apply the temperature fit
     \param[in] verbosity   Standard verbosity object.
     */
    void Extract(VectorView result,
                 ConstVectorView f_grid,
                 const Numeric& pressure,
                 const Numeric& temperature,
                 const Index& apply_tfit,
                 const Verbosity& verbosity) const;

    friend void xml_read_from_stream(std::istream& is_xml,
                                     XsecRecord& cr,
                                     bifstream *pbifs,
                                     const Verbosity& verbosity);

private:
    Index mspecies;
    Vector mcoeffs;
    Vector mrefpressure;
    Vector mreftemperature;
    ArrayOfVector mfgrids;
    ArrayOfVector mxsecs;
    ArrayOfVector mtslope;
    ArrayOfVector mtintersect;
};

typedef Array<XsecRecord> ArrayOfXsecRecord;

Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            const Index species);

std::ostream& operator<<(std::ostream& os, const XsecRecord& xd);

#endif // HITRAN_XSEC_H
