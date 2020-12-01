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

#include "array.h"
#include "arts.h"
#include "bifstream.h"
#include "gridded_fields.h"
#include "matpackI.h"
#include "messages.h"
#include "mystring.h"

class XsecRecord {
 public:
  /** Return species index */
  Species::Species Species() const { return mspecies; };

  /** Return species name */
  String SpeciesName() const;

  /** Set species name */
  void SetSpecies(const Species::Species species) { mspecies = species; };

  /** Return species index */
  Index Version() const { return mversion; };

  /** Set species name */
  void SetVersion(const Index version);

  /************ VERSION 1 *************/
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

  /** Get coefficients */
  Vector& Coeffs() { return mcoeffs; };

  /** Get reference pressures */
  Vector& RefPressure() { return mrefpressure; };

  /** Get reference temperatures */
  Vector& RefTemperature() { return mreftemperature; };

  /** Get frequency grids of cross sections */
  ArrayOfVector& Fgrids() { return mfgrids; };

  /** Get cross sections */
  ArrayOfVector& Xsecs() { return mxsecs; };

  /** Get slope of temperature fit */
  ArrayOfVector& TemperatureSlope() { return mtslope; };

  /** Get intersect of temperature fit */
  ArrayOfVector& TemperatureIntersect() { return mtintersect; };

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

  /************ VERSION 2 *************/
  /** Get mininum pressures from fit */
  const Vector& FitMinPressures() const { return mfitminpressures; };

  /** Get maximum pressures from fit */
  const Vector& FitMaxPressures() const { return mfitmaxpressures; };

  /** Get mininum temperatures from fit */
  const Vector& FitMinTemperatures() const { return mfitmintemperatures; };

  /** Get maximum temperatures */
  const Vector& FitMaxTemperatures() const { return mfitmaxtemperatures; };

  /** Get coefficients */
  const ArrayOfGriddedField2& FitCoeffs() const { return mfitcoeffs; };

  /** Get mininum pressures from fit */
  Vector& FitMinPressures() { return mfitminpressures; };

  /** Get maximum pressures from fit */
  Vector& FitMaxPressures() { return mfitmaxpressures; };

  /** Get mininum temperatures from fit */
  Vector& FitMinTemperatures() { return mfitmintemperatures; };

  /** Get maximum temperatures */
  Vector& FitMaxTemperatures() { return mfitmaxtemperatures; };

  /** Get coefficients */
  ArrayOfGriddedField2& FitCoeffs() { return mfitcoeffs; };

 private:
  void Extract1(VectorView result,
                ConstVectorView f_grid,
                const Numeric& pressure,
                const Numeric& temperature,
                const Index& apply_tfit,
                const Verbosity& verbosity) const;

  void Extract2(VectorView result,
                ConstVectorView f_grid,
                const Numeric pressure,
                const Numeric temperature,
                const Verbosity& verbosity) const;

  void CalcXsec(VectorView& xsec,
                const Index dataset,
                const Range range,
                const Numeric pressure,
                const Numeric temperature) const;

  void CalcDT(VectorView& xsec_dt,
              const Index dataset,
              const Range range,
              const Numeric pressure,
              const Numeric temperature) const;

  void CalcDP(VectorView& xsec_dp,
              const Index dataset,
              const Range range,
              const Numeric pressure,
              const Numeric temperature) const;

  static const Index P00 = 0;
  static const Index P10 = 1;
  static const Index P01 = 2;
  static const Index P20 = 3;
  static const Index P11 = 4;
  static const Index P02 = 5;

  Index mversion;
  Index mspecies;
  /* VERSION 1 */
  Vector mcoeffs;
  Vector mrefpressure;
  Vector mreftemperature;
  ArrayOfVector mfgrids;
  ArrayOfVector mxsecs;
  ArrayOfVector mtslope;
  ArrayOfVector mtintersect;
  /* VERSION 2 */
  Vector mfitminpressures;
  Vector mfitmaxpressures;
  Vector mfitmintemperatures;
  Vector mfitmaxtemperatures;
  ArrayOfGriddedField2 mfitcoeffs;
};

typedef Array<XsecRecord> ArrayOfXsecRecord;

Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            const Species::Species species);

std::ostream& operator<<(std::ostream& os, const XsecRecord& xd);

#endif  // HITRAN_XSEC_H
