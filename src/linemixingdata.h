/* Copyright (C) 2012 
Richard Larsson <ric.larsson@gmail.com>

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

/** Contains the line mixing data class definition
 * \file   linemixingdata.h
 * 
 * \author Richard Larsson
 * \date   2014-10-31
 **/

#ifndef linemixingdata_h
#define linemixingdata_h

#include <stdexcept>
#include <cmath>
#include "interpolation.h"
#include "interpolation_poly.h"
#include "matpackI.h"
#include "array.h"
#include "mystring.h"
#include "complex.h"
#include "jacobian.h"

#ifndef NEWARTSCAT

class LineMixingData
{
public:
    
    enum LM_Type : Index {
      LM_NONE,                          // Reserved for no line mixing
      LM_LBLRTM,                        // Reserved for LBLRTM line mixing
      LM_LBLRTM_O2NonResonant,          // Reserved for the non-resonant O2 line in LBLRTM
      LM_1STORDER,                      // Reserved for Tretyakov et al. 2005 1st order of line mixing
      LM_2NDORDER,                      // Reserved for Makarov et al. 2011 second order of line mixing
      LM_BYBAND                         // Reserved for Paris data of relaxation matrix line mixing for band
    };
  
    // Defining an object
    LineMixingData() : mtype(LM_NONE), mdata() {}
    
    // Use these to get the raw data from this class
    const LM_Type& Type() const {return mtype;}
    const ArrayOfVector& Data() const {return mdata;}
    
    // Use these to return data in the format required by the line shape calculator
    void GetLineMixingParams(Numeric& Y, Numeric& G, Numeric& DV, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit, const Index& order=1) const;
    void GetLineMixingParams_dT(Numeric& dY_dT, Numeric& dG_dG, Numeric& dDV_dT, const Numeric& Temperature, const Numeric& dt, const Numeric& Pressure, const Numeric& Pressure_Limit, const Index& order=1) const;
    void GetLineMixingParams_dZerothOrder(Numeric& dY0, Numeric& dG0, Numeric& dDV0, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void GetLineMixingParams_dFirstOrder(Numeric& dY1, Numeric& dG1, Numeric& dDV1, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void GetLineMixingParams_dExponent(Numeric& dYexp, Numeric& dGexp, Numeric& dDVexp, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get2ndOrder(Numeric& Y, Numeric& G, Numeric& DV, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get2ndOrder_dT(Numeric& dY_dT, Numeric& dG_dT, Numeric& dDV_dT, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get2ndOrder_dZerothOrder(Numeric& dY0, Numeric& dG0, Numeric& dDV0, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get2ndOrder_dFirstOrder(Numeric& dY1, Numeric& dG1, Numeric& dDV1, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get2ndOrder_dExponent(Numeric& dYexp, Numeric& dGexp, Numeric& dDVexp, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get1stOrder(Numeric& Y, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get1stOrder_dT(Numeric& dY_dT, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get1stOrder_dZerothOrder(Numeric& dY0, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void Get1stOrder_dExponent(Numeric& dYexp, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const;
    void GetLBLRTM(Numeric& Y, Numeric& G, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit, const Index& order) const; 
    void GetLBLRTM_dT(Numeric& dY_dT, Numeric& dG_dT, const Numeric& Temperature, const Numeric& dt, const Numeric& Pressure, const Numeric& Pressure_Limit, const Index& order) const;
    void GetLBLRTM_O2NonResonant(Numeric& G) const;
    
    // Use these to change values
    void ChangeY0(const Numeric& change, const bool relative=false);
    void ChangeY1(const Numeric& change, const bool relative=false);
    void ChangeYexp(const Numeric& change, const bool relative=false);
    void ChangeG0(const Numeric& change, const bool relative=false);
    void ChangeG1(const Numeric& change, const bool relative=false);
    void ChangeGexp(const Numeric& change, const bool relative=false);
    void ChangeDF0(const Numeric& change, const bool relative=false);
    void ChangeDF1(const Numeric& change, const bool relative=false);
    void ChangeDFexp(const Numeric& change, const bool relative=false);
    
    // Use these to insert the data in the required format from catalog readings
    void SetLBLRTMFromTheirCatalog(const Vector& t, const Vector& y, const Vector& g) 
    {
      mtype = LM_LBLRTM;
      mdata.resize(3);
      mdata[0] = t;
      mdata[1] = y;
      mdata[2] = g;
    }
    
    void SetLBLRTM_O2NonResonantFromTheirCatalog(const Vector& t, const Vector& gamma1, const Vector& gamma2) 
    {
      mtype = LM_LBLRTM_O2NonResonant;
      mdata.resize(3);
      mdata[0] = t;
      mdata[1] = gamma1;
      mdata[2] = gamma2;
    }
    
    void SetByBandType() { mtype = LM_BYBAND; }
    void Set2ndOrderType() { mtype = LM_2NDORDER; }
    
    // Use these to read data from XML-formats
    void StorageTag2SetType(const String& input);
    void SetTypeFromIndex(const Index& type) {mtype = LM_Type(type);};
    Index ExpectedVectorLengthFromType() const;
    void SetDataFromVectorWithKnownType(ConstVectorView);
    
    // Use these to read data from ARTS catalog
    void Vector2LBLRTMData(const Vector& input);
    void Vector2LBLRTM_O2NonResonantData(const Vector& input);
    void Vector2NoneData(const Vector&);
    void Vector2SecondOrderData(const Vector& input);
    void Vector2FirstOrderData(const Vector& input);
    
    // Use these to save output vector in ARTS catalog
    void GetVectorFromData(Vector& output) const;
    void LBLRTMData2Vector(Vector& output) const;
    void LBLRTM_O2NonResonantData2Vector(Vector& output) const;
    String Type2StorageTag() const;
    void SecondOrderData2Vector(Vector& output) const;
    void FirstOrderData2Vector(Vector& output) const;
    
    void SetInternalDerivatives(ComplexVector& derivatives, const ArrayOfRetrievalQuantity& ppd, const QuantumIdentifier& QI, const Numeric& temperature, const Numeric& pressure, const Numeric& pressure_limit) const;
    
private:
    // mtype identifies the type of line mixing and mdata should contain the required data
    LM_Type mtype;
    ArrayOfVector mdata;
};

#endif

#endif // linemixingdata_h
