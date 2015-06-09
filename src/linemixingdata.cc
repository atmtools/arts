/* Copyright (C) 2014 
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

/** Contains some additional functionality of the line mixing data class
   \file   linemixingdata.cc
   
   \author Richard Larsson
   \date   2014-10-31
**/

#include "linemixingdata.h"
#include "check_input.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Line mixing interactions to get cross section goes below here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void LineMixingData::GetLineMixingParams(Numeric& Y, Numeric& G, Numeric& DV, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit, const Index& order=1) const
{
    if(mtype == LM_NONE) {} // The standard case
    else if(mtype == LM_LBLRTM) // The LBLRTM case
        GetLBLRTM(Y,G,Temperature,Pressure,Pressure_Limit,order);
    else if(mtype == LM_LBLRTM_O2NonResonant) // The LBLRTM case
        GetLBLRTM_O2NonResonant(G); // Note that they only act like Y and G
    else if(mtype == LM_2NDORDER) // The 2nd order case
        Get2ndOrder(Y,G,DV,Temperature,Pressure,Pressure_Limit);
    else if(mtype == LM_1STORDER) // The 1st order case
        Get1stOrder(Y,Temperature,Pressure,Pressure_Limit);
    else
        throw std::runtime_error("You are trying to store a line mixing type that is unknown to ARTS.\n");
}



//Note that first order is used by LBLRTM on the data we have.
void LineMixingData::GetLBLRTM(Numeric& Y, Numeric& G, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit, const Index& order=1) const
{
      assert( mtype == LM_LBLRTM );
      assert(mdata.nelem() == 3);
      assert(mdata[0].nelem() == 4 && mdata[1].nelem() == 4 && mdata[2].nelem() == 4);
      
      if(Pressure>Pressure_Limit)
      {
        // Helper to understand the following interpolation
        const Vector& t = mdata[0];
        const Vector& y = mdata[1];
        const Vector& g = mdata[2];
        
        const Vector T0(1,Temperature);
        Vector tmp(1);
        
        // Interpolation variables
        ArrayOfGridPosPoly gp(1);
        
        Matrix itw;
        itw.resize(gp.nelem(),order+1);

        chk_interpolation_grids("Line mixing data temperature interpolation",
                                t,
                                T0,
                                order);

        // Interpolation variale determination
        gridpos_poly(gp, t, T0, order);
        interpweights(itw, gp);
        
        // Interpolated values
        interp(tmp, itw, y, gp);
        Y = tmp[0] * Pressure;
        interp(tmp,itw, g, gp);
        G = tmp[0] * Pressure * Pressure;
      }
      else
        Y=G=0;
}


void LineMixingData::GetLBLRTM_O2NonResonant(Numeric& G) const
{
      assert( mtype == LM_LBLRTM_O2NonResonant );
      assert(mdata.nelem() == 1);
      assert(mdata[0].nelem() == 1);
     
      G = mdata[0][0];
}


void LineMixingData::Get2ndOrder(Numeric& Y, Numeric& G, Numeric& DV, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const
{
      assert( mtype == LM_2NDORDER );
      assert(mdata.nelem() == 4);
      assert(mdata[0].nelem() == 1 && mdata[1].nelem() == 3 && mdata[2].nelem() == 3 && mdata[3].nelem() == 3);
      
      if(Pressure>Pressure_Limit)
      {
        // Helper to understand the following interpolation
        const Numeric& T0  = mdata[0][0];
        const Vector& y  = mdata[1];
        const Vector& g  = mdata[2];
        const Vector& dv = mdata[3];
        
        Y  =  ( ( y[0] + y[1] * ( T0/Temperature-1. ) ) * pow( T0/Temperature, y[2] ) ) * Pressure;
        G  =  ( ( g[0] + g[1] * ( T0/Temperature-1. ) ) * pow( T0/Temperature, g[2] ) ) * Pressure * Pressure;
        DV  = ( ( dv[0] + dv[1] * ( T0/Temperature-1. ) ) * pow( T0/Temperature, dv[2] ) ) * Pressure * Pressure;
      }
      else
        Y=DV=G=0;
}


void LineMixingData::Get1stOrder(Numeric& Y, const Numeric& Temperature, const Numeric& Pressure, const Numeric& Pressure_Limit) const
{
    assert( mtype == LM_1STORDER );
    assert(mdata.nelem() == 3);
    assert(mdata[0].nelem() == 1 && mdata[1].nelem() == 1 && mdata[2].nelem() == 1);
    
    if(Pressure>Pressure_Limit)
    {
        // Helper to understand the following interpolation
        const Numeric& T0  = mdata[0][0];
        const Numeric& y  = mdata[1][0];
        const Numeric& x  = mdata[2][0];
        
        Y  =  y * pow(T0/Temperature, x) * Pressure;
    }
    else
        Y=0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Line mixing interactions to get cross section goes above here
// Below is line mixing storage functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This will parse any Vector by the own mtype to the right settings for mdata
void LineMixingData::SetDataFromVectorWithKnownType(const Vector& input)
{
  if(mtype == LM_NONE) // The standard case
    Vector2NoneData(input);
  else if(mtype == LM_LBLRTM) // The LBLRTM case
    Vector2LBLRTMData(input);
  else if(mtype == LM_LBLRTM_O2NonResonant) // The LBLRTM case
    Vector2LBLRTM_O2NonResonantData(input);
  else if(mtype == LM_2NDORDER) // The 2nd order case
    Vector2SecondOrderData(input);
  else if(mtype == LM_1STORDER) // The 1st order case
      Vector2FirstOrderData(input);
  else
    throw std::runtime_error("You are trying to store a line mixing type that is unknown to ARTS.\n");
}


// This will be used to know how many parameters must be read from the catalog
Index LineMixingData::ExpectedVectorLengthFromType()
{
  if(mtype == LM_NONE) // The standard case
    return 0;
  else if(mtype == LM_LBLRTM) // The LBLRTM case
    return 12;
  else if(mtype == LM_LBLRTM_O2NonResonant) // Nonresonant is just a tag
    return 1;
  else if(mtype == LM_2NDORDER) // The 2nd order case
    return 10;
  else if(mtype == LM_1STORDER) // The 2nd order case
      return 3;
  else
    throw std::runtime_error("You are trying to store a line mixing type that is unknown to ARTS.\n");
  return 0;
}


// This will convert the read vector to the LBLRTM data format
void LineMixingData::Vector2LBLRTMData(const Vector& input)
{
      assert(mtype == LineMixingData::LM_LBLRTM);
      
      if(input.nelem() != ExpectedVectorLengthFromType())
      {
          throw std::runtime_error("The line mixing data vector is not of the right length for LBLRTM.\n");
      }
      
      // Then this is a three-long ArrayOfVector
      mdata.resize(3);
      
      // This is supposed to be the temperature vector
      mdata[0].resize(4);
      mdata[0][0] = input[0];
      mdata[0][1] = input[1];
      mdata[0][2] = input[2];
      mdata[0][3] = input[3];
      
      // This is supposed to be the Y-vector
      mdata[1].resize(4);
      mdata[1][0] = input[4];
      mdata[1][1] = input[5];
      mdata[1][2] = input[6];
      mdata[1][3] = input[7];
      
      // This is supposed to be the G vector
      mdata[2].resize(4);
      mdata[2][0] = input[8];
      mdata[2][1] = input[9];
      mdata[2][2] = input[10];
      mdata[2][3] = input[11];
}


// This will convert the read vector to the LBLRTM O2 non-resonant data format
void LineMixingData::Vector2LBLRTM_O2NonResonantData(const Vector& input)
{
      assert(mtype == LineMixingData::LM_LBLRTM_O2NonResonant);
      if(input.nelem() != ExpectedVectorLengthFromType())
      {
          throw std::runtime_error("The line mixing data vector is not of the right length for LBLRTM non-resonant.\n");
      }
      
      mdata.resize(1);
      mdata[0].resize(1);
      mdata[0][0] = input[0];
}


// This will convert the read vector to the none data format
void LineMixingData::Vector2NoneData(const Vector& input)
{ // Setting mdata.resize(0) is probably a better method for this...
  assert(mtype == LineMixingData::LM_NONE);
  if( input.nelem() != ExpectedVectorLengthFromType() )
    throw std::runtime_error("You are trying to set line mixing data to a none line mixed line.\n");
}


// This will convert the read vector to the 2nd order data format
void LineMixingData::Vector2SecondOrderData(const Vector& input)
{
      assert(mtype == LineMixingData::LM_2NDORDER);
      if(input.nelem() != ExpectedVectorLengthFromType())
        throw std::runtime_error("The line mixing data vector is not of the right length for 2ndOrder.\n");
      
      // Then this is a four-long ArrayOfVector
      mdata.resize(4);
      
      // This is supposed to be the temperature vector
      mdata[0].resize(1);
      mdata[0][0] = input[6];
      
      // This is supposed to be the Y components
      mdata[1].resize(3);
      mdata[1][0] = input[0];
      mdata[1][1] = input[1];
      mdata[1][2] = input[7];
      
      // This is supposed to be the G components
      mdata[2].resize(3);
      mdata[2][0] = input[2];
      mdata[2][1] = input[3];
      mdata[2][2] = input[8];
      
      // This is supposed to be the DV components
      mdata[3].resize(3);
      mdata[3][0] = input[4];
      mdata[3][1] = input[5];
      mdata[3][2] = input[9];
}


// This will convert the read vector to the 1st order data format
void LineMixingData::Vector2FirstOrderData(const Vector& input)
{
    assert(mtype == LineMixingData::LM_1STORDER);
    if(input.nelem() != ExpectedVectorLengthFromType())
        throw std::runtime_error("The line mixing data vector is not of the right length for 1stOrder.\n");
    
    // Then this is a three-long ArrayOfVector
    mdata.resize(3);
    
    // This is supposed to be the temperature vector
    mdata[0].resize(1);
    mdata[0][0] = input[0];
    
    // This is supposed to be value
    mdata[1].resize(1);
    mdata[1][0] = input[1];
    
    // This is supposed to be exponential value
    mdata[2].resize(1);
    mdata[2][0] = input[2];
    
}


// This will convert the stored two char string to LM_Type
void LineMixingData::StorageTag2SetType(const String& input)
 {
  if(input == "NA") // The standard case
    mtype = LM_NONE;
  else if(input == "LL") // The LBLRTM case
    mtype = LM_LBLRTM;
  else if(input == "NR") // The LBLRTM O2 non-resonant case
    mtype = LM_LBLRTM_O2NonResonant;
  else if(input == "L2") // The 2nd order case
    mtype = LM_2NDORDER;
  else if(input == "L1") // The 2nd order case
      mtype = LM_1STORDER;
  else
    throw std::runtime_error("You are trying to read a line mixing type that is unknown to ARTS.\n");
}
 
 
// This will convert the LBLRTM data format to a vector for storage
void LineMixingData::LBLRTMData2Vector(Vector& output) const
{
      output.resize(12);
      
      // This is the T-vector
      output[0]  = mdata[0][0];
      output[1]  = mdata[0][1];
      output[2]  = mdata[0][2];
      output[3]  = mdata[0][3];
      
      // This is the Y-vector
      output[4]  = mdata[1][0];
      output[5]  = mdata[1][1];
      output[6]  = mdata[1][2];
      output[7]  = mdata[1][3];
      
      // This is the G-vector
      output[8]  = mdata[2][0];
      output[9]  = mdata[2][1];
      output[10] = mdata[2][2];
      output[11] = mdata[2][3];
}


// This will convert the LBLRTM data format to a vector for storage
void LineMixingData::LBLRTM_O2NonResonantData2Vector(Vector& output) const
{
      output.resize(1);
      output[0]=mdata[0][0];
}


// This will convert the 2ndOrder data format to a vector for storage
void LineMixingData::SecondOrderData2Vector(Vector& output) const
{
      output.resize(10);
      
      // This is the temperature vector
      output[6] = mdata[0][0];
      
      // This is the Y components
      output[0] = mdata[1][0];
      output[1] = mdata[1][1];
      output[7] = mdata[1][2];
      
      // This is the G components
      output[2] = mdata[2][0];
      output[3] = mdata[2][1];
      output[8] = mdata[2][2];
      
      // This is the DV components
      output[4] = mdata[3][0];
      output[5] = mdata[3][1];
      output[9] = mdata[3][2];
}


// This will convert the 1stOrder data format to a vector for storage
void LineMixingData::FirstOrderData2Vector(Vector& output) const
{
    output.resize(3);
    
    // This is the temperature
    output[0] = mdata[0][0];
    
    // This is the y-value
    output[1] = mdata[1][0];
    
     // This is the exponent
    output[2] = mdata[2][0];
}


void LineMixingData::GetVectorFromData(Vector& output) const
{
    if(mtype == LM_NONE) // The standard case
        output.resize(0);
    else if(mtype == LM_LBLRTM) // The LBLRTM case
        LBLRTMData2Vector(output);
    else if(mtype == LM_LBLRTM_O2NonResonant) // The LBLRTM O2 non-resonant case
        LBLRTM_O2NonResonantData2Vector(output);
    else if(mtype == LM_2NDORDER) // The 2nd order case
        SecondOrderData2Vector(output);
    else if(mtype == LM_1STORDER) // The 1st order case
        FirstOrderData2Vector(output);
    else
        throw std::runtime_error("You are trying to store a line mixing type that is unknown to ARTS.\n");
}


// This will convert LM_Type to a two char string for storage
String LineMixingData::Type2StorageTag() const
{
  String output;
  output.resize(2); // ARTS format specify that this is the size of the tag
  
  if(mtype == LM_NONE) // The standard case
    output = "NA";
  else if(mtype == LM_LBLRTM) // The LBLRTM case
    output = "LL";
  else if(mtype == LM_LBLRTM_O2NonResonant) // The LBLRTM O2 non-resonant case
    output = "NR";
  else if(mtype == LM_2NDORDER) // The 2nd order case
    output = "L2"; 
  else if(mtype == LM_1STORDER) // The 2nd order case
      output = "L1"; 
  else
    throw std::runtime_error("You are trying to store a line mixing type that is unknown to ARTS.\n");

  return output;
}