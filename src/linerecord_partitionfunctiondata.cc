/* Copyright (C) 2015 
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

/** Contains some additional functionality of the partition function data class
   \file   linerecord_partitionfunctiondata.cc
   
   \author Richard Larsson
   \date   2015-03-27
**/

#include "linerecord_partitionfunctiondata.h"
#include <cmath>
#include "check_input.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Partition function interactions to get cross section goes below here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Index PartitionFunctionData::GetPartitionFunctionDataParams(Numeric& part, 
                                                            const Numeric& line_t0, 
                                                            const Numeric& atm_t) const
{
    if(mtype == PF_NONE) // The standard case
        return 1;
    else if(mtype == PF_Coeff) // The Coefficient case
    {
        GetCoeff(part, line_t0, atm_t);
    }
    else
        throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
    return 0;
}


//Similar to old partition functionality.
void PartitionFunctionData::GetCoeff(Numeric& part, const Numeric& line_t0, const Numeric& atm_t) const
{
  assert(mtype == PF_Coeff );
  assert(mdata.nelem() == 1);
  
  Numeric result_t0 = 0.;
  Numeric result_t  = 0.;
  Numeric exponent_t0 = 1.;
  Numeric exponent_t  = 1.;

  for (Index ii=0; ii < mnelem; ii++)
    {
      result_t    += mdata[0][ii] * exponent_t;
      exponent_t  *= atm_t;
      result_t0   += mdata[0][ii] * exponent_t0;
      exponent_t0 *= line_t0;
    }
    part = result_t0 / result_t;
}


//Similar to old partition functionality.
// void PartitionFunctionData::GetTv(Numeric& part, const Numeric& line_t0, const Numeric& atm_t, const Numeric& atm_p, const Numeric& E_low, const Numeric& E_up) const
// {    
//     extern const Numeric BOLTZMAN_CONST;
//     
//     const Numeric& Ev_low      = mdata[0][0];
//     const Numeric& Ev_up       = mdata[0][1];
//     const Vector&  p_grid      = mdata[1];
//     const Vector&  tv_low_grid = mdata[2];
//     const Vector&  tv_up_grid  = mdata[3];
//     
//     const Vector P0(1,atm_p);
//     Vector tmp1(1),tmp2(1);
//     
//     // Interpolation variables
//     ArrayOfGridPosPoly gp(1);
//     
//     Matrix itw;
//     itw.resize(gp.nelem(),2);
//     chk_interpolation_grids("Partition function NLTE pressure interpolation",
//                             p_grid,
//                             P0,
//                             1);
//     
//     // Interpolation variale determination
//     gridpos_poly(gp, p_grid, P0, 1);
//     interpweights(itw, gp);
//     
//     // Interpolated values
//     interp(tmp1, itw, tv_low_grid, gp);
//     interp(tmp2, itw, tv_up_grid, gp);
//     
//     // Set temperatures but if the temperature is negative this is a flag for setting temperature to atmospheric temperature
//     const Numeric& Tv_low = tmp1[0]<0?atm_t:tmp1[0];
//     const Numeric& Tv_up  = tmp2[0]<0?atm_t:tmp2[0];
//     
//     const Numeric denom = exp(- E_low / ( BOLTZMAN_CONST * line_t0 ) ) -
//     exp(- E_up / ( BOLTZMAN_CONST * line_t0 ) );
//     
//     const Numeric nom = exp(-Ev_low / ( BOLTZMAN_CONST * Tv_low ) ) * exp(- (E_low-Ev_low) / ( BOLTZMAN_CONST * atm_t ) ) -
//     exp(-Ev_up / ( BOLTZMAN_CONST * Tv_up ) ) * exp(- (E_up-Ev_up) / ( BOLTZMAN_CONST * atm_t ) );
//     
//     part = nom/denom;
// }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Partition function interactions to get cross section goes above here
// Below is partition function storage functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This will parse any Vector by the own mtype to the right settings for mdata
void PartitionFunctionData::SetDataFromVectorWithKnownType(const Vector& input)
{
  if(mtype == PF_NONE) // The standard case
      Vector2NoneData(input);
  else if(mtype == PF_Coeff) // The Coefficient case
      Vector2CoeffData(input);
  else
    throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
}


// This will convert the read vector to the Coefficient data format
void PartitionFunctionData::Vector2CoeffData(const Vector& input)
{
      assert(mtype == PF_Coeff); 
      
      if(input.nelem() != ExpectedVectorLengthFromType())
      {
          throw std::runtime_error("The partition function data vector is not of the right length for the coefficient case.\n");
      }
      
      mdata.resize(1);
      
      mdata[0] = input;
}

// // This will convert the read vector to the Coefficient data format
// void PartitionFunctionData::Vector2TvData(const Vector& input)
// {
//     assert(mtype == PF_Tv); 
//     
//     if(input.nelem() != ExpectedVectorLengthFromType())
//     {
//         throw std::runtime_error("The partition function data vector is not of the right length for vibrational temperature case.\n");
//     }
//     
//     mdata.resize(4);
//     
//     // The vibrational energy levels
//     mdata[0].resize(2);
//     mdata[0][0] = input[0]; // Lower vibrational energy level
//     mdata[0][1] = input[1]; // Upper vibrational energy level
//     
//     
//     mdata[1].resize((mnelem-2)/3);//Pressure Grid
//     mdata[2].resize((mnelem-2)/3);//Tv low Grid
//     mdata[3].resize((mnelem-2)/3);//Tv up Grid
//     
//     // Loop and insert values in the right spot
//     for(Index ii=0;ii<(mnelem-2)/3; ii++)
//     {
//         mdata[1][ii] = input[2+ii];
//         mdata[2][ii] = input[2+ii+(mnelem-2)/3];
//         mdata[3][ii] = input[2+ii+(mnelem-2)/3*2];
//     }
//     
// }


// This will convert the read vector to the none data format
void PartitionFunctionData::Vector2NoneData(const Vector& input)
{ // Setting mdata.resize(0) is probably a better method for this...
  assert(mtype == PartitionFunctionData::PF_NONE);
  if( input.nelem() != ExpectedVectorLengthFromType() )
    throw std::runtime_error("You are trying to set partition function data to a none line mixed line.\n");
}


// This will convert the stored two char string to PF_Type
void PartitionFunctionData::StorageTag2SetType(const String& input)
 {
  if(input == "NA") // The standard case
    mtype = PF_NONE;
  else if(input == "CN") // The Coefficient case
    mtype = PF_Coeff;
  else
    throw std::runtime_error("You are trying to read a partition function type that is unknown to ARTS.\n");
}
 
 
// This will convert the Coeff data format to a vector for storage
void PartitionFunctionData::CoeffData2Vector(Vector& output) const
{
  assert(mtype==PF_Coeff);
  assert(output.nelem()==mnelem);
  
  output = mdata[0];
}


// This will convert the Tv data format to a vector for storage
// void PartitionFunctionData::TvData2Vector(Vector& output) const
// {
//     assert(mtype==PF_Tv);
//     output.resize(ExpectedVectorLengthFromType());
//     assert(output.nelem()==mnelem);
//     
//     // First two values are for vibrational energy levels
//     output[0]=mdata[0][0];
//     output[1]=mdata[0][1];
//     
//     // Loop and insert values in the right spot
//     for(Index ii=0;ii<(mnelem-2)/3; ii++)
//     {
//         output[2+ii]                = mdata[1][ii];
//         output[2+ii+(mnelem-2)/3]   = mdata[2][ii];
//         output[2+ii+(mnelem-2)/3*2] = mdata[3][ii];
//     }
// }

void PartitionFunctionData::GetVectorFromData(Vector& output) const
{
    if(mtype == PF_NONE) // The standard case
        output.resize(0);
    else if(mtype == PF_Coeff) // The Coefficient case
        CoeffData2Vector(output);
    else
        throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
}


// This will convert PF_Type to a two char string for storage
String PartitionFunctionData::Type2StorageTag() const
{
  String output;

  if(mtype == PF_NONE) // The standard case
    output = "NA";
  else if(mtype == PF_Coeff) // The Coefficient case
    output = "CN"; 
  else
    throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
  return output;
}