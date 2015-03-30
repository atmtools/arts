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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Partition function interactions to get cross section goes below here
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool PartitionFunctionData::GetPartitionFunctionDataParams(Numeric& part, 
                                                           const Numeric& line_t0, 
                                                           const Numeric& atm_t) const
{
    if(mtype == PF_NONE) // The standard case
        return 1;
    else if(mtype == PF_Coeff) // The Coefficient case
    {
        GetCoeff(part, line_t0, atm_t);
        return 0;
    }
    else
        throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
}


//Note that first order is used by LBLRTM on the data we have.
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
          throw std::runtime_error("The partition function data vector is not of the right length for LBLRTM.\n");
      }
      
      mdata.resize(1);
      
      mdata[0] = input;
}


// This will convert the read vector to the none data format
void PartitionFunctionData::Vector2NoneData(const Vector& input)
{ // Setting mdata.resize(0) is probably a better method for this...
  assert(mtype == PartitionFunctionData::PF_NONE);
  if( input.nelem() != ExpectedVectorLengthFromType() )
    throw std::runtime_error("You are trying to set partition function data to a none line mixed line.\n");
}


// This will convert the stored two char string to LM_Type
void PartitionFunctionData::StorageTag2SetType(const String& input)
 {
  if(input == "NA") // The standard case
    mtype = PF_NONE;
  else if(input == "CN") // The Coefficient case
    mtype = PF_Coeff;
  else
    throw std::runtime_error("You are trying to read a partition function type that is unknown to ARTS.\n");
}
 
 
// This will convert the C4 data format to a vector for storage
void PartitionFunctionData::CoeffData2Vector(Vector& output) const
{
  assert(mtype==PF_Coeff);
  assert(output.nelem()==mnelem);
  
  output = mdata[0];
}


void PartitionFunctionData::GetVectorFromData(Vector& output) const
{
    if(mtype == PF_NONE) // The standard case
        output.resize(0);
    else if(mtype == PF_Coeff) // The Coefficient case
        CoeffData2Vector(output);
    else
        throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
}


// This will convert LM_Type to a two char string for storage
String PartitionFunctionData::Type2StorageTag() const
{
  String output;
  output.resize(2); // ARTS format specify that this is the size of the tag
  
  if(mtype == PF_NONE) // The standard case
    output = "NA";
  else if(mtype == PF_Coeff) // The Coefficient case
      output = "CN"; 
  else
    throw std::runtime_error("You are trying to store a partition function type that is unknown to ARTS.\n");
  return output;
}