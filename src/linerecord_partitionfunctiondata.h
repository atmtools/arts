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

/** Contains the partition function data class definition
 * \file   linerecord_partitionfunctiondata.h
 * 
 * \author Richard Larsson
 * \date   2015-03-27
 **/

#ifndef linerecord_partitionfunctiondata_h
#define linerecord_partitionfunctiondata_h

#include <stdexcept>
#include "matpackI.h"
#include "array.h"
#include "mystring.h"

class PartitionFunctionData
{
public:
    
    enum PF_Type {
      PF_NONE,                          // Reserved for no linerecord partition function data.  Defaults to abs_species partition function data
      PF_Coeff                          // Reserved for coefficients method similar to default
    };
  
    // Defining an object
    PartitionFunctionData() : mtype(PF_NONE), mdata(), mnelem(0) {}
    
    // Use these to get the raw data from this class
    const PF_Type& Type() const {return mtype;}
    const ArrayOfVector& Data() const {return mdata;}
    const Index& GetNelem() const {return mnelem;}
    
    // Use these to return data in the format required by the line shape calculator
    Index GetPartitionFunctionDataParams(Numeric& part, const Numeric& line_t0, const Numeric& atm_t) const;
    void GetCoeff(Numeric& part, const Numeric& line_t0, const Numeric& atm_t) const;
    
    // Use these to read data from XML-formats
    void StorageTag2SetType(const String& input);
    void SetNelem(const Index& input) {mnelem = input;}
    Index ExpectedVectorLengthFromType() const {return mnelem;}
    void SetDataFromVectorWithKnownType(const Vector& input);
    void Vector2CoeffData(const Vector& input);
    void Vector2NoneData(const Vector&);
    
    // Use these to save output vector in ARTS catalog
    void GetVectorFromData(Vector& output) const;
    void CoeffData2Vector(Vector& output) const;
    String Type2StorageTag() const;
    
private:
    // mtype identifies the type of partition function and mdata should contain the required data
    PF_Type mtype;
    ArrayOfVector mdata;
    Index mnelem;
};

#endif // linerecord_partitionfunctiondata_h