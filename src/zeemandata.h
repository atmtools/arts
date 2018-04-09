/* Copyright (C) 2018 Richard Larsson

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

/** Contains Zeeman Effect data class
 * \file   zeemandata.h
 * 
 * \author Richard Larsson
 * \date   2018-04-06
 **/

#ifndef zeemandata_h
#define zeemandata_h

#include "matpackI.h"
#include "mystring.h"
#include "quantum.h"

enum class ZeemanSplittingType : Index { None, ByHund, ByGData, ByPrecalc };
enum class ZeemanPolarizationType : Index { None, Pi, SigmaPlus, SigmaMinus };

class ZeemanEffectData 
{
public:
  enum class ByPrecalcPos : Index { DF=0, LEN=1 };
  enum class ByGDataPos : Index { GU=0, GL=1, LEN=2 };
  enum class ByHundPos : Index { LEN=0 };
  enum class NonePos : Index { LEN=0 };
  
  ZeemanEffectData() : msplit(ZeemanSplittingType::None), mpolar(ZeemanPolarizationType::None), mdata(Index(NonePos::LEN)) {}
  
  ZeemanSplittingType SplittingType() const { return msplit; };
  ZeemanPolarizationType PolarizationType() const { return mpolar; };
  
  bool ok() const
  {
    if(not (nelem() == mdata.nelem())) return false;
    return true;
  };
  
  Index nelem() const
  {
    switch(msplit) {
      case ZeemanSplittingType::None:        return Index(NonePos::LEN);
      case ZeemanSplittingType::ByHund:      return Index(ByHundPos::LEN);
      case ZeemanSplittingType::ByGData:     return Index(ByGDataPos::LEN);
      case ZeemanSplittingType::ByPrecalc:   return Index(ByPrecalcPos::LEN);
      default: throw std::runtime_error("cannot determine size of type");
    }
  };
  
  String splitting_tag() const noexcept
  {
    switch(msplit) {
      case ZeemanSplittingType::ByHund: return "BH";
      case ZeemanSplittingType::ByGData: return "BG";
      case ZeemanSplittingType::ByPrecalc: return "DF";
      case ZeemanSplittingType::None: default: return "NoSplit";
    }
  };
  
  String polarization_tag() const noexcept
  {
    switch(mpolar) {
      case ZeemanPolarizationType::Pi:         return "Pi";
      case ZeemanPolarizationType::SigmaPlus:  return "S+";
      case ZeemanPolarizationType::SigmaMinus: return "S-";
      case ZeemanPolarizationType::None: default: return "NoPolarization";
    }
  };
  
  ConstVectorView data() const { return mdata; };
  
  void setSplittingType(const String& tag) 
  {
    if(tag == "BH") { msplit = ZeemanSplittingType::ByHund; mdata = Vector(Index(ByHundPos::LEN)); }
    else if(tag == "BG") { msplit = ZeemanSplittingType::ByGData; mdata = Vector(Index(ByGDataPos::LEN), 0.0); }
    else if(tag == "DF") { msplit = ZeemanSplittingType::ByPrecalc; mdata = Vector(Index(ByPrecalcPos::LEN), 0.0); }
    else throw std::runtime_error("Cannot recognize provided Zeeman data type");
  };
  
  void setPolarizationType(const String& tag) 
  {
    if(tag == "NoPolarization") mpolar = ZeemanPolarizationType::None;
    else if(tag == "Pi") mpolar = ZeemanPolarizationType::Pi;
    else if(tag == "S+") mpolar = ZeemanPolarizationType::SigmaPlus;
    else if(tag == "S-") mpolar = ZeemanPolarizationType::SigmaMinus;
    else throw std::runtime_error("Cannot recognize provided Zeeman data type");
  };
  
  void setDataFromVectorWithKnownSplittingType(const ConstVectorView data)
  {
    mdata = data;
    if(not ok())
      throw std::runtime_error("Bad data for type");
  };
  
  Numeric frequency_shift_per_tesla(const QuantumNumberRecord& qnr, const Index species) const;
  
  void convertNoneToHund(const QuantumNumberRecord& qnr);
  
  void setNumericalAndPolarization(const QuantumNumberRecord& qnr, const Index species);
  
  Index dM();
  
private:
  ZeemanSplittingType    msplit;
  ZeemanPolarizationType mpolar;
  Vector                 mdata;
};

#endif /* zeemandata_h */
