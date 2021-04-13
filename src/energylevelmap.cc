/* Copyright (C) 2019
   Richard Larsson <larsson@mps.mpg.de>

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

/**
 * @file energylevelmap.h
 * @author Richard Larsson
 * @date 2019-10-28
 * 
 * @brief Class to map energy levels
 */

#include "energylevelmap.h"
#include "special_interp.h"


bool EnergyLevelMap::OK() const noexcept {
  if (not (mvalue.nbooks() == mlevels.nelem() and 
          (mvib_energy.nelem() == mlevels.nelem() or mvib_energy.nelem() == 0))) {
    return false;  // Bad dimensions, vibrational energies and IDs and data of strange size
  }
  
  if (mtype == EnergyLevelMapType::Tensor3_t) {
  } else if (mtype == EnergyLevelMapType::Vector_t) {
    if (mvalue.npages() not_eq 1 or mvalue.nrows() not_eq 1) {
      return false;  // Bad dimensions for vector type
    }
  } else if (mtype == EnergyLevelMapType::Numeric_t) {
    if (mvalue.npages() not_eq 1 or mvalue.nrows() not_eq 1 or mvalue.ncols() not_eq 1) {
      return false;  // Bad dimensions for numeric type
    }
  } else if (mtype == EnergyLevelMapType::None_t) {
    if (mvalue.npages() not_eq 0 or mvalue.nrows() not_eq 0 or mvalue.ncols() not_eq 0) {
      return false;  // Bad dimensions for none type
    }
  }
  
  for (auto& e: mvib_energy) if (e < 0) return false;  // Bad energies
  
  return true;
}


Output2 EnergyLevelMap::get_ratio_params(
  const AbsorptionLines& band, const Index& line_index) const
{
  ARTS_USER_ERROR_IF (mtype not_eq EnergyLevelMapType::Numeric_t,
                      "Must have Numeric_t, input type is bad");
    
  Output2 x{/*.r_low=*/0, /*.r_upp=*/0};
  
  bool found1=false;
  bool found2=false;
  for (size_t i=0; i<mlevels.size(); i++) {
    const Absorption::QuantumIdentifierLineTarget lt = Absorption::QuantumIdentifierLineTarget(mlevels[i], band, line_index);
    if (lt == Absorption::QuantumIdentifierLineTargetType::Level and lt.lower) {
      found1 = true;
      x.r_low = mvalue(i, 0, 0, 0);
    }
    
    if (lt == Absorption::QuantumIdentifierLineTargetType::Level and lt.upper) {
      found2 = true;
      x.r_upp = mvalue(i, 0, 0, 0);
    }
    
    if (found1 and found2)
      break;
  }
  return x;
}

Output4 EnergyLevelMap::get_vibtemp_params(
  const AbsorptionLines& band,
  const Index& line_index,
  const Numeric T) const
{
  ARTS_USER_ERROR_IF (mtype not_eq EnergyLevelMapType::Numeric_t,
                      "Must have Numeric_t, input type is bad");
  
  Output4 x{/*.E_low=*/0, /*.E_upp=*/0, /*.T_low=*/T, /*.T_upp=*/T};
  
  bool found1=false;
  bool found2=false;
  for (Index i=0; i<mlevels.nelem(); i++) {
    const Absorption::QuantumIdentifierLineTarget lt = Absorption::QuantumIdentifierLineTarget(mlevels[i], band, line_index);
    
    if (lt == Absorption::QuantumIdentifierLineTargetType::Level and lt.lower) {
      found1 = true;
      x.T_low = mvalue(i, 0, 0, 0);
      x.E_low = mvib_energy[i];
    }
    
    if (lt == Absorption::QuantumIdentifierLineTargetType::Level and lt.upper) {
      found2 = true;
      x.T_upp = mvalue(i, 0, 0, 0);
      x.E_upp = mvib_energy[i];
    }
    
    if (found1 and found2) {
      break;
    }
  }
  return x;
}

EnergyLevelMap::EnergyLevelMap(const Tensor4& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies) :
  mtype(EnergyLevelMapType::Tensor3_t),
  mlevels(levels),
  mvib_energy(energies),
  mvalue(data)
{
  ThrowIfNotOK();
}

EnergyLevelMap::EnergyLevelMap(const Matrix& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies) :
  mtype(EnergyLevelMapType::Vector_t),
  mlevels(levels),
  mvib_energy(energies),
  mvalue(data.nrows(), 1, 1, data.ncols())
{
  mvalue(joker, 0, 0, joker) = data;
  ThrowIfNotOK();
}

EnergyLevelMap::EnergyLevelMap(const Vector& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies) :
  mtype(EnergyLevelMapType::Numeric_t),
  mlevels(levels),
  mvib_energy(energies),
  mvalue(data.nelem(), 1, 1, 1)
{
  mvalue(joker, 0, 0, 0) = data;
  ThrowIfNotOK();
}

EnergyLevelMap EnergyLevelMap::InterpToGridPos(Index atmosphere_dim, const ArrayOfGridPos& p, const ArrayOfGridPos& lat, const ArrayOfGridPos& lon) const
{
  if (mtype == EnergyLevelMapType::None_t)
    return EnergyLevelMap();
  ARTS_USER_ERROR_IF (mtype not_eq EnergyLevelMapType::Tensor3_t,
                      "Must have Tensor3_t, input type is bad");
  
  EnergyLevelMap elm(EnergyLevelMapType::Vector_t, 1, 1, p.nelem(), *this);
  
  Matrix itw_field;
  interp_atmfield_gp2itw(itw_field, atmosphere_dim, p, lat, lon);
  
  const Index nnlte = mlevels.nelem();
  for (Index itnlte = 0; itnlte < nnlte; itnlte++)
    interp_atmfield_by_itw(elm.mvalue(itnlte, 0, 0, joker), atmosphere_dim,
                           mvalue(itnlte, joker, joker, joker),
                           p, lat, lon, itw_field);
  return elm;
}

EnergyLevelMap EnergyLevelMap::operator[](Index ip) const
{
  if (mtype == EnergyLevelMapType::None_t)
    return EnergyLevelMap();
  if (mtype == EnergyLevelMapType::Numeric_t)
    return *this;
  ARTS_USER_ERROR_IF (mtype not_eq EnergyLevelMapType::Vector_t
  ,"Must have Vector_t, input type is bad");
  
  ARTS_USER_ERROR_IF (ip >= mvalue.ncols() or ip < 0,
    "Bad dims for data:\n\tThe pressure dim of data contains: ",
    mvalue.ncols(), " values and you are requesting element index ", ip, "\n")
  
  EnergyLevelMap elm(EnergyLevelMapType::Numeric_t, 1, 1, 1, *this);
  elm.mvalue(joker, 0, 0, 0) = mvalue(joker, 0, 0, ip);
  return elm;
}

EnergyLevelMapType string2energylevelmaptype(const String& s) {
  if (s == "Tensor3")
    return EnergyLevelMapType::Tensor3_t;
  else if (s == "Vector")
    return EnergyLevelMapType::Vector_t;
  else if (s == "Numeric")
    return EnergyLevelMapType::Numeric_t;
  else if (s == "None")
    return EnergyLevelMapType::None_t;
  else {
    ARTS_USER_ERROR ("Only \"None\", \"Numeric\", \"Vector\", and \"Tensor3\" types accepted\n"
                     "You request to have an EnergyLevelMap of type: ", s, '\n')
  }
}

String energylevelmaptype2string(EnergyLevelMapType type)
{
  switch(type) {
    case EnergyLevelMapType::Tensor3_t:
      return "Tensor3";
    case EnergyLevelMapType::Vector_t:
      return "Vector";
    case EnergyLevelMapType::Numeric_t:
      return "Numeric";
    case EnergyLevelMapType::None_t:
      return "None";
  }
  return "";
}

std::ostream& operator<<(std::ostream& os, const EnergyLevelMap& elm) {
  return os << energylevelmaptype2string(elm.Type()) << '\n'
            << elm.Levels() << '\n'
            << elm.Data() << '\n'
            << elm.Energies() << '\n';
}

EnergyLevelMap EnergyLevelMap::operator()(Index ip, Index ilat, Index ilon) const
{
  if (mtype == EnergyLevelMapType::None_t or mtype == EnergyLevelMapType::Numeric_t)
    return *this;
  ARTS_USER_ERROR_IF (mtype not_eq EnergyLevelMapType::Tensor3_t,
                      "Must have Tensor3_t, input type is bad");
  
  auto elm = EnergyLevelMap(EnergyLevelMapType::Numeric_t, 1, 1, 1, *this);
  elm.mvalue(joker, 0, 0, 0) = mvalue(joker, ip, ilat, ilon);
  return elm;
}
