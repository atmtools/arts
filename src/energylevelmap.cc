/* Copyright (C) 2019
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
  if (not (value.nbooks() == levels.nelem() and 
          (vib_energy.nelem() == levels.nelem() or vib_energy.nelem() == 0))) {
    return false;  // Bad dimensions, vibrational energies and IDs and data of strange size
  }
  
  if (type == EnergyLevelMapType::Tensor3_t) {
  } else if (type == EnergyLevelMapType::Vector_t) {
    if (value.npages() not_eq 1 or value.nrows() not_eq 1) {
      return false;  // Bad dimensions for vector type
    }
  } else if (type == EnergyLevelMapType::Numeric_t) {
    if (value.npages() not_eq 1 or value.nrows() not_eq 1 or value.ncols() not_eq 1) {
      return false;  // Bad dimensions for numeric type
    }
  } else if (type == EnergyLevelMapType::None_t) {
    if (value.npages() not_eq 0 or value.nrows() not_eq 0 or value.ncols() not_eq 0) {
      return false;  // Bad dimensions for none type
    }
  }

  return std::all_of(vib_energy.begin(),
                     vib_energy.end(),
                     [](const auto& val) { return val >= 0; });
}


Output2 EnergyLevelMap::get_ratio_params(
  const AbsorptionLines& band,
   const Index& line_index) const
{
  ARTS_USER_ERROR_IF (type not_eq EnergyLevelMapType::Numeric_t,
                      "Must have Numeric_t, input type is bad");
    
  Output2 x{/*.r_low=*/0, /*.r_upp=*/0};
  
  bool found1=false;
  bool found2=false;
  for (size_t i=0; i<levels.size(); i++) {
    const Quantum::Number::StateMatch lt(levels[i], band.lines[line_index].localquanta, band.quantumidentity);

    if (lt == Quantum::Number::StateMatchType::Level and lt.low) {
      found1 = true;
      x.r_low = value(i, 0, 0, 0);
    }
    
    if (lt == Quantum::Number::StateMatchType::Level and lt.upp) {
      found2 = true;
      x.r_upp = value(i, 0, 0, 0);
    }
    
    if (found1 and found2)
      break;
  }
  
  return x;
}

Output4 EnergyLevelMap::get_vibtemp_params(
  const AbsorptionLines& band,
  const Numeric T) const
{
  ARTS_USER_ERROR_IF (type not_eq EnergyLevelMapType::Numeric_t,
                      "Must have Numeric_t, input type is bad");
  
  Output4 x{/*.E_low=*/0, /*.E_upp=*/0, /*.T_low=*/T, /*.T_upp=*/T};
  
  bool found1=false;
  bool found2=false;
  for (Index i=0; i<levels.nelem(); i++) {
    const Quantum::Number::StateMatch lt(levels[i], band.quantumidentity);

    if (lt == Quantum::Number::StateMatchType::Level and lt.low) {
      found1 = true;
      x.T_low = value(i, 0, 0, 0);
      x.E_low = vib_energy[i];
    }
    
    if (lt == Quantum::Number::StateMatchType::Level and lt.upp) {
      found2 = true;
      x.T_upp = value(i, 0, 0, 0);
      x.E_upp = vib_energy[i];
    }
    
    if (found1 and found2) {
      break;
    }
  }
  return x;
}

EnergyLevelMap::EnergyLevelMap(Tensor4 data, ArrayOfQuantumIdentifier levels_, Vector energies) :
  type(EnergyLevelMapType::Tensor3_t),
  levels(std::move(levels_)),
  vib_energy(std::move(energies)),
  value(std::move(data))
{
  ThrowIfNotOK();
}

EnergyLevelMap::EnergyLevelMap(const Matrix& data, ArrayOfQuantumIdentifier levels_, Vector energies) :
  type(EnergyLevelMapType::Vector_t),
  levels(std::move(levels_)),
  vib_energy(std::move(energies)),
  value(data.nrows(), 1, 1, data.ncols())
{
  value(joker, 0, 0, joker) = data;
  ThrowIfNotOK();
}

EnergyLevelMap::EnergyLevelMap(const Vector& data, ArrayOfQuantumIdentifier levels_, Vector energies) :
  type(EnergyLevelMapType::Numeric_t),
  levels(std::move(levels_)),
  vib_energy(std::move(energies)),
  value(data.nelem(), 1, 1, 1)
{
  value(joker, 0, 0, 0) = data;
  ThrowIfNotOK();
}

EnergyLevelMap EnergyLevelMap::InterpToGridPos(Index atmosphere_dim, const ArrayOfGridPos& p, const ArrayOfGridPos& lat, const ArrayOfGridPos& lon) const
{
  if (type == EnergyLevelMapType::None_t) return EnergyLevelMap{};
  ARTS_USER_ERROR_IF (type not_eq EnergyLevelMapType::Tensor3_t,
                      "Must have Tensor3_t, input type is bad");
  
  EnergyLevelMap elm(EnergyLevelMapType::Vector_t, 1, 1, p.nelem(), *this);
  
  Matrix itw_field;
  interp_atmfield_gp2itw(itw_field, atmosphere_dim, p, lat, lon);
  
  const Index nnlte = levels.nelem();
  for (Index itnlte = 0; itnlte < nnlte; itnlte++)
    interp_atmfield_by_itw(elm.value(itnlte, 0, 0, joker), atmosphere_dim,
                           value(itnlte, joker, joker, joker),
                           p, lat, lon, itw_field);
  return elm;
}

EnergyLevelMap EnergyLevelMap::operator[](Index ip) const
{
  if (type == EnergyLevelMapType::None_t) return EnergyLevelMap{};
  if (type == EnergyLevelMapType::Numeric_t)
    return *this;
  ARTS_USER_ERROR_IF (type not_eq EnergyLevelMapType::Vector_t
  ,"Must have Vector_t, input type is bad");
  
  ARTS_USER_ERROR_IF (ip >= value.ncols() or ip < 0,
    "Bad dims for data:\n\tThe pressure dim of data contains: ",
    value.ncols(), " values and you are requesting element index ", ip, "\n")
  
  EnergyLevelMap elm(EnergyLevelMapType::Numeric_t, 1, 1, 1, *this);
  elm.value(joker, 0, 0, 0) = value(joker, 0, 0, ip);
  return elm;
}

std::ostream& operator<<(std::ostream& os, const EnergyLevelMap& elm) {
  return os << elm.type << '\n'
            << elm.levels << '\n'
            << elm.value << '\n'
            << elm.vib_energy << '\n';
}

EnergyLevelMap EnergyLevelMap::operator()(Index ip, Index ilat, Index ilon) const
{
  if (type == EnergyLevelMapType::None_t or type == EnergyLevelMapType::Numeric_t)
    return *this;
  ARTS_USER_ERROR_IF (type not_eq EnergyLevelMapType::Tensor3_t,
                      "Must have Tensor3_t, input type is bad");
  
  auto elm = EnergyLevelMap(EnergyLevelMapType::Numeric_t, 1, 1, 1, *this);
  elm.value(joker, 0, 0, 0) = value(joker, ip, ilat, ilon);
  return elm;
}

EnergyLevelMapType toEnergyLevelMapTypeOrThrow(std::string_view s) {
  auto out = toEnergyLevelMapType(s);
  ARTS_USER_ERROR_IF(
      out == EnergyLevelMapType::Final_t,
      "Only \"None\", \"Numeric\", \"Vector\", and \"Tensor3\" types accepted\n"
      "You request to have an EnergyLevelMap of type: ",
      s,
      '\n')
  return out;
}

std::ostream& operator<<(std::ostream& os, EnergyLevelMapType x) {return os << toString(x);}
