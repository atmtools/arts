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

#ifndef energylevelmap_h
#define energylevelmap_h

#include "absorptionlines.h"
#include "matpackIV.h"
#include "mystring.h"
#include "quantum_numbers.h"

// Output from EnergyLevelMap
struct Output2{
  Numeric r_low;
  Numeric r_upp;
};

// Output from EnergyLevelMap
struct Output4{
  Numeric E_low;
  Numeric E_upp;
  Numeric T_low;
  Numeric T_upp;
};

enum class EnergyLevelMapType {
  Tensor3_t,
  Vector_t,
  Numeric_t,
  None_t,
  Final_t
};

constexpr EnergyLevelMapType toEnergyLevelMapType(std::string_view s) noexcept {
  if (s == "Tensor3")
    return EnergyLevelMapType::Tensor3_t;
  if (s == "Vector")
    return EnergyLevelMapType::Vector_t;
  if (s == "Numeric")
    return EnergyLevelMapType::Numeric_t;
  if (s == "None")
    return EnergyLevelMapType::None_t;
  return EnergyLevelMapType::Final_t;
}

EnergyLevelMapType toEnergyLevelMapTypeOrThrow(std::string_view s);

constexpr std::string_view toString(EnergyLevelMapType x) noexcept {
  switch(x) {
    case EnergyLevelMapType::Tensor3_t:
      return "Tensor3";
    case EnergyLevelMapType::Vector_t:
      return "Vector";
    case EnergyLevelMapType::Numeric_t:
      return "Numeric";
    case EnergyLevelMapType::None_t:
      return "None";
    case EnergyLevelMapType::Final_t:
      {/* leave last */}
  }
  return "BAD EnergyLevelMapType";
}

inline std::ostream& operator<<(std::ostream& os, EnergyLevelMapType x) {return os << toString(x);}

struct EnergyLevelMap {
  EnergyLevelMapType type{EnergyLevelMapType::None_t};
  ArrayOfQuantumIdentifier levels;
  Vector vib_energy;
  Tensor4 value;
  
  [[nodiscard]] bool OK() const noexcept;
  
  void ThrowIfNotOK() const ARTS_NOEXCEPT {ARTS_ASSERT (OK(), "Class in bad state");}
  
  EnergyLevelMap() : levels(0),
  vib_energy(0), value(0, 0, 0, 0) {}
  
  EnergyLevelMap(EnergyLevelMapType new_type, Index pages, Index rows,
                 Index cols, const EnergyLevelMap& old) : 
  type(new_type), levels(old.levels), vib_energy(old.vib_energy),
  value(old.levels.nelem(), pages, rows, cols) {ThrowIfNotOK();};
  
  // Create Tensor3_t from the raw inputs
  EnergyLevelMap(Tensor4 data, ArrayOfQuantumIdentifier levels, Vector energies=Vector(0));
  
  // Create Vector_t from the raw inputs
  EnergyLevelMap(const Matrix& data, ArrayOfQuantumIdentifier levels, Vector energies=Vector(0));
  
  // Create Numeric_t from the raw inputs
  EnergyLevelMap(const Vector& data, ArrayOfQuantumIdentifier levels, Vector energies=Vector(0));
  
  // Create Vector_t from Tensor3_t
  [[nodiscard]] EnergyLevelMap InterpToGridPos(Index atmosphere_dim, const ArrayOfGridPos& p, const ArrayOfGridPos& lat, const ArrayOfGridPos& lon) const;
  
  // Create Numeric_t from Vector_t
  [[nodiscard]] EnergyLevelMap operator[](Index ip) const;
  
  // Create Numeric_t from Tensor3_t
  [[nodiscard]] EnergyLevelMap operator()(Index ip, Index ilat, Index ilon) const;
  
  //////////////////////
  // Numeric_t access //
  //////////////////////
  
  /** Get the output required for Population::NLTE
   * 
   * @param[in] transition A line-by-line transition
   * @return Upper and lower level distributions
   */
  [[nodiscard]] Output2 get_ratio_params(const AbsorptionLines& band, const Index& line_index) const;
  
  /** Get the output required for Population::NLTE-VibrationalTemperatures
   * 
   * @param[in] transition A line-by-line transition
   * @return Upper and lower level distributions and energies
   */
  [[nodiscard]] Output4 get_vibtemp_params(const AbsorptionLines& band, const Numeric T) const;
};

std::ostream& operator<<(std::ostream& os, const EnergyLevelMap& elm);

#endif  // energylevelmap_h
