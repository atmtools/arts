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

#ifndef energylevelmap_h
#define energylevelmap_h

#include "absorptionlines.h"
#include "matpackIV.h"
#include "mystring.h"
#include "quantum.h"

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
};

EnergyLevelMapType string2energylevelmaptype(const String& s);

String energylevelmaptype2string(EnergyLevelMapType type);

class EnergyLevelMap {
private:
  EnergyLevelMapType mtype;
  ArrayOfQuantumIdentifier mlevels;
  Vector mvib_energy;
  Tensor4 mvalue;
  
public:
  bool OK() const {
    if (not (mvalue.nbooks() == mlevels.nelem() and 
      (mvib_energy.nelem() == mlevels.nelem() or mvib_energy.nelem() == 0)))
      return false;  // Bad dimensions, vibrational energies and IDs and data of strange size
    
    if (mtype == EnergyLevelMapType::Tensor3_t) {
    } else if (mtype == EnergyLevelMapType::Vector_t) {
      if (mvalue.npages() not_eq 1 or mvalue.nrows() not_eq 1)
        return false;  // Bad dimensions for vector type
    } else if (mtype == EnergyLevelMapType::Numeric_t) {
      if (mvalue.npages() not_eq 1 or mvalue.nrows() not_eq 1 or mvalue.ncols() not_eq 1)
        return false;  // Bad dimensions for numeric type
    } else if (mtype == EnergyLevelMapType::None_t) {
      if (mvalue.npages() not_eq 0 or mvalue.nrows() not_eq 0 or mvalue.ncols() not_eq 0)
        return false;  // Bad dimensions for none type
    }
    
    // FIXME: matpack does not do pointers in a clear manner
    // return not std::any_of(mvib_energy.begin(), mvib_energy.end(), [](auto& x){return x < 0;});
    
    for (auto& e: mvib_energy) if (e < 0) return false;  // Bad energies
    
    return true;
  }
  
  void ThrowIfNotOK() const ARTS_NOEXCEPT {ARTS_ASSERT (not OK(), "Class in bad state");}
  
  EnergyLevelMap() : mtype(EnergyLevelMapType::None_t), mlevels(0),
  mvib_energy(0), mvalue(0, 0, 0, 0) {ThrowIfNotOK();}
  
  EnergyLevelMap(EnergyLevelMapType new_type, Index pages, Index rows,
                 Index cols, const EnergyLevelMap& old) : 
  mtype(new_type), mlevels(old.mlevels), mvib_energy(old.mvib_energy),
  mvalue(old.mlevels.nelem(), pages, rows, cols) {ThrowIfNotOK();};
  
  // Create Tensor3_t from the raw inputs
  EnergyLevelMap(const Tensor4& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies=Vector(0));
  
  // Create Vector_t from the raw inputs
  EnergyLevelMap(const Matrix& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies=Vector(0));
  
  // Create Numeric_t from the raw inputs
  EnergyLevelMap(const Vector& data, const ArrayOfQuantumIdentifier& levels, const Vector& energies=Vector(0));
  
  // Create Vector_t from Tensor3_t
  EnergyLevelMap InterpToGridPos(Index atmosphere_dim, const ArrayOfGridPos& p, const ArrayOfGridPos& lat, const ArrayOfGridPos& lon) const;
  
  // Create Numeric_t from Vector_t
  EnergyLevelMap operator[](Index ip) const;
  
  // Create Numeric_t from Tensor3_t
  EnergyLevelMap operator()(Index ip, Index ilat, Index ilon) const;
  
  /** Energy level type */
  EnergyLevelMapType Type() const noexcept {return mtype;}
  
  /** Energy level type */
  const ArrayOfQuantumIdentifier& Levels() const noexcept {return mlevels;}
  
  /** Energy level type */
  const Vector& Energies() const noexcept {return mvib_energy;}
  
  /** Energy level type */
  const Tensor4& Data() const noexcept {return mvalue;}
  
  /** Energy level type */
  EnergyLevelMapType& Type() noexcept {return mtype;}
  
  /** Energy level type */
  ArrayOfQuantumIdentifier& Levels() noexcept {return mlevels;}
  
  /** Energy level type */
  Vector& Energies() noexcept {return mvib_energy;}
  
  /** Energy level type */
  Tensor4& Data() noexcept {return mvalue;}
  
  ////////////////////////////
  // C API interface access //
  ////////////////////////////
  
  static bool validIndexForType(Index x) noexcept
  {
    constexpr auto keys = stdarrayify(Index(EnergyLevelMapType::Tensor3_t), EnergyLevelMapType::Vector_t, EnergyLevelMapType::Numeric_t, EnergyLevelMapType::None_t);
    return std::any_of(keys.cbegin(), keys.cend(), [x](auto y){return x == y;});
  }
  
  /** Energy level type */
  void Type(EnergyLevelMapType x) noexcept {mtype = x;}
  
  static EnergyLevelMapType string2Type(const String& s) noexcept
  {
    if (s == "Tensor3")
      return EnergyLevelMapType::Tensor3_t;
    else if (s == "Vector")
      return EnergyLevelMapType::Vector_t;
    else if (s == "Numeric")
      return EnergyLevelMapType::Numeric_t;
    else if (s == "None")
      return EnergyLevelMapType::None_t;
    else {
      return EnergyLevelMapType(-1);
    }
  }
  
  //////////////////////
  // Numeric_t access //
  //////////////////////
  
  /** Get the output required for Population::NLTE
   * 
   * @param[in] transition A line-by-line transition
   * @return Upper and lower level distributions
   */
  Output2 get_ratio_params(const AbsorptionLines& band, const Index& line_index) const;
  
  /** Get the output required for Population::NLTE-VibrationalTemperatures
   * 
   * @param[in] transition A line-by-line transition
   * @return Upper and lower level distributions and energies
   */
  Output4 get_vibtemp_params(const AbsorptionLines& band, const Index& line_index, const Numeric T) const;
};

std::ostream& operator<<(std::ostream& os, const EnergyLevelMap& elm);

#endif  // energylevelmap_h
