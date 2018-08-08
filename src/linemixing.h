/* Copyright 2018, Richard Larsson
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/


// NOTE: This file contains only experimental code to test the F90-routines.
// It might evolve to replace it at some point but that is beyond present intent


#ifndef linemixing_h
#define linemixing_h

#include "rational.h"
#include "linerecord.h"
#include "absorption.h"


class AdiabaticFactor {
public:
  enum class Type {Hartmann};
  enum class HartmannPos : Index {dc, SIZE};
  
  AdiabaticFactor(const Vector& v, Type t) : mtype(t), mdata(v) {
    bool error = false;
    switch(mtype) {
      case Type::Hartmann:
        if(Index(HartmannPos::SIZE) not_eq mdata.nelem())
          error = true;
        break;
    }
    
    if(error)
      throw std::runtime_error("Bad initializaton of BasisRate, type and size disagree...");
  }
  
  
  Numeric mol_X(const int L,
                const int step,
                const Numeric& B0,
                const Numeric& T,
                const Numeric& main_mass,
                const Numeric& collider_mass,
                const bool init);
  
  Numeric get(const int L,
              const int step,
              const Numeric& B0,
              const Numeric& T,
              const Numeric& main_mass,
              const Numeric& collider_mass,
              const bool init) {
    switch(mtype) {
      case Type::Hartmann:
        return mol_X(L, step, B0, T, main_mass, collider_mass, init);
    }
    return 0;
  }
  
private:
  Type mtype;
  Numeric msave;
  Vector mdata;
};


class BasisRate {
public:
  enum class Type {Hartmann};
  enum class HartmannPos : Index {a1, a2, a3, SIZE};
  
  BasisRate(const Vector& v, Type t) : mtype(t), mdata(v) {
    bool error = false;
    switch(mtype) {
      case Type::Hartmann:
        if(Index(HartmannPos::SIZE) not_eq mdata.nelem())
          error = true;
        break;
    }
    
    if(error)
      throw std::runtime_error("Bad initializaton of BasisRate, type and size disagree...");
  }
  
  Numeric mol_X(const int L, const Numeric& B0, const Numeric& T) const;
  
  Numeric get(const int L, const Numeric& B0, const Numeric& T) const {
    switch(mtype) {
      case Type::Hartmann:
        return mol_X(L, B0, T);
    }
    return 0;
  }
  
private:
  Type mtype;
  Vector mdata;
};


namespace OffDiagonalElement {
  enum class Type {CO2_IR};
  
  tuple<Numeric, Numeric> CO2_IR(const LineRecord& j_line,
                                 const LineRecord& k_line,
                                 const Numeric& j_rho,
                                 const Numeric& k_rho,
                                 const BasisRate& br,
                                 AdiabaticFactor& af,
                                 const Numeric& T,
                                 const Numeric& B0,
                                 const Numeric& main_mass,
                                 const Numeric& collider_mass);
};


Matrix hartmann_ecs_interface(const ArrayOfLineRecord& abs_lines,
                              const SpeciesTag& main_species,
                              const ArrayOfSpeciesTag& collider_species,
                              const Vector& collider_species_vmr,
                              const SpeciesAuxData::AuxType& partition_type,
                              const ArrayOfGriddedField1& partition_data,
                              const Numeric& T,
                              const Index& size,
                              const Index type);

#endif // linemixing_h
