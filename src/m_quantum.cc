/* Copyright (C) 2019 Richard Larsson
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/**
 * @file   m_quantum.cc
 * @author Richard Larsson <larsson (at) mps.mpg.de>
 * @date   2019-10-22
 * 
 * @brief Quantum handling functions
 */

#include "energylevelmap.h"
#include "quantum_numbers.h"
#include "messages.h"

void nlte_fieldFromRaw(EnergyLevelMap& nlte_field,
                       const ArrayOfQuantumIdentifier& nlte_level_identifiers,
                       const Vector& nlte_vibrational_energies,
                       const Tensor4& data,
                       const Verbosity&)
{
  nlte_field = EnergyLevelMap(data, nlte_level_identifiers, nlte_vibrational_energies);
}

void abs_nlteFromRaw(EnergyLevelMap& abs_nlte,
                     const ArrayOfQuantumIdentifier& nlte_level_identifiers,
                     const Vector& nlte_vibrational_energies,
                     const Matrix& data,
                     const Verbosity&)
{
  abs_nlte = EnergyLevelMap(data, nlte_level_identifiers, nlte_vibrational_energies);
}

void rtp_nlteFromRaw(EnergyLevelMap& rtp_nlte,
                     const ArrayOfQuantumIdentifier& nlte_level_identifiers,
                     const Vector& nlte_vibrational_energies,
                     const Vector& data,
                     const Verbosity&)
{
  rtp_nlte = EnergyLevelMap(data, nlte_level_identifiers, nlte_vibrational_energies);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void EnergyLevelMapSet(EnergyLevelMap& x, const EnergyLevelMap& y, const Verbosity&) {
  x = y;
}
