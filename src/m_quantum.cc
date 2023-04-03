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
