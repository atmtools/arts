#pragma once

#include <atm.h>

#include "lbl_data.h"

namespace lbl::nlte {
std::unordered_map<QuantumIdentifier, AtmData> from_lte(
    const AtmField& atmospheric_field,
    const AbsorptionBands& absorption_bands,
    const Numeric& normalizing_factor);
}  // namespace lbl::nlte

using QuantumIdentifierGriddedField1Map = std::unordered_map<QuantumIdentifier, GriddedField1>;
using QuantumIdentifierVectorMap = std::unordered_map<QuantumIdentifier, Vector>;
using QuantumIdentifierNumericMap = std::unordered_map<QuantumIdentifier, Numeric>;
