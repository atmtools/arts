#pragma once

#include <atm.h>

#include "lbl_data.h"

namespace lbl::nlte {
std::unordered_map<QuantumIdentifier, AtmData> from_lte(
    const AtmField& atmospheric_field,
    const AbsorptionBands& absorption_bands,
    const Numeric& normalizing_factor);
}  // namespace lbl::nlte

//! Match collision rate to quantum identifier - as a function of temperature
using NonlteCollision = std::unordered_map<QuantumIdentifier, GriddedField1>;
