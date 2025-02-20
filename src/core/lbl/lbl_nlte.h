#pragma once

#include <atm.h>

#include "lbl_data.h"

using QuantumIdentifierGriddedField1Map =
    std::unordered_map<QuantumIdentifier, GriddedField1>;
using QuantumIdentifierVectorMap =
    std::unordered_map<QuantumIdentifier, Vector>;
using QuantumIdentifierNumericMap =
    std::unordered_map<QuantumIdentifier, Numeric>;

namespace lbl::nlte {
std::unordered_map<QuantumIdentifier, AtmData> from_lte(
    const AtmField& atmospheric_field,
    const AbsorptionBands& absorption_bands,
    const Numeric& normalizing_factor);

QuantumIdentifierNumericMap createAij(const AbsorptionBands& absorption_bands);

QuantumIdentifierNumericMap createBij(const AbsorptionBands& absorption_bands);

QuantumIdentifierNumericMap createBji(const QuantumIdentifierNumericMap& Bij,
                                      const AbsorptionBands& absorption_bands);

QuantumIdentifierVectorMap createCij(
    const AbsorptionBands& absorption_bands,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfAtmPoint& ray_path_atmospheric_point);

QuantumIdentifierVectorMap createCji(
    const QuantumIdentifierVectorMap& Cij,
    const AbsorptionBands& absorption_bands,
    const ArrayOfAtmPoint& ray_path_atmospheric_point);

Vector nlte_ratio_sum(const ArrayOfAtmPoint& ray_path_atmospheric_point,
                      const ArrayOfQuantumIdentifier& levels);

struct UppLow {
  Size upp, low;
};

Size band_level_mapUniquestIndex(
    const std::unordered_map<QuantumIdentifier, UppLow>& band_level_map);

std::unordered_map<QuantumIdentifier, UppLow> band_level_mapFromLevelKeys(
    const AbsorptionBands& absorption_bands,
    const ArrayOfQuantumIdentifier& level_keys);

Size level_count(
    const std::unordered_map<QuantumIdentifier, UppLow>& band_level_map);

Matrix statistical_equilibrium_equation(
    const QuantumIdentifierNumericMap& Aij,
    const QuantumIdentifierNumericMap& Bij,
    const QuantumIdentifierNumericMap& Bji,
    const QuantumIdentifierVectorMap& Cij,
    const QuantumIdentifierVectorMap& Cji,
    const QuantumIdentifierVectorMap& Jij,
    const std::unordered_map<QuantumIdentifier, UppLow>& level_map,
    const Size atmi,
    const Size nlevels);

Numeric set_nlte(AtmPoint& atmospheric_point,
                 const ArrayOfQuantumIdentifier& level_keys,
                 const Vector& x);
}  // namespace lbl::nlte
