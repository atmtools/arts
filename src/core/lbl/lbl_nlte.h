#pragma once

#include <atm.h>

#include "lbl_data.h"

using QuantumIdentifierGriddedField1Map =
    std::unordered_map<QuantumIdentifier, SortedGriddedField1>;
using QuantumIdentifierVectorMap =
    std::unordered_map<QuantumIdentifier, Vector>;
using QuantumIdentifierNumericMap =
    std::unordered_map<QuantumIdentifier, Numeric>;

namespace lbl::nlte {
std::unordered_map<QuantumLevelIdentifier, AtmData> from_lte(
    const AtmField& atm_field,
    const AbsorptionBands& abs_bands,
    const Numeric& normalizing_factor);

QuantumIdentifierNumericMap createAij(const AbsorptionBands& abs_bands);

QuantumIdentifierNumericMap createBij(const AbsorptionBands& abs_bands);

QuantumIdentifierNumericMap createBji(const QuantumIdentifierNumericMap& Bij,
                                      const AbsorptionBands& abs_bands);

QuantumIdentifierVectorMap createCij(
    const AbsorptionBands& abs_bands,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfAtmPoint& atm_path);

QuantumIdentifierVectorMap createCji(const QuantumIdentifierVectorMap& Cij,
                                     const AbsorptionBands& abs_bands,
                                     const ArrayOfAtmPoint& atm_path);

Vector nlte_ratio_sum(const ArrayOfAtmPoint& atm_path,
                      const ArrayOfQuantumLevelIdentifier& levels);

struct UppLow {
  Size upp, low;
};

Size band_level_mapUniquestIndex(
    const std::unordered_map<QuantumIdentifier, UppLow>& band_level_map);

std::unordered_map<QuantumIdentifier, UppLow> band_level_mapFromLevelKeys(
    const AbsorptionBands& abs_bands,
    const ArrayOfQuantumLevelIdentifier& level_keys);

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

Numeric set_nlte(AtmPoint& atm_point,
                 const ArrayOfQuantumLevelIdentifier& level_keys,
                 const Vector& x);
}  // namespace lbl::nlte
