#pragma once

#include <matpack.h>

#include <optional>
#include "rtepack.h"

namespace surface_scattering {

/// BRDF Mueller matrix over (f_grid, za_inc, aa_inc, za_scat, aa_scat, 4, 4)
using BRDFMatrix = MuelmatTensor5;

/// Emissivity vector over (f_grid, za_scat, 4)
using EmissivityVector = StokvecTensor3;

/** Bulk surface scattering properties accumulated across all surface models.
 *
 * Holds an optional BRDF Mueller matrix and an emissivity vector.
 * Multiple models can be accumulated via operator+=.
 */
struct SurfaceScatteringModelProperties {
  /// Optional BRDF matrix: dims [n_f, n_za_inc, n_aa_inc, n_za_scat, n_aa_scat, 4, 4]
  std::optional<BRDFMatrix> brdf_matrix;
  /// Emissivity vector: dims [n_f, n_za_scat, 4]
  EmissivityVector emissivity_vector;

  SurfaceScatteringModelProperties& operator+=(
      const SurfaceScatteringModelProperties& other);
};

}  // namespace surface_scattering

