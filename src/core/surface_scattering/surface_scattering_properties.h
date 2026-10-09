#pragma once

#include <matpack.h>

#include "rtepack.h"

namespace surface_scattering {

/// BRDF Mueller matrix over (f_grid, za_inc, aa_inc, za_scat, aa_scat, 4, 4)
using BRDFMatrix = MuelmatTensor5;

/// Emissivity Mueller matrix over (f_grid, za_scat, aa_scat, 4, 4)
using EmissivityVector = MuelmatTensor3;

/** Bulk surface scattering properties accumulated across all surface models.
 *
 * Holds the diffuse and specular BRDF matrices and the corresponding
 * emissivity tensors.  Multiple models can be accumulated via operator+=
 * and weighted via operator*=.
 */
struct SurfaceScatteringModelProperties {
  /// Diffuse BRDF matrix: dims [n_f, n_za_inc, n_aa_inc, n_za_scat, n_aa_scat, 4, 4]
  BRDFMatrix brdf_matrix_diffuse;
  /// Diffuse emissivity vector: dims [n_f, n_za_scat, n_aa_scat, 4, 4]
  EmissivityVector emissivity_vector_diffuse;
  /// Specular BRDF matrix: dims [n_f, n_za_inc, n_aa_inc, n_za_scat, n_aa_scat, 4, 4]
  BRDFMatrix brdf_matrix_specular;
  /// Specular emissivity vector: dims [n_f, n_za_scat, n_aa_scat, 4, 4]
  EmissivityVector emissivity_vector_specular;

  SurfaceScatteringModelProperties& operator+=(
      const SurfaceScatteringModelProperties& other);

  SurfaceScatteringModelProperties& operator*=(Numeric scalar);    

};

}  // namespace surface_scattering

