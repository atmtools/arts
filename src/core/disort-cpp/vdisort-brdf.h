#pragma once

#include <vdisort.h>

namespace vdisort::brdf {

using RawFunction = std::function<rtepack::muelmat(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth)>;

/** Polarized one-dimensional Gaussian Cox-Munk ocean reflection.
 *
 * Implements Lin et al. (2022), Eqs. (147)-(151), in the ARTS [I,Q,U,V]
 * meridian-plane convention.  The returned matrix is a physical BPrDF per
 * steradian before Fourier projection.
 */
struct CoxMunk {
  Numeric wind_speed{5.0};
  Numeric refractive_index{1.34};
  bool    shadowing{true};

  [[nodiscard]] rtepack::muelmat operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

/** Ideal polarized Fresnel reflection from a flat dielectric interface.
 *
 * Unlike Cox-Munk, this is a specular delta distribution.  Its Fourier-mode
 * factory therefore constructs the quadrature-normalized discrete boundary
 * operator directly rather than passing through fourier_modes().
 */
struct Fresnel {
  Numeric refractive_index{1.5};

  [[nodiscard]] rtepack::muelmat operator()(Numeric incident_mu) const;
};

/** Project an azimuth-resolved polarized BPrDF into VDISORT combined modes. */
[[nodiscard]] std::vector<BDRF> fourier_modes(RawFunction raw,
                                              Index       number_of_modes,
                                              Index       azimuth_quadrature_points = 100);

/** Construct VDISORT-ready combined Fourier modes for a Cox-Munk ocean. */
[[nodiscard]] std::vector<BDRF> cox_munk_fourier_modes(Numeric wind_speed,
                                                       Numeric refractive_index,
                                                       bool    shadowing,
                                                       Index   number_of_modes,
                                                       Index   azimuth_quadrature_points = 100);

/** Construct VDISORT modes for an ideal flat Fresnel surface.
 *
 * The callbacks recognize matching incoming and outgoing discrete ordinates.
 * Consequently this representation is intended for the native quadrature
 * solution, not direct-beam reflection at arbitrary user directions.
 */
[[nodiscard]] std::vector<BDRF> fresnel_fourier_modes(Numeric refractive_index, Index number_of_modes);

}  // namespace vdisort::brdf
