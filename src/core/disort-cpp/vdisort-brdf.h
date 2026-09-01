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
  Complex refractive_index{1.34, 0.0};
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
  Complex refractive_index{1.5, 0.0};

  [[nodiscard]] rtepack::muelmat operator()(Numeric incident_mu) const;
};

/** Construct exact fully depolarizing Lambertian VDISORT modes.
 *
 * The surface maps incident I to reflected I and removes Q, U, and V.  Only
 * combined cosine mode zero is nonzero; additional zero modes are returned
 * for direct composition with azimuth-dependent models.
 */
[[nodiscard]] std::vector<BDRF> lambertian_fourier_modes(Numeric albedo, Index number_of_modes);

/** Project an azimuth-resolved polarized BPrDF into VDISORT combined modes. */
[[nodiscard]] std::vector<BDRF> fourier_modes(RawFunction raw,
                                              Index       number_of_modes,
                                              Index       azimuth_quadrature_points = 100);

/** Construct VDISORT-ready combined Fourier modes for a Cox-Munk ocean. */
[[nodiscard]] std::vector<BDRF> cox_munk_fourier_modes(Numeric wind_speed,
                                                       Complex refractive_index,
                                                       bool    shadowing,
                                                       Index   number_of_modes,
                                                       Index   azimuth_quadrature_points = 100);

/** Construct VDISORT modes for an ideal flat Fresnel surface.
 *
 * The callbacks recognize matching incoming and outgoing discrete ordinates.
 * Consequently this representation is intended for the native quadrature
 * solution, not direct-beam reflection at arbitrary user directions.
 */
[[nodiscard]] std::vector<BDRF> fresnel_fourier_modes(Complex refractive_index, Index number_of_modes);

/** Form a weighted sum of two already-projected polarized surface models.
 *
 * Diffuse and direct-beam callbacks are combined independently.  A missing
 * higher mode is zero.  Weights must be finite and nonnegative, but need not
 * sum to one.
 */
[[nodiscard]] std::vector<BDRF> combine_fourier_modes(std::vector<BDRF> first,
                                                      Numeric           first_weight,
                                                      std::vector<BDRF> second,
                                                      Numeric           second_weight);

/** Construct an ideal Fresnel/Lambertian mixture.
 *
 * `fresnel_fraction` weights the ideal specular operator and its complement
 * weights the depolarizing Lambertian operator.  The ideal component remains
 * unsuitable for an off-grid direct beam; use the Cox-Munk mixture below when
 * a finite-width, beam-compatible specular lobe is required.
 */
[[nodiscard]] std::vector<BDRF> fresnel_lambertian_fourier_modes(Numeric fresnel_fraction,
                                                                 Numeric lambertian_albedo,
                                                                 Complex refractive_index,
                                                                 Index   number_of_modes);

/** Construct a rough Fresnel Cox-Munk/depolarizing Lambertian mixture. */
[[nodiscard]] std::vector<BDRF> cox_munk_lambertian_fourier_modes(Numeric cox_munk_fraction,
                                                                  Numeric lambertian_albedo,
                                                                  Numeric wind_speed,
                                                                  Complex refractive_index,
                                                                  bool    shadowing,
                                                                  Index   number_of_modes,
                                                                  Index   azimuth_quadrature_points = 100);

}  // namespace vdisort::brdf
