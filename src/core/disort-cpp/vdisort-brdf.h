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

}  // namespace vdisort::brdf
