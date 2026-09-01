#pragma once

#include <arts_constants.h>
#include <disort.h>

#include <functional>

namespace disort::brdf {

using RawFunction = std::function<Numeric(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth)>;

struct Hapke {
  Numeric opposition_amplitude{1.0};
  Numeric opposition_width{0.06};
  Numeric single_scattering_albedo{0.6};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

struct CoxMunk {
  Numeric wind_speed{12.0};
  Complex refractive_index{1.34, 0.0};
  bool    shadowing{false};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

struct RPV {
  Numeric rho0{0.027};
  Numeric kappa{0.647};
  Numeric asymmetry{-0.169};
  Numeric hotspot{0.1};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

struct RossLi {
  Numeric isotropic{0.091};
  Numeric volumetric{0.02};
  Numeric geometric{0.01};
  Numeric hotspot_angle{1.5 * Constant::pi / 180.0};

  Numeric operator()(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth) const;
};

/** Construct exact scalar Lambertian Fourier modes.
 *
 * Only mode zero is nonzero.  Additional zero modes are returned so the
 * result can be combined directly with azimuth-dependent surface models.
 */
[[nodiscard]] std::vector<BDRF> lambertian_fourier_modes(Numeric albedo, Index number_of_modes);

/** Convert an azimuth-resolved scalar BRDF into the cosine Fourier modes
 * expected by disort::main_data.
 *
 * Raw functions return a physical BRDF per steradian.  The callback
 * coefficient convention is DISORT's RHOQ quantity,
 *   R_m = (2 - delta_m0) * integral_0^pi rho(dphi) cos(m dphi) dphi.
 * Thus a physical Lambertian BRDF A/pi becomes R_0=A.  incoming_mu may be
 * signed when the resulting callbacks are evaluated; raw functions always
 * receive its positive magnitude.
 */
std::vector<BDRF> fourier_modes(RawFunction raw, Index number_of_modes, Index azimuth_quadrature_points = 100);

/** Form a weighted sum of two already-projected scalar surface models.
 *
 * A missing higher mode in either input is treated as zero.  Weights must be
 * finite and nonnegative, but need not sum to one; this permits both convex
 * mixtures and deliberate scaling of a surface operator.
 */
[[nodiscard]] std::vector<BDRF> combine_fourier_modes(std::vector<BDRF> first,
                                                      Numeric           first_weight,
                                                      std::vector<BDRF> second,
                                                      Numeric           second_weight);

/** Construct a scalar Cox-Munk/Lambertian mixture.
 *
 * `cox_munk_fraction` weights the Cox-Munk operator and its complement
 * weights a Lambertian surface with the supplied albedo.
 */
[[nodiscard]] std::vector<BDRF> cox_munk_lambertian_fourier_modes(Numeric cox_munk_fraction,
                                                                  Numeric lambertian_albedo,
                                                                  Numeric wind_speed,
                                                                  Complex refractive_index,
                                                                  bool    shadowing,
                                                                  Index   number_of_modes,
                                                                  Index   azimuth_quadrature_points = 100);

}  // namespace disort::brdf
