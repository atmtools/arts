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

/* The surface reflection of a model, evaluated on given directions (see disort::surface_reflection):
 *   mu_out, the cosines of the upward directions it reflects into (the
 *     solver's upward streams, or user directions);
 *   mu_in, the cosines of the solver's downward streams it reflects from
 *     (their magnitudes are used);
 *   mu0, the cosine of the direct beam, or NaN without a beam, in which case
 *     the beam reflection is zero and marked as computed without a beam.
 */

/** The exact scalar Lambertian surface: R_0 = A, higher modes zero.
 *
 * Only mode zero is nonzero.  The additional zero modes let the result be
 * combined directly with azimuth-dependent surface models.
 */
[[nodiscard]] surface_reflection lambertian_fourier_modes(
    Numeric albedo, Index number_of_modes, const ConstVectorView& mu_out, const ConstVectorView& mu_in, Numeric mu0);

/** The cosine Fourier modes of an azimuth-resolved scalar BRDF, as disort::main_data takes them.
 *
 * Raw functions return a physical BRDF per steradian.  The coefficient
 * convention is DISORT's RHOQ quantity,
 *   R_m = (2 - delta_m0) * integral_0^pi rho(dphi) cos(m dphi) dphi.
 * Thus a physical Lambertian BRDF A/pi becomes R_0=A.  The raw function
 * receives the magnitudes of the cosines.
 */
[[nodiscard]] surface_reflection fourier_modes(const RawFunction&     raw,
                                               Index                  number_of_modes,
                                               const ConstVectorView& mu_out,
                                               const ConstVectorView& mu_in,
                                               Numeric                mu0,
                                               Index                  azimuth_quadrature_points = 100);

/** The weighted sum of two surface reflections on the same directions.
 *
 * A missing higher mode in either input is treated as zero.  Weights must be
 * finite and nonnegative, but need not sum to one; this permits both convex
 * mixtures and deliberate scaling of a surface operator.  Both must have been
 * evaluated for the same beam (or both without one).
 */
[[nodiscard]] surface_reflection combine_fourier_modes(const surface_reflection& first,
                                                       Numeric                   first_weight,
                                                       const surface_reflection& second,
                                                       Numeric                   second_weight);

/** A scalar Cox-Munk/Lambertian mixture.
 *
 * `cox_munk_fraction` weights the Cox-Munk operator and its complement
 * weights a Lambertian surface with the supplied albedo.
 */
[[nodiscard]] surface_reflection cox_munk_lambertian_fourier_modes(Numeric                cox_munk_fraction,
                                                                   Numeric                lambertian_albedo,
                                                                   Numeric                wind_speed,
                                                                   Complex                refractive_index,
                                                                   bool                   shadowing,
                                                                   Index                  number_of_modes,
                                                                   const ConstVectorView& mu_out,
                                                                   const ConstVectorView& mu_in,
                                                                   Numeric                mu0,
                                                                   Index azimuth_quadrature_points = 100);

}  // namespace disort::brdf
