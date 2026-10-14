#pragma once

#include <vdisort.h>

namespace vdisort::brdf {

using RawFunction = std::function<rtepack::muelmat(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth)>;
using ScalarRawFunction = std::function<Numeric(Numeric outgoing_mu, Numeric incoming_mu, Numeric relative_azimuth)>;

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

/* The surface reflection of a model, evaluated on given directions (see vdisort::surface_reflection):
 *   mu_out, the cosines of the upward directions it reflects into (the
 *     solver's upward streams, or user directions);
 *   mu_in, the cosines of the solver's downward streams it reflects from
 *     (their magnitudes are used);
 *   mu0, the cosine of the direct beam, or NaN without a beam, in which case
 *     the beam reflection is zero and marked as computed without a beam.
 */

/** The exact, fully depolarizing Lambertian surface.
 *
 * The surface maps incident I to reflected I and removes Q, U, and V.  Only
 * combined cosine mode zero is nonzero; additional zero modes are returned
 * for direct composition with azimuth-dependent models.
 */
[[nodiscard]] surface_reflection lambertian_fourier_modes(Numeric                albedo,
                                                          Index                  number_of_modes,
                                                          const ConstVectorView& mu_out,
                                                          const ConstVectorView& mu_in,
                                                          Numeric                mu0);

/** The combined Fourier modes of an azimuth-resolved polarized BPrDF. */
[[nodiscard]] surface_reflection fourier_modes(const RawFunction&     raw,
                                               Index                  number_of_modes,
                                               const ConstVectorView& mu_out,
                                               const ConstVectorView& mu_in,
                                               Numeric                mu0,
                                               Index                  azimuth_quadrature_points = 100);

/** A scalar BRDF as a fully depolarizing surface.
 *
 * The scalar value is placed in M00 and every other Mueller element is zero.
 * Fourier projection then follows the same polarized normalization and
 * combined cosine/sine layout as any other VDISORT surface.
 */
[[nodiscard]] surface_reflection depolarizing_fourier_modes(const ScalarRawFunction& raw,
                                                            Index                    number_of_modes,
                                                            const ConstVectorView&   mu_out,
                                                            const ConstVectorView&   mu_in,
                                                            Numeric                  mu0,
                                                            Index                    azimuth_quadrature_points = 100);

/** A fully depolarizing Hapke surface from the shared scalar model. */
[[nodiscard]] surface_reflection hapke_fourier_modes(Numeric                opposition_amplitude,
                                                     Numeric                opposition_width,
                                                     Numeric                single_scattering_albedo,
                                                     Index                  number_of_modes,
                                                     const ConstVectorView& mu_out,
                                                     const ConstVectorView& mu_in,
                                                     Numeric                mu0,
                                                     Index                  azimuth_quadrature_points = 100);

/** A fully depolarizing RPV surface from the shared scalar model. */
[[nodiscard]] surface_reflection rpv_fourier_modes(Numeric                rho0,
                                                   Numeric                kappa,
                                                   Numeric                asymmetry,
                                                   Numeric                hotspot,
                                                   Index                  number_of_modes,
                                                   const ConstVectorView& mu_out,
                                                   const ConstVectorView& mu_in,
                                                   Numeric                mu0,
                                                   Index                  azimuth_quadrature_points = 100);

/** A fully depolarizing Ross-Li surface from the shared scalar model. */
[[nodiscard]] surface_reflection ross_li_fourier_modes(Numeric                isotropic,
                                                       Numeric                volumetric,
                                                       Numeric                geometric,
                                                       Numeric                hotspot_angle,
                                                       Index                  number_of_modes,
                                                       const ConstVectorView& mu_out,
                                                       const ConstVectorView& mu_in,
                                                       Numeric                mu0,
                                                       Index                  azimuth_quadrature_points = 100);

/** A polarized Cox-Munk ocean. */
[[nodiscard]] surface_reflection cox_munk_fourier_modes(Numeric                wind_speed,
                                                        Complex                refractive_index,
                                                        bool                   shadowing,
                                                        Index                  number_of_modes,
                                                        const ConstVectorView& mu_out,
                                                        const ConstVectorView& mu_in,
                                                        Numeric                mu0,
                                                        Index                  azimuth_quadrature_points = 100);

/** An ideal flat Fresnel surface.
 *
 * The surface reflects only specularly: the modes have no reflection kernel,
 * and the specular part is the Fresnel matrix R(mu_out), which VDISORT
 * applies on the streams and at user directions.  A direct beam (mu0 not
 * NaN) is refused: its mirror image is an upward beam.
 */
[[nodiscard]] surface_reflection fresnel_fourier_modes(Complex                refractive_index,
                                                       Index                  number_of_modes,
                                                       const ConstVectorView& mu_out,
                                                       const ConstVectorView& mu_in,
                                                       Numeric                mu0);

/** The weighted sum of two surface reflections on the same directions.
 *
 * Kernels, beam reflection and specular parts are combined independently.  A
 * missing higher mode is zero.  Weights must be finite and nonnegative, but
 * need not sum to one.  Both must have been evaluated for the same beam (or
 * both without one).
 */
[[nodiscard]] surface_reflection combine_fourier_modes(const surface_reflection& first,
                                                       Numeric                   first_weight,
                                                       const surface_reflection& second,
                                                       Numeric                   second_weight);

/** An ideal Fresnel/Lambertian mixture.
 *
 * `fresnel_fraction` weights the ideal specular operator and its complement
 * weights the depolarizing Lambertian operator.  The ideal component remains
 * unsuitable for a direct beam; use the Cox-Munk mixture below when a
 * finite-width, beam-compatible specular lobe is required.
 */
[[nodiscard]] surface_reflection fresnel_lambertian_fourier_modes(Numeric                fresnel_fraction,
                                                                  Numeric                lambertian_albedo,
                                                                  Complex                refractive_index,
                                                                  Index                  number_of_modes,
                                                                  const ConstVectorView& mu_out,
                                                                  const ConstVectorView& mu_in,
                                                                  Numeric                mu0);

/** A rough Fresnel Cox-Munk/depolarizing Lambertian mixture. */
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

}  // namespace vdisort::brdf
