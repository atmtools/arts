#pragma once

#include <atm.h>
#include <matpack.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>

#include <variant>

#include "vdisort.h"

/**
 * VDISORT inputs from ARTS-native data.
 *
 * These helpers translate ARTS scattering species, atmospheric points and
 * propagation paths into VDISORT's phase-matrix Fourier coefficients and a
 * solved vdisort::main_data.  VDISORT has a scalar extinction, so only
 * species with totally randomly oriented (TRO) data
 * (ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_gridded) can
 * be used.
 *
 * Geometry and Stokes basis (those of VDISORT's polarized tests and of RT3,
 * see rt3.h): z points up; a ray propagating at the direction cosine mu
 * (> 0 upward) and azimuth phi has k = (sin cos phi, sin sin phi, mu); the
 * Stokes vector is [I, Q, U, V] in the meridional basis h = k x z / |k x z|,
 * v = h x k, with Q = I_v - I_h and U = 2 Re(E_v E_h*).  VDISORT's direct
 * beam propagates downward toward the azimuth phi0, and its radiance at
 * phi0 + psi is RT3's at psi.
 *
 * Relation to ARTS.  ARTS's scattering data take the propagation directions
 * (za, aa), with aa clockwise seen from above.  ARTS's laboratory-frame
 * phase matrix equals the vector-geometry phase matrix built here from the
 * same scattering matrix F with mu = cos(za) and phi = -aa, element by
 * element, with the same I, Q, U and V (tested in
 * cpp.fast.vdisort-arts-test).  So F is used as ARTS stores it,
 * [F11, F12, F22, F33, F34, F44] in the scattering-plane basis with
 * Q = I_par - I_perp, and VDISORT transports ARTS's Stokes vector.
 *
 * Units: extinction and scattering per metre, heights in metres, radiances
 * in W m-2 Hz-1 sr-1 (ARTS's planck()).
 */
namespace vdisort {
/** The normalised phase-matrix Fourier coefficients of scattering species. */
struct fourier_optics {
  //! K11 of the particles per metre
  Numeric extinction{};
  //! K11 - a1, the extinction minus the absorption, per metre
  Numeric scattering{};
  //! [nfourier, mu_out.size(), mu_in.size()] ordinary cosine coefficients C^m
  rtepack::muelmat_tensor3 cosine{};
  //! [nfourier, mu_out.size(), mu_in.size()] ordinary sine coefficients S^m
  rtepack::muelmat_tensor3 sine{};
};

/** The phase-matrix Fourier coefficients of scattering species at one atmospheric point.
 *
 * mu_out, mu_in: signed direction cosines (> 0 upward), each in [-1, 0) or
 *   (0, 1], e.g. VDISORT's streams (upward first) and, for the direct beam,
 *   mu_in = {-mu0}.
 *
 * Returns the ordinary Fourier coefficients, without epsilon_m,
 *   C^m(o, i) = (1 / 2 pi) int P(mu_o, 0; mu_i, phi) cos(m phi) dphi,
 *   S^m(o, i) = (1 / 2 pi) int P(mu_o, 0; mu_i, phi) sin(m phi) dphi,
 * m = 0 .. nfourier - 1, of the laboratory-frame phase matrix
 * P = 4 pi Z / sigma in the basis above, Z = L_out^T F(Theta) L_in from
 * vector geometry, normalised to 1 over 4 pi:
 *   sigma = 2 pi int F11 dcos(Theta)
 * by an n-point Gauss-Legendre rule, n = scattering_angle_count.  F is the
 * species' TRO scattering matrix at the exact scattering angle of every pair
 * of directions.  The integral over phi is the periodic midpoint rule at
 * phi_k = (k + 1/2) 2 pi / N, N = azimuth_count.  For a regular Legendre
 * series of degree L, Z is a trigonometric polynomial of degree L in phi and
 * the rule is exact for N > L + nfourier - 1; otherwise it converges as fast
 * as the Fourier series of Z in phi.  These coefficients go to
 * vdisort::combine_phase_matrices (diffuse) and
 * vdisort::combine_beam_phase_matrices (beam column).
 *
 * sigma and the scattering coefficient K11 - a1 must agree to
 * normalisation_tolerance * K11 (as in rt3::scattering_optics), otherwise
 * this is an error.  Without particles the coefficients are zero.
 */
fourier_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                 const AtmPoint&                 atm_point,
                                 Numeric                         frequency,
                                 const Vector&                   mu_out,
                                 const Vector&                   mu_in,
                                 Index                           nfourier,
                                 Index                           azimuth_count,
                                 Index                           scattering_angle_count,
                                 Numeric                         normalisation_tolerance);

//! A depolarizing Lambertian surface: vdisort::brdf::lambertian_fourier_modes, emission [(1 - A) B, 0, 0, 0]
struct lambertian_surface {
  Numeric albedo{0.0};
};

/** A flat Fresnel surface under a medium of index 1: vdisort::brdf::fresnel_fourier_modes, emission
 *  B ([1, 0, 0, 0] - R(mu)[:, 0]) with R the Fresnel reflection matrix of vdisort::brdf::Fresnel.
 *  VDISORT reflects only between its quadrature streams; off-node user angles see the emission alone. */
struct fresnel_surface {
  Complex refractive_index{1.0, 0.0};
};

using surface = std::variant<lambertian_surface, fresnel_surface>;

//! The solver settings of main_data_from_path
struct path_settings {
  //! Number of streams (both hemispheres), even
  Index nquad{16};
  //! Number of Fourier azimuth modes
  Index nfourier{1};
  //! Azimuth samples of the Fourier coefficients (see scattering_optics)
  Index azimuth_count{64};
  //! Gauss-Legendre nodes of the phase-function normalisation (see scattering_optics)
  Index scattering_angle_count{512};
  //! See scattering_optics
  Numeric normalisation_tolerance{1e-3};
  //! Thermal emission of the layers and of the surface; the sky always emits
  bool thermal{true};
  //! Direct-beam flux on the horizontal plane at the top [W m-2 Hz-1], >= 0; 0 for no beam
  Numeric beam_flux{0.0};
  //! Cosine of the zenith angle of the beam (propagating downward), in (0, 1], not a stream
  Numeric beam_mu{0.5};
  //! VDISORT's beam azimuth phi0 [rad], in [0, 2 pi)
  Numeric beam_azimuth{0.0};
};

/** A solved VDISORT problem from an ARTS propagation path.
 *
 * The path conventions are those of rt4::problem_from_path,
 * rt3::problem_from_path and the DISORT workspace methods: ray_path,
 * atm_path and spectral_propmat_path have one entry per level, top first,
 * with strictly decreasing ray_path altitudes; the gas extinction of a layer
 * is the mean of the A elements of the unpolarized gas propagation matrix
 * (spectral_propmat_path, per metre, without particles) at its two levels,
 * and polarized gas propagation matrices are rejected; the frequency is
 * freq_grid[freq_index].
 *
 * The particle optics of a layer are the mean of its two levels'
 * scattering_optics() on VDISORT's streams: mean extinction, mean
 * scattering, and the scattering-weighted mean of the normalised Fourier
 * coefficients.  Layer l then has
 *   tau_l (layer bottom) = sum_{j <= l} (gas_j + extinction_j) dz_j,
 *   omega_l = scattering_l / (gas_l + extinction_l),
 * and every layer must have a positive optical thickness.  With thermal
 * set, the source is the Planck function at the level temperatures, linear
 * in optical depth within each layer ([c0, c1] in the global optical depth,
 * as disort_settingsLayerThermalEmissionLinearInTau), and the surface emits
 * at surface_temperature.  The sky is an isotropic unpolarized blackbody at
 * sky_temperature (0 K for none).  With beam_flux > 0 the beam has the
 * Stokes irradiance [beam_flux / beam_mu, 0, 0, 0] normal to it and the
 * beam phase matrices come from scattering_optics() with mu_in = {-beam_mu}.
 */
main_data main_data_from_path(const ArrayOfPropagationPathPoint& ray_path,
                              const ArrayOfAtmPoint&             atm_path,
                              const ArrayOfPropmatVector&        spectral_propmat_path,
                              const AscendingGrid&               freq_grid,
                              Index                              freq_index,
                              const ArrayOfScatteringSpecies&    scattering_species,
                              const path_settings&               settings,
                              const surface&                     ground,
                              Numeric                            surface_temperature,
                              Numeric                            sky_temperature);
}  // namespace vdisort
