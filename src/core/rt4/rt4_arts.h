#pragma once

#include <atm.h>
#include <matpack.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>

#include "rt4.h"

/**
 * RT4 inputs from ARTS-native data.
 *
 * These helpers translate ARTS scattering species, atmospheric points and
 * propagation paths into rt4::layer_optics and rt4::problem.  They add no
 * physics of their own: the particle optics are ARTS's bulk scattering
 * properties in the laboratory frame
 * (ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_gridded), so
 * azimuthally randomly oriented (ARO) species are supported as well as
 * totally randomly oriented (TRO) ones.
 *
 * Direction conventions:
 *   - ARTS's scattering data take the zenith angles of the propagation
 *     directions (0 deg propagating straight up) and the azimuth difference
 *     of the propagation directions.  RT4's stream (down, mu) propagates
 *     toward the surface at the zenith angle 180 deg - acos(mu); (up, mu)
 *     propagates upward at acos(mu).
 *   - ARTS's laboratory-frame phase matrix equals the vector-geometry phase
 *     matrix in RT4's meridional basis (Q = I_v - I_h, v in the plane of the
 *     ray and the zenith) with the same I and Q; the [I, Q] block of the
 *     azimuthal mean does not depend on the azimuth sense.
 *
 * Units: ARTS's extinction, absorption and phase matrix are per metre and
 * the phase matrix per steradian, with (1 / 4 pi) int Z11 dOmega times 4 pi
 * equal to the scattering coefficient.  That is RT4's convention for
 * layer_optics, with problem::height in metres.
 */
namespace rt4 {
/** The particle optics of scattering species at one atmospheric point on RT4's streams.
 *
 * mu: [nmu_total] the distinct stream cosines of one hemisphere, each in
 *   (0, 1]: RT4's quadrature nodes followed by the extra angles, as in
 *   problem.  The same values are used in both hemispheres.
 * nstokes: 1 ([I]) or 2 ([I, Q]).
 * azimuth_count: N, the number of azimuth differences of the periodic
 *   midpoint rule for the azimuthal mean, even and >= 2.
 *
 * Returns, with h the hemisphere (rt4::down or rt4::up) and mu_i = mu[i]:
 *   extinction[h, i]: the ARO extinction matrix for propagation in (h, mu_i),
 *     [[K11, K12], [K12, K11]] (K11 only for nstokes = 1); TRO species give
 *     K11 on the diagonal.
 *   absorption[h, i]: [a1, a2] (a1 only for nstokes = 1).
 *   phase[ho, hi, o, i]: the [I, Q] block of the azimuthal mean
 *     (1 / 2 pi) int Z(ho, mu_o <- hi, mu_i; dphi) ddphi of ARTS's
 *     laboratory-frame phase matrix, by the midpoint rule at
 *     dphi_k = (k + 1/2) 360 deg / N.  The samples above 180 deg mirror those
 *     below, whose [I, Q] blocks are the same for a medium with mirror
 *     symmetry (which ARTS's TRO and ARO formats have), so ARTS is only asked
 *     for the N / 2 azimuths in (0, 180) deg.
 *
 * Accuracy of the azimuthal mean.  For a scattering matrix that is a regular
 * Legendre series of degree L in cos(Theta), Z is a trigonometric polynomial
 * of degree L in dphi, and the rule is exact for N > L (Rayleigh: N >= 4).
 * Otherwise it converges as fast as the Fourier series of Z in dphi.  The
 * midpoints avoid dphi = 0 and 180 deg, where ARTS's rotation coefficients
 * snap angles within about 1e-3 rad of the principal plane.  The data are as
 * accurate as ARTS's laboratory-frame phase matrix: GasScatterer and
 * HenyeyGreensteinScatterer evaluate their closed-form scattering matrix at
 * the exact scattering angle of every direction pair; particle habits
 * interpolate linearly on their own scattering-angle grid.
 *
 * Vertical rays.  When both rays of a pair are vertical (mu = 1 in both),
 * their meridional planes, and so Q, are undefined.  ARTS then applies no
 * reference-plane rotation and the mean is F at 0 or 180 deg, where a
 * reference plane that turns with the azimuth label would average Q to 0.
 * Only extra angles can be vertical, and RT4 gives their columns no weight,
 * so the solution does not depend on these entries.
 *
 * An empty species array gives all-zero optics.  Species without a phase
 * matrix, or that do not provide laboratory-frame (ARO gridded) data, are an
 * error.
 */
layer_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                               const AtmPoint&                 atm_point,
                               Numeric                         frequency,
                               const Vector&                   mu,
                               Index                           nstokes,
                               Index                           azimuth_count);

//! The solver settings of problem_from_path
struct path_settings {
  //! 1 for [I], 2 for [I, Q]
  Index nstokes{2};
  //! Quadrature nodes per hemisphere
  Index           nmu{8};
  quadrature_type quad{quadrature_type::double_gauss};
  //! Zero-weight output angles appended after the quadrature streams, each in (0, 1]
  Vector extra_mu{};
  //! Maximum vertical optical thickness of the initial doubling sublayer, > 0
  Numeric max_delta_tau{1e-6};
  //! Azimuth differences of the azimuthal mean (see scattering_optics), even, >= 2
  Index azimuth_count{64};
};

/** An RT4 problem from an ARTS propagation path.
 *
 * The conventions are those of the DISORT workspace methods
 * (disort_settingsOpticalThicknessFromPath and
 * disort_settingsLayerThermalEmissionLinearInTau):
 *
 *   - ray_path, atm_path and spectral_propmat_path have one entry per level,
 *     top first, as for a down-looking path.  The altitudes of ray_path must
 *     decrease strictly; they are the RT4 heights [m], so layer l lies
 *     between levels l and l + 1.  Only the altitudes are used, not the
 *     lines of sight or the horizontal positions.
 *   - The level temperatures are atm_path's; RT4 makes the Planck function
 *     linear in optical depth within each layer.
 *   - spectral_propmat_path is the gas propagation matrix only (no
 *     particles), per metre, with freq_grid.size() entries per level.  The
 *     gas extinction of a layer is the mean of the A elements at its two
 *     levels.  Polarized gas propagation matrices are rejected: RT4's gas
 *     extinction is scalar.
 *   - The particle optics are scattering_optics() at every level on the
 *     streams of settings.  Each layer gets the mean of its two levels'
 *     optics, or is gas-only (layer_optics_index < 0) when that mean is all
 *     zero.
 *   - The frequency is freq_grid[freq_index].
 *
 * settings, ground, surface_temperature and sky_temperature go into the
 * problem unchanged.  The streams are RT4's own quadrature, so this needs
 * ENABLE_RT4.  RT4 still requires the optics to be mirror symmetric between
 * the hemispheres (solve() checks it).
 */
problem problem_from_path(const ArrayOfPropagationPathPoint& ray_path,
                          const ArrayOfAtmPoint&             atm_path,
                          const ArrayOfPropmatVector&        spectral_propmat_path,
                          const AscendingGrid&               freq_grid,
                          Index                              freq_index,
                          const ArrayOfScatteringSpecies&    scattering_species,
                          const path_settings&               settings,
                          const surface&                     ground,
                          Numeric                            surface_temperature,
                          Numeric                            sky_temperature);
}  // namespace rt4
