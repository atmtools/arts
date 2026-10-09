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
 * properties in the laboratory frame, as the azimuthal Fourier modes the
 * species give at the streams
 * (ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_fourier), so
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
namespace polradtran::rt4 {
/** The particle optics of scattering species at one atmospheric point on RT4's streams.
 *
 * mu: [nmu_total] the distinct stream cosines of one hemisphere, each in
 *   (0, 1]: RT4's quadrature nodes followed by the extra angles, as in
 *   problem.  The same values are used in both hemispheres.
 * nstokes: 1 ([I]) or 2 ([I, Q]).
 *
 * Returns, with h the hemisphere (rt4::down or rt4::up) and mu_i = mu[i]:
 *   extinction[h, i]: the ARO extinction matrix for propagation in (h, mu_i),
 *     [[K11, K12], [K12, K11]] (K11 only for nstokes = 1); TRO species give
 *     K11 on the diagonal.
 *   absorption[h, i]: [a1, a2] (a1 only for nstokes = 1).
 *   phase[ho, hi, o, i]: the [I, Q] block of the azimuthal mean
 *     (1 / 2 pi) int Z(ho, mu_o <- hi, mu_i; dphi) ddphi of ARTS's
 *     laboratory-frame phase matrix: the m = 0 azimuthal Fourier mode the
 *     species give at the streams
 *     (ArrayOfScatteringSpecies::get_bulk_scattering_properties_aro_fourier).
 *     GasScatterer and HenyeyGreensteinScatterer give it exactly from their
 *     closed forms, particle habits from their Legendre series (gridded TRO
 *     particle data must be converted to one first).
 *
 * Vertical rays.  A vertical ray (mu = 1, e.g. the last Lobatto node or an
 * extra angle) has the meridional plane of its azimuth label in ARTS's
 * laboratory frame, the limit along its meridian, so vertical rays need no
 * special treatment.
 *
 * What rt4::solve then requires of the optics, and rejects otherwise:
 *  - Mirror symmetry between the hemispheres (RT4's SYMMETRIC).  Azimuthal
 *    random orientation does not imply it: species whose orientation is not
 *    symmetric under reflection in the horizontal plane (tilted or
 *    asymmetric particles) give optics that rt4::solve refuses.
 *  - Energy conservation on the streams (problem::normalisation_tolerance):
 *    every incident quadrature stream must scatter K11 - a1 into the
 *    quadrature streams.  The azimuthal mean of Z sampled at the streams
 *    need not, for forward-peaked optics (large drops, ice at high
 *    frequency) whose peak falls between the streams; use more streams.
 *
 * An empty species array gives all-zero optics.  Species without a phase
 * matrix, or that cannot give the Fourier modes, are an error.
 */
layer_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                               const AtmPoint&                 atm_point,
                               Numeric                         frequency,
                               const Vector&                   mu,
                               Index                           nstokes);

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
  //! problem::normalisation_tolerance, >= 0; infinity checks nothing
  Numeric normalisation_tolerance{1e-6};
};

/** An RT4 problem from an ARTS propagation path.
 *
 * The heights [m], temperatures, gas extinction and frequency are
 * polradtran::layers_from_path's (polradtran_arts.h), with the conventions
 * of the DISORT workspace methods: ray_path, atm_path and
 * spectral_propmat_path have one entry per level, top first, so layer l
 * lies between levels l and l + 1, and the gas extinction of a layer is the
 * mean of the A elements of the unpolarized gas propagation matrix at its
 * two levels.  RT4 makes the Planck function linear in optical depth within
 * each layer.
 *
 * The particle optics are scattering_optics() at every level on the streams
 * of settings.  Each layer gets the mean of its two levels' optics, or is
 * gas-only (layer_optics_index < 0) when that mean is all zero.
 *
 * settings, ground, surface_temperature and sky_temperature go into the
 * problem unchanged.  The streams are RT4's own quadrature.  RT4 still
 * requires the optics to be mirror symmetric between
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
}  // namespace polradtran::rt4
