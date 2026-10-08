#pragma once

#include <atm.h>
#include <matpack.h>
#include <path_point.h>
#include <rtepack.h>
#include <scattering_species.h>

#include "rt3.h"

/**
 * RT3 inputs from ARTS-native data.
 *
 * These helpers translate ARTS scattering species, atmospheric points and
 * propagation paths into rt3::scattering_set and rt3::problem.  RT3 takes
 * the scattering matrix of totally randomly oriented (TRO) particles as
 * Legendre series, and gets them from the species' own Legendre series
 * (ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_spectral).
 *
 * Stokes and sign conventions.  ARTS's TRO scattering matrix is stored as
 * [F11, F12, F22, F33, F34, F44] in the scattering-plane basis with
 * Q = I_par - I_perp; RT3's columns are (F11, F12, F33, F34, F22, F44).  The
 * elements are reordered without any sign change.  That is the mapping under
 * which RT3 transports ARTS's Stokes vector: ARTS's laboratory-frame phase
 * matrix (phase_matrix.h, for the propagation directions (za, aa)) equals the
 * vector-geometry phase matrix of the same F in RT3's meridional basis
 * (Q = I_v - I_h, U = 2 Re(E_v E_h*)) with mu = cos(za) and RT3's azimuth
 * phi = -aa, element by element (tested in cpp.fast.vdisort-arts-test).
 * ARTS's azimuth runs clockwise seen from above, RT3's counterclockwise.  For
 * the same physical sphere, ARTS's Mie code gives F12 and F33 equal to, and
 * F34 of the opposite sign of, the Legendre series of Evans and Stephens
 * (1991) (runmietest; tested in cpp.fast.rt3-arts-test), so V computed by RT3
 * from ARTS data is -1 times V in the convention of that paper.
 *
 * Units: extinction and scattering per metre, with problem::height in
 * metres.
 */
namespace polradtran::rt3 {
/** The scattering set of scattering species at one atmospheric point.
 *
 * The species give their Legendre series to degree themselves
 * (ArrayOfScatteringSpecies::get_bulk_scattering_properties_tro_spectral,
 * coefficients a_l on the orthonormal Y_l0); these are the Legendre
 * coefficients c_l = a_l sqrt((2 l + 1) / 4 pi) of F = sum_l c_l P_l.
 * Species that cannot give them (gridded particle data, ARO data) are an
 * error from the species.
 *
 * Returns:
 *   extinction: K11 of the particles per metre.
 *   scattering: K11 - a1, the extinction minus the absorption.
 *   legendre: [degree + 1, 6], the series c / c_0(F11), in RT3's column order
 *     (F11, F12, F33, F34, F22, F44), so legendre[0, 0] = 1.
 *
 * The phase-function integral 4 pi c_0(F11) = 2 pi int F11 dcos(Theta) and
 * the scattering coefficient K11 - a1 must agree to
 * normalisation_tolerance * K11, otherwise this is an error: RT3 normalises
 * the series and takes the albedo from the scattering coefficient, so a
 * mismatch would silently change the scattered energy.  An infinite
 * normalisation_tolerance checks nothing.
 *
 * Without particles (an empty species array, or zero extinction and phase
 * matrix) the set has zero extinction and scattering and the isotropic,
 * depolarizing series legendre[0] = [1, 0, 0, 0, 0, 0].
 */
scattering_set scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                 const AtmPoint&                 atm_point,
                                 Numeric                         frequency,
                                 Index                           degree,
                                 Numeric                         normalisation_tolerance);

//! The solver settings of problem_from_path
struct path_settings {
  //! 1 for [I], 2 for [I, Q], 3 for [I, Q, U], 4 for [I, Q, U, V]
  Index nstokes{4};
  //! Quadrature nodes per hemisphere
  Index           nmu{8};
  quadrature_type quad{quadrature_type::gauss};
  //! Zero-weight output angles, each in (0, 1]; only with quadrature_type::gauss
  Vector extra_mu{};
  //! Highest Fourier azimuth mode, >= 0
  Index aziorder{0};
  //! Maximum vertical optical thickness of the initial doubling sublayer, > 0
  Numeric max_delta_tau{1e-6};
  //! RT3's delta-M scaling
  bool delta_m{false};
  /** Degree of the Legendre series.  Negative selects RT3's limit
   *  max_legendre_degree(nmu, quad), or, with delta_m, at least
   *  2 (nmu + extra_mu.size()), the degree at which RT3 reads the delta-M
   *  fraction. */
  Index legendre_degree{-1};
  //! See scattering_optics; infinity checks nothing
  Numeric normalisation_tolerance{1e-3};
};

/** An RT3 problem from an ARTS propagation path.
 *
 * The path conventions are those of rt4::problem_from_path and of the DISORT
 * workspace methods: ray_path, atm_path and spectral_propmat_path have one
 * entry per level, top first, with strictly decreasing ray_path altitudes
 * (the heights, in metres); the level temperatures are atm_path's; the gas
 * extinction of a layer is the mean of the A elements of the unpolarized gas
 * propagation matrix (spectral_propmat_path, per metre, without particles)
 * at its two levels, and polarized gas propagation matrices are rejected;
 * the frequency is freq_grid[freq_index].
 *
 * The scattering set of a layer is the mean of its two levels'
 * scattering_optics(): mean extinction, mean scattering, and the
 * scattering-weighted mean of the two normalised series, which is the
 * normalised series of the mean phase matrix.  A layer whose mean extinction
 * is zero is gas-only.
 *
 * The problem has thermal emission and no direct beam; set direct_flux,
 * direct_mu and thermal on the result for other sources.
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
}  // namespace polradtran::rt3
