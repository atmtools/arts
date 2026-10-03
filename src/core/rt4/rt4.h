#pragma once

#include <matpack.h>

#include <variant>
#include <vector>

/**
 * Evans' RT4 (polradtran) polarized doubling-adding solver as a reference.
 *
 * RT4 solves the thermal-only radiative transfer equation for a
 * plane-parallel, azimuthally symmetric medium for the Stokes components
 * [I] or [I, Q].  The Fortran sources are the ARTS 2.6 version in
 * 3rdparty/polradtran (with the ARTS3 changes listed in its README).  This
 * wrapper exists so that RT4 can serve as an external reference for other
 * solvers; it has no workspace layer.
 *
 * Conventions:
 *   - Stokes basis [I, Q] with the meridional-plane reference: "vertical"
 *     polarization lies in the plane of the ray and the z axis, and
 *     Q = I_v - I_h.  The same basis is used in both hemispheres.
 *   - Streams are given per hemisphere by mu = |cos(zenith)| in (0, 1],
 *     ascending, the same mu values in both hemispheres.  The first nmu are
 *     RT4's quadrature nodes, followed by the zero-weight extra_mu angles.
 *     nmu_total = nmu + extra_mu.size().
 *   - Hemisphere index down (0) is radiation propagating downward, toward
 *     increasing optical depth, which is RT4's "+"; up (1) is propagating
 *     upward, RT4's "-".
 *   - Layers and levels are ordered top-down.  Level 0 is the top of the
 *     atmosphere and level nlay is just above the surface.
 *   - Radiances are in W m-2 Hz-1 sr-1.
 *   - RT4 hard-codes its doubling as SYMMETRIC: the medium must be mirror
 *     symmetric between the hemispheres, i.e. extinction[down] ==
 *     extinction[up], absorption[down] == absorption[up],
 *     phase[down, down] == phase[up, up] and phase[down, up] ==
 *     phase[up, down].  solve() rejects optics that are not.
 *
 * None of the Fortran code is reentrant (COMMON blocks and about 40 MB of
 * static local arrays).  All calls are serialised by one global mutex.
 */
namespace rt4 {
//! Whether the optional Fortran backend is built (ENABLE_RT4=ON).
bool available();

enum class quadrature_type {
  double_gauss,  //!< RT4 'D': nmu-point Gauss-Legendre rule on [0, 1]
  gauss,         //!< RT4 'G': positive half of a 2*nmu-point Gauss-Legendre rule on [-1, 1]
  lobatto,       //!< RT4 'L': positive half of a 2*nmu-point Lobatto rule on [-1, 1]; includes mu = 1
};

//! One hemisphere's streams: ascending mu in (0, 1]; weights for the
//! integral over mu in [0, 1], summing to 1.  The 2 pi azimuth factor is not
//! included.
struct quadrature {
  Vector mu;
  Vector weights;
};

//! RT4's own quadrature routines.  nmu >= 1.
quadrature get_quadrature(Index nmu, quadrature_type type);

//! Hemisphere indices (see the conventions above).
inline constexpr Index down = 0;
inline constexpr Index up   = 1;

/** Particle optics of one homogeneous layer on the solver streams.
 *
 * All quantities are per unit length, in the inverse of the unit of
 * problem::height.  Gas extinction is added separately by RT4.
 *
 * extinction: [2 hemisphere, nmu_total, nstokes (row), nstokes (col)],
 *   the extinction matrix K for propagation in that hemisphere at that mu.
 * absorption: [2 hemisphere, nmu_total, nstokes], the absorption vector a;
 *   RT4 multiplies it by the Planck function.
 * phase: [2 out hemisphere, 2 in hemisphere, nmu_total out, nmu_total in,
 *   nstokes out, nstokes in], the azimuthal mean (1 / 2 pi) int Z dDeltaphi
 *   of the phase matrix, per unit length and per steradian.  It includes the
 *   number density and no quadrature weights.  With quadrature weights w,
 *   energy conservation on the streams reads
 *     K11(h, mu_j) = a1(h, mu_j)
 *                  + 2 pi sum_i w_i [Z(up <- h)(1, i; 1, j) + Z(down <- h)(1, i; 1, j)].
 *   Neither RT4 nor this wrapper checks or enforces it.
 *
 * RT4 chooses the number of doublings from extinction[down, 0, 0, 0] plus
 * the gas extinction, so that the initial sublayer has a vertical optical
 * thickness of at most problem::max_delta_tau; the slant thickness at small
 * mu is larger.
 */
struct layer_optics {
  Tensor4 extinction;
  Tensor3 absorption;
  Tensor6 phase;
};

//! RT4 'L'.  Reflection 2 A mu_j w_j into every stream, I to I only;
//! emission [(1 - A) B, 0].  Energy is conserved on the streams only for
//! double_gauss quadrature (sum 2 mu w = 1); gauss and lobatto are off by
//! about 3e-3 A for 8 streams.
struct lambertian_surface {
  Numeric albedo{0.0};
};

//! RT4 'F'.  Specular Fresnel reflection under a medium of index 1 with
//! R = [[R1, R2], [R2, R1]], R1 = (|r_v|^2 + |r_h|^2) / 2,
//! R2 = (|r_v|^2 - |r_h|^2) / 2, and emission [(1 - R1) B, -R2 B].
struct fresnel_surface {
  Complex refractive_index{1.0, 0.0};
};

//! RT4 'S'.  reflectivity: [nstokes, nstokes], R(out, in), applied
//! specularly to every stream; emission [(1 - R(I, I)) B, -R(Q, I) B].
struct specular_surface {
  Matrix reflectivity;
};

/** RT4 'A' (external surface).
 *
 * reflection: [nmu_total out (up), nmu_total in (down), nstokes out,
 *   nstokes in], the discrete operator
 *     I_up(i) = sum_j reflection(i, j) I_down(j) + emission(i),
 *   i.e. RT4's SURF_REFLECT including any quadrature factors (a Lambertian
 *   surface is 2 A mu_j w_j).
 * emission: [nmu_total, nstokes] in W m-2 Hz-1 sr-1.
 *
 * problem::surface_temperature is not used for this surface.
 */
struct discrete_surface {
  Tensor4 reflection;
  Matrix  emission;
};

using surface = std::variant<lambertian_surface, fresnel_surface, specular_surface, discrete_surface>;

struct problem {
  //! 1 for [I], 2 for [I, Q]
  Index nstokes{2};
  //! Quadrature nodes per hemisphere
  Index           nmu{8};
  quadrature_type quad{quadrature_type::double_gauss};
  //! Zero-weight output angles appended after the quadrature streams, each in (0, 1]
  Vector extra_mu{};
  //! Maximum vertical optical thickness of the initial doubling sublayer, > 0
  Numeric max_delta_tau{1e-6};
  //! Frequency [Hz]; RT4 is given the wavelength 1e6 c / f in micrometres
  Numeric frequency{};
  //! [nlay + 1] layer interfaces, top-down.  Only |differences| are used;
  //! the unit must be the inverse of the extinction unit.
  Vector height{};
  //! [nlay + 1] temperatures [K] at the interfaces, top-down, > 0.  The
  //! Planck function, not the temperature, is linear within each layer.
  Vector temperature{};
  //! [nlay] scalar, unpolarized gas extinction per unit length, >= 0.  It is
  //! added to every Stokes diagonal of K and to the I absorption.
  Vector gas_extinction{};
  //! Particle optics sets
  std::vector<layer_optics> optics{};
  //! [nlay] index into optics per layer, or < 0 for a gas-only layer.  Gas-only
  //! layers are solved analytically (exp(-tau / mu), Planck linear in tau).
  ArrayOfIndex layer_optics_index{};
  //! Temperature [K] of the isotropic, unpolarized blackbody incident at the top
  Numeric sky_temperature{};
  //! Surface temperature [K], used by the lambertian, fresnel and specular surfaces
  Numeric surface_temperature{};
  surface ground{};
};

/** The solution at every level.
 *
 * mu, weights: [nmu_total], the streams (extra angles have weight 0).
 * up, down: [nlay + 1 level (0 = top), nmu_total, nstokes] in
 *   W m-2 Hz-1 sr-1; up is the radiance propagating upward (seen when
 *   looking down at zenith angle pi - acos(mu)), down propagating downward.
 */
struct result {
  Vector  mu;
  Vector  weights;
  Tensor3 up;
  Tensor3 down;
};

/** Run RT4.
 *
 * Validates the shapes, every precondition on which the Fortran code would
 * STOP (nstokes <= 2, nstokes * nmu_total <= 64, nlay <= 400,
 * (nlay + 1) * (nstokes * nmu_total)^2 <= 301 * 4096), and the mirror
 * symmetry that RT4 requires, then calls RADTRANO.
 */
result solve(const problem& p);
}  // namespace rt4
