#pragma once

#include <matpack.h>
#include <polradtran.h>
#include <rtepack.h>

#include <variant>
#include <vector>

/**
 * Evans' RT3 (polradtran) polarized doubling-adding solver as a reference.
 *
 * RT3 solves the radiative transfer equation for a plane-parallel medium
 * of randomly oriented particles with a plane of symmetry, with a solar
 * beam and thermal sources, for every Fourier azimuth mode and the Stokes
 * components [I], [I, Q], [I, Q, U] or [I, Q, U, V] (K. F. Evans and
 * G. L. Stephens, 1991, JQSRT 46, 413-423).  It is the C++ port of Evans'
 * Fortran (3rdparty/polradtran keeps its licence and benchmarks):
 * rt3::solve calls rt3::radtran (radtran3.h) and rt3::ground_surface.
 * This wrapper exists so that RT3 can serve as an external reference for
 * other solvers, in particular for the solar-beam, m > 0 and U, V paths of
 * VDISORT; it has no workspace layer.
 *
 * Conventions:
 *   - Geometry: a right-handed frame with z pointing up (out of the
 *     medium).  A ray propagating at zenith cosine mu_z (> 0 upward) and
 *     azimuth phi has the direction k = (sin(theta) cos(phi),
 *     sin(theta) sin(phi), mu_z).  The direct beam propagates downward
 *     toward phi = 0.  So phi is the azimuth of the propagation direction
 *     of the ray, counted from x toward y, relative to the propagation
 *     direction of the beam; phi = 0 for upwelling radiation is forward
 *     scattering in azimuth.  This is the frame of VDISORT's polarized
 *     tests (tests/core/disort/vdisort-interface.py) with phi0 = 0.
 *   - Stokes basis [I, Q, U, V] with the meridional plane as reference:
 *     h = k x z / |k x z| is horizontal, v = h x k lies in the plane of k
 *     and z, I = I_v + I_h, Q = I_v - I_h and U = 2 Re(E_v E_h*).  The same
 *     basis is used in both hemispheres.  The sign of U is pinned by the
 *     single-scattering test against this vector construction.  RT3
 *     transports V through F34 (and the Fresnel R4), so V has the sign
 *     convention of the scattering data; see rt3_arts.h for ARTS's.
 *   - Streams are given per hemisphere by mu = |cos(zenith)| in (0, 1],
 *     ascending, the same mu values in both hemispheres.  The first nmu are
 *     RT3's quadrature nodes, followed by the zero-weight extra_mu angles
 *     (RT3's 'E' type had them with the gauss quadrature only; the port
 *     takes them with any).  nmu_total = nmu + extra_mu.size().
 *   - Hemispheres: "down" is radiation propagating downward, toward
 *     increasing optical depth (RT3's "+", printed with mu > 0 by rt3.f),
 *     "up" propagating upward (RT3's "-", printed with mu < 0).
 *   - Layers and levels are ordered top-down.  Level 0 is the top of the
 *     atmosphere and level nlay is just above the surface.
 *   - Radiances are in W m-2 Hz-1 sr-1, fluxes in W m-2 Hz-1.
 *   - Fourier series: the radiance at azimuth phi is
 *       I(phi) = sum_m c_m cos(m phi) for I and Q,
 *       I(phi) = sum_m c_m sin(m phi) for U and V,
 *     with c_m the m-th coefficient of result::up or result::down (RT3's
 *     cosine modes of I, Q and sine modes of U, V; azimuth_radiance() sums
 *     the series as rt3.f's OUTPUT_FILE does).  The m = 0 coefficients of U
 *     and V are 0.
 *
 * The solver keeps no state between calls, so concurrent calls run in
 * parallel.
 */
namespace polradtran::rt3 {
/** Highest Legendre degree RT3 keeps for a quadrature (its NLEGLIM):
 *  gauss 4 nmu - 3, double_gauss 2 nmu - 3, lobatto 4 nmu - 5, at least 1.
 *  nmu counts the quadrature nodes only, not the extra angles.  RT3
 *  truncates longer series (after delta-M scaling); solve() rejects a
 *  problem in which that would drop a non-zero coefficient.
 */
Index max_legendre_degree(Index nmu, quadrature_type type);

/** Single-scattering properties of a homogeneous particle population.
 *
 * extinction, scattering: the particle extinction and scattering
 *   coefficients per unit length, in the inverse of the unit of
 *   problem::height.  Gas extinction is added separately by RT3.
 * legendre: [nleg + 1], the scattering-plane phase matrix
 *     F(Theta) = [[F11, F12, 0, 0], [F12, F22, 0, 0],
 *                 [0, 0, F33, F34], [0, 0, -F34, F44]]
 *   as a plain Legendre series in cos(Theta),
 *     F(Theta) = sum_l legendre[l] P_l(cos(Theta)),
 *   each coefficient a CompactPlanarMuelmat (its elements by name; RT3's
 *   scattering files have them in the order F11, F12, F33, F34, F22, F44).
 *   The basis is that of the scattering plane, Q = I_par - I_perp, so
 *   Rayleigh scattering has F12 = -3/4 sin^2(Theta), i.e. legendre =
 *   [{F11 1, F12 -1/2, F22 1}, {F33 3/2, F44 3/2},
 *   {F11 1/2, F12 1/2, F22 1/2}].  The phase function is normalised to 1
 *   over 4 pi: legendre[0].F11() must be 1.  Note that the coefficients
 *   include the factor 2 l + 1, e.g. Henyey-Greenstein has
 *   legendre[l].F11() = (2 l + 1) g^l.  Trailing all-zero coefficients are
 *   ignored.
 */
struct scattering_set {
  Numeric                    extinction{};
  Numeric                    scattering{};
  CompactPlanarMuelmatVector legendre{};
};

//! The grounds of RT3, polradtran's 'L' and 'F' (polradtran.h)
using surface = std::variant<lambertian_surface, fresnel_surface>;

struct problem {
  //! 1 for [I], 2 for [I, Q], 3 for [I, Q, U], 4 for [I, Q, U, V]
  Index nstokes{4};
  //! Quadrature nodes per hemisphere
  Index           nmu{8};
  quadrature_type quad{quadrature_type::gauss};
  //! Zero-weight output angles appended after the quadrature streams, each
  //! in (0, 1]
  Vector extra_mu{};
  //! Highest Fourier azimuth mode, >= 0
  Index aziorder{0};
  //! Maximum vertical optical thickness of the initial doubling sublayer, > 0
  Numeric max_delta_tau{1e-6};
  //! Delta-M scaling of every scattering set, with M = 2 nmu_total
  bool delta_m{false};
  //! Direct (solar) beam flux on the horizontal plane at the top
  //! [W m-2 Hz-1], >= 0; 0 switches the beam off
  Numeric direct_flux{0.0};
  //! Cosine of the zenith angle of the direct beam, in (0, 1]
  Numeric direct_mu{1.0};
  //! Thermal emission of the layers and of a Lambertian surface.  The sky
  //! and a Fresnel surface always emit (RT3 includes them regardless of
  //! its source code); set their temperatures to 0 K to remove them.
  bool thermal{true};
  //! Frequency [Hz]; RT3 is given the wavelength 1e6 c / f in micrometres
  Numeric frequency{};
  //! [nlay + 1] layer interfaces, top-down.  Only |differences| are used;
  //! the unit must be the inverse of the extinction unit.
  Vector height{};
  //! [nlay + 1] temperatures [K] at the interfaces, top-down; > 0 when
  //! thermal is set, unused otherwise
  Vector temperature{};
  //! [nlay] scalar, unpolarized gas extinction per unit length, >= 0
  Vector gas_extinction{};
  //! Scattering sets
  std::vector<scattering_set> scattering_sets{};
  //! [nlay] index into scattering_sets per layer, or < 0 for a gas-only layer
  ArrayOfIndex layer_scattering_index{};
  //! Temperature [K] of the isotropic, unpolarized blackbody incident at the top
  Numeric sky_temperature{};
  //! Surface temperature [K]
  Numeric surface_temperature{};
  surface ground{};
};

/** The solution at every level.
 *
 * mu, weights: [nmu_total], the streams (extra angles have weight 0).
 * up, down: [nlay + 1 level (0 = top), aziorder + 1 mode m, nmu_total,
 *   nstokes], the Fourier coefficients in W m-2 Hz-1 sr-1 (see the
 *   Fourier convention above).  up is the radiance propagating upward
 *   (seen when looking down at nadir angle acos(mu)), down propagating
 *   downward.  down does not include the direct beam.
 * up_flux, down_flux: [nlay + 1, nstokes] in W m-2 Hz-1, RT3's
 *   2 pi sum_i w_i mu_i c_0 (only I and Q can be non-zero); down_flux
 *   includes the direct beam F_direct exp(-tau / direct_mu) in I, with the
 *   delta-M scaled tau when delta_m is set.
 */
struct result {
  Vector  mu;
  Vector  weights;
  Tensor4 up;
  Tensor4 down;
  Matrix  up_flux;
  Matrix  down_flux;
};

/** Run RT3.
 *
 * Validates the shapes and the inputs, rejects a Legendre series that RT3
 * would silently truncate, then calls RADTRAN (rt3::radtran).  Every array
 * is sized to the problem: Evans' fixed array sizes (N = nstokes *
 * nmu_total <= 64, at most 200 layers and sets, his scattering, direct-beam
 * and azimuth buffers, Legendre degree <= 1023, and an FFT of at most 512
 * azimuths) are not limits of the port.
 *
 * Requirements (beam = direct_flux > 0): per scattering set,
 * legendre[0, 0] == 1 to 1e-9 and no non-zero coefficient above
 * max_legendre_degree() (both after delta-M); with a beam, a Lambertian
 * surface.
 */
result solve(const problem& p);

/** The radiance at azimuths phi [rad] from Fourier coefficients.
 *
 * coefficients: [nlevel, nmode, nmu, nstokes], e.g. result::up or
 *   result::down.
 * Returns [nlevel, phi.size(), nmu, nstokes] with
 *   sum_m coefficients[l, m, i, s] cos(m phi) for s = 0, 1 (I, Q) and
 *   sum_m coefficients[l, m, i, s] sin(m phi) for s = 2, 3 (U, V),
 * as rt3.f's OUTPUT_FILE does (in double instead of single precision).
 */
Tensor4 azimuth_radiance(const Tensor4& coefficients, const Vector& phi);
}  // namespace polradtran::rt3
