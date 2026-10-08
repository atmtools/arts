#pragma once

#include <matpack.h>

#include <variant>
#include <vector>

/**
 * Evans' RT3 (polradtran) polarized doubling-adding solver as a reference.
 *
 * RT3 solves the radiative transfer equation for a plane-parallel medium
 * of randomly oriented particles with a plane of symmetry, with a solar
 * beam and thermal sources, for every Fourier azimuth mode and the Stokes
 * components [I], [I, Q], [I, Q, U] or [I, Q, U, V] (K. F. Evans and
 * G. L. Stephens, 1991, JQSRT 46, 413-423).  The Fortran sources are in
 * 3rdparty/polradtran (the ARTS3 changes are listed in its README.ARTS).
 * RT3 is ported to C++: rt3::solve calls rt3::radtran (radtran3.h) and
 * rt3::ground_surface, which call no Fortran; the Fortran is built as the
 * reference the port is tested against.
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
 *     (only with the gauss quadrature; RT3's 'E' type).
 *     nmu_total = nmu + extra_mu.size().
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
 * None of the Fortran code is reentrant (COMMON blocks, SAVEd FFT tables
 * and static local arrays).  rt3::solve still serialises its calls with
 * the mutex that guarded it, separate from RT4's, although rt3::radtran
 * no longer calls the Fortran.
 */
namespace rt3 {
//! Whether the optional Fortran backend is built (ENABLE_RT3=ON).
bool available();

enum class quadrature_type {
  gauss,         //!< RT3 'G': positive half of a 2*nmu-point Gauss-Legendre rule on [-1, 1]
  double_gauss,  //!< RT3 'D': nmu-point Gauss-Legendre rule on [0, 1]
  lobatto,       //!< RT3 'L': positive half of a 2*nmu-point Lobatto rule on [-1, 1]; includes mu = 1
};

//! One hemisphere's streams: ascending mu in (0, 1]; weights for the
//! integral over mu in [0, 1], summing to 1.  The 2 pi azimuth factor is not
//! included.
struct quadrature {
  Vector mu;
  Vector weights;
};

/** The streams of RT3, nmu >= 1.
 *
 * ARTS's quadratures (scattering/integration.h) in place of RT3's own: the
 * positive half of scattering::DoubleGaussQuadrature,
 * GaussLegendreQuadrature or LobattoQuadrature of degree 2 nmu.  They are
 * RT3's rules, to rounding, as RT3_DOUBLE_GAUSS_QUADRATURE,
 * RT3_GAUSS_LEGENDRE_QUADRATURE and RT3_LOBATTO_QUADRATURE compute them.
 * rt3::radtran uses the same.  Needs no Fortran.
 */
quadrature get_quadrature(Index nmu, quadrature_type type);

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
 * legendre: [nleg + 1, 6], the scattering-plane phase matrix
 *     F(Theta) = [[F11, F12, 0, 0], [F12, F22, 0, 0],
 *                 [0, 0, F33, F34], [0, 0, -F34, F44]]
 *   as plain Legendre series in cos(Theta),
 *     F_c(Theta) = sum_l legendre[l, c] P_l(cos(Theta)),
 *   with the columns in the order of RT3's scattering files,
 *     c = 0: F11, 1: F12, 2: F33, 3: F34, 4: F22, 5: F44
 *   (SUM_LEGENDRE in radscat3.f).  The basis is that of the scattering
 *   plane, Q = I_par - I_perp, so Rayleigh scattering has
 *   F12 = -3/4 sin^2(Theta), i.e. legendre = [[1, -1/2, 0, 0, 1, 0],
 *   [0, 0, 3/2, 0, 0, 3/2], [1/2, 1/2, 0, 0, 1/2, 0]].  The phase function
 *   is normalised to 1 over 4 pi: legendre[0, 0] must be 1.  Note that the
 *   coefficients include the factor 2 l + 1, e.g. Henyey-Greenstein is
 *   legendre[l, 0] = (2 l + 1) g^l.  Trailing all-zero rows are ignored.
 */
struct scattering_set {
  Numeric extinction{};
  Numeric scattering{};
  Matrix  legendre{};
};

//! RT3 'L'.  Reflection 2 A mu_j w_j of the m = 0 mode into every stream,
//! I to I only; emission (1 - A) B(surface_temperature) when thermal is
//! set; reflection of the direct beam A F_direct / pi.  Energy is conserved
//! on the streams only for double_gauss quadrature (2 sum mu w = 1); gauss
//! and lobatto are off by about 3e-3 A for 8 streams.
struct lambertian_surface {
  Numeric albedo{0.0};
};

//! RT3 'F'.  Specular Fresnel reflection under a medium of index 1, for
//! every Fourier mode: with r_v, r_h the amplitude reflection coefficients,
//! R = [[R1, R2, 0, 0], [R2, R1, 0, 0], [0, 0, R3, -R4], [0, 0, R4, R3]],
//! R1 = (|r_v|^2 + |r_h|^2) / 2, R2 = (|r_v|^2 - |r_h|^2) / 2,
//! R3 = Re(r_v r_h*), R4 = Im(r_v r_h*); emission
//! [(1 - R1) B, -R2 B, 0, 0] with B = B(surface_temperature), always
//! included (also when thermal is false).  RT3 does not allow a Fresnel
//! surface with a direct beam.
struct fresnel_surface {
  Complex refractive_index{1.0, 0.0};
};

using surface = std::variant<lambertian_surface, fresnel_surface>;

struct problem {
  //! 1 for [I], 2 for [I, Q], 3 for [I, Q, U], 4 for [I, Q, U, V]
  Index nstokes{4};
  //! Quadrature nodes per hemisphere
  Index           nmu{8};
  quadrature_type quad{quadrature_type::gauss};
  //! Zero-weight output angles appended after the quadrature streams, each
  //! in (0, 1]; only with quadrature_type::gauss (RT3's 'E' type)
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
  //! Scattering sets, at most 200
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
 * Validates the shapes and every precondition on which the Fortran code
 * would STOP or overrun a buffer (see the limits below), rejects a Legendre
 * series that RT3 would silently truncate, then calls RADTRAN.
 *
 * Limits (N = nstokes * nmu_total, A = aziorder, beam = direct_flux > 0):
 *   N <= 64; nlay <= 200; (nlay + 1) N^2 <= 101 * 4096;
 *   scattering_sets.size() <= 200;
 *   scattering_sets.size() (A + 1) 2 N^2 <= 16 * 200 * 2 * 4096;
 *   beam: (A + 1) 2 N max(nlay, scattering_sets.size()) <= 16 * 200 * 2 * 64;
 *   2 A + 1 <= 512 (beam) or 1024 (no beam);
 *   per set: nleg <= 1023; legendre[0, 0] == 1 to 1e-9 (after delta-M);
 *   no non-zero coefficient above max_legendre_degree() (after delta-M);
 *   A > 0: the degree RT3 sums, min(degree, max_legendre_degree()), <= 251;
 *   beam: Lambertian surface.
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
}  // namespace rt3
