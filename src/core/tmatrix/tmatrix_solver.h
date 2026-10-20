#pragma once

#include <enumsParticleShape.h>
#include <enumsTMatrixSizeDistribution.h>
#include <matpack.h>
#include <rtepack.h>
#include <tmatrix.h>

#include <utility>
#include <vector>

/** The T-matrix solver
 *
 * Each function implements a section of "The equations solved" in
 * doc/arts/concept.tmatrix.rst, named in its documentation.  The public
 * interface in tmatrix.h validates the inputs and assembles the results.
 */
namespace tmatrix::solver {
//! A rotationally symmetric homogeneous particle (Particle shapes)
struct Particle {
  ParticleShape shape;
  /** The horizontal-to-rotational axis ratio of a spheroid, the diameter-to-length ratio of a
   * cylinder, the deformation epsilon of a Chebyshev particle */
  Numeric aspect_ratio;
  //! The volume-equivalent radius, in the length unit of the wavelength
  Numeric radius;
  Numeric wavelength;
  //! The refractive index relative to the surrounding medium
  Complex refractive_index;
  //! The degree n of a Chebyshev particle
  Index degree{0};
  /** The b_n of a generalized Chebyshev particle, r(theta) = radius sum_n b_n cos(n theta): its shape scaled to a
   * volume-equivalent radius of 1 (volume_normalized_coefficients) */
  Vector coefficients{};
};

/** The b_n of the generalized Chebyshev shape 1 + c_0 + sum_{n >= 1} c_n cos(n theta) scaled to a volume-equivalent
 * radius of 1 (Particle shapes) */
Vector volume_normalized_coefficients(const Vector& c);

//! Whether the particle is mirror symmetric in the plane z = 0 (Particle shapes)
bool mirror_symmetric(const Particle& particle);

//! r^2 and rho = r'/r of the surface at x = cos(theta), -1 < x < 1, for the particle's radius (Particle shapes)
std::pair<Numeric, Numeric> surface(const Particle& particle, Numeric x);

//! q = r_v / r_s, the volume- over the surface-area-equivalent radius, to the accuracy of the T-matrix (Particle shapes)
Numeric volume_to_surface_radius(const Particle& particle, Numeric accuracy);

//! The largest distance between two points of the particle surface (Particle shapes)
Numeric maximum_dimension(const Particle& particle);

/** The surface quadrature and the functions at its nodes (Surface quadrature, Special functions)
 *
 * The nodes cover the hemisphere 0 < x < 1 of a mirror-symmetric particle, else -1 < x < 1.
 * Per node i: x = cos(theta), s = sin(theta), W = w r^2, rho = r'/r, xi = k r.
 * Per node and degree n = 0, ..., nmax: the Bessel functions and their tilde
 * forms, and d^n_{0m} and tau^m_n for m = 0, ..., nmax (zero for m > n).
 */
struct Surface {
  Complex refractive_index;
  Index   nmax;
  bool    mirror;
  // The largest arrays first, so that an impossible size fails before the others are filled
  Tensor3       d, dd;
  Vector        x, s, W, rho, xi;
  Matrix        j, dj, y, dy;
  ComplexMatrix J, dJ;

  //! nodes is N_G; a particle without mirror symmetry has 2 N_G nodes
  Surface(const Particle& particle, Index nmax, Index nodes);
};

/** Q and Rg Q of the azimuthal order m (The T-matrix of one azimuthal order)
 *
 * solve computes the orders with storage shared between them; this allocates its own.
 *
 * @return [2N, 2N] each, N = nmax - max(m, 1) + 1, in the block order of TMatrix::orders
 */
std::pair<ComplexMatrix, ComplexMatrix> order_q(const Surface& surface, Index m);

//! The T-matrix of the orders m = 0, ..., orders().size() - 1 at a converged truncation
struct TMatrix {
  //! n_max
  Index nmax{};
  //! N_G, the nodes per hemisphere
  Index nodes{};
  //! n_1, the truncation of random orientation; n_max for fixed orientation
  Index n1{};
  /** T = -Rg Q Q^{-1} of each order m, [2N, 2N] with N = nmax - max(m, 1) + 1: the blocks T^{pq}
   * in the order [[T11, T12], [T21, T22]], degrees ascending from max(m, 1) */
  std::vector<ComplexMatrix> orders;

  //! The block [[T^{11}, T^{12}], [T^{21}, T^{22}]] of order m at the degrees max(m, 1) <= n, n' <= nmax
  [[nodiscard]] Jones block(Index m, Index n, Index nn) const;
};

/** Converge the truncation and compute the order matrices (Truncation and convergence)
 *
 * Fails once the order-0 T-matrix is not finite, and, for random orientation, once n_1 must exceed the Wigner
 * tables.
 *
 * @param particle The particle
 * @param accuracy delta, the relative accuracy of the cross sections
 * @param nodes_per_degree N_D
 * @param random For random orientation: find n_1 and compute the orders m <= n_1; otherwise m <= n_max
 */
TMatrix solve(const Particle& particle, Numeric accuracy, Index nodes_per_degree, bool random);

//! C_ext and C_sca (Orientation-averaged cross sections)
struct CrossSections {
  Numeric extinction;
  Numeric scattering;
};

//! The orientation-averaged cross sections; random selects the n_1 truncation of the scattering sum
CrossSections cross_sections(const TMatrix& T, Numeric wavelength, bool random);

//! The laboratory-frame amplitude matrix
struct Amplitude {
  //! [[S_vv, S_vh], [S_hv, S_hh]], the Jones matrix between ARTS's meridional bases of the two rays
  Jones S{Complex{}};
  //! The sum of the magnitudes of the terms of the particle-frame sums over k: the scale of their cancellation
  Numeric magnitude{};
};

/** The amplitude matrices of a particle in a fixed orientation (Amplitude matrix in a fixed orientation)
 *
 * @param T The T-matrix
 * @param wavelength The wavelength
 * @param geometries The geometries, in degrees
 * @return One amplitude matrix per geometry
 */
std::vector<Amplitude> amplitude(const TMatrix& T, Numeric wavelength, std::span<const Geometry> geometries);

/** The scattering properties of azimuthally random orientation, on ARTS's grids (Azimuthally random orientation)
 *
 * In length units of the wavelength, squared; the matrices in ARTS's meridional bases and Stokes convention.
 */
struct AzimuthallyRandom {
  //! The phase matrix [za_inc, delta_aa, za_scat]
  MuelmatTensor3 phase;
  //! K11, K12 and K34 [za_inc]
  CompactPlanarPropmatVector extinction;
  //! I and Q [za_inc]
  CompactPlanarStokvecVector absorption;
};

/** A particle of the tilt distribution beta_k, w_k with random azimuth (Azimuthally random orientation)
 *
 * The angles in degrees: za_inc and za_scat the zenith angles of the propagation directions, delta_aa the scattered
 * minus the incident azimuth of ARTS (clockwise from north), tilts the angles beta of the symmetry axis from the
 * vertical, with weights that sum to 1.  The averages over the azimuth and the integral of the absorption vector are
 * exact for the truncated T-matrix.
 */
AzimuthallyRandom azimuthally_random(const TMatrix& T,
                                     Numeric        wavelength,
                                     const Vector&  tilts,
                                     const Vector&  weights,
                                     const Vector&  za_inc,
                                     const Vector&  delta_aa,
                                     const Vector&  za_scat);

//! The expansion coefficients alpha_1..4 and beta_1..2 for l = 0, ..., l_max
struct Expansion {
  Vector alpha1, alpha2, alpha3, alpha4, beta1, beta2;

  //! The number of coefficients of each series, L + 1
  [[nodiscard]] Size size() const { return alpha1.size(); }

  //! Adds weight times other, extending this with zeros to the longer series
  void add(const Expansion& other, Numeric weight);

  Expansion& operator/=(Numeric x);
};

/** The expansion coefficients of random orientation (Expansion coefficients in random orientation)
 *
 * All of them, l = 0, ..., 2 n_1: the exact series of the truncated T-matrix.
 *
 * Uses the wigxjpf Clebsch-Gordan coefficients; the Wigner tables must have
 * been initialized (see WignerInit) for 3j symbols with j1 + j2 + j3 up to 4 n_1.
 *
 * @param T The T-matrix with orders up to n_1
 * @param scattering C_sca of T
 * @param wavelength The wavelength
 */
Expansion expansion(const TMatrix& T, Numeric scattering, Numeric wavelength);

/** The scattering matrix at equally spaced scattering angles 0..180 degrees (Scattering matrix from the expansion)
 *
 * In ARTS's Stokes convention: F34 is the negative of the series of beta_2.
 */
CompactPlanarMuelmatVector scattering_matrix(const Expansion& expansion, Index angles);

/** The Legendre series of the scattering matrix (Scattering matrix from the expansion)
 *
 * The coefficients a_l, l = 0, ..., l_max, on Y_l0 = sqrt((2 l + 1) / 4 pi) P_l(cos(Theta)) of each element
 * (scattering/tro_legendre.h), in ARTS's Stokes convention.  Every element is a polynomial of degree l_max in
 * cos(Theta), so the l_max + 1 node Gauss-Legendre rule gives them exactly.
 */
CompactPlanarMuelmatVector legendre_series(const Expansion& expansion);

//! A size distribution and its quadrature (Size distributions)
struct SizeDistribution {
  TMatrixSizeDistribution kind;
  Numeric                 a, b, gamma;
  //! r_1 / a and r_2 / a, unused for the power law
  Numeric lower_ratio, upper_ratio;
  //! N_k
  Index nodes;
};

//! The quadrature radii and normalized weights of a size distribution, the bounds r_1 and r_2 of its radii, and its
//! effective radius and variance
struct SizeQuadrature {
  Vector  radius, weight;
  Numeric lower, upper;
  Numeric effective_radius, effective_variance;
};

SizeQuadrature size_quadrature(const SizeDistribution& distribution);

//! The size-distribution average of random orientation
struct Average {
  CrossSections cross_sections;
  Expansion     expansion;
  Numeric       effective_radius, effective_variance;
};

/** Average the cross sections and expansion coefficients over a size distribution (Size distributions)
 *
 * @param particle The particle; its radius is replaced by radius_factor times each node radius
 * @param radius_factor q for surface-area-equivalent node radii, else 1
 * @param distribution The size distribution
 * @param accuracy, nodes_per_degree As for solve
 */
Average average(Particle                particle,
                Numeric                 radius_factor,
                const SizeDistribution& distribution,
                Numeric                 accuracy,
                Index                   nodes_per_degree);
}  // namespace tmatrix::solver
