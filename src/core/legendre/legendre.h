#pragma once

#include <configtypes.h>
#include <matpack.h>

namespace Legendre {
/** Returns the Schmidt normalized Lagrange polynominal and its derivative
 *
 * The derivative is undefined if theta is 0 or pi, so the function cannot
 * be called with these values
 *
 * This function is taken from the implementation by at https://github.com/klaundal/ppigrf
 * (2024-06-13).  It implements the calculations step-by-step, updating both main and derivative
 * based on previous values.
 *
 * @param[in] theta Colatitude in radians
 * @param[in] nmax Max number of n
 * @return The pair of (nmax+1) x (nmax+1) lower triangular matrices, order: main and derivative
 */
std::pair<TriangularMatrix, TriangularMatrix> schmidt(const Numeric theta, const Index nmax);

/** As schmidt(theta, nmax), into caller-sized storage
 *
 * @param[out] P  The (nmax+1) x (nmax+1) main matrix, zero for m > n
 * @param[out] dP The derivative, of the same shape
 * @param[in] theta Colatitude in radians
 */
void schmidt(MatrixView P, MatrixView dP, const Numeric theta);

/** Computes the spherical field
 * 
 * If latitude is beyond a limit, the longitude is set to zero
 * and the latitude at the limit away from the pole is used. For these
 * position, the absolute of the horizontal component is kept around in
 * the south-facing component and the east-facing component is set to zero
 * 
 * The latitude limit is defined in the ColatitudeConversion struct.
 * 
 * The longitude is constrained to the range [-180, 180) to have consistent
 * behavior at the daytime border.
 * 
 * This function is taken from the implementation by Isabela de Oliveira Martins
 * at https://github.com/de-oliveira/IsabelaFunctions/blob/master/IsabelaFunctions/fieldmodel.py
 * (2021-05-06).
 *
 * The n = 0 coefficients, as there are no magnetic monopoles, and the h with
 * m = 0 do not affect the field and are not used, so they may be NaN.  The
 * m = 0 coefficients do not affect the longitudinal component.
 * 
 * @param[in] g A N x N lower triangular matrix of g-coefficients
 * @param[in] h A N x N lower triangular matrix of h-coefficients
 * @param[in] r0 The reference radius (spherical)
 * @param[in] pos The position [r, lat, lon] (spherical)
 * @return A spherical field {Br, Btheta, Bphi}
 */
Vector3 schmidt_fieldcalc(const ConstTriangularMatrixView& g,
                          const ConstTriangularMatrixView& h,
                          const Numeric                    r0,
                          const Vector3                    pos);

/** The derivative of the schmidt_fieldcalc function wrt g and h
 *
 * The output is 3 pairs of N x N lower triangular matrices.  The index of
 * the array is 0 for the radial component, 1 for the colatitudinal component,
 * and 2 for the longitudinal component.  The first of the pair is for the
 * g-coefficients and the second for the h-coefficients.  The indices of the
 * matrices are the n and m indices of the spherical harmonics.
 *
 * The n = 0 derivatives are zero, as schmidt_fieldcalc ignores those coefficients.
 *
 * @param[in] N The number of rows of the coefficient matrices, the maximum degree plus one
 * @param[in] r0 The reference radius (spherical)
 * @param[in] pos The position [r, lat, lon] (spherical)
 * @return The derivatives
 */
std::array<std::pair<TriangularMatrix, TriangularMatrix>, 3> dschmidt_fieldcalc(const Size    N,
                                                                                const Numeric r0,
                                                                                const Vector3 pos);

/** Computes sum (s[i] P_i(x)) for all s [first is for P_0, second is for P_1, ...]
  * 
  * @param[in] s The coefficients
  * @param[in] x The x values
  * @return The sum
  */
Numeric legendre_sum(const ConstVectorView& s, const Numeric& x);

/** Computes P_n(x)
  * 
  * @param[in] n The order
  * @param[in] x The value
  * @return The polynomial
  */
Numeric legendre(Index n, Numeric x);

/** Computes P_0(x), P_1(x), ..., P_{p.size() - 1}(x)
  * 
  * @param[out] p The polynomials
  * @param[in] x The value, in [-1, 1]
  */
void legendre_polynomials(VectorView p, Numeric x);

/** Wigner's d^l_mn(theta) of the degrees l = 0, 1, ..., d.size() - 1, for integer orders m and n
 *
 * d^l_mn(theta) = <l m| exp(-i theta J_y) |l n>, 0 for l < max(|m|, |n|).  It holds
 *   d^l_00(theta) = P_l(cos(theta)), the Legendre polynomials,
 *   d^l_0m(theta) = sqrt((l - m)! / (l + m)!) P^m_l(cos(theta)), m >= 0, with P^m_l the associated Legendre
 *     functions without the Condon-Shortley phase (the Schmidt functions of schmidt() over sqrt(2) for m > 0),
 * and the generalized spherical functions of polarized scattering (Hovenier et al. 2004) are
 *   P^l_22 = d^l_22, P^l_2-2 = d^l_2-2 and P^l_02 = -d^l_02.
 *
 * The functions are those of the three-term recurrence in l,
 *   l sqrt(((l + 1)^2 - m^2) ((l + 1)^2 - n^2)) d^(l+1)
 *     = (2 l + 1) (l (l + 1) cos(theta) - m n) d^l - (l + 1) sqrt((l^2 - m^2) (l^2 - n^2)) d^(l-1),
 * from l0 = max(|m|, |n|), where for l0 = M >= |N|
 *   d^M_MN = (-1)^(M - N) sqrt((2 M)! / ((M + N)! (M - N)!)) cos^(M + N)(theta / 2) sin^(M - N)(theta / 2),
 * and d^l_mn = (-1)^(m - n) d^l_nm = d^l_-n-m gives the other orders.  The angle, not its cosine, keeps the
 * relative accuracy of the functions that vanish at 0 or pi.
 *
 * @param[out] d The functions of the degrees 0, 1, ..., d.size() - 1
 * @param[in] m The first order
 * @param[in] n The second order
 * @param[in] theta The angle [rad]
 */
void wigner_d(VectorView d, Index m, Index n, Numeric theta);

/** Computes n!
  * 
  * @param[in] n The index
  * @return The factorial value
  */
Numeric factorial(Index n);

/** Computes P^m_l(x)
  * 
  * @param[in] l The degree
  * @param[in] m The order
  * @param[in] x The value
  * @return The associated polynomial
  */
Numeric assoc_legendre(Index l, Index m, Numeric x);

/** Computes the ratio of gamma functions
  * 
  * @param[in] x The first argument
  * @param[in] y The second argument
  * @return The ratio
  */
Numeric tgamma_ratio(Numeric x, Numeric y);

/** Computes the double Gauss Legendre quadrature
 *
 * The degree of the single Gauss Legendre quadrature is the length of x and w,
 * which must be the same.
 * 
 * @param x The coordinates
 * @param w The weights
 */
void PositiveDoubleGaussLegendre(VectorView x, VectorView w);

/** Computes the Gauss Legendre quadrature
  * 
  * @param x The coordinates
  * @param w The weights
  */
void GaussLegendre(VectorView x, VectorView w);

/** Computes the positive part of the Gauss Legendre quadrature
  *
  * The degree of the full polynomial is twice the legnth of x and w,
  * which must be the same, but only the positive parts are kept.
  * 
  * @param x The coordinates
  * @param w The weights
  */
void PositiveGaussLegendre(VectorView x, VectorView w);
}  // namespace Legendre
