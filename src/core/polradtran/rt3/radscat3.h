#pragma once

#include <matpack.h>
#include <rtepack.h>

#include "rt3_fft.h"
#include "rt3_workdata.h"

/* The subroutines of 3rdparty/polradtran/radscat3.f, ported to C++ one at a
   time.  Each follows its Fortran step by step and calls the same
   subroutines.  None of them calls Fortran or keeps static state (the FFT
   is in rt3_fft.h).

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The phase matrix in the scattering plane is rtepack's
   CompactPlanarMuelmat (its six elements by name, not in COEF's order).
   Rotated into the meridional planes, the 4 x 4 phase matrices are
   Muelmat (element (r, c) of the Fortran matrix is [r - 1, c - 1]), their
   arrays over the azimuths, the modes and the streams MuelmatVector to
   MuelmatTensor4, and the direct vectors Stokvec.  MATRIX_SYMMETRY, which
   negates the off-diagonal 2 x 2 blocks, is rtepack::mirror.  As in the Fortran, only the leading
   nstokes x nstokes of the phase matrices (and nstokes of the vectors) are
   transformed; the rest is 0.  The counts are not passed; they are the
   extents of the arrays. */
namespace polradtran::rt3 {
/** Coefficient l of a delta-M scaled Legendre series, as GET_SCAT_SET
 * scales it: the forward peak, f times the identity, is removed from the
 * normalised coefficient c / (2 l + 1) and the rest renormalised,
 * (2 l + 1) (c / (2 l + 1) - f id) / (1 - f).  So the diagonal elements
 * (F11, F22, F33, F44) lose f and F12 and F34 are only renormalised.
 */
CompactPlanarMuelmat delta_m_scaled(const CompactPlanarMuelmat& c, Index l, Numeric f);

/** GET_SCAT_SET: one scattering set as RADTRAN uses it (READ_SCAT_FILE of
 * Evans' RT3, which read it from a file).  It returns the degree of the
 * Legendre series, the extinction, the scattering coefficient and the
 * Legendre coefficients of the phase matrix (each including the factor
 * 2 l + 1).
 *
 * With delta_m, the extinction, the single scattering albedo and the
 * series are delta-M scaled (delta_m_scaled) with
 * f = coefin[M].F11() / (2 M + 1), M = 2 nummu (f = 0 for a shorter
 * series), and the series is truncated to degree M - 1.  nummu counts all
 * of RADTRAN's angles, the extra ones included.
 *
 *   coefin  [nlegin + 1]                    COEFIN(6, NLEGIN+1), the set's series
 *   coef    [max(nlegin + 1, 2 nummu) or more]
 *                                           COEF(6, *), output: the coefficients to
 *                                           max(nlegin + 1, 2 nummu) are written,
 *                                           zero beyond the series
 *
 * Throws where the scaling divides by zero: delta_m with an extinction
 * that is not positive, 1 - f = 0 or 1 - albedo f = 0.
 */
void get_scat_set(bool                                delta_m,
                  Index                               nummu,
                  CompactPlanarMuelmatConstVectorView coefin,
                  Numeric                             extin,
                  Numeric                             scatin,
                  Index&                              nlegen,
                  CompactPlanarMuelmatVectorView      coef,
                  Numeric&                            extinction,
                  Numeric&                            scatter);

/** GET_SCATTERING: azimuth mode `mode` of the scattering matrix of one
 * set, from its part of SCATBUF (of scattering).  The scattering matrix is
 * four matrices: P++, P+-, P-+ and P--, where + and - are the signs of the
 * incoming and outgoing quadrature angles; P++ and P+- come from scatbuf,
 * and SCATTER_SYMMETRY makes P-- and P-+ from them.  It includes the
 * 1 / (4 pi) factor, the quadrature weights and the Fourier basis
 * constants, as well as the phase matrix.  nstokes is the extent of
 * scatter_matrix; the leading nstokes x nstokes of each Mueller matrix is
 * read.
 *
 *   scatbuf         [aziorder + 1, 2, nummu, nummu]       the set's part of SCATBUF (see scattering)
 *   scatter_matrix  [4, nummu, nstokes, nummu, nstokes]   SCATTER_MATRIX(NSTOKES, NUMMU, NSTOKES, NUMMU, 4),
 *                                                         output
 */
void get_scattering(Index mode, MuelmatConstTensor4View scatbuf, Tensor5View scatter_matrix);

/** CHECK_NORM (RT3_CHECK_NORM): checks the normalization of the scattering
 * matrix by integrating its I-I term over the outgoing directions for each
 * incoming direction (with non-zero weight), both hemispheres.  Throws if
 * a sum differs from 1 by more than 1e-7 (or is NaN), where the Fortran
 * stopped: either the first term of the Legendre series is not 1, or the
 * series has too many terms for the number of quadrature angles.
 *
 *   quad_weights    [nummu]                               QUAD_WEIGHTS
 *   scatter_matrix  [4, nummu, nstokes, nummu, nstokes]   SCATTER_MATRIX (of get_scattering)
 */
void check_norm(ConstVectorView quad_weights, ConstTensor5View scatter_matrix);

/** SCATTER_SYMMETRY: P-- from P++ and P-+ from P+- by the symmetry of the
 * scattering matrix.  For randomly oriented particles with a plane of
 * symmetry the diagonal 2 by 2 blocks of the Stokes parameters stay the
 * same, while the off-diagonal 2 by 2 blocks change sign under negation of
 * mu1 and mu2.
 *
 *   scat  [4, nummu, nstokes, nummu, nstokes]  SCAT(NSTOKES, NUMMU, NSTOKES, NUMMU, 4);
 *                                              parts 2 and 3 are output
 */
void scatter_symmetry(Tensor5View scat);

/** GET_DIRECT: azimuth mode `mode` of the direct (solar) pseudo-source
 * vectors of one set, from its part of DIRECTBUF (of direct_scattering):
 * the leading nstokes of each Stokes vector, nstokes the extent of
 * direct_vector.
 *
 *   directbuf      [aziorder + 1, 2, nummu]  the set's part of DIRECTBUF
 *   direct_vector  [2, nummu, nstokes]       DIRECT_VECTOR(NSTOKES, NUMMU, 2), output
 */
void get_direct(Index mode, StokvecConstTensor3View directbuf, Tensor3View direct_vector);

/** SUM_LEGENDRE: the phase matrix in the scattering plane at x, the cosine
 * of the scattering angle, from its Legendre series,
 * F(x) = sum_l coef[l] P_l(x), for randomly oriented particles with a plane
 * of symmetry.  For one Stokes parameter only F11 is summed; the rest is 0.
 * This is
 * NUMBER_SUMS's choice where it matters: it skipped F34 for 2 and 3 Stokes
 * parameters, where the rotation does not mix it into the leading 3 x 3,
 * and took F22 and F44 from F11 and F33 when the series were equal.
 *
 * The Legendre polynomials are ARTS's (Legendre::legendre_polynomials, the
 * generator of every Legendre series in ARTS), made once for all the
 * series in work.legendre_p, which it sizes.  x is clamped to [-1, 1]: as
 * the cosine of a scattering angle computed from the directions, it can
 * round to just outside.
 *
 *   coef  [nlegen + 1]  COEF(6, NLEGEN+1)
 */
CompactPlanarMuelmat sum_legendre(CompactPlanarMuelmatConstVectorView coef,
                                  Numeric                             x,
                                  Index                               nstokes,
                                  rt3_workdata&                       work);

/** ROTATE_PHASE_MATRIX: the phase matrix in the meridional planes of the
 * two directions, from that in the scattering plane (of sum_legendre): its
 * polarization basis is rotated from the incident meridional plane into the
 * scattering plane and from there to the outgoing meridional plane,
 * L(2 sigma2) phase_matrix L(2 sigma1), with L rtepack::stokes_rotation,
 * by rtepack::rotated.
 * mu1 is the incoming direction, mu2 the outgoing one, delphi the azimuth
 * between them and cos_scat the cosine of the scattering angle.  In
 * forward and backward scattering (sin of the scattering angle 0) the
 * rotation is fixed.
 */
Muelmat rotate_phase_matrix(
    const CompactPlanarMuelmat& phase_matrix, Numeric mu1, Numeric mu2, Numeric delphi, Numeric cos_scat);

//! FOURIER_BASIS's DIRECTION: into azimuth (real) space (negative) or
//! into the Fourier basis (positive)
enum class fourier_direction { to_real, to_basis };

/** FOURIER_BASIS: converts a vector between the Fourier basis and azimuth
 * space.  The basis functions are 1, cos(x), cos(2x), ..., cos(Mx), sin(x),
 * sin(2x), ..., sin(Mx), M = order; the basis has 2 order + 1 elements, or
 * order + 1 for even functions (no sines).  The real space is numpts
 * azimuths x = 2 pi k / numpts.  Orders from numpts / 2 on are not resolved:
 * they are 0 in the basis and ignored going to real space.
 *
 *   basis_vector  [order + 1 or 2 order + 1]  BASIS_VECTOR, input to_real, output to_basis
 *   real_vector   [numpts]                    REAL_VECTOR, output to_real, input to_basis
 *                                             (then overwritten when order > 0)
 *
 * With order > 0 it transforms with fft1dr, so numpts must be a power of
 * two of at least 2 (throws otherwise), whose state fft keeps.
 */
void fourier_basis(
    Index order, fourier_direction direction, VectorView basis_vector, VectorView real_vector, fft_workdata& fft);

/** FOURIER_MATRIX: the azimuth Fourier modes of each element of the
 * leading nstokes x nstokes of a phase matrix sampled at numpts azimuths
 * delphi = 2 pi k / numpts, by FOURIER_BASIS: modes 1, cos(m delphi) and
 * sin(m delphi) for m = 1 to aziorder, in that order.  The other elements
 * of the modes are 0.
 *
 *   real_matrix   [numpts]            REAL_MATRIX(4, 4, NUMPTS)
 *   basis_matrix  [2 aziorder + 1]    BASIS_MATRIX(4, 4, 2*AZIORDER+1), output
 *
 * With aziorder > 0 fourier_basis transforms with fft1dr, so numpts must
 * be a power of two.  REAL_VECTOR and BASIS_VECTOR are work's real_vector
 * and basis_vector, which it sizes, and work.fft is the state of the FFT.
 */
void fourier_matrix(MuelmatConstVectorView real_matrix,
                    MuelmatVectorView      basis_matrix,
                    Index                  nstokes,
                    rt3_workdata&          work);

/** COMBINE_PHASE_MODES: azimuth mode m of a phase matrix, OUT_MATRIX, from
 * its Fourier modes (of fourier_matrix), times the quadrature factor tmp.
 * The diagonal 2 by 2 blocks are cosine modes, the off-diagonal ones sine
 * modes: for m = 0, tmp times the mean and 0; for m > 0, tmp / 2 times the
 * cos(m delphi) mode, and -tmp / 2 and +tmp / 2 times the sin(m delphi)
 * mode in the lower left and upper right blocks (SINFLAG -1 and +1).
 *
 *   basis_matrix  [2 aziorder + 1]  BASIS_MATRIX(4, 4, 2*AZIORDER+1)
 */
Muelmat combine_phase_modes(Index m, Numeric tmp, MuelmatConstVectorView basis_matrix);

/** DIRECT_SCATTERING: the direct (solar) pseudo-source vectors of one
 * scattering set for every azimuth mode.  The direct vector is the
 * integral of the phase matrix times the delta function in the direction
 * of the sun: the phase matrix with the incoming direction set to the
 * sun's (direct_mu, the direct beam going down at azimuth 0) times the
 * Stokes vector {1, 0, 0, 0}, i.e. its first column.  The phase matrices
 * are made as in SCATTERING, at numpts azimuths delphi = -2 pi k / numpts;
 * the vector of mode m holds the cos(m delphi) modes of I and Q and the
 * sin(m delphi) modes of U and V (0 for m = 0).  No quadrature weight.
 * Only the leading nstokes of each vector are made; the rest is 0.
 *
 *   mu_values      [nummu]                  MU_VALUES
 *   legendre_coef  [numlegendre + 1]        LEGENDRE_COEF(6, NUMLEGENDRE+1)
 *   directbuf      [aziorder + 1, 2, nummu] output, the set's part of DIRECTBUF: [m, l, j] for the
 *                                           outgoing mu = +-mu_values[j] (+ for l = 0)
 *
 * numpts is that of SCATTERING.  Throws for a direct_mu outside [-1, 1].
 * The Fortran's limits of 512 azimuths and 2 aziorder + 1 <= 512 (FFT1DR
 * and its buffers) are not the port's.  SCAT_MATRIX and BASIS_MATRIX are
 * work's scat_matrix and basis_matrix, which it sizes; work.fft is the
 * state of the FFT (fft1dr).
 */
void direct_scattering(ConstVectorView                     mu_values,
                       CompactPlanarMuelmatConstVectorView legendre_coef,
                       Numeric                             direct_mu,
                       Index                               nstokes,
                       StokvecTensor3View                  directbuf,
                       rt3_workdata&                       work);

/** SCATTERING: the polarized scattering matrices of one scattering set for
 * every azimuth mode.  For each pair of quadrature angles (incoming and
 * outgoing) it evaluates the single scattering phase matrix for many delta
 * phi's: it sums the Legendre series at the scattering angle and rotates
 * the polarization reference from the scattering plane to the meridional
 * planes.  A Fourier transform then turns the phase matrices from phi
 * space to azimuth mode space.  The scattering matrices include the
 * quadrature weights for integrating (w / 2) and the Fourier basis
 * constants.  Only the leading nstokes x nstokes of each matrix is made;
 * the rest is 0.
 *
 *   mu_values      [nummu]                         MU_VALUES
 *   quad_weights   [nummu]                         QUAD_WEIGHTS
 *   legendre_coef  [numlegendre + 1]               LEGENDRE_COEF(6, NUMLEGENDRE+1)
 *   scatbuf        [aziorder + 1, 2, nummu, nummu] output, the set's part of SCATBUF: [m, l, j1, j2] is
 *                                                  the matrix of record
 *                                                  M*2*NUMMU**2 + (L-1)*NUMMU**2 + (J1-1)*NUMMU + J2,
 *                                                  from incoming mu = +-mu_values[j1] (+ for l = 0) to
 *                                                  outgoing mu_values[j2]
 *
 * With aziorder > 0 the phase matrices are sampled at
 * NUMPTS = 2 * 2^int(log2(numlegendre + 4) + 1) azimuths, else at
 * 2 int((numlegendre + 1) / 2) + 4.  The Fortran's limits of 512
 * azimuths with aziorder > 0 (FFT1DR) and 1024 for NUMPTS and
 * 2 aziorder + 1 (FOURIER_MATRIX's buffers) are not the port's.
 * SCAT_MATRIX and BASIS_MATRIX are work's scat_matrix and basis_matrix,
 * which it sizes; work.fft is the state of the FFT (fft1dr).
 */
void scattering(ConstVectorView                     mu_values,
                ConstVectorView                     quad_weights,
                CompactPlanarMuelmatConstVectorView legendre_coef,
                Index                               nstokes,
                MuelmatTensor4View                  scatbuf,
                rt3_workdata&                       work);
}  // namespace polradtran::rt3
