#pragma once

#include <matpack.h>

#include "rt3_fft.h"
#include "rt3_workdata.h"

/* The subroutines of 3rdparty/polradtran/radscat3.f, ported to C++ one at a
   time.  Each follows its Fortran step by step and calls the same
   subroutines.  None of them calls Fortran or keeps static state (the FFT
   is in rt3_fft.h).

   A Fortran array A(d1, ..., dk) is the row-major matpack array
   [dk, ..., d1].  The counts are not passed; they are the extents of the
   arrays. */
namespace rt3 {
//! NUMBER_SUMS's DOSUM and SUMCASES
using IndexVector6  = matpack::cdata_t<Index, 6>;
using IndexMatrix56 = matpack::cdata_t<Index, 5, 6>;

/** GET_SCAT_SET: one scattering set as RADTRAN uses it (READ_SCAT_FILE of
 * Evans' RT3, which read it from a file).  It returns the degree of the
 * Legendre series, the extinction, the scattering coefficient and the
 * Legendre coefficients of the six unique elements of the phase matrix
 * (F11, F12, F33, F34, F22, F44, each including the factor 2 l + 1).
 *
 * With delta_m, the extinction, the single scattering albedo and the
 * series are delta-M scaled with f = coefin[M, 0] / (2 M + 1), M = 2 nummu
 * (f = 0 for a shorter series), and the series is truncated to degree
 * M - 1: the diagonal elements (F11, F33, F22, F44) become
 * (2 l + 1) (c / (2 l + 1) - f) / (1 - f) and the others
 * (2 l + 1) (c / (2 l + 1)) / (1 - f).  nummu counts all of RADTRAN's
 * angles, the extra ones included.
 *
 *   coefin  [nlegin + 1, 6]                    COEFIN(6, NLEGIN+1), the set's series
 *   coef    [max(nlegin + 1, 2 nummu) or more, 6]
 *                                              COEF(6, *), output: the rows to
 *                                              max(nlegin + 1, 2 nummu) are written,
 *                                              zero beyond the series
 *
 * Throws where the scaling divides by zero: delta_m with an extinction
 * that is not positive, 1 - f = 0 or 1 - albedo f = 0.
 */
void get_scat_set(bool            delta_m,
                  Index           nummu,
                  ConstMatrixView coefin,
                  Numeric         extin,
                  Numeric         scatin,
                  Index&          nlegen,
                  MatrixView      coef,
                  Numeric&        extinction,
                  Numeric&        scatter);

/** GET_SCATTERING: azimuth mode `mode` of the scattering matrix of one
 * set, from its part of SCATBUF (of scattering).  The scattering matrix is
 * four matrices: P++, P+-, P-+ and P--, where + and - are the signs of the
 * incoming and outgoing quadrature angles; P++ and P+- come from scatbuf,
 * and SCATTER_SYMMETRY makes P-- and P-+ from them.  It includes the
 * 1 / (4 pi) factor, the quadrature weights and the Fourier basis
 * constants, as well as the phase matrix.
 *
 *   scatbuf         [aziorder + 1, 2, nummu, nummu, nstokes, nstokes]
 *                   the set's part of SCATBUF (see scattering)
 *   scatter_matrix  [4, nummu, nstokes, nummu, nstokes]
 *                   SCATTER_MATRIX(NSTOKES, NUMMU, NSTOKES, NUMMU, 4), output
 */
void get_scattering(Index mode, ConstTensor6View scatbuf, Tensor5View scatter_matrix);

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
 * vectors of one set, from its part of DIRECTBUF (of direct_scattering).
 *
 *   directbuf      [aziorder + 1, 2, nummu, nstokes]  the set's part of DIRECTBUF
 *   direct_vector  [2, nummu, nstokes]                DIRECT_VECTOR(NSTOKES, NUMMU, 2), output
 */
void get_direct(Index mode, ConstTensor4View directbuf, Tensor3View direct_vector);

/** NUMBER_SUMS: which of the six Legendre series of coef (F11, F12, F33,
 * F34, F22, F44) SUM_LEGENDRE must sum, 1 or 0 each, for nstokes Stokes
 * parameters.  The phase matrix is Rayleigh-like (F34 = 0, F22 = F11 and
 * F44 = F33), Mie-like (F34 non-zero) or general (F22 != F11 or
 * F44 != F33), and SUM_LEGENDRE takes F22 and F44 from F11 and F33 when it
 * does not sum them:
 *
 *   nstokes 1:     F11
 *   nstokes 2, 3:  F11, F12, F33, and F22 if general
 *   nstokes 4:     F11, F12, F33, and F34 if Mie, or all six if general
 *
 *   coef  [nlegen + 1, 6]  COEF(6, NLEGEN+1)
 */
IndexVector6 number_sums(Index nstokes, ConstMatrixView coef);

/** SUM_LEGENDRE: sums the Legendre series of each element of the phase
 * matrix at x, the cosine of the scattering angle.  coef holds the six
 * independent series of randomly oriented particles with a plane of
 * symmetry (F11, F12, F33, F34, F22, F44), and dosum (of number_sums)
 * selects those to sum; the others are 0.  The series go to the phase
 * matrix elements (1,1), (1,2), (3,3), (3,4), (2,2) and (4,4); then
 * (2,1) = (1,2), (4,3) = -(3,4), and (2,2) = (1,1) and (4,4) = (3,3) when
 * F22 and F44 are not summed.  The other elements are not written.
 *
 *   coef          [nlegen + 1, 6]  COEF(6, NLEGEN+1)
 *   phase_matrix  [4, 4]           PHASE_MATRIX(4, 4), output: element (r, c)
 *                                  is phase_matrix[c - 1, r - 1]
 */
void sum_legendre(ConstMatrixView coef, Numeric x, const IndexVector6& dosum, MatrixView phase_matrix);

/** ROTATE_PHASE_MATRIX: rotates the polarization basis of the phase matrix
 * from the incident plane into the scattering plane and from the
 * scattering plane to the outgoing plane, for randomly oriented particles
 * with a plane of symmetry (six unique elements).  mu1 is the incoming
 * direction, mu2 the outgoing one, delphi the azimuth between them and
 * cos_scat the cosine of the scattering angle.  In forward and backward
 * scattering (sin of the scattering angle 0) the rotation is fixed.
 *
 *   phase_matrix1  [4, 4]              PHASE_MATRIX1(4, 4), of sum_legendre
 *   phase_matrix2  [nstokes, nstokes]  the leading part of PHASE_MATRIX2(4, 4),
 *                                      output
 *
 * Element (r, c) of either is [c - 1, r - 1].
 */
void rotate_phase_matrix(ConstMatrixView   phase_matrix1,
                         Numeric           mu1,
                         Numeric           mu2,
                         Numeric           delphi,
                         Numeric           cos_scat,
                         StridedMatrixView phase_matrix2);

/** MATRIX_SYMMETRY: a symmetry operation on a phase matrix, equivalent to
 * negating mu and mu' or negating phi' - phi: the diagonal 2 by 2 blocks
 * are copied, the off-diagonal ones negated.  matrix1 and matrix2 may be
 * the same matrix (SCATTERING does so at delphi = pi).
 *
 *   matrix1  [nstokes, nstokes]  the leading part of MATRIX1(4, 4)
 *   matrix2  [nstokes, nstokes]  the leading part of MATRIX2(4, 4), output
 */
void matrix_symmetry(StridedConstMatrixView matrix1, StridedMatrixView matrix2);

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
 * two from 2 to 512 (throws otherwise), whose state fft keeps.
 */
void fourier_basis(
    Index order, fourier_direction direction, VectorView basis_vector, VectorView real_vector, fft_workdata& fft);

/** FOURIER_MATRIX: the azimuth Fourier modes of each element of a phase
 * matrix sampled at numpts azimuths delphi = 2 pi k / numpts, by
 * FOURIER_BASIS: modes 1, cos(m delphi) and sin(m delphi) for m = 1 to
 * aziorder, in that order.
 *
 *   real_matrix   [numpts, nstokes, nstokes]           REAL_MATRIX(4, 4, NUMPTS)
 *   basis_matrix  [2 aziorder + 1, nstokes, nstokes]   BASIS_MATRIX(4, 4, 2*AZIORDER+1),
 *                                                      output
 *
 * (the leading nstokes x nstokes of the Fortran 4 x 4).  With aziorder > 0
 * fourier_basis transforms with fft1dr, so numpts must be a power of two,
 * at most 512.  REAL_VECTOR and BASIS_VECTOR are work's real_vector and
 * basis_vector, which it sizes, and work.fft is the state of the FFT; the
 * matrices must not be work's real_vector or basis_vector.
 */
void fourier_matrix(StridedConstTensor3View real_matrix, StridedTensor3View basis_matrix, rt3_workdata& work);

/** COMBINE_PHASE_MODES: azimuth mode m of a phase matrix, from its Fourier
 * modes (of fourier_matrix), times the quadrature factor tmp.  The
 * diagonal 2 by 2 blocks are cosine modes, the off-diagonal ones sine
 * modes: for m = 0, tmp times the mean and 0; for m > 0, tmp / 2 times the
 * cos(m delphi) mode, and -tmp / 2 and +tmp / 2 times the sin(m delphi)
 * mode in the upper right and lower left blocks of the Fortran matrix
 * (SINFLAG -1 and +1).
 *
 *   basis_matrix  [2 aziorder + 1, nstokes, nstokes]  BASIS_MATRIX(4, 4, 2*AZIORDER+1)
 *   out_matrix    [nstokes, nstokes]                  OUT_MATRIX, output: element
 *                                                     (r, c) is [c - 1, r - 1]
 */
void combine_phase_modes(Index m, Numeric tmp, StridedConstTensor3View basis_matrix, StridedMatrixView out_matrix);

/** DIRECT_SCATTERING: the direct (solar) pseudo-source vectors of one
 * scattering set for every azimuth mode.  The direct vector is the
 * integral of the phase matrix times the delta function in the direction
 * of the sun: the phase matrix with the incoming direction set to the
 * sun's (direct_mu, the direct beam going down at azimuth 0) times the
 * Stokes vector {1, 0, 0, 0}, i.e. its first column.  The phase matrices
 * are made as in SCATTERING, at numpts azimuths delphi = -2 pi k / numpts;
 * the vector of mode m holds the cos(m delphi) modes of I and Q and the
 * sin(m delphi) modes of U and V (0 for m = 0).  No quadrature weight.
 *
 *   mu_values      [nummu]                         MU_VALUES
 *   legendre_coef  [numlegendre + 1, 6]            LEGENDRE_COEF(6, NUMLEGENDRE+1)
 *   directbuf      [aziorder + 1, 2, nummu, nstokes]
 *                  output, the set's part of DIRECTBUF: [m, l, j] for the
 *                  outgoing mu = +-mu_values[j] (+ for l = 0)
 *
 * numpts is that of SCATTERING.  Throws for a direct_mu outside [-1, 1], and
 * where the Fortran would stop or overflow: numpts > 512 with aziorder > 0
 * (FFT1DR), numpts or 2 aziorder + 1 above 512 (the buffers of the Fortran
 * DIRECT_SCATTERING, which the port does not have; kept as RT3's limit).
 * SCAT_MATRIX and BASIS_MATRIX are work's scat_matrix and basis_matrix,
 * which it sizes; work.fft is the state of the FFT (fft1dr).
 */
void direct_scattering(ConstVectorView mu_values,
                       ConstMatrixView legendre_coef,
                       Numeric         direct_mu,
                       Tensor4View     directbuf,
                       rt3_workdata&   work);

/** SCATTERING: the polarized scattering matrices of one scattering set for
 * every azimuth mode.  For each pair of quadrature angles (incoming and
 * outgoing) it evaluates the single scattering phase matrix for many delta
 * phi's: it sums the Legendre series at the scattering angle and rotates
 * the polarization reference from the scattering plane to the meridional
 * planes.  A Fourier transform then turns the phase matrices from phi
 * space to azimuth mode space.  The scattering matrices include the
 * quadrature weights for integrating (w / 2) and the Fourier basis
 * constants.
 *
 *   mu_values      [nummu]                MU_VALUES
 *   quad_weights   [nummu]                QUAD_WEIGHTS
 *   legendre_coef  [numlegendre + 1, 6]   LEGENDRE_COEF(6, NUMLEGENDRE+1)
 *   scatbuf        [aziorder + 1, 2, nummu, nummu, nstokes, nstokes]
 *                  output, the set's part of SCATBUF: [m, l, j1, j2] is the
 *                  column-major nstokes x nstokes matrix of record
 *                  M*2*NUMMU**2 + (L-1)*NUMMU**2 + (J1-1)*NUMMU + J2, from
 *                  incoming mu = +-mu_values[j1] (+ for l = 0) to outgoing
 *                  mu_values[j2]
 *
 * With aziorder > 0 the phase matrices are sampled at
 * NUMPTS = 2 * 2^int(log2(numlegendre + 4) + 1) azimuths, else at
 * 2 int((numlegendre + 1) / 2) + 4.  Throws where the Fortran would stop or
 * overflow: NUMPTS > 512 with aziorder > 0 (FFT1DR), NUMPTS or
 * 2 aziorder + 1 above 1024 (the buffers of the Fortran FOURIER_MATRIX,
 * which the port does not have; kept as RT3's limit).  SCAT_MATRIX and
 * BASIS_MATRIX are work's scat_matrix and basis_matrix, which it sizes;
 * work.fft is the state of the FFT (fft1dr).
 */
void scattering(ConstVectorView mu_values,
                ConstVectorView quad_weights,
                ConstMatrixView legendre_coef,
                Tensor6View     scatbuf,
                rt3_workdata&   work);
}  // namespace rt3
