#include "radscat3.h"

#include <arts_constants.h>
#include <debug.h>
#include <legendre.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace rt3 {
void get_scat_set(bool            delta_m,
                  Index           nummu,
                  ConstMatrixView coefin,
                  Numeric         extin,
                  Numeric         scatin,
                  Index&          nlegen,
                  MatrixView      coef,
                  Numeric&        extinction,
                  Numeric&        scatter) {
  const Index nlegin = coefin.nrows() - 1;
  const Index nrows  = std::max(nlegin + 1, 2 * nummu);
  ARTS_USER_ERROR_IF(nummu < 1 or nlegin < 0 or coefin.ncols() != 6 or coef.ncols() != 6 or coef.nrows() < nrows,
                     "GET_SCAT_SET needs nummu >= 1, coefin [nlegin + 1, 6] and coef [max(nlegin + 1, 2 nummu) or "
                     "more, 6]; got nummu {}, coefin {:B,} and coef {:B,}",
                     nummu,
                     coefin.shape(),
                     coef.shape());

  // Copy the set where READ_SCAT_FILE read the file.  With delta_m the
  // scaling below reads coef up to l = 2 nummu even when nlegen + 1 is
  // smaller; those rows are zero.
  extinction                 = extin;
  scatter                    = scatin;
  nlegen                     = nlegin;
  coef[Range{0, nrows}]      = 0.0;
  coef[Range{0, nlegen + 1}] = coefin;

  if (delta_m) {
    const Index m = 2 * nummu;
    Numeric     f = 0.0;
    if (m + 1 <= nlegen + 1) f = coef[m, 0] / static_cast<Numeric>(2 * m + 1);
    ARTS_USER_ERROR_IF(
        not(extinction > 0.0), "Delta-M scaling divides by the extinction, which must be positive, got {}", extinction);
    ARTS_USER_ERROR_IF(1.0 - f == 0.0, "Delta-M scaling divides by 1 - f, with f = coef[{}, 0] / {} = 1", m, 2 * m + 1);
    Numeric albedo = scatter / extinction;
    ARTS_USER_ERROR_IF(
        1.0 - albedo * f == 0.0, "Delta-M scaling divides by 1 - albedo f, with albedo {} and f {}", albedo, f);
    extinction = (1.0 - albedo * f) * extinction;
    albedo     = (1.0 - f) * albedo / (1.0 - albedo * f);
    scatter    = albedo * extinction;
    nlegen     = m - 1;
    // Scale the diagonal and off-diagonal phase matrix elements differently
    for (Index l = 0; l <= nlegen; l++) {
      const auto k = static_cast<Numeric>(2 * l + 1);
      coef[l, 0]   = k * (coef[l, 0] / k - f) / (1.0 - f);
      coef[l, 1]   = k * (coef[l, 1] / k) / (1.0 - f);
      coef[l, 2]   = k * (coef[l, 2] / k - f) / (1.0 - f);
      coef[l, 3]   = k * (coef[l, 3] / k) / (1.0 - f);
      coef[l, 4]   = k * (coef[l, 4] / k - f) / (1.0 - f);
      coef[l, 5]   = k * (coef[l, 5] / k - f) / (1.0 - f);
    }
  }
}

IndexVector6 number_sums(Index nstokes, ConstMatrixView coef) {
  ARTS_USER_ERROR_IF(nstokes < 1 or nstokes > 4 or coef.nrows() < 1 or coef.ncols() != 6,
                     "NUMBER_SUMS needs 1 to 4 Stokes parameters and coef [nlegen + 1, 6]; got {} and {:B,}",
                     nstokes,
                     coef.shape());

  // SUMCASES: the series summed in each case, one case per row
  constexpr IndexMatrix56 sumcases{1, 0, 0, 0, 0, 0,  //
                                   1, 1, 1, 0, 0, 0,  //
                                   1, 1, 1, 1, 0, 0,  //
                                   1, 1, 1, 0, 1, 0,  //
                                   1, 1, 1, 1, 1, 1};

  enum class scattering_type { rayleigh, mie, general };
  scattering_type scat = scattering_type::rayleigh;
  if (stdr::any_of(coef[joker, 3], [](Numeric c) { return c != 0.0; })) scat = scattering_type::mie;
  if (not stdr::equal(coef[joker, 0], coef[joker, 4]) or not stdr::equal(coef[joker, 2], coef[joker, 5]))
    scat = scattering_type::general;

  Index sum_case;
  if (nstokes == 1) {
    sum_case = 1;
  } else if (nstokes <= 3) {
    sum_case = 2;
    if (scat == scattering_type::general) sum_case = 4;
  } else {
    sum_case = 2;
    if (scat == scattering_type::mie) sum_case = 3;
    if (scat == scattering_type::general) sum_case = 5;
  }
  return matpack::to<IndexVector6>(sumcases[sum_case - 1]);
}

void sum_legendre(
    ConstMatrixView coef, Numeric x, const IndexVector6& dosum, MatrixView phase_matrix, rt3_workdata& work) {
  const Index nlegen = coef.nrows() - 1;
  ARTS_USER_ERROR_IF(nlegen < 0 or coef.ncols() != 6 or phase_matrix.nrows() != 4 or phase_matrix.ncols() != 4,
                     "SUM_LEGENDRE needs coef [nlegen + 1, 6] and phase_matrix [4, 4]; got {:B,} and {:B,}",
                     coef.shape(),
                     phase_matrix.shape());

  // ROW and COL: the element of each series in the phase matrix, 0-based
  constexpr IndexVector6 row{0, 0, 2, 2, 1, 3}, col{0, 1, 2, 3, 1, 3};

  // The Legendre polynomials P_0(x) to P_nlegen(x), for all the series
  Vector& p = work.legendre_p.resize(nlegen + 1);
  Legendre::legendre_polynomials(p, std::clamp(x, -1.0, 1.0));

  // Sum the Legendre series
  for (Index i = 0; i < 6; i++) phase_matrix[col[i], row[i]] = dosum[i] == 1 ? dot(coef[joker, i], p) : 0.0;
  phase_matrix[0, 1] = phase_matrix[1, 0];
  phase_matrix[2, 3] = -phase_matrix[3, 2];
  if (dosum[4] == 0) phase_matrix[1, 1] = phase_matrix[0, 0];
  if (dosum[5] == 0) phase_matrix[3, 3] = phase_matrix[2, 2];
}

void rotate_phase_matrix(ConstMatrixView   phase_matrix1,
                         Numeric           mu1,
                         Numeric           mu2,
                         Numeric           delphi,
                         Numeric           cos_scat,
                         StridedMatrixView phase_matrix2) {
  const Index nstokes = phase_matrix2.nrows();
  ARTS_USER_ERROR_IF(phase_matrix1.nrows() != 4 or phase_matrix1.ncols() != 4 or nstokes < 1 or nstokes > 4 or
                         phase_matrix2.ncols() != nstokes,
                     "ROTATE_PHASE_MATRIX needs phase_matrix1 [4, 4] and phase_matrix2 [nstokes, nstokes] with 1 to 4 "
                     "Stokes parameters; got {:B,} and {:B,}",
                     phase_matrix1.shape(),
                     phase_matrix2.shape());

  // The Fortran matrices, element (r, c) at [r - 1, c - 1]
  const auto pm1 = matpack::transpose(phase_matrix1);
  auto       pm2 = matpack::transpose(phase_matrix2);

  const Numeric a1 = pm1[0, 0];
  pm2[0, 0]        = a1;
  if (nstokes == 1) return;

  const Numeric sin_scat   = std::sqrt(std::max(0.0, 1.0 - cos_scat * cos_scat));
  const Numeric sin_theta1 = std::sqrt(1.0 - mu1 * mu1);
  const Numeric sin_theta2 = std::sqrt(1.0 - mu2 * mu2);
  const Numeric sinphi     = std::sin(delphi);
  const Numeric cosphi     = std::cos(delphi);
  Numeric       sin1, sin2, cos1, cos2;
  if (sin_scat == 0.0) {
    sin1 = 0.0;
    sin2 = 0.0;
    cos1 = 1.0;
    cos2 = -1.0;
  } else {
    sin1 = sin_theta2 * sinphi / sin_scat;
    sin2 = sin_theta1 * sinphi / sin_scat;
    cos1 = (sin_theta1 * mu2 - sin_theta2 * mu1 * cosphi) / sin_scat;
    cos2 = (sin_theta2 * mu1 - sin_theta1 * mu2 * cosphi) / sin_scat;
  }
  const Numeric sin21 = 2.0 * sin1 * cos1;
  const Numeric cos21 = 1.0 - 2.0 * (sin1 * sin1);
  const Numeric sin22 = 2.0 * sin2 * cos2;
  const Numeric cos22 = 1.0 - 2.0 * (sin2 * sin2);

  // nstokes > 1
  const Numeric a2 = pm1[1, 1];
  const Numeric a3 = pm1[2, 2];
  const Numeric b1 = pm1[0, 1];
  pm2[0, 1]        = b1 * cos21;
  pm2[1, 0]        = b1 * cos22;
  pm2[1, 1]        = a2 * cos21 * cos22 - a3 * sin21 * sin22;
  if (nstokes > 2) {
    pm2[0, 2] = -(b1 * sin21);
    pm2[1, 2] = -a2 * sin21 * cos22 - a3 * cos21 * sin22;
    pm2[2, 0] = b1 * sin22;
    pm2[2, 1] = a2 * cos21 * sin22 + a3 * sin21 * cos22;
    pm2[2, 2] = -a2 * sin21 * sin22 + a3 * cos21 * cos22;
  }
  if (nstokes > 3) {
    const Numeric a4 = pm1[3, 3];
    const Numeric b2 = pm1[2, 3];
    pm2[0, 3]        = 0.0;
    pm2[1, 3]        = -(b2 * sin22);
    pm2[2, 3]        = b2 * cos22;
    pm2[3, 0]        = 0.0;
    pm2[3, 1]        = -(b2 * sin21);
    pm2[3, 2]        = -(b2 * cos21);
    pm2[3, 3]        = a4;
  }
}

void matrix_symmetry(StridedConstMatrixView matrix1, StridedMatrixView matrix2) {
  const Index nstokes = matrix2.nrows();
  ARTS_USER_ERROR_IF(nstokes < 1 or nstokes > 4 or matrix2.ncols() != nstokes or matrix1.shape() != matrix2.shape(),
                     "MATRIX_SYMMETRY needs matrix1 and matrix2 [nstokes, nstokes] with 1 to 4 Stokes parameters; got "
                     "{:B,} and {:B,}",
                     matrix1.shape(),
                     matrix2.shape());

  // Copy the diagonal 2 by 2 blocks and negate the off-diagonal ones
  matrix2 = matrix1;
  if (nstokes > 2) {
    const Range iq{0, 2}, uv{2, nstokes - 2};
    matrix2[iq, uv] *= -1.0;
    matrix2[uv, iq] *= -1.0;
  }
}

void fourier_basis(
    Index order, fourier_direction direction, VectorView basis_vector, VectorView real_vector, fft_workdata& fft) {
  const Index numbasis = basis_vector.size();
  const Index numpts   = real_vector.size();
  ARTS_USER_ERROR_IF(order < 0 or (numbasis != order + 1 and numbasis != 2 * order + 1) or numpts < 1,
                     "FOURIER_BASIS needs basis_vector [order + 1 or 2 order + 1] and real_vector [numpts]; got order "
                     "{}, {} and {}",
                     order,
                     numbasis,
                     numpts);

  // REAL_VECTOR(2*I+1) and REAL_VECTOR(2*I+2), I = 1 to BASISLEN, hold the
  // cosine and sine terms of order I
  const Index        basislen = std::min(order, numpts / 2 - 1);
  const StridedRange cosines{2, basislen, 2}, sines{3, basislen, 2};

  if (direction == fourier_direction::to_real) {
    if (order == 0) {
      real_vector = basis_vector[0];
    } else {
      real_vector    = 0.0;
      real_vector[0] = basis_vector[0];
      if (basislen > 0) {
        real_vector[cosines]  = basis_vector[Range{1, basislen}];
        real_vector[cosines] /= 2.0;
        if (numbasis > order + 1) {
          real_vector[sines]  = basis_vector[Range{order + 1, basislen}];
          real_vector[sines] /= 2.0;
        }
      }
      fft1dr(real_vector, fft_direction::inverse, fft);
    }

  } else {
    const auto n = static_cast<Numeric>(numpts);
    if (order == 0) {
      basis_vector[0] = sum(real_vector) / n;
    } else {
      fft1dr(real_vector, fft_direction::forward, fft);
      basis_vector[0] = real_vector[0] / n;
      if (basislen > 0) {
        basis_vector[Range{1, basislen}]  = real_vector[cosines];
        basis_vector[Range{1, basislen}] *= 2.0;
        basis_vector[Range{1, basislen}] /= n;
      }
      basis_vector[Range{basislen + 1, order - basislen}] = 0.0;
      if (numbasis > order + 1) {
        if (basislen > 0) {
          basis_vector[Range{order + 1, basislen}]  = real_vector[sines];
          basis_vector[Range{order + 1, basislen}] *= 2.0;
          basis_vector[Range{order + 1, basislen}] /= n;
        }
        basis_vector[Range{order + 1 + basislen, order - basislen}] = 0.0;
      }
    }
  }
}

void fourier_matrix(StridedConstTensor3View real_matrix, StridedTensor3View basis_matrix, rt3_workdata& work) {
  const Index numpts   = real_matrix.extent(0);
  const Index nstokes  = real_matrix.extent(1);
  const Index numazi   = basis_matrix.extent(0);
  const Index aziorder = (numazi - 1) / 2;
  ARTS_USER_ERROR_IF(numpts < 1 or nstokes < 1 or nstokes > 4 or real_matrix.extent(2) != nstokes or numazi % 2 != 1 or
                         basis_matrix.extent(1) != nstokes or basis_matrix.extent(2) != nstokes,
                     "FOURIER_MATRIX needs real_matrix [numpts, nstokes, nstokes] and basis_matrix "
                     "[2 aziorder + 1, nstokes, nstokes] with 1 to 4 Stokes parameters; got {:B,} and {:B,}",
                     real_matrix.shape(),
                     basis_matrix.shape());

  // REAL_VECTOR and BASIS_VECTOR; FFT1DR transforms REAL_VECTOR in place
  Vector& real_vector  = work.real_vector.resize(numpts);
  Vector& basis_vector = work.basis_vector.resize(numazi);
  for (Index i = 0; i < nstokes; i++) {
    for (Index j = 0; j < nstokes; j++) {
      real_vector = real_matrix[joker, j, i];
      fourier_basis(aziorder, fourier_direction::to_basis, basis_vector, real_vector, work.fft);
      basis_matrix[joker, j, i] = basis_vector;
    }
  }
}

void combine_phase_modes(Index m, Numeric tmp, StridedConstTensor3View basis_matrix, StridedMatrixView out_matrix) {
  const Index nstokes  = out_matrix.nrows();
  const Index aziorder = (basis_matrix.extent(0) - 1) / 2;
  ARTS_USER_ERROR_IF(nstokes < 1 or nstokes > 4 or out_matrix.ncols() != nstokes or basis_matrix.extent(0) % 2 != 1 or
                         basis_matrix.extent(1) != nstokes or basis_matrix.extent(2) != nstokes or m < 0 or
                         m > aziorder,
                     "COMBINE_PHASE_MODES needs basis_matrix [2 aziorder + 1, nstokes, nstokes] and out_matrix "
                     "[nstokes, nstokes] with 1 to 4 Stokes parameters, and 0 <= m <= aziorder; got {:B,}, {:B,} "
                     "and m {}",
                     basis_matrix.shape(),
                     out_matrix.shape(),
                     m);

  // SINFLAG: 0 in the diagonal 2 by 2 blocks (cosine modes), -1 in
  // [iq, uv] and +1 in [uv, iq] (sine modes)
  const Range iq{0, std::min<Index>(nstokes, 2)}, uv{2, nstokes - 2};

  if (m == 0) {
    out_matrix  = basis_matrix[0];
    out_matrix *= tmp;
    if (nstokes > 2) {
      out_matrix[iq, uv] = 0.0;
      out_matrix[uv, iq] = 0.0;
    }
  } else {
    // MC = M+1 and MS = M+1+AZIORDER, 0-based
    const Index   mc = m, ms = m + aziorder;
    const Numeric c = 0.5 * tmp;
    out_matrix      = basis_matrix[mc];
    if (nstokes > 2) {
      out_matrix[iq, uv]  = basis_matrix[ms, iq, uv];
      out_matrix[iq, uv] *= -1.0;
      out_matrix[uv, iq]  = basis_matrix[ms, uv, iq];
    }
    out_matrix *= c;
  }
}

void scattering(ConstVectorView mu_values,
                ConstVectorView quad_weights,
                ConstMatrixView legendre_coef,
                Tensor6View     scatbuf,
                rt3_workdata&   work) {
  using Constant::two_pi;

  const Index nummu       = mu_values.size();
  const Index aziorder    = scatbuf.extent(0) - 1;
  const Index nstokes     = scatbuf.extent(5);
  const Index numlegendre = legendre_coef.nrows() - 1;
  ARTS_USER_ERROR_IF(
      nummu < 1 or quad_weights.size() != mu_values.size() or numlegendre < 0 or legendre_coef.ncols() != 6 or
          nstokes < 1 or nstokes > 4 or aziorder < 0 or
          scatbuf.shape() != (std::array<Index, 6>{aziorder + 1, 2, nummu, nummu, nstokes, nstokes}),
      "SCATTERING needs mu_values and quad_weights [nummu], legendre_coef [numlegendre + 1, 6] and scatbuf "
      "[aziorder + 1, 2, nummu, nummu, nstokes, nstokes] with 1 to 4 Stokes parameters; got {}, {}, {:B,} and {:B,}",
      mu_values.size(),
      quad_weights.size(),
      legendre_coef.shape(),
      scatbuf.shape());

  Index numpts =
      2 * (Index{1} << static_cast<Index>(std::log(static_cast<Numeric>(numlegendre + 4)) / std::log(2.0) + 1.0));
  if (aziorder == 0) numpts = 2 * ((numlegendre + 1) / 2) + 4;
  ARTS_USER_ERROR_IF(aziorder > 0 and numpts > 512,
                     "SCATTERING samples {} azimuths for a series of degree {}; FFT1DR takes at most 512",
                     numpts,
                     numlegendre);
  ARTS_USER_ERROR_IF(numpts > 1024 or 2 * aziorder + 1 > 1024,
                     "SCATTERING samples {} azimuths and makes {} Fourier modes; FOURIER_MATRIX takes at most 1024",
                     numpts,
                     2 * aziorder + 1);

  // PHASE_MATRIX(4, 4), SCAT_MATRIX(4, 4, NUMPTS+1) and BASIS_MATRIX(4, 4,
  // 2*AZIORDER+1) in Fortran layout.  SCATTERING writes SCAT_MATRIX(1, 1,
  // NUMPTS+1), which FOURIER_MATRIX does not read.  Each mode goes straight
  // to scatbuf, which held OUT_MATRIX's copy.
  const Range stokes{0, nstokes};
  Matrix44    phase_matrix{};
  Tensor3&    scat_matrix  = work.scat_matrix.resize(numpts + 1, 4, 4);
  Tensor3&    basis_matrix = work.basis_matrix.resize(2 * aziorder + 1, 4, 4);

  // Find how many Legendre series must be summed
  const IndexVector6 dosum = number_sums(nstokes, legendre_coef);

  // MU1 is the incoming direction, and MU2 is the outgoing direction.
  for (Index j1 = 0; j1 < nummu; j1++) {
    const Numeric tmp = quad_weights[j1] / 2.0;
    for (Index j2 = 0; j2 < nummu; j2++) {
      for (Index l = 1; l <= 2; l++) {
        Numeric       mu1 = mu_values[j1];
        const Numeric mu2 = mu_values[j2];
        if (l % 2 == 0) mu1 = -mu1;
        // Only need to calculate phase matrix for half of
        // the delphi's, the rest come from symmetry.
        for (Index k = 1; k <= numpts / 2 + 1; k++) {
          const Numeric delphi   = (two_pi * static_cast<Numeric>(k - 1)) / static_cast<Numeric>(numpts);
          const Numeric cos_scat = mu1 * mu2 + std::sqrt((1.0 - mu1 * mu1) * (1.0 - mu2 * mu2)) * std::cos(delphi);
          sum_legendre(legendre_coef, cos_scat, dosum, phase_matrix, work);
          rotate_phase_matrix(phase_matrix, mu1, mu2, delphi, cos_scat, scat_matrix[k - 1, stokes, stokes]);
          // k = numpts / 2 + 1 maps onto itself
          matrix_symmetry(scat_matrix[k - 1, stokes, stokes], scat_matrix[numpts - k + 1, stokes, stokes]);
        }
        fourier_matrix(scat_matrix[Range{0, numpts}, stokes, stokes], basis_matrix[joker, stokes, stokes], work);

        for (Index m = 0; m <= aziorder; m++)
          combine_phase_modes(m, tmp, basis_matrix[joker, stokes, stokes], scatbuf[m, l - 1, j1, j2]);
      }
    }
  }
}

void direct_scattering(ConstVectorView mu_values,
                       ConstMatrixView legendre_coef,
                       Numeric         direct_mu,
                       Tensor4View     directbuf,
                       rt3_workdata&   work) {
  using Constant::two_pi;

  const Index nummu       = mu_values.size();
  const Index aziorder    = directbuf.extent(0) - 1;
  const Index nstokes     = directbuf.extent(3);
  const Index numlegendre = legendre_coef.nrows() - 1;
  ARTS_USER_ERROR_IF(nummu < 1 or numlegendre < 0 or legendre_coef.ncols() != 6 or nstokes < 1 or nstokes > 4 or
                         aziorder < 0 or directbuf.shape() != (std::array<Index, 4>{aziorder + 1, 2, nummu, nstokes}),
                     "DIRECT_SCATTERING needs mu_values [nummu], legendre_coef [numlegendre + 1, 6] and directbuf "
                     "[aziorder + 1, 2, nummu, nstokes] with 1 to 4 Stokes parameters; got {}, {:B,} and {:B,}",
                     mu_values.size(),
                     legendre_coef.shape(),
                     directbuf.shape());
  ARTS_USER_ERROR_IF(
      not(direct_mu >= -1.0 and direct_mu <= 1.0), "DIRECT_SCATTERING needs direct_mu in [-1, 1], got {}", direct_mu);

  Index numpts =
      2 * (Index{1} << static_cast<Index>(std::log(static_cast<Numeric>(numlegendre + 4)) / std::log(2.0) + 1.0));
  if (aziorder == 0) numpts = 2 * ((numlegendre + 1) / 2) + 4;
  ARTS_USER_ERROR_IF(aziorder > 0 and numpts > 512,
                     "DIRECT_SCATTERING samples {} azimuths for a series of degree {}; FFT1DR takes at most 512",
                     numpts,
                     numlegendre);
  ARTS_USER_ERROR_IF(numpts > 512 or 2 * aziorder + 1 > 512,
                     "DIRECT_SCATTERING samples {} azimuths and makes {} Fourier modes; RT3 takes at most 512",
                     numpts,
                     2 * aziorder + 1);

  // PHASE_MATRIX(4, 4), SCAT_MATRIX(4, 4, NUMPTS) and BASIS_MATRIX(4, 4,
  // 2*AZIORDER+1) in Fortran layout.  The vector of a mode holds the I and
  // Q (iq) of its cosine and the U and V (uv) of its sine.
  const Range stokes{0, nstokes}, iq{0, std::min<Index>(nstokes, 2)}, uv{2, nstokes - 2};
  Matrix44    phase_matrix{};
  Tensor3&    scat_matrix  = work.scat_matrix.resize(numpts, 4, 4);
  Tensor3&    basis_matrix = work.basis_matrix.resize(2 * aziorder + 1, 4, 4);

  // Find how many Legendre series must be summed
  const IndexVector6 dosum = number_sums(nstokes, legendre_coef);

  for (Index j = 0; j < nummu; j++) {
    for (Index l = 0; l < 2; l++) {
      const Numeric mu2 = l == 0 ? mu_values[j] : -mu_values[j];
      for (Index k = 0; k < numpts; k++) {
        const Numeric delphi = -(two_pi * static_cast<Numeric>(k)) / static_cast<Numeric>(numpts);
        const Numeric cos_scat =
            mu2 * direct_mu + std::sqrt((1.0 - mu2 * mu2) * (1.0 - direct_mu * direct_mu)) * std::cos(delphi);
        sum_legendre(legendre_coef, cos_scat, dosum, phase_matrix, work);
        rotate_phase_matrix(phase_matrix, direct_mu, mu2, delphi, cos_scat, scat_matrix[k, stokes, stokes]);
      }
      fourier_matrix(scat_matrix[joker, stokes, stokes], basis_matrix[joker, stokes, stokes], work);

      // Store away the first column of the combined mode phase matrix
      directbuf[0, l, j, iq] = basis_matrix[0, 0, iq];
      if (nstokes > 2) directbuf[0, l, j, uv] = 0.0;
      for (Index m = 1; m <= aziorder; m++) {
        directbuf[m, l, j, iq] = basis_matrix[m, 0, iq];
        if (nstokes > 2) directbuf[m, l, j, uv] = basis_matrix[m + aziorder, 0, uv];
      }
    }
  }
}

void get_scattering(Index mode, ConstTensor6View scatbuf, Tensor5View scatter_matrix) {
  const Index aziorder = scatbuf.extent(0) - 1;
  const Index nummu    = scatbuf.extent(2);
  const Index nstokes  = scatbuf.extent(4);
  ARTS_USER_ERROR_IF(aziorder < 0 or nummu < 1 or nstokes < 1 or nstokes > 4 or
                         scatbuf.shape() != (std::array<Index, 6>{aziorder + 1, 2, nummu, nummu, nstokes, nstokes}) or
                         scatter_matrix.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}) or
                         mode < 0 or mode > aziorder,
                     "GET_SCATTERING needs scatbuf [aziorder + 1, 2, nummu, nummu, nstokes, nstokes], scatter_matrix "
                     "[4, nummu, nstokes, nummu, nstokes] and 0 <= mode <= aziorder; got {:B,}, {:B,} and mode {}",
                     scatbuf.shape(),
                     scatter_matrix.shape(),
                     mode);

  // Read in the phase matrix for the azimuth mode: SCATBUF's record of
  // (L, J1, J2) is the column-major (I2, I1) matrix,
  // SCATTER_MATRIX(I2, J2, I1, J1, L) is [L, J1, I1, J2, I2]
  for (Index l = 0; l < 2; l++)
    for (Index j1 = 0; j1 < nummu; j1++)
      for (Index j2 = 0; j2 < nummu; j2++) scatter_matrix[l, j1, joker, j2, joker] = scatbuf[mode, l, j1, j2];

  // Use the symmetry of the scattering matrix to get
  // P-- from P++, and P-+ from P+-.
  scatter_symmetry(scatter_matrix);
}

void scatter_symmetry(Tensor5View scat) {
  const Index nummu   = scat.extent(1);
  const Index nstokes = scat.extent(2);
  ARTS_USER_ERROR_IF(nummu < 1 or nstokes < 1 or nstokes > 4 or
                         scat.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}),
                     "SCATTER_SYMMETRY needs scat [4, nummu, nstokes, nummu, nstokes] with 1 to 4 Stokes "
                     "parameters; got {:B,}",
                     scat.shape());

  // P-- is P++ and P-+ is P+-, with the off-diagonal 2 by 2 blocks negated
  scat[3] = scat[0];
  scat[2] = scat[1];
  if (nstokes > 2) {
    const Range iq{0, 2}, uv{2, nstokes - 2};
    for (Index l : {2, 3}) {
      scat[l, joker, iq, joker, uv] *= -1.0;
      scat[l, joker, uv, joker, iq] *= -1.0;
    }
  }
}

void check_norm(ConstVectorView quad_weights, ConstTensor5View scatter_matrix) {
  const Index nummu   = quad_weights.size();
  const Index nstokes = scatter_matrix.extent(2);
  ARTS_USER_ERROR_IF(
      nummu < 1 or nstokes < 1 or scatter_matrix.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}),
      "CHECK_NORM needs quad_weights [nummu] and scatter_matrix [4, nummu, nstokes, nummu, nstokes]; "
      "got {} and {:B,}",
      quad_weights.size(),
      scatter_matrix.shape());

  // The largest deviation from 1, NaN if any is
  Numeric maxsum = 0.0;
  for (Index j1 = 0; j1 < nummu; j1++) {
    const Numeric tmp = quad_weights[j1];
    if (tmp != 0.0) {
      for (Index l = 0; l < 2; l++) {
        Numeric sum = -1.0;
        for (Index j2 = 0; j2 < nummu; j2++)
          sum += quad_weights[j2] / tmp * (scatter_matrix[l, j1, 0, j2, 0] + scatter_matrix[l + 2, j1, 0, j2, 0]);
        if (not(std::abs(sum) <= maxsum)) maxsum = std::abs(sum);
      }
    }
  }
  ARTS_USER_ERROR_IF(not(maxsum <= 1.0e-7),
                     "Phase function not normalized: its I-I term integrates to 1 within {}, more than 1e-7.  Either "
                     "the first term of its Legendre series is not 1, or the series has too many terms for the "
                     "number of quadrature angles.",
                     maxsum);
}

void get_direct(Index mode, ConstTensor4View directbuf, Tensor3View direct_vector) {
  const Index aziorder = directbuf.extent(0) - 1;
  ARTS_USER_ERROR_IF(aziorder < 0 or directbuf.extent(1) != 2 or direct_vector.shape() != directbuf[0].shape() or
                         mode < 0 or mode > aziorder,
                     "GET_DIRECT needs directbuf [aziorder + 1, 2, nummu, nstokes], direct_vector "
                     "[2, nummu, nstokes] and 0 <= mode <= aziorder; got {:B,}, {:B,} and mode {}",
                     directbuf.shape(),
                     direct_vector.shape(),
                     mode);

  direct_vector = directbuf[mode];
}
}  // namespace rt3
