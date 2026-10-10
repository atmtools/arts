#include "radscat3.h"

#include <arts_constants.h>
#include <debug.h>
#include <legendre.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace polradtran::rt3 {
CompactPlanarMuelmat delta_m_scaled(const CompactPlanarMuelmat& c, Index l, Numeric f) {
  const auto k = static_cast<Numeric>(2 * l + 1);
  return k * (c / k - f * CompactPlanarMuelmat::id()) / (1.0 - f);
}

void get_scat_set(bool                                delta_m,
                  Index                               nummu,
                  CompactPlanarMuelmatConstVectorView coefin,
                  Numeric                             extin,
                  Numeric                             scatin,
                  Index&                              nlegen,
                  CompactPlanarMuelmatVectorView      coef,
                  Numeric&                            extinction,
                  Numeric&                            scatter) {
  const Index nlegin = static_cast<Index>(coefin.size()) - 1;
  const Index ncoef  = std::max(nlegin + 1, 2 * nummu);
  ARTS_USER_ERROR_IF(nummu < 1 or nlegin < 0 or static_cast<Index>(coef.size()) < ncoef,
                     "GET_SCAT_SET needs nummu >= 1, coefin [nlegin + 1] and coef [max(nlegin + 1, 2 nummu) or "
                     "more]; got nummu {}, coefin [{}] and coef [{}]",
                     nummu,
                     coefin.size(),
                     coef.size());

  // Copy the set where READ_SCAT_FILE read the file.  With delta_m the
  // scaling below reads coef up to l = 2 nummu even when nlegen + 1 is
  // smaller; those coefficients are zero.
  extinction                 = extin;
  scatter                    = scatin;
  nlegen                     = nlegin;
  coef[Range{0, ncoef}]      = CompactPlanarMuelmat{};
  coef[Range{0, nlegen + 1}] = coefin;

  if (delta_m) {
    const Index m = 2 * nummu;
    Numeric     f = 0.0;
    if (m + 1 <= nlegen + 1) f = coef[m].F11() / static_cast<Numeric>(2 * m + 1);
    ARTS_USER_ERROR_IF(
        not(extinction > 0.0), "Delta-M scaling divides by the extinction, which must be positive, got {}", extinction);
    ARTS_USER_ERROR_IF(
        1.0 - f == 0.0, "Delta-M scaling divides by 1 - f, with f = coef[{}].F11() / {} = 1", m, 2 * m + 1);
    Numeric albedo = scatter / extinction;
    ARTS_USER_ERROR_IF(
        1.0 - albedo * f == 0.0, "Delta-M scaling divides by 1 - albedo f, with albedo {} and f {}", albedo, f);
    extinction = (1.0 - albedo * f) * extinction;
    albedo     = (1.0 - f) * albedo / (1.0 - albedo * f);
    scatter    = albedo * extinction;
    nlegen     = m - 1;
    for (Index l = 0; l <= nlegen; l++) coef[l] = delta_m_scaled(coef[l], l, f);
  }
}

CompactPlanarMuelmat sum_legendre(CompactPlanarMuelmatConstVectorView coef,
                                  Numeric                             x,
                                  Index                               nstokes,
                                  rt3_workdata&                       work) {
  const Index nlegen = static_cast<Index>(coef.size()) - 1;
  ARTS_USER_ERROR_IF(nlegen < 0 or nstokes < 1 or nstokes > 4,
                     "SUM_LEGENDRE needs coef [nlegen + 1] and 1 to 4 Stokes parameters; got [{}] and {}",
                     coef.size(),
                     nstokes);

  // The Legendre polynomials P_0(x) to P_nlegen(x), for all the series
  Vector& p = work.legendre_p.resize(nlegen + 1);
  Legendre::legendre_polynomials(p, std::clamp(x, -1.0, 1.0));

  // Sum the Legendre series, or that of F11 alone for the intensity
  CompactPlanarMuelmat f{};
  if (nstokes == 1) {
    for (Index l = 0; l <= nlegen; l++) f.F11() += p[l] * coef[l].F11();
  } else {
    for (Index l = 0; l <= nlegen; l++) f += p[l] * coef[l];
  }
  return f;
}

Muelmat rotate_phase_matrix(
    const CompactPlanarMuelmat& phase_matrix, Numeric mu1, Numeric mu2, Numeric delphi, Numeric cos_scat) {
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

  // From the incident meridional plane into the scattering plane, and from
  // it to the outgoing meridional plane
  return rotated(phase_matrix, cos21, -sin21, cos22, -sin22);
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

void fourier_matrix(MuelmatConstVectorView real_matrix,
                    MuelmatVectorView      basis_matrix,
                    Index                  nstokes,
                    rt3_workdata&          work) {
  const Index numpts   = real_matrix.size();
  const Index numazi   = basis_matrix.size();
  const Index aziorder = (numazi - 1) / 2;
  ARTS_USER_ERROR_IF(numpts < 1 or nstokes < 1 or nstokes > 4 or numazi % 2 != 1,
                     "FOURIER_MATRIX needs real_matrix [numpts], basis_matrix [2 aziorder + 1] and 1 to 4 Stokes "
                     "parameters; got {}, {} and {}",
                     numpts,
                     numazi,
                     nstokes);

  // REAL_VECTOR and BASIS_VECTOR; FFT1DR transforms REAL_VECTOR in place
  Vector& real_vector  = work.real_vector.resize(numpts);
  Vector& basis_vector = work.basis_vector.resize(numazi);
  basis_matrix         = Muelmat{0.0};
  for (Index i = 0; i < nstokes; i++) {
    for (Index j = 0; j < nstokes; j++) {
      for (Index k = 0; k < numpts; k++) real_vector[k] = real_matrix[k][i, j];
      fourier_basis(aziorder, fourier_direction::to_basis, basis_vector, real_vector, work.fft);
      for (Index k = 0; k < numazi; k++) basis_matrix[k][i, j] = basis_vector[k];
    }
  }
}

Muelmat combine_phase_modes(Index m, Numeric tmp, MuelmatConstVectorView basis_matrix) {
  const Index aziorder = (basis_matrix.size() - 1) / 2;
  ARTS_USER_ERROR_IF(basis_matrix.size() % 2 != 1 or m < 0 or m > aziorder,
                     "COMBINE_PHASE_MODES needs basis_matrix [2 aziorder + 1] and 0 <= m <= aziorder; got {} and m {}",
                     basis_matrix.size(),
                     m);

  // SINFLAG: 0 in the diagonal 2 by 2 blocks (cosine modes), -1 in the
  // lower left [uv, iq] and +1 in the upper right [iq, uv] (sine modes)
  const Range iq{0, 2}, uv{2, 2};
  if (m == 0) {
    Muelmat out  = basis_matrix[0];
    out         *= tmp;
    out[iq, uv]  = 0.0;
    out[uv, iq]  = 0.0;
    return out;
  }

  // MC = M+1 and MS = M+1+AZIORDER, 0-based
  const Index    mc = m, ms = m + aziorder;
  const Muelmat& sine  = basis_matrix[ms];
  Muelmat        out   = basis_matrix[mc];
  out[uv, iq]          = sine[uv, iq];
  out[uv, iq]         *= -1.0;
  out[iq, uv]          = sine[iq, uv];
  out                 *= 0.5 * tmp;
  return out;
}

void scattering(ConstVectorView                     mu_values,
                ConstVectorView                     quad_weights,
                CompactPlanarMuelmatConstVectorView legendre_coef,
                Index                               nstokes,
                MuelmatTensor4View                  scatbuf,
                rt3_workdata&                       work) {
  using Constant::two_pi;

  const Index nummu       = mu_values.size();
  const Index aziorder    = scatbuf.extent(0) - 1;
  const Index numlegendre = static_cast<Index>(legendre_coef.size()) - 1;
  ARTS_USER_ERROR_IF(nummu < 1 or quad_weights.size() != mu_values.size() or numlegendre < 0 or nstokes < 1 or
                         nstokes > 4 or aziorder < 0 or
                         scatbuf.shape() != (std::array<Index, 4>{aziorder + 1, 2, nummu, nummu}),
                     "SCATTERING needs mu_values and quad_weights [nummu], legendre_coef [numlegendre + 1], 1 to 4 "
                     "Stokes parameters and scatbuf [aziorder + 1, 2, nummu, nummu]; got {}, {}, [{}], {} and {:B,}",
                     mu_values.size(),
                     quad_weights.size(),
                     legendre_coef.size(),
                     nstokes,
                     scatbuf.shape());

  Index numpts =
      2 * (Index{1} << static_cast<Index>(std::log(static_cast<Numeric>(numlegendre + 4)) / std::log(2.0) + 1.0));
  if (aziorder == 0) numpts = 2 * ((numlegendre + 1) / 2) + 4;

  // SCAT_MATRIX(4, 4, NUMPTS+1) and BASIS_MATRIX(4, 4, 2*AZIORDER+1).
  // SCATTERING writes SCAT_MATRIX(1, 1, NUMPTS+1), which FOURIER_MATRIX does
  // not read.  Each mode goes straight to scatbuf, which held OUT_MATRIX's
  // copy.
  MuelmatVector& scat_matrix  = work.scat_matrix.resize(numpts + 1);
  MuelmatVector& basis_matrix = work.basis_matrix.resize(2 * aziorder + 1);

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
          scat_matrix[k - 1] =
              rotate_phase_matrix(sum_legendre(legendre_coef, cos_scat, nstokes, work), mu1, mu2, delphi, cos_scat);
          // MATRIX_SYMMETRY; k = numpts / 2 + 1 maps onto itself
          scat_matrix[numpts - k + 1] = mirror(scat_matrix[k - 1]);
        }
        fourier_matrix(scat_matrix[Range{0, numpts}], basis_matrix, nstokes, work);

        for (Index m = 0; m <= aziorder; m++) scatbuf[m, l - 1, j1, j2] = combine_phase_modes(m, tmp, basis_matrix);
      }
    }
  }
}

void direct_scattering(ConstVectorView                     mu_values,
                       CompactPlanarMuelmatConstVectorView legendre_coef,
                       Numeric                             direct_mu,
                       Index                               nstokes,
                       StokvecTensor3View                  directbuf,
                       rt3_workdata&                       work) {
  using Constant::two_pi;

  const Index nummu       = mu_values.size();
  const Index aziorder    = directbuf.extent(0) - 1;
  const Index numlegendre = static_cast<Index>(legendre_coef.size()) - 1;
  ARTS_USER_ERROR_IF(nummu < 1 or numlegendre < 0 or nstokes < 1 or nstokes > 4 or aziorder < 0 or
                         directbuf.shape() != (std::array<Index, 3>{aziorder + 1, 2, nummu}),
                     "DIRECT_SCATTERING needs mu_values [nummu], legendre_coef [numlegendre + 1], 1 to 4 Stokes "
                     "parameters and directbuf [aziorder + 1, 2, nummu]; got {}, [{}], {} and {:B,}",
                     mu_values.size(),
                     legendre_coef.size(),
                     nstokes,
                     directbuf.shape());
  ARTS_USER_ERROR_IF(
      not(direct_mu >= -1.0 and direct_mu <= 1.0), "DIRECT_SCATTERING needs direct_mu in [-1, 1], got {}", direct_mu);

  Index numpts =
      2 * (Index{1} << static_cast<Index>(std::log(static_cast<Numeric>(numlegendre + 4)) / std::log(2.0) + 1.0));
  if (aziorder == 0) numpts = 2 * ((numlegendre + 1) / 2) + 4;

  // SCAT_MATRIX(4, 4, NUMPTS) and BASIS_MATRIX(4, 4, 2*AZIORDER+1)
  MuelmatVector& scat_matrix  = work.scat_matrix.resize(numpts);
  MuelmatVector& basis_matrix = work.basis_matrix.resize(2 * aziorder + 1);

  for (Index j = 0; j < nummu; j++) {
    for (Index l = 0; l < 2; l++) {
      const Numeric mu2 = l == 0 ? mu_values[j] : -mu_values[j];
      for (Index k = 0; k < numpts; k++) {
        const Numeric delphi = -(two_pi * static_cast<Numeric>(k)) / static_cast<Numeric>(numpts);
        const Numeric cos_scat =
            mu2 * direct_mu + std::sqrt((1.0 - mu2 * mu2) * (1.0 - direct_mu * direct_mu)) * std::cos(delphi);
        scat_matrix[k] =
            rotate_phase_matrix(sum_legendre(legendre_coef, cos_scat, nstokes, work), direct_mu, mu2, delphi, cos_scat);
      }
      fourier_matrix(scat_matrix, basis_matrix, nstokes, work);

      // Store away the first column of the combined mode phase matrix: the
      // I and Q of its cosine mode and the U and V of its sine mode
      const Muelmat& mean = basis_matrix[0];
      directbuf[0, l, j]  = Stokvec{mean[0, 0], mean[1, 0], 0.0, 0.0};
      for (Index m = 1; m <= aziorder; m++) {
        const Muelmat &cosine = basis_matrix[m], &sine = basis_matrix[m + aziorder];
        directbuf[m, l, j] = Stokvec{cosine[0, 0], cosine[1, 0], sine[2, 0], sine[3, 0]};
      }
    }
  }
}

void get_scattering(Index mode, MuelmatConstTensor4View scatbuf, Tensor5View scatter_matrix) {
  const Index aziorder = scatbuf.extent(0) - 1;
  const Index nummu    = scatbuf.extent(2);
  const Index nstokes  = scatter_matrix.extent(2);
  ARTS_USER_ERROR_IF(aziorder < 0 or nummu < 1 or nstokes < 1 or nstokes > 4 or
                         scatbuf.shape() != (std::array<Index, 4>{aziorder + 1, 2, nummu, nummu}) or
                         scatter_matrix.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}) or
                         mode < 0 or mode > aziorder,
                     "GET_SCATTERING needs scatbuf [aziorder + 1, 2, nummu, nummu], scatter_matrix "
                     "[4, nummu, nstokes, nummu, nstokes] with 1 to 4 Stokes parameters and 0 <= mode <= aziorder; "
                     "got {:B,}, {:B,} and mode {}",
                     scatbuf.shape(),
                     scatter_matrix.shape(),
                     mode);

  // Read in the phase matrix for the azimuth mode: the leading nstokes x
  // nstokes of SCATBUF's matrix of (L, J1, J2), from Stokes parameter I1 to
  // I2, is SCATTER_MATRIX(I2, J2, I1, J1, L), which is [L, J1, I1, J2, I2]
  const Range stokes{0, nstokes};
  for (Index l = 0; l < 2; l++)
    for (Index j1 = 0; j1 < nummu; j1++)
      for (Index j2 = 0; j2 < nummu; j2++)
        scatter_matrix[l, j1, joker, j2, joker] = transpose(scatbuf[mode, l, j1, j2][stokes, stokes]);

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

void get_direct(Index mode, StokvecConstTensor3View directbuf, Tensor3View direct_vector) {
  const Index aziorder = directbuf.extent(0) - 1;
  const Index nummu    = directbuf.extent(2);
  const Index nstokes  = direct_vector.extent(2);
  ARTS_USER_ERROR_IF(aziorder < 0 or directbuf.extent(1) != 2 or nstokes < 1 or nstokes > 4 or
                         direct_vector.shape() != (std::array<Index, 3>{2, nummu, nstokes}) or mode < 0 or
                         mode > aziorder,
                     "GET_DIRECT needs directbuf [aziorder + 1, 2, nummu], direct_vector [2, nummu, nstokes] with 1 "
                     "to 4 Stokes parameters and 0 <= mode <= aziorder; got {:B,}, {:B,} and mode {}",
                     directbuf.shape(),
                     direct_vector.shape(),
                     mode);

  // The leading nstokes Stokes parameters of each vector
  const Range stokes{0, nstokes};
  for (Index l = 0; l < 2; l++)
    for (Index j = 0; j < nummu; j++) direct_vector[l, j] = directbuf[mode, l, j][stokes];
}
}  // namespace polradtran::rt3
