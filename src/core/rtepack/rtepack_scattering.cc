#include "rtepack_scattering.h"

#include <arts_conversions.h>
#include <debug.h>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
Array<muelmat_vector> bulk_backscatter(const ConstTensor5View &Pe,
                                       const ConstMatrixView &pnd) {
  assert(Pe.ncols() == 4 and Pe.nrows() == 4);
  const Index nv = Pe.npages();
  const Index np = Pe.nbooks();
  const Index ne = Pe.nshelves();

  Array<muelmat_vector> aotm(np, muelmat_vector(nv));
  if (ne == 0) return aotm;

  for (Index ip = 0; ip < np; ip++) {
    for (Index iv = 0; iv < nv; iv++) {
      auto m       = Pe[0, ip, iv, joker, joker];
      aotm[ip][iv] = pnd[0, ip] * muelmat(m[0, 0],
                                          m[0, 1],
                                          m[0, 2],
                                          m[0, 3],
                                          m[1, 0],
                                          m[1, 1],
                                          m[1, 2],
                                          m[1, 3],
                                          m[2, 0],
                                          m[2, 1],
                                          m[2, 2],
                                          m[2, 3],
                                          m[3, 0],
                                          m[3, 1],
                                          m[3, 2],
                                          m[3, 3]);
    }
  }

  for (Index ie = 1; ie < ne; ie++) {
    for (Index ip = 0; ip < np; ip++) {
      for (Index iv = 0; iv < nv; iv++) {
        auto m        = Pe[ie, ip, iv, joker, joker];
        aotm[ip][iv] += pnd[ie, ip] * muelmat(m[0, 0],
                                              m[0, 1],
                                              m[0, 2],
                                              m[0, 3],
                                              m[1, 0],
                                              m[1, 1],
                                              m[1, 2],
                                              m[1, 3],
                                              m[2, 0],
                                              m[2, 1],
                                              m[2, 2],
                                              m[2, 3],
                                              m[3, 0],
                                              m[3, 1],
                                              m[3, 2],
                                              m[3, 3]);
      }
    }
  }
  return aotm;
}

Array<muelmat_matrix> bulk_backscatter_derivative(
    const ConstTensor5View &Pe, const ArrayOfMatrix &dpnd_dx) {
  assert(Pe.ncols() == 4 and Pe.nrows() == 4);

  const Index nv = Pe.npages();
  const Index np = Pe.nbooks();
  const Index ne = Pe.nshelves();
  const Size nq  = dpnd_dx.size();

  Array<muelmat_matrix> aoaotm(np, muelmat_matrix(nq, nv));

  for (Index ip = 0; ip < np; ip++) {
    for (Size iq = 0; iq < nq; iq++) {
      aoaotm[ip][iq] = 0.0;
      for (Index iv = 0; iv < nv; iv++) {
        for (Index ie = 0; ie < ne; ie++) {
          auto m              = Pe[ie, ip, iv, joker, joker];
          aoaotm[ip][iq, iv] += dpnd_dx[iq][ie, ip] * muelmat(m[0, 0],
                                                              m[0, 1],
                                                              m[0, 2],
                                                              m[0, 3],
                                                              m[1, 0],
                                                              m[1, 1],
                                                              m[1, 2],
                                                              m[1, 3],
                                                              m[2, 0],
                                                              m[2, 1],
                                                              m[2, 2],
                                                              m[2, 3],
                                                              m[3, 0],
                                                              m[3, 1],
                                                              m[3, 2],
                                                              m[3, 3]);
        }
      }
    }
  }
  return aoaotm;
}

namespace {
void setBackscatterTransmission(stokvec_vector &out,
                                const stokvec_vector &I0,
                                const muelmat_vector &Tr,
                                const muelmat_vector &Tf,
                                const muelmat_vector &Z) {
  for (Size i = 0; i < out.size(); i++) out[i] = Tr[i] * Z[i] * Tf[i] * I0[i];
}

void setBackscatterTransmissionDerivative(stokvec_matrix &out,
                                          const stokvec_vector &I0,
                                          const muelmat_vector &Tr,
                                          const muelmat_vector &Tf,
                                          const muelmat_matrix &dZ) {
  for (Index j = 0; j < dZ.nrows(); j++)
    for (Index i = 0; i < dZ.ncols(); i++)
      out[j, i] += Tr[i] * dZ[j, i] * Tf[i] * I0[i];
}
}  // namespace

void bulk_backscatter_commutative_transmission_rte(
    Array<stokvec_vector> &I,
    Array<Array<stokvec_matrix>> &dI,
    const stokvec_vector &I_incoming,
    const Array<muelmat_vector> &T,
    const Array<muelmat_vector> &PiTf,
    const Array<muelmat_vector> &PiTr,
    const Array<muelmat_vector> &Z,
    const Array<muelmat_matrix> &dT1,
    const Array<muelmat_matrix> &dT2,
    const Array<muelmat_matrix> &dZ) {
  const Size np  = dT1.size();
  const Index nv = np ? dT1.front().ncols() : 0;
  const Index nq = np ? dT1.front().nrows() : 0;

  // For all transmission, the I-vector is the same
  for (Size ip = 0; ip < np; ip++)
    setBackscatterTransmission(I[ip], I_incoming, PiTr[ip], PiTf[ip], Z[ip]);

  for (Size ip = 0; ip < np; ip++) {
    setBackscatterTransmissionDerivative(
        dI[ip][ip], I_incoming, PiTr[ip], PiTf[ip], dZ[ip]);
  }

  for (Size ip = 0; ip < np; ip++) {
    for (Size j = ip; j < np; j++) {
      for (Index iq = 0; iq < nq; iq++) {
        for (Index iv = 0; iv < nv; iv++) {
          dI[ip][j][iq, iv] +=
              T[ip][iv] * ((dT1[ip][iq, iv] + dT2[ip][iq, iv]) * I[j][iv]);

          if (j < np - 1 and j > ip)
            dI[ip][j][iq, iv] +=
                inv(T[ip + 1][iv]) *
                ((dT1[ip][iq, iv] + dT2[ip][iq, iv]) * I[j][iv]);
        }
      }
    }
  }
}

namespace {
Numeric cos_scat_angle(const Vector2 &los_in, const Vector2 &los_out) {
  using Conversion::cosd;
  using Conversion::sind;

  const auto &[za_sca, aa_sca] = los_out;
  const auto &[za_inc, aa_inc] = los_in;

  const Numeric theta = cosd(za_sca) * cosd(za_inc) +
                        sind(za_sca) * sind(za_inc) * cosd(aa_sca - aa_inc);

  // Fix potential overflows by clamping numerical errors
  return std::clamp(theta, -1.0, 1.0);
}
}  // namespace

muelmat rayleigh_scattering(const Vector2 &los_in,
                            const Vector2 &los_out,
                            const Numeric depolarization_factor) try {
  using Constant::pi;
  using Math::pow2;

  const auto cos_theta        = cos_scat_angle(los_in, los_out);
  const auto theta_rad        = std::acos(cos_theta);
  const auto sin_theta        = std::sin(theta_rad);
  const auto [za_sca, aa_sca] = los_out;
  const auto [za_inc, aa_inc] = los_in;

  const Numeric delta =
      (1 - depolarization_factor) / (1 + 0.5 * depolarization_factor);
  const Numeric delta_prime =
      (1 - 2 * depolarization_factor) / (1 - depolarization_factor);

  const auto F11 = 0.75 * delta * (1 + pow2(cos_theta)) + 1 - delta;
  const auto F12 = -0.75 * delta * pow2(sin_theta);
  const auto F22 = 0.75 * delta * (1 + pow2(cos_theta));
  const auto F33 = 1.5 * delta * cos_theta;
  const auto F44 = 1.5 * delta * delta_prime * cos_theta;

  muelmat pha;
  pha[0, 0] = F11;

  const Numeric za_sca_rad = Conversion::deg2rad(za_sca);
  const Numeric za_inc_rad = Conversion::deg2rad(za_inc);
  const Numeric aa_sca_rad = Conversion::deg2rad(aa_sca);
  const Numeric aa_inc_rad = Conversion::deg2rad(aa_inc);

  constexpr Numeric ANGTOL_RAD = 1e-6;  //CPD: this constant is used to adjust
  //zenith angles close to 0 and PI.  This is
  //also used to avoid float == float statements.

  //
  // Several cases have to be considered:
  //

  if ((std::abs(theta_rad) < ANGTOL_RAD)          // forward scattering
      || (std::abs(theta_rad - pi) < ANGTOL_RAD)  // backward scattering
      || (std::abs(aa_inc_rad - aa_sca_rad) <
          ANGTOL_RAD)  // inc and sca on meridian
      || (std::abs(std::abs(aa_inc_rad - aa_sca_rad) - Constant::two_pi) <
          ANGTOL_RAD)                                                     //   "
      || (std::abs(std::abs(aa_inc_rad - aa_sca_rad) - pi) < ANGTOL_RAD)  //   "
  ) {
    pha[0, 1] = F12;
    pha[1, 0] = F12;
    pha[1, 1] = F22;
    pha[0, 2] = 0;
    pha[1, 2] = 0;
    pha[2, 0] = 0;
    pha[2, 1] = 0;
    pha[2, 2] = F33;
    pha[0, 3] = 0;
    pha[1, 3] = 0;
    pha[2, 3] = 0;
    pha[3, 0] = 0;
    pha[3, 1] = 0;
    pha[3, 2] = 0;
    pha[3, 3] = F44;
  } else {
    Numeric sigma1;
    Numeric sigma2;

    Numeric s1, s2;

    // In these cases we have to take limiting values.

    if (za_inc_rad < ANGTOL_RAD) {
      sigma1 = pi + aa_sca_rad - aa_inc_rad;
      sigma2 = 0;
    } else if (za_inc_rad > pi - ANGTOL_RAD) {
      sigma1 = aa_sca_rad - aa_inc_rad;
      sigma2 = pi;
    } else if (za_sca_rad < ANGTOL_RAD) {
      sigma1 = 0;
      sigma2 = pi + aa_sca_rad - aa_inc_rad;
    } else if (za_sca_rad > pi - ANGTOL_RAD) {
      sigma1 = pi;
      sigma2 = aa_sca_rad - aa_inc_rad;
    } else {
      s1 = (std::cos(za_sca_rad) - std::cos(za_inc_rad) * cos_theta) /
           (std::sin(za_inc_rad) * sin_theta);
      s2 = (std::cos(za_inc_rad) - std::cos(za_sca_rad) * cos_theta) /
           (std::sin(za_sca_rad) * sin_theta);

      sigma1 = std::acos(std::clamp(s1, -1., 1.));
      sigma2 = std::acos(std::clamp(s2, -1., 1.));
    }

    const Numeric C1 = std::cos(2 * sigma1);
    const Numeric C2 = std::cos(2 * sigma2);

    const Numeric S1 = std::sin(2 * sigma1);
    const Numeric S2 = std::sin(2 * sigma2);

    pha[0, 1] = C1 * F12;
    pha[1, 0] = C2 * F12;
    pha[1, 1] = C1 * C2 * F22 - S1 * S2 * F33;

    ARTS_USER_ERROR_IF(
        std::isnan(pha[0, 1]) || std::isnan(pha[1, 0]) || std::isnan(pha[1, 1]),
        "NaN value(s) detected in *phaCalc* (0/1,1). Could the "
        "input data contain NaNs? Please check with *scat_dataCheck*. If "
        "input data are OK  and you critically need the ongoing calculations, "
        "try to change the observation LOS slightly. If you can reproduce "
        "this error, please contact Patrick in order to help tracking down "
        "the reason to this problem. If you see this message occasionally "
        "when doing MC calculations, it should not be critical. This path "
        "sampling will be rejected and replaced with a new one.");

    /*CPD: For skokes_dim > 2 some of the transformation formula
            for each element have a different sign depending on whether or
            not 0<aa_scat-aa_inc<180.  For details see pages 94 and 95 of
            Mishchenkos chapter in :
            Mishchenko, M. I., and L. D. Travis, 2003: Electromagnetic
            scattering by nonspherical particles. In Exploring the Atmosphere
            by Remote Sensing Techniques (R. Guzzi, Ed.), Springer-Verlag,
            Berlin, pp. 77-127.
            This is available at http://www.giss.nasa.gov/~crmim/publications/ */
    const Numeric delta_aa = aa_sca - aa_inc + (aa_sca - aa_inc < -180) * 360 -
                             (aa_sca - aa_inc > 180) * 360;
    if (delta_aa >= 0) {
      pha[0, 2] = S1 * F12;
      pha[1, 2] = S1 * C2 * F22 + C1 * S2 * F33;
      pha[2, 0] = -S2 * F12;
      pha[2, 1] = -C1 * S2 * F22 - S1 * C2 * F33;
    } else {
      pha[0, 2] = -S1 * F12;
      pha[1, 2] = -S1 * C2 * F22 - C1 * S2 * F33;
      pha[2, 0] = S2 * F12;
      pha[2, 1] = C1 * S2 * F22 + S1 * C2 * F33;
    }
    pha[2, 2] = -S1 * S2 * F22 + C1 * C2 * F33;

    pha[1, 3] = 0;
    pha[3, 1] = 0;
    pha[0, 3] = 0;
    pha[2, 3] = 0;
    pha[3, 0] = 0;
    pha[3, 2] = 0;
    pha[3, 3] = F44;
  }

  return pha;
}
ARTS_METHOD_ERROR_CATCH
}  // namespace rtepack