#pragma once

/* The laboratory-frame phase matrix of randomly oriented particles with a
   plane of symmetry, from rtepack's polarization geometry
   (rtepack_polarization.h, tested against electric fields).  Shared by the
   VDISORT comparisons with RT4 (vdisort-rt4-comparison.cpp), RT3
   (vdisort-rt3-comparison.cpp) and ARTS's scattering data
   (vdisort-arts-test.cpp) and the delta-M test (vdisort-delta-m-test.cpp).
   Neither solver's rotation formulas are used.

   Frame: z points up.  A ray propagating at zenith cosine mu (> 0 upward)
   and azimuth phi has the direction n = (sin cos phi, sin sin phi, mu) and
   the meridional basis h = n x z / |n x z|, v = h x n, at mu = +-1 that of
   the limit along the meridian of phi, h = (sin phi, -cos phi, 0).
   Stokes: I = I_v + I_h, Q = I_v - I_h, U = 2 Re(E_v E_h*). */

#include <arts_constants.h>
#include <matpack.h>
#include <rtepack.h>

#include <algorithm>
#include <cmath>
#include <functional>

namespace vdisort_test {
//! The scattering matrix in the scattering-plane basis (e_par, e_perp), Q = I_par - I_perp, as a function of
//! cos(Theta)
using tro_matrix = std::function<CompactPlanarMuelmat(Numeric)>;

inline Vector3 direction(Numeric cos_theta, Numeric phi) {
  const Numeric st = std::sqrt(std::max(0.0, 1 - cos_theta * cos_theta));
  return {st * std::cos(phi), st * std::sin(phi), cos_theta};
}

//! The meridional basis of the ray direction(cos_theta, phi), at the poles that of the azimuth phi
inline rtepack::polarization_basis basis(Numeric cos_theta, Numeric phi) {
  return rtepack::meridional_basis(direction(cos_theta, phi), Vector3{std::sin(phi), -std::cos(phi), 0.0});
}

/* The lab-frame phase matrix Z(out <- in) for in = (ci, phii) and
   out = (co, phio), in the meridional bases.  Exactly forward or backward
   the plane through n_in and h(n_in) is used: for |mu| < 1 a forward pair
   has the same and a backward pair the opposite azimuth, both rotations
   are then by 0 or pi, and Z = F, as in RT3 (ROTATE_PHASE_MATRIX for
   SIN_SCAT = 0).  At mu = +-1 the meridional basis depends on the azimuth
   label, and Z is F rotated between the two bases.  For a regular F
   (F12 = F34 = 0 and |F22| = |F33| there) Z does not depend on the choice
   of plane. */
inline rtepack::muelmat lab_frame(const tro_matrix& F, Numeric ci, Numeric phii, Numeric co, Numeric phio) {
  return rtepack::lab_frame_phase_matrix(
      F(std::clamp(dot(direction(ci, phii), direction(co, phio)), -1.0, 1.0)), basis(ci, phii), basis(co, phio));
}

//! F(cos Theta) of a Legendre series of scattering-plane matrices, F = sum_l c_l P_l(cos Theta)
inline tro_matrix legendre_matrix(const CompactPlanarMuelmatVector& coef) {
  return [coef](Numeric x) {
    // The sum along the recursion l P_l = (2 l - 1) x P_l-1 - (l - 1) P_l-2, without storing the P_l
    const auto           n = static_cast<Index>(coef.size());
    CompactPlanarMuelmat sum{};
    Numeric              p_prev = 0.0, p = 1.0;
    for (Index l = 0; l < n; l++) {
      if (l > 0) {
        const auto    k    = static_cast<Numeric>(l);
        const Numeric next = ((2 * k - 1) * x * p - (k - 1) * p_prev) / k;
        p_prev             = p;
        p                  = next;
      }
      sum += coef[l] * p;
    }
    return sum;
  };
}

struct ordinary_modes {
  rtepack::muelmat_tensor3 C;  // [m, out, in]
  rtepack::muelmat_tensor3 S;
};

/* The ordinary azimuthal Fourier coefficients, without eps_m,
     C^m(mu_o, mu_i) = (1 / 2 pi) int Z(mu_o, 0; mu_i, phi') cos(m phi') dphi',
   and S^m with sin(m phi'), on signed cosines (> 0 upward), by the
   trapezoid rule at the nphi equidistant phi' = 2 pi k / nphi,
   k = 0 .. nphi - 1, for the modes m < min(nfourier, nphi / 2).  Rows and
   columns from nstokes on are zeroed, as for a solver that transports the
   leading nstokes x nstokes block.  midpoints shifts the samples by half a
   step, and azimuth_sign = -1 takes Z at -phi' (a deliberate mistake). */
inline ordinary_modes fourier_modes(const tro_matrix& F,
                                    const Vector&     mu_out,
                                    const Vector&     mu_in,
                                    Index             nfourier,
                                    Index             nphi,
                                    Index             nstokes      = 4,
                                    bool              midpoints    = false,
                                    Numeric           azimuth_sign = 1.0) {
  const auto     no = static_cast<Index>(mu_out.size()), ni = static_cast<Index>(mu_in.size());
  ordinary_modes r{.C = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0}),
                   .S = rtepack::muelmat_tensor3(nfourier, no, ni, rtepack::muelmat{0.0})};
  const Index    nmode = std::min(nfourier, nphi / 2);
  const auto     n     = static_cast<Numeric>(nphi);

  // The samples and their trapezoid weights cos(m phi) / nphi and sin(m phi) / nphi
  Vector phi(nphi);
  Matrix cosine(nmode, nphi), sine(nmode, nphi);
  for (Index k = 0; k < nphi; k++) {
    phi[k] = (Constant::two_pi * (static_cast<Numeric>(k) + (midpoints ? 0.5 : 0.0))) / n;
    for (Index m = 0; m < nmode; m++) {
      cosine[m, k] = std::cos(static_cast<Numeric>(m) * phi[k]) / n;
      sine[m, k]   = std::sin(static_cast<Numeric>(m) * phi[k]) / n;
    }
  }

  for (Index io = 0; io < no; io++) {
    for (Index ii = 0; ii < ni; ii++) {
      for (Index k = 0; k < nphi; k++) {
        auto Z = lab_frame(F, mu_in[ii], azimuth_sign * phi[k], mu_out[io], 0.0);
        for (Index a = 0; a < 4; a++)
          for (Index b = 0; b < 4; b++)
            if (a >= nstokes or b >= nstokes) Z[a, b] = 0.0;
        for (Index m = 0; m < nmode; m++) {
          r.C[m, io, ii] += cosine[m, k] * Z;
          r.S[m, io, ii] += sine[m, k] * Z;
        }
      }
    }
  }
  return r;
}
}  // namespace vdisort_test
