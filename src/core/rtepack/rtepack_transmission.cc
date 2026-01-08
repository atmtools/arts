#include "rtepack_transmission.h"

#include <arts_constants.h>
#include <arts_constexpr_math.h>
#include <arts_omp.h>

#include <Faddeeva.hh>
#include <algorithm>
#include <cmath>
#include <limits>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_spectral_matrix.h"

namespace rtepack {
static constexpr Numeric too_small = 1e-4;

tran::tran(const propmat &k1, const propmat &k2, const Numeric r)
    : a{-0.5 * r * (k1.A() + k2.A())},
      exp_a{std::exp(a)},
      polarized(k1.is_polarized() or k2.is_polarized()) {
  if (not polarized) return;

  b = -0.5 * r * (k1.B() + k2.B());
  c = -0.5 * r * (k1.C() + k2.C());
  d = -0.5 * r * (k1.D() + k2.D());
  u = -0.5 * r * (k1.U() + k2.U());
  v = -0.5 * r * (k1.V() + k2.V());
  w = -0.5 * r * (k1.W() + k2.W());

  b2 = b * b;
  c2 = c * c;
  d2 = d * d;
  u2 = u * u;
  v2 = v * v;
  w2 = w * w;

  /** Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
  B = u2 + v2 + w2 - b2 - c2 - d2;
  C = -Math::pow2(d * u - c * v + b * w);
  S = std::sqrt(B * B - 4 * C);

  /**
        We define: 
            x2 and y2 are the squares of x and y
            x and y are real and positive
            x is from the real part of the Eigenvalues
            y is from the imag part of the Eigenvalues
        Notes:
             S  >=  0
            |S| >= |B|
            S-B >=  0
            S+B >=  0
            The y2 sqrt is without the minus sign to avoid complex numbers
    */
  x2 = std::sqrt(0.5 * (S - B));
  y2 = std::sqrt(0.5 * (S + B));
  x  = std::sqrt(x2);
  y  = std::sqrt(y2);

  cy = std::cos(y);
  sy = std::sin(y);
  cx = std::cosh(x);
  sx = std::sinh(x);

  x_zero      = x < too_small;
  y_zero      = y < too_small;
  both_zero   = y_zero and x_zero;
  either_zero = y_zero or x_zero;

  /* Using:
     *    lim x→0 [({cosh(x),cos(x)} - 1) / x^2] → 1/2
     *    lim x→0 [{sinh(x),sin(x)} / x]  → 1
     *    inv_x2 := 1 for x == 0,
     *    -i sin(ix) → sinh(x)
     *    cos(ix) → cosh(x)
     *    C0, C1, C2 ∝ [1/x^2]
     */
  ix = x_zero ? 0.0 : 1.0 / x;
  iy = y_zero ? 0.0 : 1.0 / y;

  // The first "1.0" is the trick for above limits
  inv_x2y2 = both_zero ? 1.0 : 1.0 / (x2 + y2);

  /* FIXME:
    C0 can expand as (1.0 + x2 * y2 / 24.0)
    C1 can expand as (1.0 + x2 * y2 / 120.0)
    C2 can expand as (0.5 + (x2 - y2) / 24.0)
    C3 can expand as (1.0/6.0 + (x2 - y2) / 120.0)
    to avoid numerical issues when x2 and y2 are small.

    The derivatives are not as simple and may require separate expansions for numerical stability.
  */

  C0 = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
  C1 = either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
  C2 = both_zero ? 0.5 : (cx - cy) * inv_x2y2;
  C3 = both_zero
           ? 1.0 / 6.0
           : ((x_zero ? 1.0 : sx * ix) - (y_zero ? 1.0 : sy * iy)) * inv_x2y2;
}

muelmat tran::operator()() const noexcept {
  if (not polarized) return exp_a;

  // Do the calculation exp(a) * (C0 * I + C1 * K + C2 * K^2 + C3 * K^3)

  const Numeric C2b = C2 * (c * u + d * v);
  const Numeric C2c = C2 * (b * u - d * w);
  const Numeric C2d = C2 * (b * v + c * w);
  const Numeric C2u = C2 * (b * c - v * w);
  const Numeric C2v = C2 * (b * d + u * w);
  const Numeric C2w = C2 * (c * d - u * v);

  // B = u2 + v2 + w2 - b2 - c2 - d2
  const Numeric C3b = C3 * (b * (B - w2) + w * (c * v - d * u));
  const Numeric C3c = C3 * (c * (v2 - B) - v * (d * u + b * w));
  const Numeric C3d = C3 * (d * (u2 - B) - u * (c * v - b * w));
  const Numeric C3u = C3 * (d * (c * v - b * w) - u * (B + d2));
  const Numeric C3v = C3 * (c * (d * u + b * w) - v * (B + c2));
  const Numeric C3w = C3 * (b * (c * v - d * u) - w * (B + b2));

  const Numeric M00 = C0 + C2 * (b2 + c2 + d2);
  const Numeric M11 = C0 + C2 * (b2 - u2 - v2);
  const Numeric M22 = C0 + C2 * (c2 - u2 - w2);
  const Numeric M33 = C0 + C2 * (d2 - v2 - w2);

  // clang-format off
  return exp_a * muelmat{
    M00,                  C1 * b - C2b - C3b,   C1 * c + C2c + C3c, C1 * d + C2d + C3d,
    C1 * b + C2b - C3b,   M11,                  C1 * u + C2u + C3u, C1 * v + C2v + C3v,
    C1 * c - C2c + C3c, - C1 * u + C2u - C3u,   M22,                C1 * w + C2w + C3w,
    C1 * d - C2d + C3d, - C1 * v + C2v - C3v, - C1 * w + C2w - C3w, M33};
  // clang-format on
}

muelmat tran::expm1() const noexcept {
  if (not polarized) return {std::expm1(a)};

  // 1. Compute (C0 - 1) robustly
  Numeric C0_m1{};
  if (both_zero) {
    // Taylor expansion: C0 approx 1 + x^2*y^2/24
    C0_m1 = x2 * y2 / 24.0;
  } else {
    // Robust formula using half-angles to avoid cancellation in (cosh(x)-1) and (cos(y)-1)
    const Numeric sh_x2 = std::sinh(0.5 * x);
    const Numeric si_y2 = std::sin(0.5 * y);
    C0_m1 = 2.0 * (y2 * sh_x2 * sh_x2 - x2 * si_y2 * si_y2) * inv_x2y2;
  }

  // 2. Compute the diagonal offset: e^a * C0 - 1
  //    = (expm1(a) + 1) * (C0_m1 + 1) - 1
  //    = expm1(a) * (C0_m1 + 1) + C0_m1
  const Numeric em1_a       = std::expm1(a);
  const Numeric diag_offset = em1_a * (C0_m1 + 1.0) + C0_m1;

  // 3. Compute the rest of the matrix terms
  //    The result is: diag_offset * I + exp_a * (C1*K + C2*K^2 + C3*K^3)
  //    Note that K and K^3 have zero diagonals, but K^2 has diagonal terms.

  const Numeric C2b = C2 * (c * u + d * v);
  const Numeric C2c = C2 * (b * u - d * w);
  const Numeric C2d = C2 * (b * v + c * w);
  const Numeric C2u = C2 * (b * c - v * w);
  const Numeric C2v = C2 * (b * d + u * w);
  const Numeric C2w = C2 * (c * d - u * v);

  const Numeric C3b = C3 * (b * (B - w2) + w * (c * v - d * u));
  const Numeric C3c = C3 * (c * (v2 - B) - v * (d * u + b * w));
  const Numeric C3d = C3 * (d * (u2 - B) - u * (c * v - b * w));
  const Numeric C3u = C3 * (d * (c * v - b * w) - u * (B + d2));
  const Numeric C3v = C3 * (c * (d * u + b * w) - v * (B + c2));
  const Numeric C3w = C3 * (b * (c * v - d * u) - w * (B + b2));

  // Diagonal contributions from K^2 term
  const Numeric M00_rest = C2 * (b2 + c2 + d2);
  const Numeric M11_rest = C2 * (b2 - u2 - v2);
  const Numeric M22_rest = C2 * (c2 - u2 - w2);
  const Numeric M33_rest = C2 * (d2 - v2 - w2);

  // clang-format off
  return muelmat{
    diag_offset + exp_a * M00_rest,       exp_a * (C1 * b - C2b - C3b),         exp_a * (C1 * c + C2c + C3c),       exp_a * (C1 * d + C2d + C3d),
    exp_a * (C1 * b + C2b - C3b),         diag_offset + exp_a * M11_rest,       exp_a * (C1 * u + C2u + C3u),       exp_a * (C1 * v + C2v + C3v),
    exp_a * (C1 * c - C2c + C3c),         exp_a * (-C1 * u + C2u - C3u),        diag_offset + exp_a * M22_rest,     exp_a * (C1 * w + C2w + C3w),
    exp_a * (C1 * d - C2d + C3d),         exp_a * (-C1 * v + C2v - C3v),        exp_a * (-C1 * w + C2w - C3w),      diag_offset + exp_a * M33_rest
  };
  // clang-format on
}

muelmat tran::linsrc() const noexcept {
  // K = - k, expm1 = exp(-k) - 1, so "-" cancels out,
  // since we want k^-1 (1 - exp(-k)) = K^-1 * expm1

  if (not polarized) return std::expm1(a) / a;

  const propmat K{a, b, c, d, u, v, w};

  return inv(K) * expm1();
}

muelmat tran::linsrc_deriv(const muelmat &l,
                           const propmat &dk,
                           const muelmat &dt,
                           const Numeric r,
                           const Numeric dr) const {
  if (not polarized) {
    const Numeric da = (dr / r) * a - 0.5 * r * dk.A();
    return muelmat{(dt[0, 0] - l[0, 0] * da) / a};
  }

  const propmat K{a, b, c, d, u, v, w};

  return inv(K) * ((0.5 * r * dk - (dr / r) * K) * l + dt);
}

muelmat tran::linsrc_linprop(const muelmat &t,
                             const propmat &k1,
                             const propmat &k2,
                             const Numeric r) const noexcept {
  using Faddeeva::Dawson;

  const propmat alpha2 = (k2 - k1) / (2.0 * r);

  // Ignore when the gradient is negative
  if (alpha2.A() < 1e-8) return linsrc();

  if (not polarized) {
    const Numeric alpha = std::sqrt(alpha2.A());
    const Numeric u0    = k1.A() / (2.0 * alpha);
    const Numeric u1    = k2.A() / (2.0 * alpha);

    return (Dawson(u1) - t[0, 0] * Dawson(u0)) / (r * alpha);
  }

  const specmat alpha     = sqrt(alpha2);
  const specmat alpha_inv = inv(alpha);
  const specmat u0        = alpha_inv * (k1 / 2.0);
  const specmat u1        = alpha_inv * (k2 / 2.0);

  return real(alpha_inv * (dawson(u1) - t * dawson(u0))) / r;
}

muelmat tran::linsrc_linprop_deriv(const muelmat &lambda,
                                   const muelmat &t,
                                   const propmat &k1,
                                   const propmat &k2,
                                   const propmat &dk_in,
                                   const muelmat &dt,
                                   const Numeric r,
                                   const Numeric dr,
                                   bool k1_deriv) const {
  using Faddeeva::Dawson;

  const propmat alpha2 = (k2 - k1) / (2.0 * r);

  if (alpha2.A() < 1e-8) return linsrc_deriv(lambda, dk_in, dt, r, dr);

  if (not polarized) {
    const Numeric k1a = k1.A();
    const Numeric k2a = k2.A();
    const Numeric dk  = dk_in.A();

    const Numeric delta_k = k2a - k1a;
    const Numeric denom   = 2.0 * r;
    const Numeric alpha   = std::sqrt(std::max(0.0, delta_k / denom));

    const Numeric u0 = k1a / (2.0 * alpha);
    const Numeric u1 = k2a / (2.0 * alpha);

    const Numeric D0  = Dawson(u0);
    const Numeric D1  = Dawson(u1);
    const Numeric dD0 = 1.0 - 2.0 * u0 * D0;
    const Numeric dD1 = 1.0 - 2.0 * u1 * D1;

    const Numeric t00  = t[0, 0];
    const Numeric dt00 = dt[0, 0];

    Numeric d_alpha = 0.0, d_u0 = 0.0, d_u1 = 0.0;
    if (k1_deriv) {
      d_alpha = -0.5 * dk / (denom * alpha);
      d_u0 = (dk * 2.0 * alpha - k1a * 2.0 * d_alpha) / (4.0 * alpha * alpha);
      d_u1 = -k2a * d_alpha / (2.0 * alpha * alpha);
    } else {
      d_alpha = 0.5 * dk / (denom * alpha);
      d_u0    = -k1a * d_alpha / (2.0 * alpha * alpha);
      d_u1 = (dk * 2.0 * alpha - k2a * 2.0 * d_alpha) / (4.0 * alpha * alpha);
    }

    const Numeric d_num =
        dD1 * d_u1 - dt00 * D0 - t00 * dD0 * d_u0 - t00 * D0 * d_u0;

    const Numeric denom_val   = r * alpha;
    const Numeric d_denom_val = dr * alpha + r * d_alpha;
    const Numeric result = (d_num * denom_val - (D1 - t00 * D0) * d_denom_val) /
                           (denom_val * denom_val);

    return muelmat{result};
  }

  // These derivaties don't work so we use perturbations...

  constexpr Numeric eps = 1e-6;

  if (k1_deriv) {
    const propmat k1p = k1 + dk_in * eps;
    const Numeric rp  = r + dr * eps;

    const tran tran_p{k1p, k2, rp};
    const muelmat tp      = tran_p();
    const muelmat lambdap = tran_p.linsrc_linprop(tp, k1p, k2, rp);

    return (lambdap - lambda) / eps;
  }

  const propmat k2p = k2 + dk_in * eps;
  const Numeric rp  = r + dr * eps;

  const tran tran_p{k1, k2p, rp};
  const muelmat tp      = tran_p();
  const muelmat lambdap = tran_p.linsrc_linprop(tp, k1, k2p, rp);

  return (lambdap - lambda) / eps;
}

muelmat tran::deriv(const muelmat &t,
                    const propmat &k1,
                    const propmat &k2,
                    const propmat &dk,
                    const Numeric r,
                    const Numeric dr) const {
  const Numeric da = -0.5 * (r * dk.A() + dr * (k1.A() + k2.A()));
  if (not polarized) return {da * exp_a};

  const Numeric db = -0.5 * (r * dk.B() + dr * (k1.B() + k2.B()));
  const Numeric dc = -0.5 * (r * dk.C() + dr * (k1.C() + k2.C()));
  const Numeric dd = -0.5 * (r * dk.D() + dr * (k1.D() + k2.D()));
  const Numeric du = -0.5 * (r * dk.U() + dr * (k1.U() + k2.U()));
  const Numeric dv = -0.5 * (r * dk.V() + dr * (k1.V() + k2.V()));
  const Numeric dw = -0.5 * (r * dk.W() + dr * (k1.W() + k2.W()));

  const Numeric db2 = 2 * db * b;
  const Numeric dc2 = 2 * dc * c;
  const Numeric dd2 = 2 * dd * d;
  const Numeric du2 = 2 * du * u;
  const Numeric dv2 = 2 * dv * v;
  const Numeric dw2 = 2 * dw * w;

  /* Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
  const Numeric dB = du2 + dv2 + dw2 - db2 - dc2 - dd2;
  const Numeric dC = -2 * (d * u - c * v + b * w) *
                     (dd * u + d * du - dc * v - c * dv + db * w + b * dw);
  const Numeric dS = (B * dB - 2 * dC) / S;

  const Numeric dx2 = 0.25 * (dS - dB) / x2;
  const Numeric dy2 = 0.25 * (dS + dB) / y2;
  const Numeric dx  = 0.5 * dx2 / x;
  const Numeric dy  = 0.5 * dy2 / y;

  const Numeric dcy    = -sy * dy;
  const Numeric dsy    = cy * dy;
  const Numeric dcx    = sx * dx;
  const Numeric dsx    = cx * dx;
  const Numeric dix    = -dx * ix * ix;
  const Numeric diy    = -dy * iy * iy;
  const Numeric dx2dy2 = dx2 + dy2;

  const Numeric dC0 =
      either_zero ? 0.0
                  : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0 * dx2dy2) *
                        inv_x2y2;
  const Numeric dC1 =
      either_zero
          ? 0.0
          : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy + dsx * y2 * ix +
             sx * dy2 * ix + sx * y2 * dix - C1 * dx2dy2) *
                inv_x2y2;
  const Numeric dC2 = both_zero ? 0.0
                                : ((x_zero ? 0.0 : (dcx - C2 * dx2)) -
                                   (y_zero ? 0.0 : (dcy + C2 * dy2))) *
                                      inv_x2y2;
  const Numeric dC3 =
      both_zero ? 0.0
                : ((x_zero ? 0.0 : (dsx * ix + sx * dix - C3 * dx2)) -
                   (y_zero ? 0.0 : (dsy * iy + sy * diy + C3 * dy2))) *
                      inv_x2y2;

  const Numeric dC2b =
      dC2 * (c * u + d * v) + C2 * (dc * u + c * du + dd * v + d * dv);
  const Numeric dC2c =
      dC2 * (b * u - d * w) + C2 * (db * u + b * du - dd * w - d * dw);
  const Numeric dC2d =
      dC2 * (b * v + c * w) + C2 * (db * v + b * dv + dc * w + c * dw);
  const Numeric dC2u =
      dC2 * (b * c - v * w) + C2 * (db * c + b * dc - dv * w - v * dw);
  const Numeric dC2v =
      dC2 * (b * d + u * w) + C2 * (db * d + b * dd + du * w + u * dw);
  const Numeric dC2w =
      dC2 * (c * d - u * v) + C2 * (dc * d + c * dd - du * v - u * dv);

  const Numeric dC3b =
      dC3 * (b * (B - w2) + w * (c * v - d * u)) +
      C3 * (db * (B - w2) + b * (dB - dw2) + dw * (c * v - d * u) +
            w * (dc * v + c * dv - dd * u - d * du));
  const Numeric dC3c =
      dC3 * (c * (v2 - B) - v * (d * u + b * w)) +
      C3 * (dc * (v2 - B) + c * (dv2 - dB) - dv * (d * u + b * w) -
            v * (dd * u + d * du + db * w + b * dw));
  const Numeric dC3d =
      dC3 * (d * (u2 - B) - u * (c * v - b * w)) +
      C3 * (dd * (u2 - B) + d * (du2 - dB) - du * (c * v - b * w) -
            u * (dc * v + c * dv - db * w - b * dw));
  const Numeric dC3u =
      dC3 * (d * (c * v - b * w) - u * (B + d2)) +
      C3 * (dd * (c * v - b * w) + d * (dc * v + c * dv - db * w - b * dw) -
            du * (B + d2) - u * (dB + dd2));
  const Numeric dC3v =
      dC3 * (c * (d * u + b * w) - v * (B + c2)) +
      C3 * (dc * (d * u + b * w) + c * (dd * u + d * du + db * w + b * dw) -
            dv * (B + c2) - v * (dB + dc2));
  const Numeric dC3w =
      dC3 * (b * (c * v - d * u) - w * (B + b2)) +
      C3 * (db * (c * v - d * u) + b * (dc * v + c * dv - dd * u - d * du) -
            dw * (B + b2) - w * (dB + db2));

  const Numeric dM00 = dC0 + dC2 * (b2 + c2 + d2) + C2 * (db2 + dc2 + dd2);
  const Numeric dM11 = dC0 + dC2 * (b2 - u2 - v2) + C2 * (db2 - du2 - dv2);
  const Numeric dM22 = dC0 + dC2 * (c2 - u2 - w2) + C2 * (dc2 - du2 - dw2);
  const Numeric dM33 = dC0 + dC2 * (d2 - v2 - w2) + C2 * (dd2 - dv2 - dw2);

  // clang-format off
  return da * t + exp_a * muelmat{
    dM00,                              dC1 * b + C1 * db - dC2b - dC3b,   dC1 * c + C1 * dc + dC2c + dC3c, dC1 * d + C1 * dd + dC2d + dC3d,
    dC1 * b + C1 * db + dC2b - dC3b,   dM11,                              dC1 * u + C1 * du + dC2u + dC3u, dC1 * v + C1 * dv + dC2v + dC3v,
    dC1 * c + C1 * dc - dC2c + dC3c, - dC1 * u - C1 * du + dC2u - dC3u,   dM22,                            dC1 * w + C1 * dw + dC2w + dC3w,
    dC1 * d + C1 * dd - dC2d + dC3d, - dC1 * v - C1 * dv + dC2v - dC3v, - dC1 * w - C1 * dw + dC2w - dC3w, dM33};
  // clang-format on
}

void two_level_exp(muelmat &t,
                   muelmat_vector_view dt1,
                   muelmat_vector_view dt2,
                   const propmat &k1,
                   const propmat &k2,
                   const propmat_vector_const_view &dk1,
                   const propmat_vector_const_view &dk2,
                   const Numeric r,
                   const ConstVectorView &dr1,
                   const ConstVectorView &dr2) {
  assert(dk1.size() == dk2.size() and dk1.size() == dr1.size() and
         dk1.size() == dr2.size());

  const tran tran_state{k1, k2, r};
  t = tran_state();

  const auto deriv = [&](const propmat &dk, const Numeric &dr) -> muelmat {
    return tran_state.deriv(t, k1, k2, dk, r, dr);
  };

  std::transform(dk1.begin(), dk1.end(), dr1.begin(), dt1.begin(), deriv);
  std::transform(dk2.begin(), dk2.end(), dr2.begin(), dt2.begin(), deriv);
}

namespace {
void two_level_exp(muelmat_vector_view tv,
                   muelmat_matrix_view dt1v,
                   muelmat_matrix_view dt2v,
                   const propmat_vector_const_view &k1v,
                   const propmat_vector_const_view &k2v,
                   const propmat_matrix_const_view &dk1v,
                   const propmat_matrix_const_view &dk2v,
                   const Numeric rv,
                   const ConstVectorView &dr1v,
                   const ConstVectorView &dr2v) {
  const Size nf = tv.size();
  const Size nq = dr1v.size();

  assert(nf == k1v.size());
  assert(nf == k2v.size());
  assert(nf == static_cast<Size>(dk1v.ncols()));
  assert(nf == static_cast<Size>(dk2v.ncols()));
  assert(nq == static_cast<Size>(dk1v.nrows()));
  assert(nq == static_cast<Size>(dk2v.nrows()));
  assert(nf == static_cast<Size>(dt1v.ncols()));
  assert(nf == static_cast<Size>(dt2v.ncols()));
  assert(nq == static_cast<Size>(dt1v.nrows()));
  assert(nq == static_cast<Size>(dt2v.nrows()));
  assert(nq == dr2v.size());

  for (Size i = 0; i < nf; ++i) {
    const tran tran_state{k1v[i], k2v[i], rv};
    tv[i] = tran_state();

    for (Size j = 0; j < nq; j++) {
      dt1v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk1v[j, i], rv, dr1v[j]);
      dt2v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk2v[j, i], rv, dr2v[j]);
    }
  }
}

void two_level_exp_linsrc(muelmat_vector_view tv,
                          muelmat_vector_view lv,
                          muelmat_matrix_view dt1v,
                          muelmat_matrix_view dt2v,
                          muelmat_matrix_view dl1v,
                          muelmat_matrix_view dl2v,
                          const propmat_vector_const_view &k1v,
                          const propmat_vector_const_view &k2v,
                          const propmat_matrix_const_view &dk1v,
                          const propmat_matrix_const_view &dk2v,
                          const Numeric rv,
                          const ConstVectorView &dr1v,
                          const ConstVectorView &dr2v) {
  const Size nf = tv.size();
  const Size nq = dr1v.size();

  assert(nf == k1v.size());
  assert(nf == k2v.size());
  assert(nf == lv.size());
  assert(nf == static_cast<Size>(dk1v.ncols()));
  assert(nf == static_cast<Size>(dk2v.ncols()));
  assert(nq == static_cast<Size>(dk1v.nrows()));
  assert(nq == static_cast<Size>(dk2v.nrows()));
  assert(nf == static_cast<Size>(dt1v.ncols()));
  assert(nf == static_cast<Size>(dt2v.ncols()));
  assert(nq == static_cast<Size>(dt1v.nrows()));
  assert(nq == static_cast<Size>(dt2v.nrows()));
  assert(nf == static_cast<Size>(dl1v.ncols()));
  assert(nf == static_cast<Size>(dl2v.ncols()));
  assert(nq == static_cast<Size>(dl1v.nrows()));
  assert(nq == static_cast<Size>(dl2v.nrows()));
  assert(nq == dr2v.size());

  for (Size i = 0; i < nf; ++i) {
    const tran tran_state{k1v[i], k2v[i], rv};
    tv[i] = tran_state();
    lv[i] = tran_state.linsrc();

    for (Size j = 0; j < nq; j++) {
      dt1v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk1v[j, i], rv, dr1v[j]);
      dt2v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk2v[j, i], rv, dr2v[j]);
      dl1v[j, i] =
          tran_state.linsrc_deriv(lv[i], dk1v[j, i], dt1v[j, i], rv, dr1v[j]);
      dl2v[j, i] =
          tran_state.linsrc_deriv(lv[i], dk2v[j, i], dt2v[j, i], rv, dr2v[j]);
    }
  }
}

void two_level_exp_linsrc_linprop(muelmat_vector_view tv,
                                  muelmat_vector_view lv,
                                  muelmat_matrix_view dt1v,
                                  muelmat_matrix_view dt2v,
                                  muelmat_matrix_view dl1v,
                                  muelmat_matrix_view dl2v,
                                  const propmat_vector_const_view &k1v,
                                  const propmat_vector_const_view &k2v,
                                  const propmat_matrix_const_view &dk1v,
                                  const propmat_matrix_const_view &dk2v,
                                  const Numeric rv,
                                  const ConstVectorView &dr1v,
                                  const ConstVectorView &dr2v) {
  const Size nf = tv.size();
  const Size nq = dr1v.size();

  assert(nf == k1v.size());
  assert(nf == k2v.size());
  assert(nf == lv.size());
  assert(nf == lv.size());
  assert(nf == static_cast<Size>(dk1v.ncols()));
  assert(nf == static_cast<Size>(dk2v.ncols()));
  assert(nq == static_cast<Size>(dk1v.nrows()));
  assert(nq == static_cast<Size>(dk2v.nrows()));
  assert(nf == static_cast<Size>(dt1v.ncols()));
  assert(nf == static_cast<Size>(dt2v.ncols()));
  assert(nq == static_cast<Size>(dt1v.nrows()));
  assert(nq == static_cast<Size>(dt2v.nrows()));
  assert(nf == static_cast<Size>(dl1v.ncols()));
  assert(nf == static_cast<Size>(dl2v.ncols()));
  assert(nq == static_cast<Size>(dl1v.nrows()));
  assert(nq == static_cast<Size>(dl2v.nrows()));
  assert(nq == dr2v.size());

  for (Size i = 0; i < nf; ++i) {
    const tran tran_state{k1v[i], k2v[i], rv};
    tv[i] = tran_state();
    lv[i] = tran_state.linsrc_linprop(tv[i], k1v[i], k2v[i], rv);

    for (Size j = 0; j < nq; j++) {
      dt1v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk1v[j, i], rv, dr1v[j]);
      dt2v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk2v[j, i], rv, dr2v[j]);
      dl1v[j, i] = tran_state.linsrc_linprop_deriv(lv[i],
                                                   tv[i],
                                                   k1v[i],
                                                   k2v[i],
                                                   dk1v[j, i],
                                                   dt1v[j, i],
                                                   rv,
                                                   dr1v[j],
                                                   true);
      dl2v[j, i] = tran_state.linsrc_linprop_deriv(lv[i],
                                                   tv[i],
                                                   k1v[i],
                                                   k2v[i],
                                                   dk2v[j, i],
                                                   dt2v[j, i],
                                                   rv,
                                                   dr2v[j],
                                                   false);
    }
  }
}
}  // namespace

muelmat exp(propmat k, Numeric r) { return tran(k, k, r)(); }

void two_level_exp(muelmat_vector_view tv,
                   const propmat_vector_const_view &k1v,
                   const propmat_vector_const_view &k2v,
                   const Numeric rv) {
  assert(k2v.size() == k1v.size());
  assert(tv.size() == k1v.size());

  std::transform(
      k1v.begin(),
      k1v.end(),
      k2v.begin(),
      tv.begin(),
      [rv](const propmat &a, const propmat &b) { return tran(a, b, rv)(); });
}

void two_level_exp(std::vector<muelmat_vector> &T,
                   std::vector<muelmat_tensor3> &dT,
                   const std::vector<propmat_vector> &K,
                   const std::vector<propmat_matrix> &dK,
                   const Vector &r,
                   const Tensor3 &dr) {
  const Size N = K.size();

  ARTS_USER_ERROR_IF(
      N != dK.size(), "Must have same number of levels ({}) in K and dK", N);

  ARTS_USER_ERROR_IF(
      (N - 1) != static_cast<Size>(r.size()),
      "Must have one fewer layer distances ({}) than levels ({}) in K",
      (N - 1),
      N);

  ARTS_USER_ERROR_IF(
      (N - 1) != static_cast<Size>(dr.nrows()),
      "Must have one fewer layer distance derivatives ({}) than levels ({}) in K",
      N - 1,
      N);

  T.resize(N);

  dT.resize(N);

  if (N == 0) return;

  const Size nv  = K[0].size();
  const Index nq = dr.ncols();

  for (auto &x : T) {
    x.resize(nv);
    x = 1.0;
  }

  for (auto &x : dT) {
    x.resize(2, nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      stdr::any_of(K, Cmp::ne(nv), [](auto &x) { return x.size(); }),
      "Must have same number of frequency elements ({}) in all K:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      stdr::any_of(dK, Cmp::ne(static_cast<Index>(nv)), &propmat_matrix::ncols),
      "Must have same number of frequency elements ({}) in all dK:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      stdr::any_of(dK, Cmp::ne(nq), &propmat_matrix::nrows),
      "Must have same number of derivative elements ({}) in all dK:s as in dr",
      nq);

  ARTS_USER_ERROR_IF(
      dr.npages() != 2,
      "Must have 2 as first dimension in dr (upper and lower level distance derivatives), got {}",
      dr.npages());

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    two_level_exp(T[i],
                  dT[i - 1][0],
                  dT[i][1],
                  K[i - 1],
                  K[i],
                  dK[i - 1],
                  dK[i],
                  r[i - 1],
                  dr[0, i - 1],
                  dr[1, i - 1]);
  }
}

void two_level_exp_linsrc(std::vector<muelmat_vector> &T,
                          std::vector<muelmat_vector> &L,
                          std::vector<muelmat_tensor3> &dT,
                          std::vector<muelmat_tensor3> &dL,
                          const std::vector<propmat_vector> &K,
                          const std::vector<propmat_matrix> &dK,
                          const Vector &r,
                          const Tensor3 &dr) {
  const Size N = K.size();

  ARTS_USER_ERROR_IF(
      N != dK.size(), "Must have same number of levels ({}) in K and dK", N);

  ARTS_USER_ERROR_IF(
      (N - 1) != static_cast<Size>(r.size()),
      "Must have one fewer layer distances ({}) than levels ({}) in K",
      (N - 1),
      N);

  ARTS_USER_ERROR_IF(
      (N - 1) != static_cast<Size>(dr.nrows()),
      "Must have one fewer layer distance derivatives ({}) than levels ({}) in K",
      N - 1,
      N);

  T.resize(N);
  L.resize(N);

  dT.resize(N);
  dL.resize(N);

  if (N == 0) return;

  const Size nv  = K[0].size();
  const Index nq = dr.ncols();

  for (auto &x : T) {
    x.resize(nv);
    x = 1.0;
  }

  for (auto &x : L) {
    x.resize(nv);
    x = 1.0;
  }

  for (auto &x : dT) {
    x.resize(2, nq, nv);
    x = 0.0;
  }

  for (auto &x : dL) {
    x.resize(2, nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      stdr::any_of(K, Cmp::ne(nv), [](auto &x) { return x.size(); }),
      "Must have same number of frequency elements ({}) in all K:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      stdr::any_of(dK, Cmp::ne(static_cast<Index>(nv)), &propmat_matrix::ncols),
      "Must have same number of frequency elements ({}) in all dK:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      stdr::any_of(dK, Cmp::ne(nq), &propmat_matrix::nrows),
      "Must have same number of derivative elements ({}) in all dK:s as in dr",
      nq);

  ARTS_USER_ERROR_IF(
      dr.npages() != 2,
      "Must have 2 as first dimension in dr (upper and lower level distance derivatives), got {}",
      dr.npages());

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    two_level_exp_linsrc(T[i],
                         L[i],
                         dT[i - 1][0],
                         dT[i][1],
                         dL[i - 1][0],
                         dL[i][1],
                         K[i - 1],
                         K[i],
                         dK[i - 1],
                         dK[i],
                         r[i - 1],
                         dr[0, i - 1],
                         dr[1, i - 1]);
  }
}

void two_level_exp_linsrc_linprop(std::vector<muelmat_vector> &T,
                                  std::vector<muelmat_vector> &L,
                                  std::vector<muelmat_tensor3> &dT,
                                  std::vector<muelmat_tensor3> &dL,
                                  const std::vector<propmat_vector> &K,
                                  const std::vector<propmat_matrix> &dK,
                                  const Vector &r,
                                  const Tensor3 &dr) {
  const Size N = K.size();

  ARTS_USER_ERROR_IF(
      N != dK.size(), "Must have same number of levels ({}) in K and dK", N);

  ARTS_USER_ERROR_IF(
      (N - 1) != static_cast<Size>(r.size()),
      "Must have one fewer layer distances ({}) than levels ({}) in K",
      (N - 1),
      N);

  ARTS_USER_ERROR_IF(
      (N - 1) != static_cast<Size>(dr.nrows()),
      "Must have one fewer layer distance derivatives ({}) than levels ({}) in K",
      N - 1,
      N);

  T.resize(N);
  L.resize(N);

  dT.resize(N);
  dL.resize(N);

  if (N == 0) return;

  const Size nv  = K[0].size();
  const Index nq = dr.ncols();

  for (auto &x : T) {
    x.resize(nv);
    x = 1.0;
  }

  for (auto &x : L) {
    x.resize(nv);
    x = 1.0;
  }

  for (auto &x : dT) {
    x.resize(2, nq, nv);
    x = 0.0;
  }

  for (auto &x : dL) {
    x.resize(2, nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      stdr::any_of(K, Cmp::ne(nv), [](auto &x) { return x.size(); }),
      "Must have same number of frequency elements ({}) in all K:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      stdr::any_of(dK, Cmp::ne(static_cast<Index>(nv)), &propmat_matrix::ncols),
      "Must have same number of frequency elements ({}) in all dK:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      stdr::any_of(dK, Cmp::ne(nq), &propmat_matrix::nrows),
      "Must have same number of derivative elements ({}) in all dK:s as in dr",
      nq);

  ARTS_USER_ERROR_IF(
      dr.npages() != 2,
      "Must have 2 as first dimension in dr (upper and lower level distance derivatives), got {}",
      dr.npages());

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 1; i < N; i++) {
    two_level_exp_linsrc_linprop(T[i],
                                 L[i],
                                 dT[i - 1][0],
                                 dT[i][1],
                                 dL[i - 1][0],
                                 dL[i][1],
                                 K[i - 1],
                                 K[i],
                                 dK[i - 1],
                                 dK[i],
                                 r[i - 1],
                                 dr[0, i - 1],
                                 dr[1, i - 1]);
  }
}

propmat logK(const muelmat &m) {
  if (not m.is_polarized()) return std::log(midtr(m));

  /**
    The code tries to retrieve K from exp(K) = M input as muelmat m,
    where K is a propmat.

    K is unknown but, as a reminder, it looks like:
      K = [ a    b    c    d
            b    a    u    v
            c   -u    a    w
            d   -v   -w    a ],

    where a, b, c, d, u, v, w are real numbers.

    Since det(e^A) = e^tr(A), we can find a from the determinant of m,
    using tr(K) = 4 a.
   */

  // Small fix: in cases this matter, logK is not useful
  const Numeric det_m = det(m);
  const Numeric a     = (det_m > std::numeric_limits<Numeric>::min())
                            ? 0.25 * std::log(det_m)
                            : 0.0;

  /**
    We use this knowing
        m = exp(K)
          = exp(a) * (C0 I + C1 K + C2 K^2 + C3 K^3)

    where C0, C1, C2, C3 are the Cayley Hamilton coefficients derived in tran::tran.
    
    This means that we can extract
        T = m * exp(-a)
          = C0 I + C1 K + C2 K^2 + C3 K^3
    using the same logic and notation.
   */

  const muelmat T = m * std::exp(-a);

  /**
    These coefficients C0, C1, C2, C3 are found from the eigenvalues of K, which
    have the form: lambda = exp(x), exp(-x), exp(iy), exp(-iy), where x and y are real.

    Since the trace of a matrix is the sum of its eigenvalues, and the eigenvalues
    of a squared matrix are the squares of the eigenvalues of the original matrix,
    we can form these two relations from our known T matrix:
      S1 = tr(T)
         = e^x + e^-x + e^iy + e^-iy
         = 2 cosh(x) + 2 cos(y)
    and
      S2 = tr(T^2)
         = (e^x + e^-x)^2 + (e^iy + e^-iy)^2
         = 2 cosh(2x) + 2 cos(2y)
         = 2 (2 cosh^2(x) - 1) + 2 (2 cos^2(y) - 1)

    Letting u = cosh(x), v = cos(y), we can write:
        S1 = 2 u + 2 v  ->  
            u + v = S1 / 2
    and
        S2 = 2 (2 u^2 - 1) + 2 (2 v^2 - 1)
           = 4 u^2 + 4 v^2 - 4 ->
            u^2 + v^2 = (S2 + 4) / 4

    Since
        (u - v)^2 = 2 (u^2 + v^2) - (u + v)^2
    we can solve for u and v:
        (u - v)^2 = 2 * (S2 + 4) / 4 - (S1 / 2)^2
                  = (S2 + 8 - S1^2) / 4
    and write
        D = (u - v)^2
          = (S2 + 8 - S1^2) / 4
   */

  const Numeric S1 = tr(T);
  const Numeric D  = 0.25 * (2 * tr(T * T) + 8 - S1 * S1);

  /**
    We can use D to find u and v:
        u = (u + v + u - v) / 2
          = (S1 / 2 + sqrt(D)) / 2
    and
        v = (u + v - (u - v)) / 2
          = (S1 / 2 - sqrt(D)) / 2
    Since u - v = sqrt(D)
   */

  const Numeric sqrtD = D > 0 ? std::sqrt(D) : 0.0;
  const Numeric u     = std::max((S1 * 0.5 + sqrtD) * 0.5, 1.0);
  const Numeric v     = std::clamp((S1 * 0.5 - sqrtD) * 0.5, -1.0, 1.0);

  /**
    From the definition of u and v we can now define x and y:
        x = acosh(u)
        y = acos(v)
   */

  const Numeric x = std::acosh(u);
  const Numeric y = std::acos(v);

  // From here, we simply compute the coefficients C0, C1, C2, C3 as in tran::tran
  const Numeric x2       = x * x;
  const Numeric y2       = y * y;
  const bool x_zero      = x < too_small;
  const bool y_zero      = y < too_small;
  const bool both_zero   = x_zero && y_zero;
  const bool either_zero = x_zero || y_zero;
  const Numeric cy       = v;
  const Numeric sy       = std::sin(y);
  const Numeric cx       = u;
  const Numeric sx       = std::sinh(x);
  const Numeric ix       = x_zero ? 0.0 : 1.0 / x;
  const Numeric iy       = y_zero ? 0.0 : 1.0 / y;
  const Numeric inv_x2y2 = both_zero ? 1.0 : 1.0 / (x2 + y2);
  const Numeric C0       = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
  const Numeric C1 =
      either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
  // Skipping C2 and C3

  /**
    Now we are just one step away from reconstructing log(T) and thus K from exp(K).

    We have:
      T = C0 I + C1 K + C2 K^2 + C3 K^3

    We write:
        T02 = C0 I + C2 K^2  ->
                K^2 = (T02 - C0 I) / C2
    and
        T13 = C1 K + C3 K^3  ->
                T13 = K (C1 I + C3 K^2)

    So that:
        K = T13 * inv(C1 I + C3 K^2)
          = T13 * inv(C1 I + C3 (T02 - C0 I) / C2)
    
    OK.  We just have to be able to construct T02 and T13 from T.
    This requires some separation of terms in T.

    For T02, 
      [   T[0,0]             (T[0,1]-T[1,0])/2   (T[0,2]-T[2,0])/2   (T[0,3]-T[3,0])/2
        -(T[1,0]-T[0,1])/2    T[1,1]             (T[1,2]+T[2,1])/2   (T[1,3]+T[3,1])/2
        -(T[2,0]-T[0,2])/2   (T[2,1]+T[1,2])/2    T[2,2]             (T[2,3]+T[3,2])/2
        -(T[3,0]-T[0,3])/2   (T[3,1]+T[1,3])/2   (T[3,2]+T[2,3])/2    T[3,3] ]
    For T13,
      [  0                   (T[0,1]+T[1,0])/2    (T[0,2]+T[2,0])/2   (T[0,3]+T[3,0])/2
        (T[1,0]+T[0,1])/2     0                   (T[1,2]-T[2,1])/2   (T[1,3]-T[3,1])/2
        (T[2,0]+T[0,2])/2   -(T[2,1]-T[1,2])/2     0                  (T[2,3]-T[3,2])/2
        (T[3,0]+T[0,3])/2   -(T[3,1]-T[1,3])/2   -(T[3,2]-T[2,3])/2    0 ]
    
    (To make the above symmetries and anti-symmetries clear, the return value of
     tran::operator() has been formatted to demonstrate them.)
   */
  const Numeric factor =
      (both_zero or sqrtD <= 0)
          ? 1.0 / 3.0
          : ((x_zero ? 1.0 : sx * ix) - (y_zero ? 1.0 : sy * iy)) / sqrtD;
  const Numeric z01 = factor * std::midpoint(T[0, 1], -T[1, 0]);
  const Numeric z02 = factor * std::midpoint(T[0, 2], -T[2, 0]);
  const Numeric z03 = factor * std::midpoint(T[0, 3], -T[3, 0]);
  const Numeric z12 = factor * std::midpoint(T[1, 2], T[2, 1]);
  const Numeric z13 = factor * std::midpoint(T[1, 3], T[3, 1]);
  const Numeric z23 = factor * std::midpoint(T[2, 3], T[3, 2]);
  const Numeric z00 = factor * (T[0, 0] - C0) + C1;
  const Numeric z11 = factor * (T[1, 1] - C0) + C1;
  const Numeric z22 = factor * (T[2, 2] - C0) + C1;
  const Numeric z33 = factor * (T[3, 3] - C0) + C1;
  const Numeric b01 = std::midpoint(T[0, 1], T[1, 0]);
  const Numeric b02 = std::midpoint(T[0, 2], T[2, 0]);
  const Numeric b03 = std::midpoint(T[0, 3], T[3, 0]);
  const Numeric b12 = std::midpoint(T[1, 2], -T[2, 1]);
  const Numeric b13 = std::midpoint(T[1, 3], -T[3, 1]);
  const Numeric b23 = std::midpoint(T[2, 3], -T[3, 2]);

  // clang-format off
  const muelmat Kp = muelmat{ 0.0,  b01,  b02, b03,
                              b01,  0.0,  b12, b13,
                              b02, -b12,  0.0, b23,
                              b03, -b13, -b23, 0.0} *
              inv(muelmat{ z00,  z01,  z02, z03,
                             -z01,  z11,  z12, z13,
                             -z02,  z12,  z22, z23,
                             -z03,  z13,  z23, z33});
  // clang-format on

  return {a,
          std::midpoint(Kp[0, 1], Kp[1, 0]),
          std::midpoint(Kp[0, 2], Kp[2, 0]),
          std::midpoint(Kp[0, 3], Kp[3, 0]),
          std::midpoint(Kp[1, 2], -Kp[2, 1]),
          std::midpoint(Kp[1, 3], -Kp[3, 1]),
          std::midpoint(Kp[2, 3], -Kp[3, 2])};
}

specmat sqrt(const propmat &pm) {
  const Numeric a      = pm.A();
  const Complex sqrt_a = std::sqrt(Complex(a));

  if (not pm.is_polarized()) return sqrt_a;

  const Numeric b = pm.B();
  const Numeric c = pm.C();
  const Numeric d = pm.D();
  const Numeric u = pm.U();
  const Numeric v = pm.V();
  const Numeric w = pm.W();

  const Numeric b2 = b * b;
  const Numeric c2 = c * c;
  const Numeric d2 = d * d;
  const Numeric u2 = u * u;
  const Numeric v2 = v * v;
  const Numeric w2 = w * w;

  Complex d0c{}, d1c{}, d2c{}, d3c{};

  constexpr auto is_small = [](auto... v) {
    return ((v <= std::numeric_limits<Numeric>::epsilon()) and ...);
  };

  if (pm.is_rotational()) {
    const Numeric rho = std::hypot(u, v, w);
    if (is_small(rho)) return {0.0};

    const Numeric r = std::sqrt(2.0 * rho);
    d0c             = 0.0;
    d1c             = 1.0 / r;
    d2c             = -1.0 / (rho * r);
    d3c             = 0.0;
  } else {
    /** Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
    const Numeric B = u2 + v2 + w2 - b2 - c2 - d2;
    const Numeric C = -Math::pow2(d * u - c * v + b * w);
    const Numeric S = std::sqrt(B * B - 4 * C);

    const Numeric x2     = std::max(0.0, 0.5 * (S - B));
    const Numeric abs_y2 = std::max(0.0, 0.5 * (S + B));
    const Numeric x      = std::sqrt(x2);
    const Complex y      = Complex(0, std::sqrt(abs_y2));
    const Complex sx     = std::sqrt(Complex(a + x));
    const Complex dx     = std::sqrt(Complex(a - x));
    const Complex sy     = std::sqrt(a + y);
    const Complex dy     = std::sqrt(a - y);
    const Complex Sx     = sx + dx;
    const Complex Dx     = sx - dx;
    const Complex Sy     = sy + dy;
    const Complex Dy     = sy - dy;

    if (is_small(x2 + abs_y2)) {
      d0c = sqrt_a;
      if (is_small(a)) {
        d1c = 0.5 / sqrt_a;
        d2c = 0.125 / (a * sqrt_a);
        d3c = 0.0625 / (a * a * sqrt_a);
      }
    } else {
      const Numeric inv_sum_sq = 1.0 / (x2 + abs_y2);

      d0c = (abs_y2 * Sx + x2 * Sy) * (0.5 * inv_sum_sq);
      d2c = (Sx - Sy) * (0.5 * inv_sum_sq);

      const Complex term1 = is_small(x, a) ? 0.0
                            : is_small(x)  ? 0.5 / sqrt_a
                                           : 0.5 * Dx / x;
      const Complex term2 = is_small(abs_y2, a) ? 0.0
                            : is_small(abs_y2)  ? 0.5 / sqrt_a
                                                : 0.5 * Dy / y;

      d1c = (abs_y2 * term1 + x2 * term2) * inv_sum_sq;
      d3c = (term1 - term2) * inv_sum_sq;
    }
  }

  // Terms for K0^2
  const Numeric k2_00 = b2 + c2 + d2;
  const Numeric k2_11 = b2 - u2 - v2;
  const Numeric k2_22 = c2 - u2 - w2;
  const Numeric k2_33 = d2 - v2 - w2;

  const Numeric k2_01_val = -(c * u + d * v);
  const Numeric k2_02_val = b * u - d * w;
  const Numeric k2_03_val = b * v + c * w;
  const Numeric k2_12_val = b * c - v * w;
  const Numeric k2_13_val = b * d + u * w;
  const Numeric k2_23_val = c * d - u * v;

  // Terms for K0^3
  const Numeric k3_01_val = b * k2_00 - u * k2_02_val - v * k2_03_val;
  const Numeric k3_02_val = c * k2_00 + u * k2_01_val - w * k2_03_val;
  const Numeric k3_03_val = d * k2_00 + v * k2_01_val + w * k2_02_val;
  const Numeric k3_12_val = -c * k2_01_val + u * k2_11 - w * k2_13_val;
  const Numeric k3_13_val = -d * k2_01_val + v * k2_11 + w * k2_12_val;
  const Numeric k3_23_val = -d * k2_02_val + v * k2_12_val + w * k2_22;

  specmat K;

  K[0, 0] = d0c + d2c * k2_00;
  K[1, 1] = d0c + d2c * k2_11;
  K[2, 2] = d0c + d2c * k2_22;
  K[3, 3] = d0c + d2c * k2_33;

  K[0, 1] = d1c * b + d2c * k2_01_val + d3c * k3_01_val;
  K[1, 0] = d1c * b - d2c * k2_01_val + d3c * k3_01_val;

  K[0, 2] = d1c * c + d2c * k2_02_val + d3c * k3_02_val;
  K[2, 0] = d1c * c - d2c * k2_02_val + d3c * k3_02_val;

  K[0, 3] = d1c * d + d2c * k2_03_val + d3c * k3_03_val;
  K[3, 0] = d1c * d - d2c * k2_03_val + d3c * k3_03_val;

  K[1, 2] = d1c * u + d2c * k2_12_val + d3c * k3_12_val;
  K[2, 1] = -d1c * u + d2c * k2_12_val - d3c * k3_12_val;

  K[1, 3] = d1c * v + d2c * k2_13_val + d3c * k3_13_val;
  K[3, 1] = -d1c * v + d2c * k2_13_val - d3c * k3_13_val;

  K[2, 3] = d1c * w + d2c * k2_23_val + d3c * k3_23_val;
  K[3, 2] = -d1c * w + d2c * k2_23_val - d3c * k3_23_val;

  return K;
}
}  // namespace rtepack
