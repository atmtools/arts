#include "rtepack_transmission.h"

#include <algorithm>
#include <functional>

#include "arts_constants.h"
#include "matpack_complex.h"
#include "matpack_concepts.h"
#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"

namespace rtepack {
static constexpr Numeric lower_is_considered_zero_for_sinc_likes = 1e-4;
static constexpr Numeric sqrt_05 = Constant::inv_sqrt_2;

struct tran {
  Numeric a{}, b{}, c{}, d{}, u{}, v{}, w{};   // To not repeat input
  Numeric exp_a{};                             // To not repeat exp(a)
  Numeric b2{}, c2{}, d2{}, u2{}, v2{}, w2{};  // To shorten expressions
  Numeric B, C, S;  // From L^4 + BL^2 + C = 0; S = sqrt(B^2 - 4C)
  Numeric x2{}, y2{}, x{}, y{}, cy{}, sy{}, cx{},
      sx{};  // Eigenvalues and their used trigonometric functions
  Numeric ix{}, iy{}, inv_x2y2{};  // Computational helpers
  Numeric C0{}, C1{}, C2{}, C3{};  // The Cayley-Hamilton coefficients
  bool unpolarized{}, x_zero{}, y_zero{}, both_zero{}, either_zero{};

  constexpr tran() = default;

  tran(const propmat &k1, const propmat &k2, const Numeric r) {
    a = -0.5 * r * (k1.A() + k2.A());
    b = -0.5 * r * (k1.B() + k2.B());
    c = -0.5 * r * (k1.C() + k2.C());
    d = -0.5 * r * (k1.D() + k2.D());
    u = -0.5 * r * (k1.U() + k2.U());
    v = -0.5 * r * (k1.V() + k2.V());
    w = -0.5 * r * (k1.W() + k2.W());
    exp_a = std::exp(a);

    unpolarized =
        b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.;
    if (unpolarized) {
      return;
    }

    b2 = b * b;
    c2 = c * c;
    d2 = d * d;
    u2 = u * u;
    v2 = v * v;
    w2 = w * w;

    /* Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
    B = u2 + v2 + w2 - b2 - c2 - d2;
    C = -Math::pow2(d * u - c * v + b * w);
    S = std::sqrt(B * B - 4 * C);

    // One will be complex, so we need to take the square root of a real-valued complex
    y2 = std::sqrt(0.5 * (S + B));
    x2 = std::sqrt(0.5 * (S - B));

    /*
        We have define: 
            x2 and y2 are the squares of x and y
            x and y are the real and positive Eigenvalues
            x is from the real part of the Eigenvalues
            y is from the imag part of the Eigenvalues
    */
    x = std::sqrt(x2);
    y = std::sqrt(y2);

    cy = std::cos(y);
    sy = std::sin(y);
    cx = std::cosh(x);
    sx = std::sinh(x);

    x_zero = x < lower_is_considered_zero_for_sinc_likes;
    y_zero = y < lower_is_considered_zero_for_sinc_likes;
    both_zero = y_zero and x_zero;
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

    C0 = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
    C1 = either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
    C2 = both_zero ? 0.5 : (cx - cy) * inv_x2y2;
    C3 = both_zero ? 1.0 / 6.0
                   : (x_zero   ? 1.0 - sy * iy
                      : y_zero ? sx * ix - 1.0
                               : sx * ix - sy * iy) *
                         inv_x2y2;
  }

  constexpr muelmat operator()() const noexcept {
    return unpolarized ? muelmat{exp_a}
                       : muelmat{exp_a * (C0 + C2 * b2 + C2 * c2 + C2 * d2),
                                 -exp_a * (-C1 * b + C2 * c * u + C2 * d * v -
                                           C3 * b * (b2 + c2 + d2) +
                                           C3 * u * (b * u - d * w) +
                                           C3 * v * (b * v + c * w)),
                                 exp_a * (C1 * c + C2 * b * u - C2 * d * w +
                                          C3 * c * (b2 + c2 + d2) -
                                          C3 * u * (c * u + d * v) -
                                          C3 * w * (b * v + c * w)),
                                 exp_a * (C1 * d + C2 * b * v + C2 * c * w +
                                          C3 * d * (b2 + c2 + d2) -
                                          C3 * v * (c * u + d * v) +
                                          C3 * w * (b * u - d * w)),
                                 exp_a * (C1 * b + C2 * c * u + C2 * d * v -
                                          C3 * b * (-b2 + u2 + v2) +
                                          C3 * c * (b * c - v * w) +
                                          C3 * d * (b * d + u * w)),
                                 exp_a * (C0 + C2 * b2 - C2 * u2 - C2 * v2),
                                 exp_a * (C1 * u + C2 * b * c - C2 * v * w +
                                          C3 * c * (c * u + d * v) -
                                          C3 * u * (-b2 + u2 + v2) -
                                          C3 * w * (b * d + u * w)),
                                 exp_a * (C1 * v + C2 * b * d + C2 * u * w +
                                          C3 * d * (c * u + d * v) -
                                          C3 * v * (-b2 + u2 + v2) +
                                          C3 * w * (b * c - v * w)),
                                 exp_a * (C1 * c - C2 * b * u + C2 * d * w +
                                          C3 * b * (b * c - v * w) -
                                          C3 * c * (-c2 + u2 + w2) +
                                          C3 * d * (c * d - u * v)),
                                 -exp_a * (C1 * u - C2 * b * c + C2 * v * w +
                                           C3 * b * (b * u - d * w) -
                                           C3 * u * (-c2 + u2 + w2) +
                                           C3 * v * (c * d - u * v)),
                                 exp_a * (C0 + C2 * c2 - C2 * u2 - C2 * w2),
                                 exp_a * (C1 * w + C2 * c * d - C2 * u * v -
                                          C3 * d * (b * u - d * w) +
                                          C3 * v * (b * c - v * w) -
                                          C3 * w * (-c2 + u2 + w2)),
                                 exp_a * (C1 * d - C2 * b * v - C2 * c * w +
                                          C3 * b * (b * d + u * w) +
                                          C3 * c * (c * d - u * v) -
                                          C3 * d * (-d2 + v2 + w2)),
                                 -exp_a * (C1 * v - C2 * b * d - C2 * u * w +
                                           C3 * b * (b * v + c * w) +
                                           C3 * u * (c * d - u * v) -
                                           C3 * v * (-d2 + v2 + w2)),
                                 -exp_a * (C1 * w - C2 * c * d + C2 * u * v +
                                           C3 * c * (b * v + c * w) -
                                           C3 * u * (b * d + u * w) -
                                           C3 * w * (-d2 + v2 + w2)),
                                 exp_a * (C0 + C2 * d2 - C2 * v2 - C2 * w2)};
  }

  [[nodiscard]] muelmat deriv(const muelmat &t,
                              const propmat &k1,
                              const propmat &k2,
                              const propmat &dk,
                              const Numeric r,
                              const Numeric dr) const {
    const Numeric da = -0.5 * (r * dk.A() + dr * (k1.A() + k2.A()));
    if (unpolarized) {
      return {da * exp_a};
    }

    const Numeric db = -0.5 * (r * dk.B() + dr * (k1.B() + k2.B())),
                  dc = -0.5 * (r * dk.C() + dr * (k1.C() + k2.C())),
                  dd = -0.5 * (r * dk.D() + dr * (k1.D() + k2.D())),
                  du = -0.5 * (r * dk.U() + dr * (k1.U() + k2.U())),
                  dv = -0.5 * (r * dk.V() + dr * (k1.V() + k2.V())),
                  dw = -0.5 * (r * dk.W() + dr * (k1.W() + k2.W()));
    const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                  du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;

    /* Solve: 
        0 = L^4 + B L^2 + C
        B = U^2+V^2+W^2-B^2-C^2-D^2
        C = - (DU - CV + BW)^2
    */
    const Numeric dB = du2 + dv2 + dw2 - db2 - dc2 - dd2;
    const Numeric dC = -2 * (b * w - c * v + d * u) *
                       (b * dw - c * dv + d * du + db * w - dc * v + dd * u);
    const Numeric dS = (2 * B * dB - 4 * dC) * 0.5 / S;

    const Numeric dy2 = 0.25 * (dS + dB) / y2;
    const Numeric dx2 = 0.25 * (dS - dB) / x2;
    const Numeric dx = 0.5 * dx2 / x;
    const Numeric dy = 0.5 * dy2 / y;

    const Numeric dcy = -sy * dy;
    const Numeric dsy = cy * dy;
    const Numeric dcx = sx * dx;
    const Numeric dsx = cx * dx;
    const Numeric dix = -dx * ix * ix;
    const Numeric diy = -dy * iy * iy;
    const Numeric dx2dy2 = dx2 + dy2;
    const Numeric dC0 =
        either_zero
            ? 0.0
            : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0 * dx2dy2) *
                  inv_x2y2;
    const Numeric dC1 =
        either_zero
            ? 0.0
            : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy + dsx * y2 * ix +
               sx * dy2 * ix + sx * y2 * dix - C1 * dx2dy2) *
                  inv_x2y2;
    const Numeric dC2 = both_zero ? 0.0 : (dcx - dcy - C2 * dx2dy2) * inv_x2y2;
    const Numeric dC3 =
        both_zero ? 0.0
                  : ((x_zero   ? -dsy * iy - sy * diy
                      : y_zero ? dsx * ix + sx * dix
                               : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                     C3 * dx2dy2) *
                        inv_x2y2;

    return {t(0, 0) * da + exp_a * dC0 + dC2 * (b2 + c2 + d2) +
                C2 * (db2 + dc2 + dd2),
            t(0, 1) * da + exp_a * db * C1 + b * dC1 + dC2 * (-c * u - d * v) +
                C2 * (-dc * u - dd * v - c * du - d * dv) +
                dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                       v * (b * v + c * w)) +
                C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                      dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                      u * (db * u - dd * w) - v * (db * v + dc * w) -
                      u * (b * du - d * dw) - v * (b * dv + c * dw)),
            t(0, 2) * da + exp_a * dC1 * c + C1 * dc + dC2 * (b * u - d * w) +
                C2 * (db * u - dd * w + b * du - d * dw) +
                dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                       w * (b * v + c * w)) +
                C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                      dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                      u * (dc * u + dd * v) - w * (db * v + dc * w) -
                      u * (c * du + d * dv) - w * (b * dv + c * dw)),
            t(0, 3) * da + exp_a * dC1 * d + C1 * dd + dC2 * (b * v + c * w) +
                C2 * (db * v + dc * w + b * dv + c * dw) +
                dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                       w * (b * u - d * w)) +
                C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                      dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                      v * (dc * u + dd * v) + w * (db * u - dd * w) -
                      v * (c * du + d * dv) + w * (b * du - d * dw)),

            t(1, 0) * da + exp_a * db * C1 + b * dC1 + dC2 * (c * u + d * v) +
                C2 * (dc * u + dd * v + c * du + d * dv) +
                dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                       d * (b * d + u * w)) +
                C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                      dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                      c * (db * c - dv * w) + d * (db * d + du * w) +
                      c * (b * dc - v * dw) + d * (b * dd + u * dw)),
            t(1, 1) * da + exp_a * dC0 + dC2 * (b2 - u2 - v2) +
                C2 * (db2 - du2 - dv2),
            t(1, 2) * da + exp_a * dC2 * (b * c - v * w) +
                C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u + C1 * du +
                dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                       w * (b * d + u * w)) +
                C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                      dw * (b * d + u * w) + c * (dc * u + dd * v) -
                      u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                      c * (c * du + d * dv) - w * (b * dd + u * dw)),
            t(1, 3) * da + exp_a * dC2 * (b * d + u * w) +
                C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v + C1 * dv +
                dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                       w * (b * c - v * w)) +
                C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                      dw * (b * c - v * w) + d * (dc * u + dd * v) -
                      v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                      d * (c * du + d * dv) + w * (b * dc - v * dw)),

            t(2, 0) * da + exp_a * dC1 * c + C1 * dc + dC2 * (-b * u + d * w) +
                C2 * (-db * u + dd * w - b * du + d * dw) +
                dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                       d * (c * d - u * v)) +
                C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                      dd * (c * d - u * v) + b * (db * c - dv * w) -
                      c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                      b * (b * dc - v * dw) + d * (c * dd - u * dv)),
            t(2, 1) * da + exp_a * dC2 * (b * c - v * w) +
                C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u - C1 * du +
                dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                       v * (c * d - u * v)) +
                C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                      dv * (c * d - u * v) - b * (db * u - dd * w) +
                      u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                      b * (b * du - d * dw) - v * (c * dd - u * dv)),
            t(2, 2) * da + exp_a * dC0 + dC2 * (c2 - u2 - w2) +
                C2 * (dc2 - du2 - dw2),
            t(2, 3) * da + exp_a * dC2 * (c * d - u * v) +
                C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w + C1 * dw +
                dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                       w * (-c2 + u2 + w2)) +
                C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                      dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                      v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                      d * (b * du - d * dw) + v * (b * dc - v * dw)),

            t(3, 0) * da + exp_a * dC1 * d + C1 * dd + dC2 * (-b * v - c * w) +
                C2 * (-db * v - dc * w - b * dv - c * dw) +
                dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                       d * (-d2 + v2 + w2)) +
                C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                      dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                      c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                      b * (b * dd + u * dw) + c * (c * dd - u * dv)),
            t(3, 1) * da + exp_a * dC2 * (b * d + u * w) +
                C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v - C1 * dv +
                dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                       v * (-d2 + v2 + w2)) +
                C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                      dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                      u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                      b * (b * dv + c * dw) - u * (c * dd - u * dv)),
            t(3, 2) * da + exp_a * dC2 * (c * d - u * v) +
                C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w - C1 * dw +
                dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                       w * (-d2 + v2 + w2)) +
                C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                      dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                      u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                      c * (b * dv + c * dw) + u * (b * dd + u * dw)),
            t(3, 3) * da + exp_a * dC0 + dC2 * (d2 - v2 - w2) +
                C2 * (dd2 - dv2 - dw2)};
  }
};

void two_level_exp_test(muelmat &t,
                        muelmat_vector_view &dt1,
                        muelmat_vector_view &dt2,
                        const propmat &k1,
                        const propmat &k2,
                        const propmat_vector_const_view &dk1,
                        const propmat_vector_const_view &dk2,
                        const Numeric r,
                        const ExhaustiveConstVectorView &dr1,
                        const ExhaustiveConstVectorView &dr2) {
  ARTS_ASSERT(dk1.size() == dk2.size() and dk1.size() == dr1.size() and
              dk1.size() == dr2.size())

  const tran tran_state{k1, k2, r};
  t = tran_state();

  const auto deriv = [&](const propmat &dk, const Numeric &dr) -> muelmat {
    return tran_state.deriv(t, k1, k2, dk, r, dr);
  };

  std::transform(dk1.begin(), dk1.end(), dr1.begin(), dt1.begin(), deriv);
  std::transform(dk2.begin(), dk2.end(), dr2.begin(), dt2.begin(), deriv);
}

void two_level_exp(muelmat &t,
                   muelmat_vector_view &dt1,
                   muelmat_vector_view &dt2,
                   const propmat &k1,
                   const propmat &k2,
                   const propmat_vector_const_view &dk1,
                   const propmat_vector_const_view &dk2,
                   const Numeric r,
                   const ExhaustiveConstVectorView &dr1,
                   const ExhaustiveConstVectorView &dr2) {
  const Index N = dk1.size();
  ARTS_ASSERT(N == dk2.size() and N == dr1.size() and N == dr2.size())

  const Numeric a = -0.5 * r * (k1.A() + k2.A()),
                b = -0.5 * r * (k1.B() + k2.B()),
                c = -0.5 * r * (k1.C() + k2.C()),
                d = -0.5 * r * (k1.D() + k2.D()),
                u = -0.5 * r * (k1.U() + k2.U()),
                v = -0.5 * r * (k1.V() + k2.V()),
                w = -0.5 * r * (k1.W() + k2.W());
  const Numeric exp_a = std::exp(a);

  if (b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.) {
    t = muelmat{exp_a};

    for (Index j = 0; j < N; j++) {
      dt1[j] =
          muelmat{-0.5 * (r * dk1[j].A() + dr1[j] * (k1.A() + k2.A())) * exp_a};
      dt2[j] =
          muelmat{-0.5 * (r * dk2[j].A() + dr2[j] * (k1.A() + k2.A())) * exp_a};
    }
  } else {
    const Numeric b2 = b * b, c2 = c * c, d2 = d * d, u2 = u * u, v2 = v * v,
                  w2 = w * w;
    const Numeric tmp =
        w2 * w2 + 2 * (b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                       c2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                       d2 * (d2 * 0.5 + u2 - v2 - w2) +
                       u2 * (u2 * 0.5 + v2 + w2) + v2 * (v2 * 0.5 + w2) +
                       4 * (b * d * u * w - b * c * v * w - c * d * u * v));
    const Complex Const1 = std::sqrt(Complex(tmp, 0));
    const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
    const Complex tmp_x_sqrt = std::sqrt(Const2 + Const1);
    const Complex tmp_y_sqrt = std::sqrt(Const2 - Const1);
    const Complex x = tmp_x_sqrt * sqrt_05;
    const Complex y = tmp_y_sqrt * sqrt_05 * Complex(0, 1);
    const Complex x2 = x * x;
    const Complex y2 = y * y;
    const Complex cy = std::cos(y);
    const Complex sy = std::sin(y);
    const Complex cx = std::cosh(x);
    const Complex sx = std::sinh(x);

    const bool x_zero = std::abs(x) < lower_is_considered_zero_for_sinc_likes;
    const bool y_zero = std::abs(y) < lower_is_considered_zero_for_sinc_likes;
    const bool both_zero = y_zero and x_zero;
    const bool either_zero = y_zero or x_zero;

    /* Using:
     *    lim x→0 [({cosh(x),cos(x)} - 1) / x^2] → 1/2
     *    lim x→0 [{sinh(x),sin(x)} / x]  → 1
     *    inv_x2 := 1 for x == 0,
     *    -i sin(ix) → sinh(x)
     *    cos(ix) → cosh(x)
     *    C0, C1, C2 ∝ [1/x^2]
     */
    const Complex ix = x_zero ? 0.0 : 1.0 / x;
    const Complex iy = y_zero ? 0.0 : 1.0 / y;
    const Complex inv_x2y2 =
        both_zero
            ? 1.0
            : 1.0 / (x2 + y2);  // The first "1.0" is the trick for above limits
    const Complex C0c = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
    const Complex C1c =
        either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
    const Complex C2c = both_zero ? 0.5 : (cx - cy) * inv_x2y2;
    const Complex C3c = both_zero ? 1.0 / 6.0
                                  : (x_zero   ? 1.0 - sy * iy
                                     : y_zero ? sx * ix - 1.0
                                              : sx * ix - sy * iy) *
                                        inv_x2y2;

    const Numeric &C0 = real_val(C0c);
    const Numeric &C1 = real_val(C1c);
    const Numeric &C2 = real_val(C2c);
    const Numeric &C3 = real_val(C3c);
    t = muelmat{exp_a * C0 + C2 * (b2 + c2 + d2),
                exp_a * C1 * b + C2 * (-c * u - d * v) +
                    C3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                          v * (b * v + c * w)),
                exp_a * C1 * c + C2 * (b * u - d * w) +
                    C3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                          w * (b * v + c * w)),
                exp_a * C1 * d + C2 * (b * v + c * w) +
                    C3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                          w * (b * u - d * w)),

                exp_a * C1 * b + C2 * (c * u + d * v) +
                    C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                          d * (b * d + u * w)),
                exp_a * C0 + C2 * (b2 - u2 - v2),
                exp_a * C2 * (b * c - v * w) + C1 * u +
                    C3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                          w * (b * d + u * w)),
                exp_a * C2 * (b * d + u * w) + C1 * v +
                    C3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                          w * (b * c - v * w)),

                exp_a * C1 * c + C2 * (-b * u + d * w) +
                    C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                          d * (c * d - u * v)),
                exp_a * C2 * (b * c - v * w) - C1 * u +
                    C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                          v * (c * d - u * v)),
                exp_a * C0 + C2 * (c2 - u2 - w2),
                exp_a * C2 * (c * d - u * v) + C1 * w +
                    C3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                          w * (-c2 + u2 + w2)),

                exp_a * C1 * d + C2 * (-b * v - c * w) +
                    C3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                          d * (-d2 + v2 + w2)),
                exp_a * C2 * (b * d + u * w) - C1 * v +
                    C3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                          v * (-d2 + v2 + w2)),
                exp_a * C2 * (c * d - u * v) - C1 * w +
                    C3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                          w * (-d2 + v2 + w2)),
                exp_a * C0 + C2 * (d2 - v2 - w2)};

    for (Index j = 0; j < dk1.nelem(); j++) {
      {
        const Numeric da = -0.5 * (r * dk1[j].A() + dr1[j] * (k1.A() + k2.A())),
                      db = -0.5 * (r * dk1[j].B() + dr1[j] * (k1.B() + k2.B())),
                      dc = -0.5 * (r * dk1[j].C() + dr1[j] * (k1.C() + k2.C())),
                      dd = -0.5 * (r * dk1[j].D() + dr1[j] * (k1.D() + k2.D())),
                      du = -0.5 * (r * dk1[j].U() + dr1[j] * (k1.U() + k2.U())),
                      dv = -0.5 * (r * dk1[j].V() + dr1[j] * (k1.V() + k2.V())),
                      dw = -0.5 * (r * dk1[j].W() + dr1[j] * (k1.W() + k2.W()));
        const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                      du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;
        const Numeric dtmp =
            2 * w2 * dw2 +
            2 * (db2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                 b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2) +
                 dc2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                 c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2) +
                 dd2 * (d2 * 0.5 + u2 - v2 - w2) +
                 d2 * (dd2 * 0.5 + du2 - dv2 - dw2) +
                 du2 * (u2 * 0.5 + v2 + w2) + u2 * (du2 * 0.5 + dv2 + dw2) +
                 dv2 * (v2 * 0.5 + w2) + v2 * (dv2 * 0.5 + dw2) +
                 4 * (db * d * u * w - db * c * v * w - dc * d * u * v +
                      b * dd * u * w - b * dc * v * w - c * dd * u * v +
                      b * d * du * w - b * c * dv * w - c * d * du * v +
                      b * d * u * dw - b * c * v * dw - c * d * u * dv));
        const Complex dConst1 = 0.5 * dtmp / Const1;
        const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
        const Complex dx =
            x_zero ? 0 : (0.5 * (dConst2 + dConst1) / tmp_x_sqrt) * sqrt_05;
        const Complex dy = y_zero ? 0
                                  : (0.5 * (dConst2 - dConst1) / tmp_y_sqrt) *
                                        sqrt_05 * Complex(0, 1);
        const Complex dx2 = 2.0 * x * dx;
        const Complex dy2 = 2.0 * y * dy;
        const Complex dcy = -sy * dy;
        const Complex dsy = cy * dy;
        const Complex dcx = sx * dx;
        const Complex dsx = cx * dx;
        const Complex dix = -dx * ix * ix;
        const Complex diy = -dy * iy * iy;
        const Complex dx2dy2 = dx2 + dy2;
        const Complex dC0c =
            either_zero
                ? 0.0
                : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0c * dx2dy2) *
                      inv_x2y2;
        const Complex dC1c =
            either_zero ? 0.0
                        : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy +
                           dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix -
                           C1c * dx2dy2) *
                              inv_x2y2;
        const Complex dC2c =
            both_zero ? 0.0 : (dcx - dcy - C2c * dx2dy2) * inv_x2y2;
        const Complex dC3c =
            both_zero
                ? 0.0
                : ((x_zero   ? -dsy * iy - sy * diy
                    : y_zero ? dsx * ix + sx * dix
                             : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                   C3c * dx2dy2) *
                      inv_x2y2;

        const Numeric &dC0 = real_val(dC0c);
        const Numeric &dC1 = real_val(dC1c);
        const Numeric &dC2 = real_val(dC2c);
        const Numeric &dC3 = real_val(dC3c);
        dt1[j] = muelmat{
            t(0, 0) * da + exp_a * dC0 + dC2 * (b2 + c2 + d2) +
                C2 * (db2 + dc2 + dd2),
            t(0, 1) * da + exp_a * db * C1 + b * dC1 + dC2 * (-c * u - d * v) +
                C2 * (-dc * u - dd * v - c * du - d * dv) +
                dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                       v * (b * v + c * w)) +
                C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                      dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                      u * (db * u - dd * w) - v * (db * v + dc * w) -
                      u * (b * du - d * dw) - v * (b * dv + c * dw)),
            t(0, 2) * da + exp_a * dC1 * c + C1 * dc + dC2 * (b * u - d * w) +
                C2 * (db * u - dd * w + b * du - d * dw) +
                dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                       w * (b * v + c * w)) +
                C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                      dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                      u * (dc * u + dd * v) - w * (db * v + dc * w) -
                      u * (c * du + d * dv) - w * (b * dv + c * dw)),
            t(0, 3) * da + exp_a * dC1 * d + C1 * dd + dC2 * (b * v + c * w) +
                C2 * (db * v + dc * w + b * dv + c * dw) +
                dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                       w * (b * u - d * w)) +
                C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                      dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                      v * (dc * u + dd * v) + w * (db * u - dd * w) -
                      v * (c * du + d * dv) + w * (b * du - d * dw)),

            t(1, 0) * da + exp_a * db * C1 + b * dC1 + dC2 * (c * u + d * v) +
                C2 * (dc * u + dd * v + c * du + d * dv) +
                dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                       d * (b * d + u * w)) +
                C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                      dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                      c * (db * c - dv * w) + d * (db * d + du * w) +
                      c * (b * dc - v * dw) + d * (b * dd + u * dw)),
            t(1, 1) * da + exp_a * dC0 + dC2 * (b2 - u2 - v2) +
                C2 * (db2 - du2 - dv2),
            t(1, 2) * da + exp_a * dC2 * (b * c - v * w) +
                C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u + C1 * du +
                dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                       w * (b * d + u * w)) +
                C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                      dw * (b * d + u * w) + c * (dc * u + dd * v) -
                      u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                      c * (c * du + d * dv) - w * (b * dd + u * dw)),
            t(1, 3) * da + exp_a * dC2 * (b * d + u * w) +
                C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v + C1 * dv +
                dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                       w * (b * c - v * w)) +
                C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                      dw * (b * c - v * w) + d * (dc * u + dd * v) -
                      v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                      d * (c * du + d * dv) + w * (b * dc - v * dw)),

            t(2, 0) * da + exp_a * dC1 * c + C1 * dc + dC2 * (-b * u + d * w) +
                C2 * (-db * u + dd * w - b * du + d * dw) +
                dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                       d * (c * d - u * v)) +
                C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                      dd * (c * d - u * v) + b * (db * c - dv * w) -
                      c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                      b * (b * dc - v * dw) + d * (c * dd - u * dv)),
            t(2, 1) * da + exp_a * dC2 * (b * c - v * w) +
                C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u - C1 * du +
                dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                       v * (c * d - u * v)) +
                C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                      dv * (c * d - u * v) - b * (db * u - dd * w) +
                      u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                      b * (b * du - d * dw) - v * (c * dd - u * dv)),
            t(2, 2) * da + exp_a * dC0 + dC2 * (c2 - u2 - w2) +
                C2 * (dc2 - du2 - dw2),
            t(2, 3) * da + exp_a * dC2 * (c * d - u * v) +
                C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w + C1 * dw +
                dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                       w * (-c2 + u2 + w2)) +
                C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                      dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                      v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                      d * (b * du - d * dw) + v * (b * dc - v * dw)),

            t(3, 0) * da + exp_a * dC1 * d + C1 * dd + dC2 * (-b * v - c * w) +
                C2 * (-db * v - dc * w - b * dv - c * dw) +
                dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                       d * (-d2 + v2 + w2)) +
                C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                      dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                      c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                      b * (b * dd + u * dw) + c * (c * dd - u * dv)),
            t(3, 1) * da + exp_a * dC2 * (b * d + u * w) +
                C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v - C1 * dv +
                dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                       v * (-d2 + v2 + w2)) +
                C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                      dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                      u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                      b * (b * dv + c * dw) - u * (c * dd - u * dv)),
            t(3, 2) * da + exp_a * dC2 * (c * d - u * v) +
                C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w - C1 * dw +
                dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                       w * (-d2 + v2 + w2)) +
                C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                      dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                      u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                      c * (b * dv + c * dw) + u * (b * dd + u * dw)),
            t(3, 3) * da + exp_a * dC0 + dC2 * (d2 - v2 - w2) +
                C2 * (dd2 - dv2 - dw2)};
      }

      {
        const Numeric da = -0.5 * (r * dk2[j].A() + dr2[j] * (k1.A() + k2.A())),
                      db = -0.5 * (r * dk2[j].B() + dr2[j] * (k1.B() + k2.B())),
                      dc = -0.5 * (r * dk2[j].C() + dr2[j] * (k1.C() + k2.C())),
                      dd = -0.5 * (r * dk2[j].D() + dr2[j] * (k1.D() + k2.D())),
                      du = -0.5 * (r * dk2[j].U() + dr2[j] * (k1.U() + k2.U())),
                      dv = -0.5 * (r * dk2[j].V() + dr2[j] * (k1.V() + k2.V())),
                      dw = -0.5 * (r * dk2[j].W() + dr2[j] * (k1.W() + k2.W()));
        const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                      du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;
        const Numeric dtmp =
            2 * w2 * dw2 +
            2 * (db2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                 b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2) +
                 dc2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                 c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2) +
                 dd2 * (d2 * 0.5 + u2 - v2 - w2) +
                 d2 * (dd2 * 0.5 + du2 - dv2 - dw2) +
                 du2 * (u2 * 0.5 + v2 + w2) + u2 * (du2 * 0.5 + dv2 + dw2) +
                 dv2 * (v2 * 0.5 + w2) + v2 * (dv2 * 0.5 + dw2) +
                 4 * (db * d * u * w - db * c * v * w - dc * d * u * v +
                      b * dd * u * w - b * dc * v * w - c * dd * u * v +
                      b * d * du * w - b * c * dv * w - c * d * du * v +
                      b * d * u * dw - b * c * v * dw - c * d * u * dv));
        const Complex dConst1 = 0.5 * dtmp / Const1;
        const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
        const Complex dx =
            x_zero ? 0 : (0.5 * (dConst2 + dConst1) / tmp_x_sqrt) * sqrt_05;
        const Complex dy = y_zero ? 0
                                  : (0.5 * (dConst2 - dConst1) / tmp_y_sqrt) *
                                        sqrt_05 * Complex(0, 1);
        const Complex dx2 = 2.0 * x * dx;
        const Complex dy2 = 2.0 * y * dy;
        const Complex dcy = -sy * dy;
        const Complex dsy = cy * dy;
        const Complex dcx = sx * dx;
        const Complex dsx = cx * dx;
        const Complex dix = -dx * ix * ix;
        const Complex diy = -dy * iy * iy;
        const Complex dx2dy2 = dx2 + dy2;
        const Complex dC0c =
            either_zero
                ? 0.0
                : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0c * dx2dy2) *
                      inv_x2y2;
        const Complex dC1c =
            either_zero ? 0.0
                        : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy +
                           dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix -
                           C1c * dx2dy2) *
                              inv_x2y2;
        const Complex dC2c =
            both_zero ? 0.0 : (dcx - dcy - C2c * dx2dy2) * inv_x2y2;
        const Complex dC3c =
            both_zero
                ? 0.0
                : ((x_zero   ? -dsy * iy - sy * diy
                    : y_zero ? dsx * ix + sx * dix
                             : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                   C3c * dx2dy2) *
                      inv_x2y2;

        const Numeric &dC0 = real_val(dC0c);
        const Numeric &dC1 = real_val(dC1c);
        const Numeric &dC2 = real_val(dC2c);
        const Numeric &dC3 = real_val(dC3c);
        dt2[j] = muelmat{
            t(0, 0) * da + exp_a * dC0 + dC2 * (b2 + c2 + d2) +
                C2 * (db2 + dc2 + dd2),
            t(0, 1) * da + exp_a * db * C1 + b * dC1 + dC2 * (-c * u - d * v) +
                C2 * (-dc * u - dd * v - c * du - d * dv) +
                dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                       v * (b * v + c * w)) +
                C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                      dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                      u * (db * u - dd * w) - v * (db * v + dc * w) -
                      u * (b * du - d * dw) - v * (b * dv + c * dw)),
            t(0, 2) * da + exp_a * dC1 * c + C1 * dc + dC2 * (b * u - d * w) +
                C2 * (db * u - dd * w + b * du - d * dw) +
                dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                       w * (b * v + c * w)) +
                C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                      dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                      u * (dc * u + dd * v) - w * (db * v + dc * w) -
                      u * (c * du + d * dv) - w * (b * dv + c * dw)),
            t(0, 3) * da + exp_a * dC1 * d + C1 * dd + dC2 * (b * v + c * w) +
                C2 * (db * v + dc * w + b * dv + c * dw) +
                dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                       w * (b * u - d * w)) +
                C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                      dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                      v * (dc * u + dd * v) + w * (db * u - dd * w) -
                      v * (c * du + d * dv) + w * (b * du - d * dw)),

            t(1, 0) * da + exp_a * db * C1 + b * dC1 + dC2 * (c * u + d * v) +
                C2 * (dc * u + dd * v + c * du + d * dv) +
                dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                       d * (b * d + u * w)) +
                C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                      dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                      c * (db * c - dv * w) + d * (db * d + du * w) +
                      c * (b * dc - v * dw) + d * (b * dd + u * dw)),
            t(1, 1) * da + exp_a * dC0 + dC2 * (b2 - u2 - v2) +
                C2 * (db2 - du2 - dv2),
            t(1, 2) * da + exp_a * dC2 * (b * c - v * w) +
                C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u + C1 * du +
                dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                       w * (b * d + u * w)) +
                C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                      dw * (b * d + u * w) + c * (dc * u + dd * v) -
                      u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                      c * (c * du + d * dv) - w * (b * dd + u * dw)),
            t(1, 3) * da + exp_a * dC2 * (b * d + u * w) +
                C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v + C1 * dv +
                dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                       w * (b * c - v * w)) +
                C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                      dw * (b * c - v * w) + d * (dc * u + dd * v) -
                      v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                      d * (c * du + d * dv) + w * (b * dc - v * dw)),

            t(2, 0) * da + exp_a * dC1 * c + C1 * dc + dC2 * (-b * u + d * w) +
                C2 * (-db * u + dd * w - b * du + d * dw) +
                dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                       d * (c * d - u * v)) +
                C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                      dd * (c * d - u * v) + b * (db * c - dv * w) -
                      c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                      b * (b * dc - v * dw) + d * (c * dd - u * dv)),
            t(2, 1) * da + exp_a * dC2 * (b * c - v * w) +
                C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u - C1 * du +
                dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                       v * (c * d - u * v)) +
                C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                      dv * (c * d - u * v) - b * (db * u - dd * w) +
                      u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                      b * (b * du - d * dw) - v * (c * dd - u * dv)),
            t(2, 2) * da + exp_a * dC0 + dC2 * (c2 - u2 - w2) +
                C2 * (dc2 - du2 - dw2),
            t(2, 3) * da + exp_a * dC2 * (c * d - u * v) +
                C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w + C1 * dw +
                dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                       w * (-c2 + u2 + w2)) +
                C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                      dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                      v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                      d * (b * du - d * dw) + v * (b * dc - v * dw)),

            t(3, 0) * da + exp_a * dC1 * d + C1 * dd + dC2 * (-b * v - c * w) +
                C2 * (-db * v - dc * w - b * dv - c * dw) +
                dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                       d * (-d2 + v2 + w2)) +
                C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                      dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                      c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                      b * (b * dd + u * dw) + c * (c * dd - u * dv)),
            t(3, 1) * da + exp_a * dC2 * (b * d + u * w) +
                C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v - C1 * dv +
                dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                       v * (-d2 + v2 + w2)) +
                C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                      dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                      u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                      b * (b * dv + c * dw) - u * (c * dd - u * dv)),
            t(3, 2) * da + exp_a * dC2 * (c * d - u * v) +
                C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w - C1 * dw +
                dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                       w * (-d2 + v2 + w2)) +
                C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                      dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                      u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                      c * (b * dv + c * dw) + u * (b * dd + u * dw)),
            t(3, 3) * da + exp_a * dC0 + dC2 * (d2 - v2 - w2) +
                C2 * (dd2 - dv2 - dw2)};
      }
    }
  }
}

void two_level_exp(muelmat_vector_view &tv,
                   muelmat_matrix_view &dt1v,
                   muelmat_matrix_view &dt2v,
                   const propmat_vector_const_view &k1v,
                   const propmat_vector_const_view &k2v,
                   const propmat_matrix_const_view &dk1v,
                   const propmat_matrix_const_view &dk2v,
                   const Numeric rv,
                   const ExhaustiveConstVectorView &dr1v,
                   const ExhaustiveConstVectorView &dr2v) {
  const Index nf = tv.nelem();
  const Index nq = dr1v.nelem();

  ARTS_ASSERT(nf == k1v.nelem());
  ARTS_ASSERT(nf == k2v.nelem());
  ARTS_ASSERT(nf == dk1v.ncols());
  ARTS_ASSERT(nf == dk2v.ncols());
  ARTS_ASSERT(nq == dk1v.nrows());
  ARTS_ASSERT(nq == dk2v.nrows());
  ARTS_ASSERT(nf == dt1v.ncols());
  ARTS_ASSERT(nf == dt2v.ncols());
  ARTS_ASSERT(nq == dt1v.nrows());
  ARTS_ASSERT(nq == dt2v.nrows());
  ARTS_ASSERT(nq == dr2v.nelem());

  for (Index i = 0; i < nf; ++i) {
    const Numeric a = -0.5 * rv * (k1v[i].A() + k2v[i].A()),
                  b = -0.5 * rv * (k1v[i].B() + k2v[i].B()),
                  c = -0.5 * rv * (k1v[i].C() + k2v[i].C()),
                  d = -0.5 * rv * (k1v[i].D() + k2v[i].D()),
                  u = -0.5 * rv * (k1v[i].U() + k2v[i].U()),
                  v = -0.5 * rv * (k1v[i].V() + k2v[i].V()),
                  w = -0.5 * rv * (k1v[i].W() + k2v[i].W());
    const Numeric exp_a = std::exp(a);

    if (b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.) {
      tv[i] = muelmat{exp_a};

      for (Index j = 0; j < nq; j++) {
        dt1v(j, i) = muelmat{
            -0.5 * (rv * dk1v(j, i).A() + dr1v[j] * (k1v[i].A() + k2v[i].A())) *
            exp_a};
        dt2v(j, i) = muelmat{
            -0.5 * (rv * dk2v(j, i).A() + dr2v[j] * (k1v[i].A() + k2v[i].A())) *
            exp_a};
      }
    } else {
      const Numeric b2 = b * b, c2 = c * c, d2 = d * d, u2 = u * u, v2 = v * v,
                    w2 = w * w;
      const Numeric tmp =
          w2 * w2 + 2 * (b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                         c2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                         d2 * (d2 * 0.5 + u2 - v2 - w2) +
                         u2 * (u2 * 0.5 + v2 + w2) + v2 * (v2 * 0.5 + w2) +
                         4 * (b * d * u * w - b * c * v * w - c * d * u * v));
      const Complex Const1 = std::sqrt(Complex(tmp, 0));
      const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
      const Complex tmp_x_sqrt = std::sqrt(Const2 + Const1);
      const Complex tmp_y_sqrt = std::sqrt(Const2 - Const1);
      const Complex x = tmp_x_sqrt * sqrt_05;
      const Complex y = tmp_y_sqrt * sqrt_05 * Complex(0, 1);
      const Complex x2 = x * x;
      const Complex y2 = y * y;
      const Complex cy = std::cos(y);
      const Complex sy = std::sin(y);
      const Complex cx = std::cosh(x);
      const Complex sx = std::sinh(x);

      const bool x_zero = std::abs(x) < lower_is_considered_zero_for_sinc_likes;
      const bool y_zero = std::abs(y) < lower_is_considered_zero_for_sinc_likes;
      const bool both_zero = y_zero and x_zero;
      const bool either_zero = y_zero or x_zero;

      /* Using:
       *    lim x→0 [({cosh(x),cos(x)} - 1) / x^2] → 1/2
       *    lim x→0 [{sinh(x),sin(x)} / x]  → 1
       *    inv_x2 := 1 for x == 0,
       *    -i sin(ix) → sinh(x)
       *    cos(ix) → cosh(x)
       *    C0, C1, C2 ∝ [1/x^2]
       */
      const Complex ix = x_zero ? 0.0 : 1.0 / x;
      const Complex iy = y_zero ? 0.0 : 1.0 / y;
      const Complex inv_x2y2 =
          both_zero
              ? 1.0
              : 1.0 /
                    (x2 + y2);  // The first "1.0" is the trick for above limits
      const Complex C0c = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
      const Complex C1c =
          either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
      const Complex C2c = both_zero ? 0.5 : (cx - cy) * inv_x2y2;
      const Complex C3c = both_zero ? 1.0 / 6.0
                                    : (x_zero   ? 1.0 - sy * iy
                                       : y_zero ? sx * ix - 1.0
                                                : sx * ix - sy * iy) *
                                          inv_x2y2;

      const Numeric &C0 = real_val(C0c);
      const Numeric &C1 = real_val(C1c);
      const Numeric &C2 = real_val(C2c);
      const Numeric &C3 = real_val(C3c);
      tv[i] = muelmat{exp_a * C0 + C2 * (b2 + c2 + d2),
                      exp_a * C1 * b + C2 * (-c * u - d * v) +
                          C3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                                v * (b * v + c * w)),
                      exp_a * C1 * c + C2 * (b * u - d * w) +
                          C3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                                w * (b * v + c * w)),
                      exp_a * C1 * d + C2 * (b * v + c * w) +
                          C3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                                w * (b * u - d * w)),

                      exp_a * C1 * b + C2 * (c * u + d * v) +
                          C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                                d * (b * d + u * w)),
                      exp_a * C0 + C2 * (b2 - u2 - v2),
                      exp_a * C2 * (b * c - v * w) + C1 * u +
                          C3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                                w * (b * d + u * w)),
                      exp_a * C2 * (b * d + u * w) + C1 * v +
                          C3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                                w * (b * c - v * w)),

                      exp_a * C1 * c + C2 * (-b * u + d * w) +
                          C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                                d * (c * d - u * v)),
                      exp_a * C2 * (b * c - v * w) - C1 * u +
                          C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                                v * (c * d - u * v)),
                      exp_a * C0 + C2 * (c2 - u2 - w2),
                      exp_a * C2 * (c * d - u * v) + C1 * w +
                          C3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                                w * (-c2 + u2 + w2)),

                      exp_a * C1 * d + C2 * (-b * v - c * w) +
                          C3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                                d * (-d2 + v2 + w2)),
                      exp_a * C2 * (b * d + u * w) - C1 * v +
                          C3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                                v * (-d2 + v2 + w2)),
                      exp_a * C2 * (c * d - u * v) - C1 * w +
                          C3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                                w * (-d2 + v2 + w2)),
                      exp_a * C0 + C2 * (d2 - v2 - w2)};

      for (Index j = 0; j < nq; j++) {
        {
          const Numeric da = -0.5 * (rv * dk1v(j, i).A() +
                                     dr1v[j] * (k1v[i].A() + k2v[i].A())),
                        db = -0.5 * (rv * dk1v(j, i).B() +
                                     dr1v[j] * (k1v[i].B() + k2v[i].B())),
                        dc = -0.5 * (rv * dk1v(j, i).C() +
                                     dr1v[j] * (k1v[i].C() + k2v[i].C())),
                        dd = -0.5 * (rv * dk1v(j, i).D() +
                                     dr1v[j] * (k1v[i].D() + k2v[i].D())),
                        du = -0.5 * (rv * dk1v(j, i).U() +
                                     dr1v[j] * (k1v[i].U() + k2v[i].U())),
                        dv = -0.5 * (rv * dk1v(j, i).V() +
                                     dr1v[j] * (k1v[i].V() + k2v[i].V())),
                        dw = -0.5 * (rv * dk1v(j, i).W() +
                                     dr1v[j] * (k1v[i].W() + k2v[i].W()));
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                        du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          const Numeric dtmp =
              2 * w2 * dw2 +
              2 * (db2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                   b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2) +
                   dc2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                   c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2) +
                   dd2 * (d2 * 0.5 + u2 - v2 - w2) +
                   d2 * (dd2 * 0.5 + du2 - dv2 - dw2) +
                   du2 * (u2 * 0.5 + v2 + w2) + u2 * (du2 * 0.5 + dv2 + dw2) +
                   dv2 * (v2 * 0.5 + w2) + v2 * (dv2 * 0.5 + dw2) +
                   4 * (db * d * u * w - db * c * v * w - dc * d * u * v +
                        b * dd * u * w - b * dc * v * w - c * dd * u * v +
                        b * d * du * w - b * c * dv * w - c * d * du * v +
                        b * d * u * dw - b * c * v * dw - c * d * u * dv));
          const Complex dConst1 = 0.5 * dtmp / Const1;
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          const Complex dx =
              x_zero ? 0 : (0.5 * (dConst2 + dConst1) / tmp_x_sqrt) * sqrt_05;
          const Complex dy = y_zero ? 0
                                    : (0.5 * (dConst2 - dConst1) / tmp_y_sqrt) *
                                          sqrt_05 * Complex(0, 1);
          const Complex dx2 = 2.0 * x * dx;
          const Complex dy2 = 2.0 * y * dy;
          const Complex dcy = -sy * dy;
          const Complex dsy = cy * dy;
          const Complex dcx = sx * dx;
          const Complex dsx = cx * dx;
          const Complex dix = -dx * ix * ix;
          const Complex diy = -dy * iy * iy;
          const Complex dx2dy2 = dx2 + dy2;
          const Complex dC0c =
              either_zero
                  ? 0.0
                  : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0c * dx2dy2) *
                        inv_x2y2;
          const Complex dC1c =
              either_zero ? 0.0
                          : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy +
                             dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix -
                             C1c * dx2dy2) *
                                inv_x2y2;
          const Complex dC2c =
              both_zero ? 0.0 : (dcx - dcy - C2c * dx2dy2) * inv_x2y2;
          const Complex dC3c =
              both_zero
                  ? 0.0
                  : ((x_zero   ? -dsy * iy - sy * diy
                      : y_zero ? dsx * ix + sx * dix
                               : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                     C3c * dx2dy2) *
                        inv_x2y2;

          const Numeric &dC0 = real_val(dC0c);
          const Numeric &dC1 = real_val(dC1c);
          const Numeric &dC2 = real_val(dC2c);
          const Numeric &dC3 = real_val(dC3c);
          dt1v(j, i) = muelmat{
              tv[i](0, 0) * da + exp_a * dC0 + dC2 * (b2 + c2 + d2) +
                  C2 * (db2 + dc2 + dd2),
              tv[i](0, 1) * da + exp_a * db * C1 + b * dC1 +
                  dC2 * (-c * u - d * v) +
                  C2 * (-dc * u - dd * v - c * du - d * dv) +
                  dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                         v * (b * v + c * w)) +
                  C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                        dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                        u * (db * u - dd * w) - v * (db * v + dc * w) -
                        u * (b * du - d * dw) - v * (b * dv + c * dw)),
              tv[i](0, 2) * da + exp_a * dC1 * c + C1 * dc +
                  dC2 * (b * u - d * w) +
                  C2 * (db * u - dd * w + b * du - d * dw) +
                  dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                         w * (b * v + c * w)) +
                  C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                        dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                        u * (dc * u + dd * v) - w * (db * v + dc * w) -
                        u * (c * du + d * dv) - w * (b * dv + c * dw)),
              tv[i](0, 3) * da + exp_a * dC1 * d + C1 * dd +
                  dC2 * (b * v + c * w) +
                  C2 * (db * v + dc * w + b * dv + c * dw) +
                  dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                         w * (b * u - d * w)) +
                  C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                        dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                        v * (dc * u + dd * v) + w * (db * u - dd * w) -
                        v * (c * du + d * dv) + w * (b * du - d * dw)),

              tv[i](1, 0) * da + exp_a * db * C1 + b * dC1 +
                  dC2 * (c * u + d * v) +
                  C2 * (dc * u + dd * v + c * du + d * dv) +
                  dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                         d * (b * d + u * w)) +
                  C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                        dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                        c * (db * c - dv * w) + d * (db * d + du * w) +
                        c * (b * dc - v * dw) + d * (b * dd + u * dw)),
              tv[i](1, 1) * da + exp_a * dC0 + dC2 * (b2 - u2 - v2) +
                  C2 * (db2 - du2 - dv2),
              tv[i](1, 2) * da + exp_a * dC2 * (b * c - v * w) +
                  C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u + C1 * du +
                  dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                         w * (b * d + u * w)) +
                  C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                        dw * (b * d + u * w) + c * (dc * u + dd * v) -
                        u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                        c * (c * du + d * dv) - w * (b * dd + u * dw)),
              tv[i](1, 3) * da + exp_a * dC2 * (b * d + u * w) +
                  C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v + C1 * dv +
                  dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                         w * (b * c - v * w)) +
                  C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                        dw * (b * c - v * w) + d * (dc * u + dd * v) -
                        v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                        d * (c * du + d * dv) + w * (b * dc - v * dw)),

              tv[i](2, 0) * da + exp_a * dC1 * c + C1 * dc +
                  dC2 * (-b * u + d * w) +
                  C2 * (-db * u + dd * w - b * du + d * dw) +
                  dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                         d * (c * d - u * v)) +
                  C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                        dd * (c * d - u * v) + b * (db * c - dv * w) -
                        c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                        b * (b * dc - v * dw) + d * (c * dd - u * dv)),
              tv[i](2, 1) * da + exp_a * dC2 * (b * c - v * w) +
                  C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u - C1 * du +
                  dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                         v * (c * d - u * v)) +
                  C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                        dv * (c * d - u * v) - b * (db * u - dd * w) +
                        u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                        b * (b * du - d * dw) - v * (c * dd - u * dv)),
              tv[i](2, 2) * da + exp_a * dC0 + dC2 * (c2 - u2 - w2) +
                  C2 * (dc2 - du2 - dw2),
              tv[i](2, 3) * da + exp_a * dC2 * (c * d - u * v) +
                  C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w + C1 * dw +
                  dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                         w * (-c2 + u2 + w2)) +
                  C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                        dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                        v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                        d * (b * du - d * dw) + v * (b * dc - v * dw)),

              tv[i](3, 0) * da + exp_a * dC1 * d + C1 * dd +
                  dC2 * (-b * v - c * w) +
                  C2 * (-db * v - dc * w - b * dv - c * dw) +
                  dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                         d * (-d2 + v2 + w2)) +
                  C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                        dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                        c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                        b * (b * dd + u * dw) + c * (c * dd - u * dv)),
              tv[i](3, 1) * da + exp_a * dC2 * (b * d + u * w) +
                  C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v - C1 * dv +
                  dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                         v * (-d2 + v2 + w2)) +
                  C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                        dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                        u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                        b * (b * dv + c * dw) - u * (c * dd - u * dv)),
              tv[i](3, 2) * da + exp_a * dC2 * (c * d - u * v) +
                  C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w - C1 * dw +
                  dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                         w * (-d2 + v2 + w2)) +
                  C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                        dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                        u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                        c * (b * dv + c * dw) + u * (b * dd + u * dw)),
              tv[i](3, 3) * da + exp_a * dC0 + dC2 * (d2 - v2 - w2) +
                  C2 * (dd2 - dv2 - dw2)};
        }

        {
          const Numeric da = -0.5 * (rv * dk2v(j, i).A() +
                                     dr2v[j] * (k1v[i].A() + k2v[i].A())),
                        db = -0.5 * (rv * dk2v(j, i).B() +
                                     dr2v[j] * (k1v[i].B() + k2v[i].B())),
                        dc = -0.5 * (rv * dk2v(j, i).C() +
                                     dr2v[j] * (k1v[i].C() + k2v[i].C())),
                        dd = -0.5 * (rv * dk2v(j, i).D() +
                                     dr2v[j] * (k1v[i].D() + k2v[i].D())),
                        du = -0.5 * (rv * dk2v(j, i).U() +
                                     dr2v[j] * (k1v[i].U() + k2v[i].U())),
                        dv = -0.5 * (rv * dk2v(j, i).V() +
                                     dr2v[j] * (k1v[i].V() + k2v[i].V())),
                        dw = -0.5 * (rv * dk2v(j, i).W() +
                                     dr2v[j] * (k1v[i].W() + k2v[i].W()));
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                        du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          const Numeric dtmp =
              2 * w2 * dw2 +
              2 * (db2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                   b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2) +
                   dc2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                   c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2) +
                   dd2 * (d2 * 0.5 + u2 - v2 - w2) +
                   d2 * (dd2 * 0.5 + du2 - dv2 - dw2) +
                   du2 * (u2 * 0.5 + v2 + w2) + u2 * (du2 * 0.5 + dv2 + dw2) +
                   dv2 * (v2 * 0.5 + w2) + v2 * (dv2 * 0.5 + dw2) +
                   4 * (db * d * u * w - db * c * v * w - dc * d * u * v +
                        b * dd * u * w - b * dc * v * w - c * dd * u * v +
                        b * d * du * w - b * c * dv * w - c * d * du * v +
                        b * d * u * dw - b * c * v * dw - c * d * u * dv));
          const Complex dConst1 = 0.5 * dtmp / Const1;
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          const Complex dx =
              x_zero ? 0 : (0.5 * (dConst2 + dConst1) / tmp_x_sqrt) * sqrt_05;
          const Complex dy = y_zero ? 0
                                    : (0.5 * (dConst2 - dConst1) / tmp_y_sqrt) *
                                          sqrt_05 * Complex(0, 1);
          const Complex dx2 = 2.0 * x * dx;
          const Complex dy2 = 2.0 * y * dy;
          const Complex dcy = -sy * dy;
          const Complex dsy = cy * dy;
          const Complex dcx = sx * dx;
          const Complex dsx = cx * dx;
          const Complex dix = -dx * ix * ix;
          const Complex diy = -dy * iy * iy;
          const Complex dx2dy2 = dx2 + dy2;
          const Complex dC0c =
              either_zero
                  ? 0.0
                  : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0c * dx2dy2) *
                        inv_x2y2;
          const Complex dC1c =
              either_zero ? 0.0
                          : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy +
                             dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix -
                             C1c * dx2dy2) *
                                inv_x2y2;
          const Complex dC2c =
              both_zero ? 0.0 : (dcx - dcy - C2c * dx2dy2) * inv_x2y2;
          const Complex dC3c =
              both_zero
                  ? 0.0
                  : ((x_zero   ? -dsy * iy - sy * diy
                      : y_zero ? dsx * ix + sx * dix
                               : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                     C3c * dx2dy2) *
                        inv_x2y2;

          const Numeric &dC0 = real_val(dC0c);
          const Numeric &dC1 = real_val(dC1c);
          const Numeric &dC2 = real_val(dC2c);
          const Numeric &dC3 = real_val(dC3c);
          dt2v(j, i) = muelmat{
              tv[i](0, 0) * da + exp_a * dC0 + dC2 * (b2 + c2 + d2) +
                  C2 * (db2 + dc2 + dd2),
              tv[i](0, 1) * da + exp_a * db * C1 + b * dC1 +
                  dC2 * (-c * u - d * v) +
                  C2 * (-dc * u - dd * v - c * du - d * dv) +
                  dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                         v * (b * v + c * w)) +
                  C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                        dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                        u * (db * u - dd * w) - v * (db * v + dc * w) -
                        u * (b * du - d * dw) - v * (b * dv + c * dw)),
              tv[i](0, 2) * da + exp_a * dC1 * c + C1 * dc +
                  dC2 * (b * u - d * w) +
                  C2 * (db * u - dd * w + b * du - d * dw) +
                  dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                         w * (b * v + c * w)) +
                  C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                        dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                        u * (dc * u + dd * v) - w * (db * v + dc * w) -
                        u * (c * du + d * dv) - w * (b * dv + c * dw)),
              tv[i](0, 3) * da + exp_a * dC1 * d + C1 * dd +
                  dC2 * (b * v + c * w) +
                  C2 * (db * v + dc * w + b * dv + c * dw) +
                  dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                         w * (b * u - d * w)) +
                  C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                        dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                        v * (dc * u + dd * v) + w * (db * u - dd * w) -
                        v * (c * du + d * dv) + w * (b * du - d * dw)),

              tv[i](1, 0) * da + exp_a * db * C1 + b * dC1 +
                  dC2 * (c * u + d * v) +
                  C2 * (dc * u + dd * v + c * du + d * dv) +
                  dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                         d * (b * d + u * w)) +
                  C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                        dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                        c * (db * c - dv * w) + d * (db * d + du * w) +
                        c * (b * dc - v * dw) + d * (b * dd + u * dw)),
              tv[i](1, 1) * da + exp_a * dC0 + dC2 * (b2 - u2 - v2) +
                  C2 * (db2 - du2 - dv2),
              tv[i](1, 2) * da + exp_a * dC2 * (b * c - v * w) +
                  C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u + C1 * du +
                  dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                         w * (b * d + u * w)) +
                  C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                        dw * (b * d + u * w) + c * (dc * u + dd * v) -
                        u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                        c * (c * du + d * dv) - w * (b * dd + u * dw)),
              tv[i](1, 3) * da + exp_a * dC2 * (b * d + u * w) +
                  C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v + C1 * dv +
                  dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                         w * (b * c - v * w)) +
                  C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                        dw * (b * c - v * w) + d * (dc * u + dd * v) -
                        v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                        d * (c * du + d * dv) + w * (b * dc - v * dw)),

              tv[i](2, 0) * da + exp_a * dC1 * c + C1 * dc +
                  dC2 * (-b * u + d * w) +
                  C2 * (-db * u + dd * w - b * du + d * dw) +
                  dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                         d * (c * d - u * v)) +
                  C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                        dd * (c * d - u * v) + b * (db * c - dv * w) -
                        c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                        b * (b * dc - v * dw) + d * (c * dd - u * dv)),
              tv[i](2, 1) * da + exp_a * dC2 * (b * c - v * w) +
                  C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u - C1 * du +
                  dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                         v * (c * d - u * v)) +
                  C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                        dv * (c * d - u * v) - b * (db * u - dd * w) +
                        u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                        b * (b * du - d * dw) - v * (c * dd - u * dv)),
              tv[i](2, 2) * da + exp_a * dC0 + dC2 * (c2 - u2 - w2) +
                  C2 * (dc2 - du2 - dw2),
              tv[i](2, 3) * da + exp_a * dC2 * (c * d - u * v) +
                  C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w + C1 * dw +
                  dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                         w * (-c2 + u2 + w2)) +
                  C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                        dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                        v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                        d * (b * du - d * dw) + v * (b * dc - v * dw)),

              tv[i](3, 0) * da + exp_a * dC1 * d + C1 * dd +
                  dC2 * (-b * v - c * w) +
                  C2 * (-db * v - dc * w - b * dv - c * dw) +
                  dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                         d * (-d2 + v2 + w2)) +
                  C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                        dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                        c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                        b * (b * dd + u * dw) + c * (c * dd - u * dv)),
              tv[i](3, 1) * da + exp_a * dC2 * (b * d + u * w) +
                  C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v - C1 * dv +
                  dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                         v * (-d2 + v2 + w2)) +
                  C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                        dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                        u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                        b * (b * dv + c * dw) - u * (c * dd - u * dv)),
              tv[i](3, 2) * da + exp_a * dC2 * (c * d - u * v) +
                  C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w - C1 * dw +
                  dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                         w * (-d2 + v2 + w2)) +
                  C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                        dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                        u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                        c * (b * dv + c * dw) + u * (b * dd + u * dw)),
              tv[i](3, 3) * da + exp_a * dC0 + dC2 * (d2 - v2 - w2) +
                  C2 * (dd2 - dv2 - dw2)};
        }
      }
    }
  }
}

void two_level_exp_nopolar(muelmat_vector_view &tv,
                           const propmat_vector_const_view &k1v,
                           const propmat_vector_const_view &k2v,
                           const Numeric rv) {
  ARTS_ASSERT(k2v.nelem() == k1v.nelem());
  ARTS_ASSERT(tv.nelem() == k1v.nelem());

  std::transform(k1v.begin(),
                 k1v.end(),
                 k2v.begin(),
                 tv.begin(),
                 [rv](const propmat &a, const propmat &b) {
                   return exp(-0.5 * rv * (a + b));
                 });
}
}  // namespace rtepack
