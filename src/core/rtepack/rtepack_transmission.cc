#include "rtepack_transmission.h"

#include <algorithm>

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"

namespace rtepack {
static constexpr Numeric lower_is_considered_zero_for_sinc_likes = 1e-4;

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
    a     = -0.5 * r * (k1.A() + k2.A());
    b     = -0.5 * r * (k1.B() + k2.B());
    c     = -0.5 * r * (k1.C() + k2.C());
    d     = -0.5 * r * (k1.D() + k2.D());
    u     = -0.5 * r * (k1.U() + k2.U());
    v     = -0.5 * r * (k1.V() + k2.V());
    w     = -0.5 * r * (k1.W() + k2.W());
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

    /*
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

    x_zero      = x < lower_is_considered_zero_for_sinc_likes;
    y_zero      = y < lower_is_considered_zero_for_sinc_likes;
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
    // Do the calculation exp(a) * (C0 * I + C1 * K + C2 * K^2 + C3 * K^3)
    return unpolarized
               ? muelmat{exp_a}
               : muelmat{
                     exp_a * (C0 + C2 * (b2 + c2 + d2)),
                     -exp_a * (-C1 * b + C2 * (c * u + d * v) +
                               C3 * (u * (b * u - d * w) - b * (b2 + c2 + d2) +
                                     v * (b * v + c * w))),
                     exp_a * (C1 * c + C2 * (b * u - d * w) +
                              C3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                                    w * (b * v + c * w))),
                     exp_a * (C1 * d + C2 * (b * v + c * w) +
                              C3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                                    w * (b * u - d * w))),
                     exp_a * (C1 * b + C2 * (c * u + d * v) +
                              C3 * (c * (b * c - v * w) - b * (-b2 + u2 + v2) +
                                    d * (b * d + u * w))),
                     exp_a * (C0 + C2 * (b2 - u2 - v2)),
                     exp_a * (C1 * u + C2 * (b * c - v * w) +
                              C3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                                    w * (b * d + u * w))),
                     exp_a * (C1 * v + C2 * (b * d + u * w) +
                              C3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                                    w * (b * c - v * w))),
                     exp_a * (C1 * c + C2 * (d * w - b * u) +
                              C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                                    d * (c * d - u * v))),
                     -exp_a * (C1 * u + C2 * (v * w - b * c) +
                               C3 * (b * (b * u - d * w) - u * (-c2 + u2 + w2) +
                                     v * (c * d - u * v))),
                     exp_a * (C0 + C2 * (c2 - u2 - w2)),
                     exp_a * (C1 * w + C2 * (c * d - u * v) +
                              C3 * (v * (b * c - v * w) - d * (b * u - d * w) -
                                    w * (-c2 + u2 + w2))),
                     exp_a * (C1 * d - C2 * (b * v + c * w) +
                              C3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                                    d * (-d2 + v2 + w2))),
                     -exp_a * (C1 * v - C2 * (b * d + u * w) +
                               C3 * (b * (b * v + c * w) + u * (c * d - u * v) -
                                     v * (-d2 + v2 + w2))),
                     -exp_a * (C1 * w + C2 * (u * v - c * d) +
                               C3 * (c * (b * v + c * w) - u * (b * d + u * w) -
                                     w * (-d2 + v2 + w2))),
                     exp_a * (C0 + C2 * (d2 - v2 - w2))};
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
    const Numeric dC = -2 * (b * w - c * v + d * u) *
                       (b * dw - c * dv + d * du + u * dd - v * dc + w * db);
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
    const Numeric dC2 = both_zero ? 0.0
                                  : ((x_zero ? 0.0 : (dcx - C2 * dx2)) -
                                     (y_zero ? 0.0 : (dcy + C2 * dy2))) *
                                        inv_x2y2;
    const Numeric dC3 =
        both_zero ? 0.0
                  : ((x_zero ? 0.0 : (dsx * ix + sx * dix - C3 * dx2)) -
                     (y_zero ? 0.0 : (dsy * iy + sy * diy + C3 * dy2))) *
                        inv_x2y2;

    return {
        (dC0 + dC2 * (b2 + c2 + d2) + C2 * (db2 + dc2 + dd2)) * exp_a +
            da * t[0, 0],
        (db * C1 + b * dC1 + dC2 * (-c * u - d * v) +
         C2 * (-dc * u - dd * v - c * du - d * dv) +
         dC3 *
             (b * (b2 + c2 + d2) - u * (b * u - d * w) - v * (b * v + c * w)) +
         C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
               dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
               u * (db * u - dd * w) - v * (db * v + dc * w) -
               u * (b * du - d * dw) - v * (b * dv + c * dw))) *
                exp_a +
            da * t[0, 1],
        (dC1 * c + C1 * dc + dC2 * (b * u - d * w) +
         C2 * (db * u - dd * w + b * du - d * dw) +
         dC3 *
             (c * (b2 + c2 + d2) - u * (c * u + d * v) - w * (b * v + c * w)) +
         C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
               dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
               u * (dc * u + dd * v) - w * (db * v + dc * w) -
               u * (c * du + d * dv) - w * (b * dv + c * dw))) *
                exp_a +
            da * t[0, 2],
        (dC1 * d + C1 * dd + dC2 * (b * v + c * w) +
         C2 * (db * v + dc * w + b * dv + c * dw) +
         dC3 *
             (d * (b2 + c2 + d2) - v * (c * u + d * v) + w * (b * u - d * w)) +
         C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
               dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
               v * (dc * u + dd * v) + w * (db * u - dd * w) -
               v * (c * du + d * dv) + w * (b * du - d * dw))) *
                exp_a +
            da * t[0, 3],

        (db * C1 + b * dC1 + dC2 * (c * u + d * v) +
         C2 * (dc * u + dd * v + c * du + d * dv) +
         dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                d * (b * d + u * w)) +
         C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
               dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
               c * (db * c - dv * w) + d * (db * d + du * w) +
               c * (b * dc - v * dw) + d * (b * dd + u * dw))) *
                exp_a +
            da * t[1, 0],
        (dC0 + dC2 * (b2 - u2 - v2) + C2 * (db2 - du2 - dv2)) * exp_a +
            da * t[1, 1],
        (dC2 * (b * c - v * w) + C2 * (db * c + b * dc - dv * w - v * dw) +
         dC1 * u + C1 * du +
         dC3 *
             (c * (c * u + d * v) - u * (-b2 + u2 + v2) - w * (b * d + u * w)) +
         C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
               dw * (b * d + u * w) + c * (dc * u + dd * v) -
               u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
               c * (c * du + d * dv) - w * (b * dd + u * dw))) *
                exp_a +
            da * t[1, 2],
        (dC2 * (b * d + u * w) + C2 * (db * d + b * dd + du * w + u * dw) +
         dC1 * v + C1 * dv +
         dC3 *
             (d * (c * u + d * v) - v * (-b2 + u2 + v2) + w * (b * c - v * w)) +
         C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
               dw * (b * c - v * w) + d * (dc * u + dd * v) -
               v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
               d * (c * du + d * dv) + w * (b * dc - v * dw))) *
                exp_a +
            da * t[1, 3],

        (dC1 * c + C1 * dc + dC2 * (-b * u + d * w) +
         C2 * (-db * u + dd * w - b * du + d * dw) +
         dC3 *
             (b * (b * c - v * w) - c * (-c2 + u2 + w2) + d * (c * d - u * v)) +
         C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
               dd * (c * d - u * v) + b * (db * c - dv * w) -
               c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
               b * (b * dc - v * dw) + d * (c * dd - u * dv))) *
                exp_a +
            da * t[2, 0],
        (dC2 * (b * c - v * w) + C2 * (db * c + b * dc - dv * w - v * dw) -
         dC1 * u - C1 * du +
         dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                v * (c * d - u * v)) +
         C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
               dv * (c * d - u * v) - b * (db * u - dd * w) +
               u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
               b * (b * du - d * dw) - v * (c * dd - u * dv))) *
                exp_a +
            da * t[2, 1],
        (dC0 + dC2 * (c2 - u2 - w2) + C2 * (dc2 - du2 - dw2)) * exp_a +
            da * t[2, 2],
        (dC2 * (c * d - u * v) + C2 * (dc * d + c * dd - du * v - u * dv) +
         dC1 * w + C1 * dw +
         dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                w * (-c2 + u2 + w2)) +
         C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
               dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
               v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
               d * (b * du - d * dw) + v * (b * dc - v * dw))) *
                exp_a +
            da * t[2, 3],

        (dC1 * d + C1 * dd + dC2 * (-b * v - c * w) +
         C2 * (-db * v - dc * w - b * dv - c * dw) +
         dC3 *
             (b * (b * d + u * w) + c * (c * d - u * v) - d * (-d2 + v2 + w2)) +
         C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
               dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
               c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
               b * (b * dd + u * dw) + c * (c * dd - u * dv))) *
                exp_a +
            da * t[3, 0],
        (dC2 * (b * d + u * w) + C2 * (db * d + b * dd + du * w + u * dw) -
         dC1 * v - C1 * dv +
         dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                v * (-d2 + v2 + w2)) +
         C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
               dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
               u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
               b * (b * dv + c * dw) - u * (c * dd - u * dv))) *
                exp_a +
            da * t[3, 1],
        (dC2 * (c * d - u * v) + C2 * (dc * d + c * dd - du * v - u * dv) -
         dC1 * w - C1 * dw +
         dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                w * (-d2 + v2 + w2)) +
         C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
               dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
               u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
               c * (b * dv + c * dw) + u * (b * dd + u * dw))) *
                exp_a +
            da * t[3, 2],
        (C2 * (dd2 - dv2 - dw2) + dC2 * (d2 - v2 - w2) + dC0) * exp_a +
            da * t[3, 3]};
  }
};

void two_level_exp(muelmat &t,
                   muelmat_vector_view dt1,
                   muelmat_vector_view dt2,
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

void two_level_exp(muelmat_vector_view tv,
                   muelmat_matrix_view dt1v,
                   muelmat_matrix_view dt2v,
                   const propmat_vector_const_view &k1v,
                   const propmat_vector_const_view &k2v,
                   const propmat_matrix_const_view &dk1v,
                   const propmat_matrix_const_view &dk2v,
                   const Numeric rv,
                   const ExhaustiveConstVectorView &dr1v,
                   const ExhaustiveConstVectorView &dr2v) {
  const Index nf = tv.size();
  const Index nq = dr1v.size();

  ARTS_ASSERT(nf == k1v.size());
  ARTS_ASSERT(nf == k2v.size());
  ARTS_ASSERT(nf == dk1v.ncols());
  ARTS_ASSERT(nf == dk2v.ncols());
  ARTS_ASSERT(nq == dk1v.nrows());
  ARTS_ASSERT(nq == dk2v.nrows());
  ARTS_ASSERT(nf == dt1v.ncols());
  ARTS_ASSERT(nf == dt2v.ncols());
  ARTS_ASSERT(nq == dt1v.nrows());
  ARTS_ASSERT(nq == dt2v.nrows());
  ARTS_ASSERT(nq == dr2v.size());

  for (Index i = 0; i < nf; ++i) {
    const tran tran_state{k1v[i], k2v[i], rv};
    tv[i] = tran_state();

    for (Index j = 0; j < nq; j++) {
      dt1v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk1v[j, i], rv, dr1v[j]);
      dt2v[j, i] =
          tran_state.deriv(tv[i], k1v[i], k2v[i], dk2v[j, i], rv, dr2v[j]);
    }
  }
}

muelmat exp(propmat k, Numeric r) { return tran(k, k, r)(); }

void two_level_exp(muelmat_vector_view tv,
                   const propmat_vector_const_view &k1v,
                   const propmat_vector_const_view &k2v,
                   const Numeric rv) {
  ARTS_ASSERT(k2v.size() == k1v.size());
  ARTS_ASSERT(tv.size() == k1v.size());

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

  ARTS_USER_ERROR_IF(N != static_cast<Size>(r.size()),
                     "Must have same number of levels ({}) in K and r",
                     N);

  ARTS_USER_ERROR_IF(N != static_cast<Size>(dr.nrows()),
                     "Must have same number of levels ({}) in K and dr",
                     N);

  T.resize(N);

  dT.resize(N);

  if (N == 0) return;

  const Index nv = K[0].size();
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
      std::ranges::any_of(K, Cmp::ne(nv), &propmat_vector::size),
      "Must have same number of frequency elements ({}) in all K:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dK, Cmp::ne(nv), &propmat_matrix::ncols),
      "Must have same number of frequency elements ({}) in all dK:s as in K[0]",
      nv);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dK, Cmp::ne(nq), &propmat_matrix::nrows),
      "Must have same number of derivative elements ({}) in all dK:s as in dr",
      nq);

  ARTS_USER_ERROR_IF(
      dr.npages() != 2,
      "Must have 2 as first dimension in dr (upper and lower level distance derivatives), got {}",
      dr.npages());

  for (Size i = 1; i < N; i++) {
    two_level_exp(T[i],
                  dT[i - 1][0],
                  dT[i][1],
                  K[i - 1],
                  K[i],
                  dK[i - 1],
                  dK[i],
                  r[i],
                  dr[0][i],
                  dr[1][i]);
  }
}
}  // namespace rtepack
