#include "arts_constants.h"

#include "rtepack_multitype.h"

namespace rtepack {
constexpr Numeric lower_is_considered_zero_for_sinc_likes = 1e-4;

muelmat exp(const propmat &k, const Numeric r) {
  static constexpr Numeric sqrt_05 = Constant::inv_sqrt_2;
  const auto [a, b, c, d, u, v, w] = propmat{-r * k};
  const Numeric exp_a = std::exp(a);
  if (b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.)
    return muelmat{exp_a};

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
        : 1.0 / (x2 + y2); // The first "1.0" is the trick for above limits
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
return muelmat{exp_a * C0 + C2 * (b2 + c2 + d2),
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
}
} // namespace rtepack
