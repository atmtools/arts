#pragma once

#include "matpack_mdspan.h"

namespace matpack {
Vector uniform_grid(Numeric x0, Index N, Numeric dx);
ComplexVector uniform_grid(Complex x0, Index N, Complex dx);

template <ranked_md<1> VEC1, ranked_md<1> VEC2>
constexpr Vector3 cross(const VEC1& a, const VEC2& b) noexcept {
  assert(a.size() == 3 and b.size() == 3);
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

template <mut_ranked_md<1> VEC1, ranked_md<1> VEC2, ranked_md<1> VEC3>
constexpr void cross3(VEC1&& x, const VEC2& a, const VEC3& b) noexcept {
  assert(x.size() == 3 and a.size() == 3 and b.size() == 3);
  x[0] = a[1] * b[2] - a[2] * b[1];
  x[1] = a[2] * b[0] - a[0] * b[2];
  x[2] = a[0] * b[1] - a[1] * b[0];
}

void gaussian_grid(Vector &y,
                   const Vector &x,
                   const Numeric &x0,
                   const Numeric &si,
                   const Numeric &fwhm);

Vector gaussian_grid(Vector x,
                     const Numeric &x0,
                     const Numeric &si,
                     const Numeric &fwhm);
}  // namespace matpack
