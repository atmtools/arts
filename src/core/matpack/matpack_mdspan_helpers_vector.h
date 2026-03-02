#pragma once

#include "matpack_mdspan_cdata_t.h"

namespace matpack {
Vector uniform_grid(Numeric x0, Index N, Numeric dx);
ComplexVector uniform_grid(Complex x0, Index N, Complex dx);

constexpr Vector3 cross(const Vector3 &a, const Vector3 &b) noexcept {
  assert(a.size() == 3 and b.size() == 3);
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
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
