#include "matpack_mdspan_helpers.h"

#include <arts_constants.h>

namespace matpack {
Vector uniform_grid(Numeric x0, Index N, Numeric dx) {
  Vector out(N);
  std::generate(out.begin(), out.end(), [x = x0, dx]() mutable {
    auto xd  = x;
    x       += dx;
    return xd;
  });
  return out;
}

ComplexVector uniform_grid(Complex x0, Index N, Complex dx) {
  ComplexVector out(N);
  std::generate(out.begin(), out.end(), [x = x0, dx]() mutable {
    auto xd  = x;
    x       += dx;
    return xd;
  });
  return out;
}

void gaussian_grid(Vector &y,
                   const Vector &x,
                   const Numeric &x0,
                   const Numeric &si,
                   const Numeric &fwhm) {
  ARTS_USER_ERROR_IF(
      (si <= 0 && fwhm <= 0) || (si > 0 && fwhm > 0),
      "One of the GINs *si* and *fwhm* shall be >0, but just one.");

  const Index n = x.size();

  // Note that y and x can be the same vector
  if (&y != &x) {
    y.resize(n);
  }

  const Numeric si2use =
      si > 0 ? si : fwhm / (2 * std::sqrt(2 * Constant::ln_2));
  const Numeric fac = 1 / (std::sqrt(2 * Constant::pi) * si2use);
  for (Index i = 0; i < n; ++i) {
    y[i] = fac * std::exp(-0.5 * Math::pow2((x[i] - x0) / si2use));
  }
}

Vector gaussian_grid(Vector x,
                     const Numeric &x0,
                     const Numeric &si,
                     const Numeric &fwhm) {
  gaussian_grid(x, x, x0, si, fwhm);
  return x;
}
}  // namespace matpack
