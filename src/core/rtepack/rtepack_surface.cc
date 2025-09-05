#include "rtepack_surface.h"

#include <arts_constexpr_math.h>

namespace rtepack {
muelmat fresnel_reflectance(Complex Rv, Complex Rh) {
  const Numeric rv    = std::norm(Rv);
  const Numeric rh    = std::norm(Rh);
  const Numeric rmean = 0.5 * (rv + rh);
  const Numeric rdiff = 0.5 * (rv - rh);
  const Complex a     = Rh * std::conj(Rv);
  const Complex b     = Rv * std::conj(Rh);
  const Numeric c     = 0.5 * std::real(a + b);
  const Numeric d     = 0.5 * std::imag(a - b);

  muelmat out{};

  out[0, 0] = rmean;
  out[1, 0] = rdiff;
  out[0, 1] = rdiff;
  out[1, 1] = rmean;
  out[2, 2] = c;
  out[2, 3] = d;
  out[3, 2] = -d;
  out[3, 3] = c;

  return out;
}

stokvec flat_scalar_reflection(stokvec I, const Numeric R, const stokvec B) {
  I      = I * R;
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  I     += B * (1.0 - R);
  return I;
}

stokvec dflat_scalar_reflection_dr(stokvec I, const stokvec B) {
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  I     -= B;
  return I;
}

stokvec reflection(stokvec I, const muelmat R, const stokvec B) {
  I      = R * I;
  I.V() *= -1.0;  //! FIXME: Is this correct?
  I     += (1.0 - R) * B;
  return I;
}

stokvec dreflection(stokvec I, const muelmat dR, const stokvec B) {
  I      = dR * I;
  I.V() *= -1.0;  //! FIXME: Is this correct?
  I     -= dR * B;
  return I;
}
}  // namespace rtepack
