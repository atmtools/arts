#include "rtepack_surface.h"

#include <arts_constexpr_math.h>

namespace rtepack {
muelmat fresnel_reflectance(Complex Rv, Complex Rh) {
  const Numeric rv    = Math::pow2(std::abs(Rv));
  const Numeric rh    = Math::pow2(std::abs(Rh));
  const Complex a     = Rh * std::conj(Rv);
  const Complex b     = Rv * std::conj(Rh);
  const Numeric rmean = (rv + rh) / 2.0;
  const Numeric rdiff = (rv - rh) / 2.0;
  const Numeric c     = std::real(a + b) / 2.0;
  const Numeric d     = std::imag(a - b) / 2.0;

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

stokvec dreflection_dn2(stokvec I, const muelmat dR, const stokvec B) {
  I      = dR * I;
  I.V() *= -1.0;  //! FIXME: Is this correct?
  I     -= dR * B;
  return I;
}
}  // namespace rtepack
