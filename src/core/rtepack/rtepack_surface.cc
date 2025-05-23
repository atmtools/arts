#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
stokvec flat_scalar_reflection(stokvec I, const Numeric R, const Numeric B) {
  assert(R >= 0.0 and R <= 1.0);

  I      = I * R;
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  I.I() += B * (1.0 - R);
  return I;
}

stokvec flat_scalar_reflection_pure_reflect(stokvec I, const Numeric R) {
  assert(R >= 0.0 and R <= 1.0);

  I      = I * R;
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  return I;
}

stokvec dflat_scalar_reflection_dr(stokvec I, const Numeric dRdx, const Numeric B) {
  I      = I * dRdx;
  I.V() *= -1.0;  //! NOTE: The elementwise multiplication is [R, R, R, -R]
  I.I() -= B * dRdx;
  return I;
}

stokvec dflat_scalar_reflection_db(const Numeric R, const Numeric dBdx) {
  assert(R >= 0.0 and R <= 1.0);

  return {dBdx * (1.0 - R), 0.0, 0.0, 0.0};
}
}  // namespace rtepack
