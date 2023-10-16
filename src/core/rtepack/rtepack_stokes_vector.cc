#include "rtepack_stokes_vector.h"

namespace rtepack {
stokvec_vector operator*(Numeric x, const stokvec_vector_const_view& y) {
  stokvec_vector z(y.size());
  for (Index i = 0; i < y.size(); ++i) {
    z[i] = x * y[i];
  }
  return z;
}
}  // namespace rtepack
