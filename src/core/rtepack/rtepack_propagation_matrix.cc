#include "rtepack_propagation_matrix.h"

namespace rtepack {
propmat_vector operator*(Numeric x, const propmat_vector_const_view& y) {
  propmat_vector z(y.size());
  for (Size i = 0; i < y.size(); ++i) {
    z[i] = x * y[i];
  }
  return z;
}
}  // namespace rtepack