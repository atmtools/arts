#include "rtepack_propagation_matrix.h"

#include <Faddeeva.hh>

namespace rtepack {
propmat_vector operator*(Numeric x, const propmat_vector_const_view &y) {
  propmat_vector z(y.size());
  for (Size i = 0; i < y.size(); ++i) {
    z[i] = x * y[i];
  }
  return z;
}

propmat dawson(const propmat &pm) {
  return {Faddeeva::Dawson(pm.A()),
          Faddeeva::Dawson(pm.B()),
          Faddeeva::Dawson(pm.C()),
          Faddeeva::Dawson(pm.D()),
          Faddeeva::Dawson(pm.U()),
          Faddeeva::Dawson(pm.V()),
          Faddeeva::Dawson(pm.W())};
}
}  // namespace rtepack
