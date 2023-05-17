#include "rtepack_mueller_matrix.h"

namespace rtepack {
muelmat_matrix reverse_cumulative_transmission(const Array<muelmat_vector> &T) {
  const Index n = T.nelem();
  const Index nf = n ? T.front().nelem() : 0;

  muelmat_matrix PiT(n, nf);
  for (Index i = 1; i < n; i++) {
    for (Index j = 0; j < nf; j++) {
      PiT(i, j) = T[i][j] * PiT(i - 1, j);
    }
  }
  return PiT;
}

muelmat_matrix forward_cumulative_transmission(const Array<muelmat_vector> &T) {
  const Index n = T.nelem();
  const Index nf = n ? T.front().nelem() : 0;

  muelmat_matrix PiT(n, nf);
  for (Index i = 1; i < n; i++) {
    for (Index j = 0; j < nf; j++) {
      PiT(j, j) = PiT(i - 1, j) * T[i][j];
    }
  }
  return PiT;
}

} // namespace rtepack