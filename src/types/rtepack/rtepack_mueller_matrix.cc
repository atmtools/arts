#include "rtepack_mueller_matrix.h"

namespace rtepack {
Array<muelmat_vector> reverse_cumulative_transmission(const Array<muelmat_vector> &T) {
  const Size n = T.size();
  const Index nf = n ? T.front().nelem() : 0;

  Array<muelmat_vector> PiT(n, muelmat_vector(nf, muelmat::id()));
  for (Size i = 1; i < n; i++) {
    for (Index j = 0; j < nf; j++) {
      PiT[i][j] = T[i][j] * PiT[i - 1][j];
    }
  }
  return PiT;
}

Array<muelmat_vector> forward_cumulative_transmission(const Array<muelmat_vector> &T) {
  const Size n = T.size();
  const Index nf = n ? T.front().nelem() : 0;

  Array<muelmat_vector> PiT(n, muelmat_vector(nf, muelmat::id()));
  for (Size i = 1; i < n; i++) {
    for (Index j = 0; j < nf; j++) {
      PiT[i][j] = PiT[i - 1][j] * T[i][j];
    }
  }
  return PiT;
}
} // namespace rtepack
