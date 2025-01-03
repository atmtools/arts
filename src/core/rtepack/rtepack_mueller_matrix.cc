#include "rtepack_mueller_matrix.h"

#include <configtypes.h>

namespace rtepack {
void forward_cumulative_transmission(Array<muelmat_vector> &Pi,
                                     const Array<muelmat_vector> &T) {
  const Size N = T.size();

  Pi.resize(N);

  if (N == 0) return;

  const Index nv = T.front().size();
  for (auto &p : Pi) p.resize(nv);

  Pi.front() = T.front();

  for (Size i = 1; i < N; i++) {
    for (auto &&[Pi1, Pi0, T] :
         std::ranges::views::zip(Pi[i], Pi[i - 1], T[i])) {
      Pi1 = Pi0 * T;
    }
  }
}

Array<muelmat_vector> forward_cumulative_transmission(
    const Array<muelmat_vector> &T) {
  Array<muelmat_vector> out;
  forward_cumulative_transmission(out, T);
  return out;
}
}  // namespace rtepack
