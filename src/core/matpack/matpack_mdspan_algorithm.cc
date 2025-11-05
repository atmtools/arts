#include "matpack_mdspan_algorithm.h"

namespace matpack {
std::vector<Range> omp_offset_count(const Index N, const Index n) {
  std::vector<Range> result(n, {0, 0});
  const Index dn       = N / n;
  result.front().nelem = dn;

  for (Index i = 1; i < n - 1; i++) {
    result[i].offset = result[i - 1].offset + dn;
    result[i].nelem  = dn;
  }

  result.back().offset = result[n - 2].offset + dn;
  result.back().nelem  = N - result.back().offset;

  return result;
}
}  // namespace matpack
