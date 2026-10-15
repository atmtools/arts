#include "rtepack_jones_matrix.h"

#include <debug.h>

#include <cmath>

#include "rtepack_multitype.h"

namespace rtepack {
/* M = A (J x J*) A^-1, where A takes the coherencies (E_v E_v*, E_v E_h*,
   E_h E_v*, E_h E_h*) to the Stokes vector, so K = A^-1 M A is J x J* with
   K[2 i + j, 2 k + l] = J[i, k] conj(J[j, l]).  Its element
   K[3 a, 3 b] = |J[a, b]|^2 picks the largest element of J, taken real and
   positive, and K[2 i + a, 2 k + b] = J[i, k] J[a, b] gives the others. */
jones::jones(const muelmat &M) : jones(Complex{}) {
  constexpr Complex     i{0.0, 1.0};
  const ComplexMatrix44 A{1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 0, -i, i, 0};
  const ComplexMatrix44 inv_A{0.5, 0.5, 0, 0, 0, 0, 0.5, 0.5 * i, 0, 0, 0.5, -0.5 * i, 0.5, -0.5, 0, 0};

  ComplexMatrix44 MA{}, K{};
  for (Index r = 0; r < 4; r++)
    for (Index c = 0; c < 4; c++)
      for (Index n = 0; n < 4; n++) MA[r, c] += M[r, n] * A[n, c];
  for (Index r = 0; r < 4; r++)
    for (Index c = 0; c < 4; c++)
      for (Index n = 0; n < 4; n++) K[r, c] += inv_A[r, n] * MA[n, c];

  Index a = 0, b = 0;
  for (Index p = 0; p < 2; p++)
    for (Index q = 0; q < 2; q++)
      if (K[3 * p, 3 * q].real() > K[3 * a, 3 * b].real()) {
        a = p;
        b = q;
      }

  const Numeric largest = K[3 * a, 3 * b].real();
  if (largest > 0.0) {
    const Numeric J_ab = std::sqrt(largest);
    for (Index r = 0; r < 2; r++)
      for (Index c = 0; c < 2; c++) (*this)[r, c] = K[2 * r + a, 2 * c + b] / J_ab;
  }

  const muelmat back      = to_muelmat(*this);
  Numeric       deviation = 0.0, scale = 0.0;
  for (Size n = 0; n < 16; n++) {
    deviation = std::max(deviation, std::abs(back.data[n] - M.data[n]));
    scale     = std::max(scale, std::abs(M.data[n]));
  }
  ARTS_USER_ERROR_IF(deviation > 1e-9 * scale,
                     R"(The Mueller matrix is not that of a Jones matrix (it depolarizes):

{}

The closest Jones matrix gives a Mueller matrix that differs by {} (relative to its largest element).
)",
                     M,
                     deviation / scale);
}
}  // namespace rtepack
