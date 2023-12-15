#include "rtepack_multitype.h"

#include <algorithm>

#include "debug.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
stokvec_vector absvec(const propmat_vector_const_view &k) {
  stokvec_vector out(k.size());
  std::transform(k.begin(), k.end(), out.begin(), [](const propmat &v) {
    return absvec(v);
  });
  return out;
}

Tensor3 to_tensor3(const muelmat_vector_const_view &m) {
  Tensor3 out(m.size(), 4, 4);

  std::transform(m.begin(), m.end(), out.begin(), [](const muelmat &v) {
    return Matrix{v};
  });

  return out;
}

stokvec_vector to_stokvec_vector(const ExhaustiveConstMatrixView &v) {
  ARTS_ASSERT(v.ncols() == 4)

  stokvec_vector out(v.nrows());
  std::transform(
      v.begin(), v.end(), out.begin(), [](const ExhaustiveConstVectorView &a) {
        return stokvec{a[0], a[1], a[2], a[3]};
      });
  return out;
}

Matrix to_matrix(const stokvec_vector_const_view &v) {
  Matrix out(v.size(), 4);
  std::transform(v.begin(), v.end(), out.begin(), [](const stokvec &s) {
    return Vector{s[0], s[1], s[2], s[3]};
  });
  return out;
}

Matrix to_matrix(const propmat &v) {
  return Vector{v.A(),
                v.B(),
                v.C(),
                v.D(),
                v.B(),
                v.A(),
                v.U(),
                v.V(),
                v.C(),
                -v.U(),
                v.A(),
                v.W(),
                v.D(),
                -v.V(),
                -v.W(),
                v.A()}
      .reshape(4, 4);
}

Vector to_vector(const stokvec &v) { return {v.I(), v.Q(), v.U(), v.V()}; }

stokvec to_stokvec(const ConstVectorView &a) {
  ARTS_ASSERT(a.size() == 4, "Must be 4-long")

  return {a[0], a[1], a[2], a[3]};
}

propmat to_propmat(const ConstMatrixView &a) {
  ARTS_ASSERT(a.ncols() == 4 and a.nrows() == 4, "Must be 4x4")

  ARTS_ASSERT(a(0, 0) == a(1, 1) and a(1, 1) == a(2, 2) and a(2, 2) == a(3, 3),
              "Must be 00 - 11 - 22 - 33 symmetric")
  ARTS_ASSERT(a(0, 1) == a(1, 0), "Must be 10 - 01 symmetric")
  ARTS_ASSERT(a(0, 2) == a(2, 0), "Must be 20 - 02 symmetric")
  ARTS_ASSERT(a(0, 3) == a(3, 0), "Must be 30 - 03 symmetric")
  ARTS_ASSERT(a(1, 2) == -a(2, 1), "Must be 12 - 21 asymmetric")
  ARTS_ASSERT(a(1, 3) == -a(3, 1), "Must be 13 - 31 asymmetric")
  ARTS_ASSERT(a(2, 3) == -a(3, 2), "Must be 23 - 32 asymmetric")

  return {a[0][0], a[0][1], a[0][2], a[0][3], a[1][2], a[1][3], a[2][3]};
}

muelmat to_muelmat(const ConstMatrixView &a) {
  ARTS_ASSERT(a.ncols() == 4)
  ARTS_ASSERT(a.nrows() == 4)

  return {a[0][0],
          a[0][1],
          a[0][2],
          a[0][3],
          a[1][0],
          a[1][1],
          a[1][2],
          a[1][3],
          a[2][0],
          a[2][1],
          a[2][2],
          a[2][3],
          a[3][0],
          a[3][1],
          a[3][2],
          a[3][3]};
}
}  // namespace rtepack
