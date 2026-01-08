#include "rtepack_multitype.h"

#include <debug.h>

#include <Eigen/Dense>
#include <algorithm>

#include "rtepack_propagation_matrix.h"
#include "rtepack_spectral_matrix.h"
#include "rtepack_stokes_vector.h"
#include "rtepack_transmission.h"

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

stokvec_vector to_stokvec_vector(const ConstMatrixView &v) {
  assert(v.ncols() == 4);

  stokvec_vector out(v.nrows());
  std::transform(v.begin(), v.end(), out.begin(), [](const ConstVectorView &a) {
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

stokvec to_stokvec(const propmat &a) { return {a[0], a[1], a[2], a[3]}; }

stokvec to_stokvec(const ConstVectorView &a) {
  assert(a.size() == 4);

  return {a[0], a[1], a[2], a[3]};
}

propmat to_propmat(const ConstMatrixView &a) {
  assert(a.ncols() == 4 and a.nrows() == 4);

  assert((a[0, 0] == a[1, 1] and a[1, 1] == a[2, 2] and a[2, 2] == a[3, 3]));
  assert((a[0, 1] == a[1, 0]));
  assert((a[0, 2] == a[2, 0]));
  assert((a[0, 3] == a[3, 0]));
  assert((a[1, 2] == -a[2, 1]));
  assert((a[1, 3] == -a[3, 1]));
  assert((a[2, 3] == -a[3, 2]));

  return {a[0][0], a[0][1], a[0][2], a[0][3], a[1][2], a[1][3], a[2][3]};
}

muelmat to_muelmat(const ConstMatrixView &a) {
  assert(a.ncols() == 4);
  assert(a.nrows() == 4);

  return {a[0, 0],
          a[0, 1],
          a[0, 2],
          a[0, 3],
          a[1, 0],
          a[1, 1],
          a[1, 2],
          a[1, 3],
          a[2, 0],
          a[2, 1],
          a[2, 2],
          a[2, 3],
          a[3, 0],
          a[3, 1],
          a[3, 2],
          a[3, 3]};
}

muelmat real(const specmat &A) {
  return muelmat{std::real(A.data[0]),
                 std::real(A.data[1]),
                 std::real(A.data[2]),
                 std::real(A.data[3]),
                 std::real(A.data[4]),
                 std::real(A.data[5]),
                 std::real(A.data[6]),
                 std::real(A.data[7]),
                 std::real(A.data[8]),
                 std::real(A.data[9]),
                 std::real(A.data[10]),
                 std::real(A.data[11]),
                 std::real(A.data[12]),
                 std::real(A.data[13]),
                 std::real(A.data[14]),
                 std::real(A.data[15])};
}

specmat frechet_sqrt(const propmat &X, const propmat &E) {
  specmat A  = sqrt(X);
  muelmat Em = to_muelmat(E);

  // For small matrices, use eigen-decomposition
  Eigen::Matrix4cd a_mat;
  Eigen::Matrix4d e_mat;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      a_mat(i, j) = Complex{A[i, j]};
      e_mat(i, j) = Em[i, j];
    }
  }

  Eigen::ComplexEigenSolver<Eigen::Matrix4cd> es(a_mat);
  const Eigen::Matrix4cd &V = es.eigenvectors();
  Eigen::Matrix4cd D        = es.eigenvalues().asDiagonal();

  // Transform E into eigenbasis
  const auto E_tilde = V.inverse() * e_mat * V;

  // Solve for S_tilde: (d_i + d_j) S_ij = E_ij
  Eigen::Matrix4cd S_tilde = Eigen::Matrix4cd::Zero();
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      S_tilde(i, j) = E_tilde(i, j) / (D(i, i) + D(j, j));
    }
  }

  // Transform back
  const auto S = V * S_tilde * V.inverse();

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      A[i, j] = S(i, j);
    }
  }

  return A;  // Convert back to specmat
}
}  // namespace rtepack
