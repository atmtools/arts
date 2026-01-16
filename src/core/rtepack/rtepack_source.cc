#include "rtepack_source.h"

#include <array_algo.h>
#include <arts_omp.h>
#include <debug.h>
#include <physics_funcs.h>

#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
void level_nlte(stokvec_vector_view J,
                stokvec_matrix_view dJ,
                const propmat_vector_const_view &K,
                const stokvec_vector_const_view &S,
                const propmat_matrix_const_view &dK,
                const stokvec_matrix_const_view &dS,
                const ConstVectorView &f,
                const Numeric &t,
                const Index &it) {
  const Size N = J.size();
  assert(N == static_cast<Size>(dJ.ncols()));
  assert(N == K.size());
  assert(N == S.size());
  assert(N == static_cast<Size>(dK.ncols()));
  assert(N == static_cast<Size>(dS.ncols()));
  assert(N == f.size());

  const Index M = dJ.nrows();
  assert(M == dK.nrows());
  assert(M == dS.nrows());
  assert(M > it);

  for (Size i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i]         = 0.0;
      dJ[joker, i] = 0.0;
    } else {
      const auto b    = stokvec{planck(f[i], t), 0, 0, 0};
      const auto invK = inv(K[i]);
      J[i]            = b + invK * S[i];

      for (Index j = 0; j < M; j++) {
        dJ[j, i]  = {it == j ? dplanck_dt(f[i], t) : 0, 0, 0, 0};
        dJ[j, i] -= invK * (dK[j, i] * invK * S[i] - dS[j, i]);
      }
    }
  }
}

void SourceVector::init(const std::span<const propmat_vector> &K,
                        const std::span<const propmat_matrix> &dK,
                        const std::span<const stokvec_vector> &nlte,
                        const std::span<const stokvec_matrix> &dnlte,
                        const std::span<const AscendingGrid> &freq_grid,
                        const std::span<const Numeric> &ts,
                        const Size &it) {
  ARTS_USER_ERROR_IF(not arr::same_size(K, dK, nlte, dnlte, freq_grid, ts),
                     "All input must have the same size.");

  const Size np = dK.size();

  if (np == 0) {
    J.resize(J.nrows(), 0);
    dJ.resize(dJ.npages(), 0, dJ.ncols());
    return;
  }

  const Size nq = dK.front().nrows();
  const Size nf = dK.front().ncols();

  J.resize(nf, np);
  dJ.resize(nf, np, nq);

  ARTS_USER_ERROR_IF(not all_same_shape({nq, nf}, dK, dnlte),
                     "All derivative parameters must have same shape ({}, {}).",
                     nq,
                     nf);

  ARTS_USER_ERROR_IF(not all_same_shape({nf}, K, nlte, freq_grid),
                     "All forward parameters must have same shape ({}).",
                     nf);

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 0; i < np; i++) {
    for (Size j = 0; j < nf; j++) {
      if (K[i][j].is_rotational()) {
        J[j, i]  = 0.0;
        dJ[j, i] = 0.0;
      } else {
        const auto b    = planck(freq_grid[i][j], ts[i]);
        const auto invK = inv(K[i][j]);
        const auto n    = invK * nlte[i][j];
        J[j, i]         = b + n;

        for (Size k = 0; k < nq; k++) {
          dJ[j, i, k] =
              stokvec{
                  it == k ? dplanck_dt(freq_grid[i][j], ts[i]) : 0, 0, 0, 0} -
              invK * (dK[i][k, j] * n - dnlte[i][k, j]);
        }
      }
    }
  }
}

void SourceVector::init(const std::span<const propmat> &K,
                        const std::span<const propmat_vector> &dK,
                        const std::span<const stokvec> &nlte,
                        const std::span<const stokvec_vector> &dnlte,
                        const std::span<const Numeric> &freq,
                        const std::span<const Numeric> &ts,
                        const Size &it) {
  ARTS_USER_ERROR_IF(not arr::same_size(K, dK, nlte, dnlte, freq, ts),
                     "All input must have the same size.");

  constexpr Size nf = 1;
  const Size np     = dK.size();

  J.resize(nf, np);
  auto &&J_ = J[0];

  if (np == 0) {
    dJ.resize(nf, np, dJ.nrows());
    return;
  }

  const Size nq = dK.front().nrows();

  dJ.resize(nf, np, nq);
  auto &&dJ_ = dJ[0];

  ARTS_USER_ERROR_IF(not all_same_shape({nq}, dK, dnlte),
                     "All derivative parameters must have same shape ({}).",
                     nq);

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size i = 0; i < np; i++) {
    if (K[i].is_rotational()) {
      J_[i]  = 0.0;
      dJ_[i] = 0.0;
    } else {
      const auto b    = planck(freq[i], ts[i]);
      const auto invK = inv(K[i]);
      const auto n    = invK * nlte[i];
      J_[i]           = b + n;

      for (Size k = 0; k < nq; k++) {
        dJ_[i, k] = stokvec{it == k ? dplanck_dt(freq[i], ts[i]) : 0, 0, 0, 0} -
                    invK * (dK[i][k] * n - dnlte[i][k]);
      }
    }
  }
}

void SourceVector::check(Size np,
                         Size nq,
                         Size nf,
                         const std::string_view caller) const {
  ARTS_USER_ERROR_IF(
      not same_shape({nf, np}, J),
      R"(Shape mismatch in source called from {}: J : {:B,} expected, got [{}, {}])",
      caller,
      J.shape(),
      nf,
      np);

  ARTS_USER_ERROR_IF(
      not same_shape({nf, np, nq}, dJ),
      R"(Shape mismatch in derivative called from {}: dJ : {:B,} expected, got [{}, {}, {}])",
      caller,
      dJ.shape(),
      nf,
      np,
      nq);
}

std::array<Size, 3> SourceVector::shape() const noexcept {
  return {static_cast<Size>(dJ.npages()),
          static_cast<Size>(dJ.nrows()),
          static_cast<Size>(dJ.ncols())};
}
}  // namespace rtepack
