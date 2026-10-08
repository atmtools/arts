#include "radintg.h"

#include <debug.h>
#include <lin_alg.h>

#include <array>

namespace polradtran {
void combine_layers(ConstTensor3View reflect1,
                    ConstTensor3View trans1,
                    ConstMatrixView  source1,
                    ConstTensor3View reflect2,
                    ConstTensor3View trans2,
                    ConstMatrixView  source2,
                    Tensor3View      out_reflect,
                    Tensor3View      out_trans,
                    MatrixView       out_source,
                    workdata&        work) {
  const Index n = reflect1.extent(1);
  for (auto shape :
       {reflect1.shape(), trans1.shape(), reflect2.shape(), trans2.shape(), out_reflect.shape(), out_trans.shape()})
    ARTS_USER_ERROR_IF(shape != (std::array<Index, 3>{2, n, n}),
                       "COMBINE_LAYERS needs the reflections and transmissions [2, n, n] (n = {}); got {:B,}",
                       n,
                       shape);
  for (auto shape : {source1.shape(), source2.shape(), out_source.shape()})
    ARTS_USER_ERROR_IF(
        shape != (std::array<Index, 2>{2, n}), "COMBINE_LAYERS needs the sources [2, n] (n = {}); got {:B,}", n, shape);

  /* The plus (1) and minus (2) halves of the layers.  A row-major matpack
     matrix holds the transpose of the Fortran matrix: MMULT's C = A B is
     mult(C, B, A), and y = A x is mult(y, transpose(A), x). */
  const auto r1p = reflect1[0], r1m = reflect1[1], t1p = trans1[0], t1m = trans1[1];
  const auto r2p = reflect2[0], r2m = reflect2[1], t2p = trans2[0], t2m = trans2[1];
  const auto s1p = source1[0], s1m = source1[1], s2p = source2[0], s2m = source2[1];

  // X, Y and GAMMA (COMMON /SCRATCH1/ and /SCRATCH2/, RT3's /RT3_SCRATCH1/
  // and /RT3_SCRATCH2/), and the vectors X and Y take as well: the work
  // data's
  ARTS_USER_ERROR_IF(
      not work.scratch_sized(n), "COMBINE_LAYERS needs a workdata sized for {} streams (workdata::resize)", n);
  Matrix&       x     = work.x;
  Matrix&       y     = work.y;
  Matrix&       gamma = work.gamma;
  Vector&       xv    = work.xv;
  Vector&       yv    = work.yv;
  inv_workdata& wo    = work.inv;

  // GAMMAp = inv[1 - R1p * R2m]     (p for +,  m for -)
  mult(identity(gamma), r2m, r1p, -1.0, 1.0);
  inv_inplace(gamma, wo);

  // RTp = R2p + T2p * GAMMAp * R1p * T2m
  mult(x, t2m, r1p);
  mult(y, x, gamma);
  out_reflect[0] = r2p;
  mult(out_reflect[0], y, t2p, 1.0, 1.0);

  // TTp = T2p * GAMMAp * T1p
  mult(x, t1p, gamma);
  mult(out_trans[0], x, t2p);

  // STp = S2p + T2p * GAMMAp * (S1p + R1p * S2m)
  yv = s1p;
  mult(yv, transpose(r1p), s2m, 1.0, 1.0);
  mult(xv, transpose(gamma), yv);
  out_source[0] = s2p;
  mult(out_source[0], transpose(t2p), xv, 1.0, 1.0);

  // GAMMAm = inv[1 - R2m * R1p]
  mult(identity(gamma), r1p, r2m, -1.0, 1.0);
  inv_inplace(gamma, wo);

  // RTm = R1m + T1m * GAMMAm * R2m * T1p
  mult(x, t1p, r2m);
  mult(y, x, gamma);
  out_reflect[1] = r1m;
  mult(out_reflect[1], y, t1m, 1.0, 1.0);

  // TTm = T1m * GAMMAm * T2m
  mult(x, t2m, gamma);
  mult(out_trans[1], x, t1m);

  // STm = S1m + T1m * GAMMAm * (S2m + R2m * S1p)
  yv = s2m;
  mult(yv, transpose(r2m), s1p, 1.0, 1.0);
  mult(xv, transpose(gamma), yv);
  out_source[1] = s1m;
  mult(out_source[1], transpose(t1m), xv, 1.0, 1.0);
}

void internal_radiance(ConstTensor3View upreflect,
                       ConstTensor3View uptrans,
                       ConstMatrixView  upsource,
                       ConstTensor3View downreflect,
                       ConstTensor3View downtrans,
                       ConstMatrixView  downsource,
                       ConstVectorView  intoprad,
                       ConstVectorView  inbottomrad,
                       VectorView       uprad,
                       VectorView       downrad,
                       workdata&        work) {
  const Index n = upreflect.extent(1);
  for (auto shape : {upreflect.shape(), uptrans.shape(), downreflect.shape(), downtrans.shape()})
    ARTS_USER_ERROR_IF(shape != (std::array<Index, 3>{2, n, n}),
                       "INTERNAL_RADIANCE needs the reflections and transmissions [2, n, n] (n = {}); got {:B,}",
                       n,
                       shape);
  for (auto shape : {upsource.shape(), downsource.shape()})
    ARTS_USER_ERROR_IF(shape != (std::array<Index, 2>{2, n}),
                       "INTERNAL_RADIANCE needs the sources [2, n] (n = {}); got {:B,}",
                       n,
                       shape);
  for (Index size : {intoprad.extent(0), inbottomrad.extent(0), uprad.extent(0), downrad.extent(0)})
    ARTS_USER_ERROR_IF(size != n, "INTERNAL_RADIANCE needs the radiances [n] (n = {}); got {}", n, size);

  /* The plus side of the atmosphere above and the minus side of the one
     below.  A row-major matpack matrix holds the transpose of the Fortran
     matrix: MMULT's C = A B is mult(C, B, A), and y = A x is
     mult(y, transpose(A), x). */
  const auto rup = upreflect[0], tup = uptrans[0];
  const auto rdm = downreflect[1], tdm = downtrans[1];
  const auto sup = upsource[0], sdm = downsource[1];

  // X (COMMON /SCRATCH1/, RT3's /RT3_SCRATCH1/), here the inverse, and the
  // vectors S and V: the work data's gamma, xv and yv
  ARTS_USER_ERROR_IF(
      not work.scratch_sized(n), "INTERNAL_RADIANCE needs a workdata sized for {} streams (workdata::resize)", n);
  Matrix&       x  = work.gamma;
  Vector&       s  = work.xv;
  Vector&       v  = work.yv;
  inv_workdata& wo = work.inv;

  // Compute gamma plus: inv[1 - UPREFLECT(+) DOWNREFLECT(-)]
  mult(identity(x), rdm, rup, -1.0, 1.0);
  inv_inplace(x, wo);
  // Calculate the internal downwelling (plus) radiance vector
  mult(v, transpose(tdm), inbottomrad);
  mult(s, transpose(rup), v);
  mult(s, transpose(tup), intoprad, 1.0, 1.0);
  mult(s, transpose(rup), sdm, 1.0, 1.0);
  s += sup;
  mult(downrad, transpose(x), s);

  // Compute gamma minus: inv[1 - DOWNREFLECT(-) UPREFLECT(+)]
  mult(identity(x), rup, rdm, -1.0, 1.0);
  inv_inplace(x, wo);
  // Calculate the internal upwelling (minus) radiance vector
  mult(v, transpose(tup), intoprad);
  mult(s, transpose(rdm), v);
  mult(s, transpose(tdm), inbottomrad, 1.0, 1.0);
  mult(s, transpose(rdm), sup, 1.0, 1.0);
  s += sdm;
  mult(uprad, transpose(x), s);
}
}  // namespace polradtran
