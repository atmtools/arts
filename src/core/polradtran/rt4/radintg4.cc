#include "radintg4.h"

#include <arts_constants.h>
#include <debug.h>
#include <lin_alg.h>

#include <array>
#include <cmath>

namespace polradtran::rt4 {
void doubling_integration(Index         num_doubles,
                          bool          symmetric,
                          Tensor3View   reflect,
                          Tensor3View   trans,
                          MatrixView    lin_source,
                          Numeric       linfactor,
                          Tensor3View   t_reflect,
                          Tensor3View   t_trans,
                          MatrixView    t_source,
                          rt4_workdata& work) {
  const Index n = reflect.extent(1);
  ARTS_USER_ERROR_IF(
      reflect.shape() != (std::array<Index, 3>{2, n, n}) or trans.shape() != (std::array<Index, 3>{2, n, n}) or
          lin_source.shape() != (std::array<Index, 2>{2, n}) or t_reflect.shape() != (std::array<Index, 3>{2, n, n}) or
          t_trans.shape() != (std::array<Index, 3>{2, n, n}) or t_source.shape() != (std::array<Index, 2>{2, n}),
      "DOUBLING_INTEGRATION needs reflect, trans, t_reflect and t_trans [2, n, n] and lin_source and "
      "t_source [2, n]; got {:B,}, {:B,}, {:B,}, {:B,}, {:B,} and {:B,}",
      reflect.shape(),
      trans.shape(),
      t_reflect.shape(),
      t_trans.shape(),
      lin_source.shape(),
      t_source.shape());

  /* The plus (1) and minus (2) halves of REFLECT, TRANS and LIN_SOURCE,
     which are updated in place.  The n x n matrices are the Fortran's
     column-major arrays, so a row-major matpack matrix holds the transpose
     of the Fortran matrix: MMULT's C = A B (DGEMM) is mult(C, B, A), and
     the matrix-vector y = A x is mult(y, transpose(A), x) (DGEMV).  The
     MIDENTITY, MSUB and MADD around a product are DGEMM's and DGEMV's
     alpha and beta. */
  const auto rp = reflect[0], rm = reflect[1];
  const auto tp = trans[0], tm = trans[1];
  const auto sp = lin_source[0], sm = lin_source[1];

  // X, Y and GAMMA (COMMON /SCRATCH1/ and /SCRATCH2/), the vectors X and Y
  // take as well, T_LIN, CONST and T_CONST: the work data's
  ARTS_USER_ERROR_IF(not work.scratch_sized(n),
                     "DOUBLING_INTEGRATION needs an rt4_workdata sized for {} streams (rt4_workdata::resize)",
                     n);
  Matrix&       x       = work.x;
  Matrix&       y       = work.y;
  Matrix&       gamma   = work.gamma;
  Vector&       xv      = work.xv;
  Vector&       yv      = work.yv;
  Matrix&       t_lin   = work.t_lin;
  Matrix&       t_const = work.t_const;
  Matrix&       cnst    = work.cnst;
  inv_workdata& wo      = work.inv;

  Numeric linfac = linfactor;
  cnst           = lin_source;
  const auto cp = cnst[0], cm = cnst[1];

  for (Index i = 0; i < num_doubles; i++) {
    // Make gamma plus matrix: GAMMA = inv[1 - Rp*Rm]
    mult(identity(gamma), rm, rp, -1.0, 1.0);
    inv_inplace(gamma, wo);

    // Rp(2N) = Rp + Tp * GAMMA * Rp * Tm
    mult(x, tm, rp);
    mult(y, x, gamma);
    t_reflect[0] = rp;
    mult(t_reflect[0], y, tp, 1.0, 1.0);

    // Tp(2N) = Tp * GAMMA * Tp
    mult(x, tp, gamma);
    mult(t_trans[0], x, tp);

    // Linear source doubling
    //   Sp(2N) = (Sp+f*Cp) + Tp * GAMMA * (Sp + Rp * (Sm+f*Cm))
    xv  = cm;
    xv *= linfac;
    yv  = sm;
    yv += xv;
    xv  = sp;
    mult(xv, transpose(rp), yv, 1.0, 1.0);
    mult(yv, transpose(gamma), xv);
    t_lin[0] = sp;
    mult(t_lin[0], transpose(tp), yv, 1.0, 1.0);
    yv        = cp;
    yv       *= linfac;
    t_lin[0] += yv;
    //   Cp(2N) = Cp + Tp * GAMMA * (Cp + Rp * Cm)
    yv = cp;
    mult(yv, transpose(rp), cm, 1.0, 1.0);
    mult(xv, transpose(gamma), yv);
    t_const[0] = cp;
    mult(t_const[0], transpose(tp), xv, 1.0, 1.0);

    if (symmetric) {
      t_reflect[1] = t_reflect[0];
      t_trans[1]   = t_trans[0];
    } else {
      // Make gamma minus matrix: GAMMA = inv[1 - Rm*Rp]
      mult(identity(gamma), rp, rm, -1.0, 1.0);
      inv_inplace(gamma, wo);

      // Rm(2N) = Rm + Tm * GAMMA * Rm * Tp
      mult(x, tp, rm);
      mult(y, x, gamma);
      t_reflect[1] = rm;
      mult(t_reflect[1], y, tm, 1.0, 1.0);

      // Tm(2N) = Tm * GAMMA * Tm
      mult(x, tm, gamma);
      mult(t_trans[1], x, tm);
    }

    // Linear source doubling
    //   Sm(2N) = Sm + Tm * GAMMA * (Sm+f*Cm + Rm * Sp)
    xv  = cm;
    xv *= linfac;
    yv  = sm;
    yv += xv;
    mult(yv, transpose(rm), sp, 1.0, 1.0);
    mult(xv, transpose(gamma), yv);
    t_lin[1] = sm;
    mult(t_lin[1], transpose(tm), xv, 1.0, 1.0);
    //   Cm(2N) = Cm + Tm * GAMMA * (Cm + Rm * Cp)
    yv = cm;
    mult(yv, transpose(rm), cp, 1.0, 1.0);
    mult(xv, transpose(gamma), yv);
    t_const[1] = cm;
    mult(t_const[1], transpose(tm), xv, 1.0, 1.0);

    lin_source  = t_lin;
    cnst        = t_const;
    linfac     *= 2.0;

    reflect = t_reflect;
    trans   = t_trans;
  }

  if (num_doubles <= 0) {
    t_reflect = reflect;
    t_trans   = trans;
  }

  t_source = lin_source;
}

void combine_layers(ConstTensor3View reflect1,
                    ConstTensor3View trans1,
                    ConstMatrixView  source1,
                    ConstTensor3View reflect2,
                    ConstTensor3View trans2,
                    ConstMatrixView  source2,
                    Tensor3View      out_reflect,
                    Tensor3View      out_trans,
                    MatrixView       out_source,
                    rt4_workdata&    work) {
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

  /* The plus (1) and minus (2) halves of the layers.  As in
     DOUBLING_INTEGRATION, a row-major matpack matrix holds the transpose of
     the Fortran matrix: MMULT's C = A B is mult(C, B, A), and y = A x is
     mult(y, transpose(A), x). */
  const auto r1p = reflect1[0], r1m = reflect1[1], t1p = trans1[0], t1m = trans1[1];
  const auto r2p = reflect2[0], r2m = reflect2[1], t2p = trans2[0], t2m = trans2[1];
  const auto s1p = source1[0], s1m = source1[1], s2p = source2[0], s2m = source2[1];

  // X, Y and GAMMA (COMMON /SCRATCH1/ and /SCRATCH2/), and the vectors X
  // and Y take as well: the work data's
  ARTS_USER_ERROR_IF(
      not work.scratch_sized(n), "COMBINE_LAYERS needs an rt4_workdata sized for {} streams (rt4_workdata::resize)", n);
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
                       rt4_workdata&    work) {
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
     below.  As in DOUBLING_INTEGRATION, a row-major matpack matrix holds
     the transpose of the Fortran matrix: MMULT's C = A B is mult(C, B, A),
     and y = A x is mult(y, transpose(A), x). */
  const auto rup = upreflect[0], tup = uptrans[0];
  const auto rdm = downreflect[1], tdm = downtrans[1];
  const auto sup = upsource[0], sdm = downsource[1];

  // X (COMMON /SCRATCH1/), here the inverse, and the vectors S and V: the
  // work data's gamma, xv and yv
  ARTS_USER_ERROR_IF(not work.scratch_sized(n),
                     "INTERNAL_RADIANCE needs an rt4_workdata sized for {} streams (rt4_workdata::resize)",
                     n);
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

void initialize(Numeric          delta_z,
                ConstVectorView  mu_values,
                ConstVectorView  quad_weights,
                Numeric          gas_extinct,
                ConstTensor4View extinct_matrix,
                ConstTensor5View scatter_matrix,
                Tensor5View      reflect,
                Tensor5View      trans) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = reflect.extent(4);
  ARTS_USER_ERROR_IF(quad_weights.extent(0) != nummu or
                         extinct_matrix.shape() != (std::array<Index, 4>{2, nummu, nstokes, nstokes}) or
                         scatter_matrix.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}) or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}),
                     "INITIALIZE with {} mu_values needs quad_weights [nummu], extinct_matrix [2, nummu, nstokes, "
                     "nstokes], scatter_matrix [4, nummu, nstokes, nummu, nstokes] and reflect and trans [2, nummu, "
                     "nstokes, nummu, nstokes]; got {}, {:B,}, {:B,}, {:B,} and {:B,}",
                     nummu,
                     quad_weights.extent(0),
                     extinct_matrix.shape(),
                     scatter_matrix.shape(),
                     reflect.shape(),
                     trans.shape());

  const Numeric c = Constant::two_pi;

  for (Index i2 = 0; i2 < nstokes; i2++) {
    for (Index j2 = 0; j2 < nummu; j2++) {
      const Numeric tmp = delta_z / mu_values[j2];
      for (Index i1 = 0; i1 < nstokes; i1++) {
        Numeric gext = 0.0;
        if (i1 == i2) gext = gas_extinct;
        for (Index j1 = 0; j1 < nummu; j1++) {
          reflect[0, j1, i1, j2, i2] = c * tmp * quad_weights[j1] * scatter_matrix[1, j1, i1, j2, i2];
          reflect[1, j1, i1, j2, i2] = c * tmp * quad_weights[j1] * scatter_matrix[2, j1, i1, j2, i2];
          Numeric diag               = 0.0;
          if (i1 == i2 and j1 == j2) diag = 1.0;
          Numeric ext = 0.0;
          if (j1 == j2) ext = extinct_matrix[0, j2, i1, i2] + gext;
          trans[0, j1, i1, j2, i2] = diag - tmp * (ext - c * quad_weights[j1] * scatter_matrix[0, j1, i1, j2, i2]);
          if (j1 == j2) ext = extinct_matrix[1, j2, i1, i2] + gext;
          trans[1, j1, i1, j2, i2] = diag - tmp * (ext - c * quad_weights[j1] * scatter_matrix[3, j1, i1, j2, i2]);
        }
      }
    }
  }
}

void initial_source(Numeric          delta_z,
                    ConstVectorView  mu_values,
                    Numeric          planck,
                    ConstTensor3View emis_vector,
                    Numeric          gas_extinct,
                    Tensor3View      source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(emis_vector.shape() != (std::array<Index, 3>{2, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "INITIAL_SOURCE with {} mu_values needs emis_vector and source [2, nummu, nstokes]; got "
                     "{:B,} and {:B,}",
                     nummu,
                     emis_vector.shape(),
                     source.shape());

  source = 0.0;

  for (Index i = 0; i < nstokes; i++) {
    for (Index j = 0; j < nummu; j++) {
      Numeric ext = 0.0;
      if (i == 0) ext = gas_extinct;
      const Numeric tmp = planck * delta_z / mu_values[j];
      source[0, j, i]   = tmp * (emis_vector[0, j, i] + ext);
      source[1, j, i]   = tmp * (emis_vector[1, j, i] + ext);
    }
  }
}

void nonscatter_layer(Numeric         deltatau,
                      ConstVectorView mu_values,
                      Numeric         planck0,
                      Numeric         planck1,
                      Tensor5View     reflect,
                      Tensor5View     trans,
                      Tensor3View     source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "NONSCATTER_LAYER with {} mu_values needs reflect and trans [2, nummu, nstokes, nummu, "
                     "nstokes] and source [2, nummu, nstokes]; got {:B,}, {:B,} and {:B,}",
                     nummu,
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;

  trans = 0.0;
  for (Index j = 0; j < nummu; j++) {
    const Numeric factor = std::exp(-deltatau / mu_values[j]);
    identity(trans[0, j, joker, j, joker], factor);
    identity(trans[1, j, joker, j, joker], factor);
  }

  source = 0.0;
  if (deltatau > 0.0) {
    for (Index j = 0; j < nummu; j++) {
      const Numeric path  = deltatau / mu_values[j];
      const Numeric slope = (planck1 - planck0) / path;
      source[0, j, 0]     = planck1 - slope - (planck1 - slope * (1.0 + path)) * std::exp(-path);
      source[1, j, 0]     = planck0 + slope - (planck0 + slope * (1.0 + path)) * std::exp(-path);
    }
  }
}
}  // namespace polradtran::rt4
