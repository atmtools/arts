#include "radintg3.h"

#include <debug.h>
#include <lin_alg.h>

#include <array>
#include <cmath>

namespace rt3 {
void doubling_integration(Index       num_doubles,
                          Index       src_code,
                          bool        symmetric,
                          Tensor3View reflect,
                          Tensor3View trans,
                          MatrixView  exp_source,
                          Numeric     expfactor,
                          MatrixView  lin_source,
                          Numeric     linfactor,
                          Tensor3View t_reflect,
                          Tensor3View t_trans,
                          MatrixView  t_source) {
  const Index n = reflect.extent(1);
  ARTS_USER_ERROR_IF(
      n < 1 or reflect.shape() != (std::array<Index, 3>{2, n, n}) or trans.shape() != reflect.shape() or
          t_reflect.shape() != reflect.shape() or t_trans.shape() != reflect.shape() or
          exp_source.shape() != (std::array<Index, 2>{2, n}) or lin_source.shape() != exp_source.shape() or
          t_source.shape() != exp_source.shape(),
      "DOUBLING_INTEGRATION needs reflect, trans, t_reflect and t_trans [2, n, n] and exp_source, lin_source and "
      "t_source [2, n]; got {:B,}, {:B,}, {:B,}, {:B,}, {:B,}, {:B,} and {:B,}",
      reflect.shape(),
      trans.shape(),
      t_reflect.shape(),
      t_trans.shape(),
      exp_source.shape(),
      lin_source.shape(),
      t_source.shape());
  ARTS_USER_ERROR_IF(src_code < 0 or src_code > 3, "DOUBLING_INTEGRATION needs src_code 0 to 3, got {}", src_code);
  const bool solar   = src_code == 1 or src_code == 3;
  const bool thermal = src_code >= 2;

  /* The plus (1) and minus (2) halves of REFLECT, TRANS, EXP_SOURCE and
     LIN_SOURCE, which are updated in place.  The n x n matrices are the
     Fortran's column-major arrays, so a row-major matpack matrix holds the
     transpose of the Fortran matrix: MMULT's C = A B (DGEMM) is
     mult(C, B, A), and the matrix-vector y = A x is mult(y, transpose(A), x)
     (DGEMV).  The MIDENTITY, MSUB and MADD around a product are DGEMM's and
     DGEMV's alpha and beta. */
  const auto rp = reflect[0], rm = reflect[1];
  const auto tp = trans[0], tm = trans[1];
  const auto ep = exp_source[0], em = exp_source[1];
  const auto sp = lin_source[0], sm = lin_source[1];

  // X, Y and GAMMA (COMMON /RT3_SCRATCH1/ and /RT3_SCRATCH2/), the vectors
  // X and Y take as well, T_EXP, T_LIN, CONST and T_CONST
  Matrix     x(n, n), y(n, n), gamma(n, n);
  Vector     xv(n), yv(n);
  Matrix     t_exp(2, n), t_lin(2, n), t_const(2, n);
  Matrix     cnst{lin_source};
  const auto cp = cnst[0], cm = cnst[1];

  Numeric expfac = expfactor;
  Numeric linfac = linfactor;
  for (Index i = 0; i < num_doubles; i++) {
    // Make gamma plus matrix: GAMMA = inv[1 - Rp*Rm]
    mult(identity(gamma), rm, rp, -1.0, 1.0);
    inv_inplace(gamma);

    // Rp(2N) = Rp + Tp * GAMMA * Rp * Tm
    mult(x, tm, rp);
    mult(y, x, gamma);
    t_reflect[0] = rp;
    mult(t_reflect[0], y, tp, 1.0, 1.0);

    // Tp(2N) = Tp * GAMMA * Tp
    mult(x, tp, gamma);
    mult(t_trans[0], x, tp);

    // Exponential source doubling
    if (solar) {
      //   Sp(2N) = e*Sp + Tp * GAMMA * (Sp + Rp * e*Sm)
      yv  = em;
      yv *= expfac;
      xv  = ep;
      mult(xv, transpose(rp), yv, 1.0, 1.0);
      mult(yv, transpose(gamma), xv);
      t_exp[0]  = ep;
      t_exp[0] *= expfac;
      mult(t_exp[0], transpose(tp), yv, 1.0, 1.0);
    }

    // Linear source doubling
    if (thermal) {
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
    }

    if (symmetric) {
      t_reflect[1] = t_reflect[0];
      t_trans[1]   = t_trans[0];
    } else {
      // Make gamma minus matrix: GAMMA = inv[1 - Rm*Rp]
      mult(identity(gamma), rp, rm, -1.0, 1.0);
      inv_inplace(gamma);

      // Rm(2N) = Rm + Tm * GAMMA * Rm * Tp
      mult(x, tp, rm);
      mult(y, x, gamma);
      t_reflect[1] = rm;
      mult(t_reflect[1], y, tm, 1.0, 1.0);

      // Tm(2N) = Tm * GAMMA * Tm
      mult(x, tm, gamma);
      mult(t_trans[1], x, tm);
    }

    // Exponential source doubling
    if (solar) {
      //   Sm(2N) = Sm + Tm * GAMMA * (e*Sm + Rm * Sp)
      yv  = em;
      yv *= expfac;
      mult(yv, transpose(rm), ep, 1.0, 1.0);
      mult(xv, transpose(gamma), yv);
      t_exp[1] = em;
      mult(t_exp[1], transpose(tm), xv, 1.0, 1.0);

      exp_source  = t_exp;
      expfac     *= expfac;
    }

    // Linear source doubling
    if (thermal) {
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
    }

    reflect = t_reflect;
    trans   = t_trans;
  }

  if (num_doubles <= 0) {
    t_reflect = reflect;
    t_trans   = trans;
  }

  if (src_code == 3) {
    t_source  = exp_source;
    t_source += lin_source;
  } else if (src_code == 2) {
    t_source = lin_source;
  } else if (src_code == 1) {
    t_source = exp_source;
  } else {
    t_source = 0.0;
  }
}

void initialize(Numeric          delta_z,
                ConstVectorView  mu_values,
                Numeric          extinction,
                Numeric          albedo,
                ConstTensor5View phase_function,
                Tensor5View      reflect,
                Tensor5View      trans) {
  const Index nummu   = mu_values.size();
  const Index nstokes = reflect.extent(4);
  ARTS_USER_ERROR_IF(nummu < 1 or nstokes < 1 or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != reflect.shape() or
                         phase_function.shape() != (std::array<Index, 5>{4, nummu, nstokes, nummu, nstokes}),
                     "INITIALIZE with {} mu_values needs phase_function [4, nummu, nstokes, nummu, nstokes] and "
                     "reflect and trans [2, nummu, nstokes, nummu, nstokes]; got {:B,}, {:B,} and {:B,}",
                     nummu,
                     phase_function.shape(),
                     reflect.shape(),
                     trans.shape());

  // The rows of angle j1, [l, joker, joker, j1, joker]:
  // REFLECT = f ALBEDO PHASE_FUNCTION(2 or 3) and
  // TRANS = DIAG - f (DIAG - ALBEDO PHASE_FUNCTION(1 or 4)).  The diagonal
  // of TRANS is kept in that form, 1 minus a small number rounded once: it
  // carries the thin layer's extinction, which the doubling amplifies.
  for (Index j1 = 0; j1 < nummu; j1++) {
    const Numeric f = delta_z / mu_values[j1] * extinction;
    for (Index l = 0; l < 2; l++) {
      reflect[l, joker, joker, j1, joker]  = phase_function[l + 1, joker, joker, j1, joker];
      reflect[l, joker, joker, j1, joker] *= f * albedo;
      trans[l, joker, joker, j1, joker]    = phase_function[3 * l, joker, joker, j1, joker];
      trans[l, joker, joker, j1, joker]   *= f * albedo;
      for (Index i = 0; i < nstokes; i++)
        trans[l, j1, i, j1, i] = 1.0 - f * (1.0 - albedo * phase_function[3 * l, j1, i, j1, i]);
    }
  }
}

void initial_source(Numeric          delta_z,
                    ConstVectorView  mu_values,
                    Numeric          extinction,
                    ConstTensor3View source_vector,
                    Tensor3View      source) {
  const Index nummu   = mu_values.size();
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(nummu < 1 or nstokes < 1 or source.shape() != (std::array<Index, 3>{2, nummu, nstokes}) or
                         source_vector.shape() != source.shape(),
                     "INITIAL_SOURCE with {} mu_values needs source_vector and source [2, nummu, nstokes]; got "
                     "{:B,} and {:B,}",
                     nummu,
                     source_vector.shape(),
                     source.shape());

  for (Index j = 0; j < nummu; j++) {
    const Numeric tmp        = delta_z / mu_values[j];
    source[joker, j, joker]  = source_vector[joker, j, joker];
    source[joker, j, joker] *= tmp * extinction;
  }
}

void nonscatter_layer(Index           mode,
                      Numeric         deltatau,
                      ConstVectorView mu_values,
                      Numeric         planck0,
                      Numeric         planck1,
                      Tensor5View     reflect,
                      Tensor5View     trans,
                      Tensor3View     source) {
  const Index nummu   = mu_values.size();
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(
      nummu < 1 or nstokes < 1 or reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
          trans.shape() != reflect.shape() or source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
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
  if (mode == 0 and deltatau > 0.0) {
    for (Index j = 0; j < nummu; j++) {
      const Numeric path  = deltatau / mu_values[j];
      const Numeric slope = (planck1 - planck0) / path;
      source[0, j, 0]     = planck1 - slope - (planck1 - slope * (1.0 + path)) * std::exp(-path);
      source[1, j, 0]     = planck0 + slope - (planck0 + slope * (1.0 + path)) * std::exp(-path);
    }
  }
}

void combine_layers(ConstTensor3View reflect1,
                    ConstTensor3View trans1,
                    ConstMatrixView  source1,
                    ConstTensor3View reflect2,
                    ConstTensor3View trans2,
                    ConstMatrixView  source2,
                    Tensor3View      out_reflect,
                    Tensor3View      out_trans,
                    MatrixView       out_source) {
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
     mult(y, transpose(A), x); the MIDENTITY, MSUB and MADD around a
     product are its alpha and beta. */
  const auto r1p = reflect1[0], r1m = reflect1[1], t1p = trans1[0], t1m = trans1[1];
  const auto r2p = reflect2[0], r2m = reflect2[1], t2p = trans2[0], t2m = trans2[1];
  const auto s1p = source1[0], s1m = source1[1], s2p = source2[0], s2m = source2[1];

  // X, Y and GAMMA (COMMON /RT3_SCRATCH1/ and /RT3_SCRATCH2/), and the
  // vectors X and Y take as well
  Matrix x(n, n), y(n, n), gamma(n, n);
  Vector xv(n), yv(n);

  // GAMMAp = inv[1 - R1p * R2m]     (p for +,  m for -)
  mult(identity(gamma), r2m, r1p, -1.0, 1.0);
  inv_inplace(gamma);

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
  inv_inplace(gamma);

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
                       VectorView       downrad) {
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
     and y = A x is mult(y, transpose(A), x); the MIDENTITY, MSUB and MADD
     around a product are its alpha and beta. */
  const auto rup = upreflect[0], tup = uptrans[0];
  const auto rdm = downreflect[1], tdm = downtrans[1];
  const auto sup = upsource[0], sdm = downsource[1];

  // X (COMMON /RT3_SCRATCH1/), here the inverse, and the vectors S and V
  Matrix x(n, n);
  Vector s(n), v(n);

  // Compute gamma plus: inv[1 - UPREFLECT(+) DOWNREFLECT(-)]
  mult(identity(x), rdm, rup, -1.0, 1.0);
  inv_inplace(x);
  // Calculate the internal downwelling (plus) radiance vector
  mult(v, transpose(tdm), inbottomrad);
  mult(s, transpose(rup), v);
  mult(s, transpose(tup), intoprad, 1.0, 1.0);
  mult(s, transpose(rup), sdm, 1.0, 1.0);
  s += sup;
  mult(downrad, transpose(x), s);

  // Compute gamma minus: inv[1 - DOWNREFLECT(-) UPREFLECT(+)]
  mult(identity(x), rup, rdm, -1.0, 1.0);
  inv_inplace(x);
  // Calculate the internal upwelling (minus) radiance vector
  mult(v, transpose(tup), intoprad);
  mult(s, transpose(rdm), v);
  mult(s, transpose(tdm), inbottomrad, 1.0, 1.0);
  mult(s, transpose(rdm), sup, 1.0, 1.0);
  s += sdm;
  mult(uprad, transpose(x), s);
}
}  // namespace rt3
