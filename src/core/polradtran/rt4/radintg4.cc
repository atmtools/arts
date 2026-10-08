#include "radintg4.h"

#include <arts_constants.h>
#include <debug.h>
#include <lin_alg.h>

#include <array>
#include <cmath>

namespace polradtran::rt4 {
void doubling_integration(Index       num_doubles,
                          bool        symmetric,
                          Tensor3View reflect,
                          Tensor3View trans,
                          MatrixView  lin_source,
                          Numeric     linfactor,
                          Tensor3View t_reflect,
                          Tensor3View t_trans,
                          MatrixView  t_source,
                          workdata&   work) {
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
  ARTS_USER_ERROR_IF(
      not work.scratch_sized(n), "DOUBLING_INTEGRATION needs a workdata sized for {} streams (workdata::resize)", n);
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
