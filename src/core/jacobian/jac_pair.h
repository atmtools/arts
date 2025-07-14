#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurface.h>
#include <surf.h>

#include "jacobian.h"

template <class T>
struct pairfwd {
  using VectorFunc = std::function<Vector(ConstVectorView, const T&)>;

  VectorFunc vfunc_1;
  VectorFunc vfunc_2;

  Vector operator()(ConstVectorView x, const T& y) const {
    return vfunc_2(vfunc_1(x, y), y);
  }
};

template <class T>
struct pairinv {
  using VectorFunc = std::function<Vector(ConstVectorView, const T&)>;
  using MatrixFunc =
      std::function<Matrix(ConstMatrixView, ConstVectorView, const T&)>;

  VectorFunc vfunc_1;
  VectorFunc vfunc_2;

  MatrixFunc mfunc_1;
  MatrixFunc mfunc_2;

  Vector operator()(ConstVectorView x, const T& y) const {
    return vfunc_2(vfunc_1(x, y), y);
  }

  Matrix operator()(ConstMatrixView dy, ConstVectorView x, const T& y) const {
    return mfunc_2(mfunc_1(dy, x, y), x, y);
  }
};

void make_logrelfit(Jacobian::AtmTarget&, const AtmField&);
void make_logrelfit(Jacobian::SurfaceTarget&, const SurfaceField&);
void make_logrelfit(Jacobian::SubsurfaceTarget&, const SubsurfaceField&);
