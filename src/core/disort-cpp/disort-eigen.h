#pragma once

#include <matpack.h>

struct real_diagonalize_workdata {
#if ARTS_LGPL
  Vector imag;
  diagonalize_workdata workdata;

  real_diagonalize_workdata(Size N) : imag(N), workdata(N) {}
#else
  Vector reals;
  std::vector<Index> ints;

  real_diagonalize_workdata(Size N) : reals(2 * N), ints(N) {}
#endif

  real_diagonalize_workdata()                                     = default;
  real_diagonalize_workdata(const real_diagonalize_workdata&)     = default;
  real_diagonalize_workdata(real_diagonalize_workdata&&) noexcept = default;
  real_diagonalize_workdata& operator=(const real_diagonalize_workdata&) =
      default;
  real_diagonalize_workdata& operator=(real_diagonalize_workdata&&) noexcept =
      default;
};

Index diagonalize_inplace(MatrixView P,
                          VectorView W,
                          MatrixView A,
                          real_diagonalize_workdata& workdata);
