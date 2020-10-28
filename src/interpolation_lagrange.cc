#include "interpolation_lagrange.h"

namespace Interpolation {

std::vector<Lagrange> LagrangeVector(const ConstVectorView xs,
                                     const ConstVectorView xi,
                                     const Index polyorder,
                                     const Numeric extrapol) {
  std::vector<Lagrange> out;
  out.reserve(xs.nelem());
  bool has_one=false;
  for (auto x : xs) {
    if (has_one) {
      out.emplace_back(out.back()).update(x, xi, extrapol);
    } else {
      out.push_back(Lagrange(x, xi, polyorder, extrapol));
      has_one = true;
    }
  }
  return out;
}

////////////////////////////////////////////////////////
///////////////////////////////////////////////// Vector
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

const Vector& interpweights(const Lagrange& dim0) { return dim0.lx; }

Grid<Vector, 1> interpweights(const std::vector<Lagrange>& dim0) {
  Grid<Vector, 1> out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++) out(i) = interpweights(dim0[i]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

const Vector& dinterpweights(const Lagrange& dim0) { return dim0.dlx; }

Grid<Vector, 1> dinterpweights(const std::vector<Lagrange>& dim0) {
  Grid<Vector, 1> out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    out(i) = dinterpweights(dim0[i]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstVectorView yi, const ConstVectorView iw,
               const Lagrange& dim0) {
  Numeric out(0.0);
  for (Index i = 0; i < dim0.size(); i++) out += iw[i] * yi[i + dim0.pos];
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Vector reinterp(const ConstVectorView iy, const Grid<Vector, 1>& iw,
                const std::vector<Lagrange>& dim0) {
  Vector out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    out[i] = interp(iy, iw(i), dim0[i]);
  return out;
}

////////////////////////////////////////////////////////
///////////////////////////////////////////////// Matrix
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

Matrix interpweights(const Lagrange& dim0, const Lagrange& dim1) {
  Matrix out(dim0.size(), dim1.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++) out(i, j) = dim0.lx[i] * dim1.lx[j];
  return out;
}

Grid<Matrix, 2> interpweights(const std::vector<Lagrange>& dim0,
                              const std::vector<Lagrange>& dim1) {
  Grid<Matrix, 2> out(dim0.size(), dim1.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = interpweights(dim0[i], dim1[j]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Matrix dinterpweights(const Lagrange& dim0, const Lagrange& dim1, Index dim) {
  Matrix out(dim0.size(), dim1.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      out(i, j) = (dim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                  (dim == 1 ? dim1.dlx[j] : dim1.lx[j]);
  return out;
}

Grid<Matrix, 2> dinterpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1, Index dim) {
  Grid<Matrix, 2> out(dim0.size(), dim1.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = dinterpweights(dim0[i], dim1[j], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstMatrixView yi, const ConstMatrixView iw,
               const Lagrange& dim0, const Lagrange& dim1) {
  Numeric out(0.0);

  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      out += iw(i, j) * yi(i + dim0.pos, j + dim1.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Matrix reinterp(const ConstMatrixView iy, const Grid<Matrix, 2>& iw,
                const std::vector<Lagrange>& dim0,
                const std::vector<Lagrange>& dim1) {
  Matrix out(dim0.size(), dim1.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = interp(iy, iw(i, j), dim0[i], dim1[j]);
  return out;
}

////////////////////////////////////////////////////////
//////////////////////////////////////////////// Tensor3
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

Tensor3 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2) {
  Tensor3 out(dim0.size(), dim1.size(), dim2.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out(i, j, k) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k];
  return out;
}

Grid<Tensor3, 3> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2) {
  Grid<Tensor3, 3> out(dim0.size(), dim1.size(), dim2.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        out(i, j, k) = interpweights(dim0[i], dim1[j], dim2[k]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Tensor3 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, Index dim) {
  Tensor3 out(dim0.size(), dim1.size(), dim2.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out(i, j, k) = (dim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                       (dim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                       (dim == 2 ? dim2.dlx[k] : dim2.lx[k]);
  return out;
}

Grid<Tensor3, 3> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2, Index dim) {
  Grid<Tensor3, 3> out(dim0.size(), dim1.size(), dim2.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        out(i, j, k) = dinterpweights(dim0[i], dim1[j], dim2[k], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstTensor3View yi, const ConstTensor3View iw,
               const Lagrange& dim0, const Lagrange& dim1,
               const Lagrange& dim2) {
  Numeric out(0.0);

  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out += iw(i, j, k) * yi(i + dim0.pos, j + dim1.pos, k + dim2.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Tensor3 reinterp(const ConstTensor3View iy, const Grid<Tensor3, 3>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2) {
  Tensor3 out(dim0.size(), dim1.size(), dim2.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        out(i, j, k) = interp(iy, iw(i, j, k), dim0[i], dim1[j], dim2[k]);
  return out;
}

////////////////////////////////////////////////////////
//////////////////////////////////////////////// Tensor4
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

Tensor4 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3) {
  Tensor4 out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          out(i, j, k, l) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l];
  return out;
}

Grid<Tensor4, 4> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3) {
  Grid<Tensor4, 4> out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          out(i, j, k, l) = interpweights(dim0[i], dim1[j], dim2[k], dim3[l]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Tensor4 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3, Index dim) {
  Tensor4 out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          out(i, j, k, l) = (dim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                            (dim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                            (dim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                            (dim == 3 ? dim3.dlx[l] : dim3.lx[l]);
  return out;
}

Grid<Tensor4, 4> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3, Index dim) {
  Grid<Tensor4, 4> out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          out(i, j, k, l) =
              dinterpweights(dim0[i], dim1[j], dim2[k], dim3[l], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstTensor4View yi, const ConstTensor4View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3) {
  Numeric out(0.0);

  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          out += iw(i, j, k, l) *
                 yi(i + dim0.pos, j + dim1.pos, k + dim2.pos, l + dim3.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Tensor4 reinterp(const ConstTensor4View iy, const Grid<Tensor4, 4>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3) {
  Tensor4 out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          out(i, j, k, l) =
              interp(iy, iw(i, j, k, l), dim0[i], dim1[j], dim2[k], dim3[l]);
  return out;
}

////////////////////////////////////////////////////////
//////////////////////////////////////////////// Tensor5
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

Tensor5 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3,
                      const Lagrange& dim4) {
  Tensor5 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) =
                dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l] * dim4.lx[m];
  return out;
}

Grid<Tensor5, 5> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3,
                               const std::vector<Lagrange>& dim4) {
  Grid<Tensor5, 5> out(dim0.size(), dim1.size(), dim2.size(), dim3.size(),
                       dim4.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) =
                interpweights(dim0[i], dim1[j], dim2[k], dim3[l], dim4[m]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Tensor5 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3,
                       const Lagrange& dim4, Index dim) {
  Tensor5 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) = (dim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                                 (dim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                                 (dim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                                 (dim == 3 ? dim3.dlx[l] : dim3.lx[l]) *
                                 (dim == 4 ? dim4.dlx[m] : dim4.lx[m]);
  return out;
}

Grid<Tensor5, 5> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3,
                                const std::vector<Lagrange>& dim4, Index dim) {
  Grid<Tensor5, 5> out(dim0.size(), dim1.size(), dim2.size(), dim3.size(),
                       dim4.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) = dinterpweights(dim0[i], dim1[j], dim2[k],
                                                dim3[l], dim4[m], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstTensor5View yi, const ConstTensor5View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4) {
  Numeric out(0.0);
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            out +=
                iw(i, j, k, l, m) * yi(i + dim0.pos, j + dim1.pos, k + dim2.pos,
                                       l + dim3.pos, m + dim4.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Tensor5 reinterp(const ConstTensor5View iy, const Grid<Tensor5, 5>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4) {
  Tensor5 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) = interp(iy, iw(i, j, k, l, m), dim0[i], dim1[j],
                                        dim2[k], dim3[l], dim4[m]);
  return out;
}

////////////////////////////////////////////////////////
//////////////////////////////////////////////// Tensor6
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

Tensor6 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3,
                      const Lagrange& dim4, const Lagrange& dim5) {
  Tensor6 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
              dim5.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] *
                                      dim3.lx[l] * dim4.lx[m] * dim5.lx[n];
  return out;
}

Grid<Tensor6, 6> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3,
                               const std::vector<Lagrange>& dim4,
                               const std::vector<Lagrange>& dim5) {
  Grid<Tensor6, 6> out(dim0.size(), dim1.size(), dim2.size(), dim3.size(),
                       dim4.size(), dim5.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) = interpweights(dim0[i], dim1[j], dim2[k],
                                                    dim3[l], dim4[m], dim5[n]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Tensor6 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3,
                       const Lagrange& dim4, const Lagrange& dim5, Index dim) {
  Tensor6 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
              dim5.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) = (dim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                                      (dim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                                      (dim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                                      (dim == 3 ? dim3.dlx[l] : dim3.lx[l]) *
                                      (dim == 4 ? dim4.dlx[m] : dim4.lx[m]) *
                                      (dim == 5 ? dim5.dlx[n] : dim5.lx[n]);
  return out;
}

Grid<Tensor6, 6> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3,
                                const std::vector<Lagrange>& dim4,
                                const std::vector<Lagrange>& dim5, Index dim) {
  Grid<Tensor6, 6> out(dim0.size(), dim1.size(), dim2.size(), dim3.size(),
                       dim4.size(), dim5.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) = dinterpweights(
                  dim0[i], dim1[j], dim2[k], dim3[l], dim4[m], dim5[n], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstTensor6View yi, const ConstTensor6View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4,
               const Lagrange& dim5) {
  Numeric out(0.0);
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              out += iw(i, j, k, l, m, n) * yi(i + dim0.pos, j + dim1.pos,
                                               k + dim2.pos, l + dim3.pos,
                                               m + dim4.pos, n + dim5.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Tensor6 reinterp(const ConstTensor6View iy, const Grid<Tensor6, 6>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4,
                 const std::vector<Lagrange>& dim5) {
  Tensor6 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
              dim5.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) =
                  interp(iy, iw(i, j, k, l, m, n), dim0[i], dim1[j], dim2[k],
                         dim3[l], dim4[m], dim5[n]);
  return out;
}

////////////////////////////////////////////////////////
//////////////////////////////////////////////// Tensor7
////////////////////////////////////////////////////////

////////////////////////////////////////////////
////////////////////////// Interpolation Weights
////////////////////////////////////////////////

Tensor7 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3,
                      const Lagrange& dim4, const Lagrange& dim5,
                      const Lagrange& dim6) {
  Tensor7 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
              dim5.size(), dim6.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              for (Index o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) = dim0.lx[i] * dim1.lx[j] *
                                           dim2.lx[k] * dim3.lx[l] *
                                           dim4.lx[m] * dim5.lx[n] * dim6.lx[o];
  return out;
}

Grid<Tensor7, 7> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3,
                               const std::vector<Lagrange>& dim4,
                               const std::vector<Lagrange>& dim5,
                               const std::vector<Lagrange>& dim6) {
  Grid<Tensor7, 7> out(dim0.size(), dim1.size(), dim2.size(), dim3.size(),
                       dim4.size(), dim5.size(), dim6.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              for (std::size_t o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) =
                    interpweights(dim0[i], dim1[j], dim2[k], dim3[l], dim4[m],
                                  dim5[n], dim6[o]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Tensor7 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3,
                       const Lagrange& dim4, const Lagrange& dim5,
                       const Lagrange& dim6, Index dim) {
  Tensor7 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
              dim5.size(), dim6.size());
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              for (Index o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) =
                    (dim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                    (dim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                    (dim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                    (dim == 3 ? dim3.dlx[l] : dim3.lx[l]) *
                    (dim == 4 ? dim4.dlx[m] : dim4.lx[m]) *
                    (dim == 5 ? dim5.dlx[n] : dim5.lx[n]) *
                    (dim == 6 ? dim6.dlx[o] : dim6.lx[o]);
  return out;
}

Grid<Tensor7, 7> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3,
                                const std::vector<Lagrange>& dim4,
                                const std::vector<Lagrange>& dim5,
                                const std::vector<Lagrange>& dim6, Index dim) {
  Grid<Tensor7, 7> out(dim0.size(), dim1.size(), dim2.size(), dim3.size(),
                       dim4.size(), dim5.size(), dim6.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              for (std::size_t o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) =
                    dinterpweights(dim0[i], dim1[j], dim2[k], dim3[l], dim4[m],
                                   dim5[n], dim6[o], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstTensor7View yi, const ConstTensor7View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5,
               const Lagrange& dim6) {
  Numeric out(0.0);
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              for (Index o = 0; o < dim6.size(); o++)
                out += iw(i, j, k, l, m, n, o) * yi(i + dim0.pos, j + dim1.pos,
                                                    k + dim2.pos, l + dim3.pos,
                                                    m + dim4.pos, n + dim5.pos,
                                                    o + dim6.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

Tensor7 reinterp(const ConstTensor7View iy, const Grid<Tensor7, 7>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4,
                 const std::vector<Lagrange>& dim5,
                 const std::vector<Lagrange>& dim6) {
  Tensor7 out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
              dim5.size(), dim6.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              for (std::size_t o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) =
                    interp(iy, iw(i, j, k, l, m, n, o), dim0[i], dim1[j],
                           dim2[k], dim3[l], dim4[m], dim5[n], dim6[o]);
  return out;
}
}  // namespace Interpolation
