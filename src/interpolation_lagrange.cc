#include "interpolation_lagrange.h"
#include "matpack_math.h"

namespace Interpolation {

Array<Lagrange> LagrangeVector(
    const ConstVectorView& xs, const ConstVectorView& xi, const Index polyorder,
    const Numeric extrapol, const bool do_derivs, const GridType type, const std::pair<Numeric, Numeric> cycle) {
  if (xs.size() == 1) {
    check_lagrange_interpolation(xi, polyorder, xs[0], extrapol, type, cycle);
  } else {
    check_lagrange_interpolation(xi, polyorder, {min(xs), max(xs)}, extrapol, type, cycle);
  }
  Array<Lagrange> out;
  out.reserve(xs.size());
  bool has_one = false;
  for (auto x : xs) {
    if (has_one) {
      out.emplace_back(out.back().pos, x, xi, polyorder, do_derivs, type, cycle);
    } else {
      out.emplace_back(start_pos_finder(x, xi), x, xi, polyorder, do_derivs, type, cycle);
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

Vector interpweights(const Lagrange& dim0) { return dim0.lx; }

Grid<Vector, 1> interpweights(const Array<Lagrange>& dim0) {
  Grid<Vector, 1> out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++) out(i) = interpweights(dim0[i]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

Vector dinterpweights(const Lagrange& dim0) { return dim0.dlx; }

Grid<Vector, 1> dinterpweights(const Array<Lagrange>& dim0) {
  Grid<Vector, 1> out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    out(i) = dinterpweights(dim0[i]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstVectorView& yi, const ConstVectorView& iw,
               const Lagrange& dim0) {
  Numeric out(0.0);
  const Index I = yi.size();
  for (Index i = 0; i < dim0.size(); i++)
    out += iw[i] * yi[cycler(i + dim0.pos, I)];
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(VectorView out, const ConstVectorView& iy,
              const Grid<Vector, 1>& iw, const Array<Lagrange>& dim0) {
  for (std::size_t i = 0; i < dim0.size(); i++)
    out[i] = interp(iy, iw(i), dim0[i]);
}

Vector reinterp(const ConstVectorView& iy, const Grid<Vector, 1>& iw,
                const Array<Lagrange>& dim0) {
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

Grid<Matrix, 2> interpweights(const Array<Lagrange>& dim0,
                              const Array<Lagrange>& dim1) {
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

Grid<Matrix, 2> dinterpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1, Index dim) {
  Grid<Matrix, 2> out(dim0.size(), dim1.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = dinterpweights(dim0[i], dim1[j], dim);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

Numeric interp(const ConstMatrixView& yi, const ConstMatrixView& iw,
               const Lagrange& dim0, const Lagrange& dim1) {
  Numeric out(0.0);
  const Index I = yi.nrows();
  const Index J = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      out += iw(i, j) * yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(MatrixView out, const ConstMatrixView& iy,
              const Grid<Matrix, 2>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1) {
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = interp(iy, iw(i, j), dim0[i], dim1[j]);
}

Matrix reinterp(const ConstMatrixView& iy, const Grid<Matrix, 2>& iw,
                const Array<Lagrange>& dim0,
                const Array<Lagrange>& dim1) {
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

Grid<Tensor3, 3> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2) {
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

Grid<Tensor3, 3> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2, Index dim) {
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

Numeric interp(const ConstTensor3View& yi, const ConstTensor3View& iw,
               const Lagrange& dim0, const Lagrange& dim1,
               const Lagrange& dim2) {
  Numeric out(0.0);
  const Index I = yi.npages();
  const Index J = yi.nrows();
  const Index K = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out += iw(i, j, k) * yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J), cycler(k + dim2.pos, K));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(Tensor3View out, const ConstTensor3View& iy,
              const Grid<Tensor3, 3>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2) {
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        out(i, j, k) = interp(iy, iw(i, j, k), dim0[i], dim1[j], dim2[k]);
}

Tensor3 reinterp(const ConstTensor3View& iy, const Grid<Tensor3, 3>& iw,
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2) {
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

Grid<Tensor4, 4> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3) {
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

Grid<Tensor4, 4> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3, Index dim) {
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

Numeric interp(const ConstTensor4View& yi, const ConstTensor4View& iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3) {
  Numeric out(0.0);
  const Index I = yi.nbooks();
  const Index J = yi.npages();
  const Index K = yi.nrows();
  const Index L = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          out += iw(i, j, k, l) *
                 yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J),
                    cycler(k + dim2.pos, K), cycler(l + dim3.pos, L));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(Tensor4View out, const ConstTensor4View& iy,
              const Grid<Tensor4, 4>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3) {
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          out(i, j, k, l) =
              interp(iy, iw(i, j, k, l), dim0[i], dim1[j], dim2[k], dim3[l]);
}

Tensor4 reinterp(const ConstTensor4View& iy, const Grid<Tensor4, 4>& iw,
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3) {
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

Grid<Tensor5, 5> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3,
                               const Array<Lagrange>& dim4) {
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

Grid<Tensor5, 5> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3,
                                const Array<Lagrange>& dim4, Index dim) {
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

Numeric interp(const ConstTensor5View& yi, const ConstTensor5View& iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4) {
  Numeric out(0.0);
  const Index I = yi.nshelves();
  const Index J = yi.nbooks();
  const Index K = yi.npages();
  const Index L = yi.nrows();
  const Index M = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            out += iw(i, j, k, l, m) *
                   yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J),
                      cycler(k + dim2.pos, K), cycler(l + dim3.pos, L),
                      cycler(m + dim4.pos, M));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(Tensor5View out, const ConstTensor5View& iy,
              const Grid<Tensor5, 5>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3,
              const Array<Lagrange>& dim4) {
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) = interp(iy, iw(i, j, k, l, m), dim0[i], dim1[j],
                                        dim2[k], dim3[l], dim4[m]);
}

Tensor5 reinterp(const ConstTensor5View& iy, const Grid<Tensor5, 5>& iw,
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3,
                 const Array<Lagrange>& dim4) {
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

Grid<Tensor6, 6> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3,
                               const Array<Lagrange>& dim4,
                               const Array<Lagrange>& dim5) {
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

Grid<Tensor6, 6> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3,
                                const Array<Lagrange>& dim4,
                                const Array<Lagrange>& dim5, Index dim) {
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

Numeric interp(const ConstTensor6View& yi, const ConstTensor6View& iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4,
               const Lagrange& dim5) {
  Numeric out(0.0);
  const Index I = yi.nvitrines();
  const Index J = yi.nshelves();
  const Index K = yi.nbooks();
  const Index L = yi.npages();
  const Index M = yi.nrows();
  const Index N = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              out += iw(i, j, k, l, m, n) *
                     yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J),
                        cycler(k + dim2.pos, K), cycler(l + dim3.pos, L),
                        cycler(m + dim4.pos, M), cycler(n + dim5.pos, N));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(Tensor6View out, const ConstTensor6View& iy,
              const Grid<Tensor6, 6>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3,
              const Array<Lagrange>& dim4,
              const Array<Lagrange>& dim5) {
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) =
                  interp(iy, iw(i, j, k, l, m, n), dim0[i], dim1[j], dim2[k],
                         dim3[l], dim4[m], dim5[n]);
}

Tensor6 reinterp(const ConstTensor6View& iy, const Grid<Tensor6, 6>& iw,
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3,
                 const Array<Lagrange>& dim4,
                 const Array<Lagrange>& dim5) {
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

Grid<Tensor7, 7> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3,
                               const Array<Lagrange>& dim4,
                               const Array<Lagrange>& dim5,
                               const Array<Lagrange>& dim6) {
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

Grid<Tensor7, 7> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3,
                                const Array<Lagrange>& dim4,
                                const Array<Lagrange>& dim5,
                                const Array<Lagrange>& dim6, Index dim) {
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

Numeric interp(const ConstTensor7View& yi, const ConstTensor7View& iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5,
               const Lagrange& dim6) {
  Numeric out(0.0);
  const Index I = yi.nlibraries();
  const Index J = yi.nvitrines();
  const Index K = yi.nshelves();
  const Index L = yi.nbooks();
  const Index M = yi.npages();
  const Index N = yi.nrows();
  const Index O = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              for (Index o = 0; o < dim6.size(); o++)
                out += iw(i, j, k, l, m, n, o) *
                       yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J),
                          cycler(k + dim2.pos, K), cycler(l + dim3.pos, L),
                          cycler(m + dim4.pos, M), cycler(n + dim5.pos, N),
                          cycler(o + dim6.pos, O));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

void reinterp(Tensor7View out, const ConstTensor7View& iy,
              const Grid<Tensor7, 7>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3,
              const Array<Lagrange>& dim4,
              const Array<Lagrange>& dim5,
              const Array<Lagrange>& dim6) {
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
}

Tensor7 reinterp(const ConstTensor7View& iy, const Grid<Tensor7, 7>& iw,
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3,
                 const Array<Lagrange>& dim4,
                 const Array<Lagrange>& dim5,
                 const Array<Lagrange>& dim6) {
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
