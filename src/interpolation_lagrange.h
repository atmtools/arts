#ifndef lag_interp_h
#define lag_interp_h

#include <array>
#include <numeric>

#include "matpackVII.h"

namespace Interpolation {

/*! A Lagrange interpolation computer */
struct Lagrange {
  Index pos;
  Vector lx;
  Vector dlx;

  /*! Update the weights for a new x */
  void update(const Numeric x, const ConstVectorView xi,
              const Numeric extrapol = 0.5) {
    const Index n = xi.nelem();
    const Index p = lx.nelem();

    if (p >= n) {
      throw std::runtime_error(
          "Requesting greater interpolation order than possible with given "
          "input grid\n");
    } else if (extrapol >= 0 and
               (x < (xi[0] - extrapol * (xi[1] - xi[0])) or
                x > (xi[n - 1] + extrapol * (xi[n - 1] - xi[n - 2])))) {
      std::ostringstream os;
      os << "Extrapolation factor too small at: " << extrapol
         << ", for grid: " << xi << '\n';
      throw std::runtime_error(os.str());
    } else {
      // Find first larger x
      for (pos = 0; pos < n - p; pos++)
        if (xi[pos] > x) break;

      // Adjust back so x is between the two Xs (if possible)
      if (pos > 0) pos--;

      // Set weights
      for (Index j = 0; j < p; j++) {
        lx[j] = l(x, xi, j, p);
        dlx[j] = dl(x, xi, j, p);
      }
    }
  }

  /*! Standard and only initializer, assumes sorted xi and atleast 2 of them */
  Lagrange(const Numeric x, const ConstVectorView xi, const Index polyorder = 1,
           const Numeric extrapol = 0.5)
      : pos(0), lx(polyorder + 1), dlx(polyorder + 1) {
    update(x, xi, extrapol);
  }

  friend std::ostream& operator<<(std::ostream& os, const Lagrange& l) {
    return os << l.pos << ' ' << l.lx << ' ' << l.dlx;
  }

 private:
  /*! Computes the weights for a given coefficient */
  Numeric l(Numeric x, const ConstVectorView xi, Index j, Index n) {
    Numeric val = 1.0;
    for (Index m = 0; m < n; m++) {
      if (m not_eq j) {
        val *= (x - xi[m + pos]) / (xi[j + pos] - xi[m + pos]);
      }
    }
    return val;
  }

  /*! Computes the derivatives of the weights for a given coefficient */
  Numeric dl(Numeric x, const ConstVectorView xi, Index j, Index n) {
    Numeric dval = 0.0;
    for (Index i = 0; i < n; i++) {
      if (i not_eq j) {
        Numeric val = 1.0;
        for (Index m = 0; m < n; m++) {
          if (m not_eq j and m not_eq i) {
            val *= (x - xi[m + pos]) / (xi[j + pos] - xi[m + pos]);
          }
        }
        dval += val / (xi[j + pos] - xi[i + pos]);
      }
    }
    return dval;
  }
};

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
Vector interpweights(const Lagrange& dim0);

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
Vector dinterpweights(const Lagrange& dim0);

/*! Interpolation weights for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @return Matrix - interpweights
 */
Matrix interpweights(const Lagrange& dim0, const Lagrange& dim1);

/*! Interpolation weights derivative for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Matrix - interpweights derivative along dim dimension
 */
Matrix dinterpweights(const Lagrange& dim0, const Lagrange& dim1, Index dim);

/*! Interpolation weights for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @return Tensor3 - interpweights
 */
Tensor3 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2);

/*! Interpolation weights derivative for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor3 - interpweights derivative along dim dimension
 */
Tensor3 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, Index dim);

/*! Interpolation weights for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @return Tensor4 - interpweights
 */
Tensor4 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3);

/*! Interpolation weights derivative for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor4 - interpweights derivative along dim dimension
 */
Tensor4 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3, Index dim);

/*! Interpolation weights for a 5D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @return Tensor5 - interpweights
 */
Tensor5 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3,
                      const Lagrange& dim4);

/*! Interpolation weights derivative for a 5D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor5 - interpweights derivative along dim dimension
 */
Tensor5 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3,
                       const Lagrange& dim4, Index dim);

/*! Interpolation weights for a 6D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @param[in] dim5 - Interpolation along dimension 5
 * @return Tensor6 - interpweights
 */
Tensor6 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3,
                      const Lagrange& dim4, const Lagrange& dim5);

/*! Interpolation weights derivative for a 6D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @param[in] dim5 - Interpolation along dimension 5
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor6 - interpweights derivative along dim dimension
 */
Tensor6 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3,
                       const Lagrange& dim4, const Lagrange& dim5, Index dim);

/*! Interpolation weights for a 7D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @param[in] dim5 - Interpolation along dimension 5
 * @param[in] dim6 - Interpolation along dimension 6
 * @return Tensor7 - interpweights
 */
Tensor7 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2, const Lagrange& dim3,
                      const Lagrange& dim4, const Lagrange& dim5,
                      const Lagrange& dim6);

/*! Interpolation weights derivative for a 7D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @param[in] dim5 - Interpolation along dimension 5
 * @param[in] dim6 - Interpolation along dimension 6
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor7 - interpweights derivative along dim dimension
 */
Tensor7 dinterpweights(const Lagrange& dim0, const Lagrange& dim1,
                       const Lagrange& dim2, const Lagrange& dim3,
                       const Lagrange& dim4, const Lagrange& dim5,
                       const Lagrange& dim6, Index dim);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstVectorView yi, const ConstVectorView iw,
               const Lagrange& dim0);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstMatrixView yi, const ConstMatrixView iw,
               const Lagrange& dim0, const Lagrange& dim1);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstTensor3View yi, const ConstTensor3View iw,
               const Lagrange& dim0, const Lagrange& dim1,
               const Lagrange& dim2);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstTensor4View yi, const ConstTensor4View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] dim4 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstTensor5View yi, const ConstTensor5View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] dim4 - Lagrange weigths along the dimension
 * @param[in] dim5 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstTensor6View yi, const ConstTensor6View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4,
               const Lagrange& dim5);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] dim4 - Lagrange weigths along the dimension
 * @param[in] dim5 - Lagrange weigths along the dimension
 * @param[in] dim6 - Lagrange weigths along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstTensor7View yi, const ConstTensor7View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5,
               const Lagrange& dim6);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Vector with selected axis of original data squashed
 */
Vector interp(const ConstMatrixView yi, const ConstVectorView iw,
              const Lagrange& dim0, const std::array<Index, 1>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Vector with selected axis of original data squashed
 */
Vector interp(const ConstTensor3View yi, const ConstMatrixView iw,
              const Lagrange& dim0, const Lagrange& dim1,
              const std::array<Index, 2>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Vector with selected axis of original data squashed
 */
Vector interp(const ConstTensor4View yi, const ConstTensor3View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const std::array<Index, 3>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Vector with selected axis of original data squashed
 */
Vector interp(const ConstTensor5View yi, const ConstTensor4View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const Lagrange& dim3, const std::array<Index, 4>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] dim4 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Vector with selected axis of original data squashed
 */
Vector interp(const ConstTensor6View yi, const ConstTensor5View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const Lagrange& dim3, const Lagrange& dim4,
              const std::array<Index, 5>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] dim4 - Lagrange weigths along the dimension
 * @param[in] dim5 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Vector with selected axis of original data squashed
 */
Vector interp(const ConstTensor7View yi, const ConstTensor6View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5,
              const std::array<Index, 6>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Matrix with selected axis of original data squashed
 */
Matrix interp(const ConstTensor3View yi, const ConstVectorView iw,
              const Lagrange& dim0, const std::array<Index, 1>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Matrix with selected axis of original data squashed
 */
Matrix interp(const ConstTensor4View yi, const ConstMatrixView iw,
              const Lagrange& dim0, const Lagrange& dim1,
              const std::array<Index, 2>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Matrix with selected axis of original data squashed
 */
Matrix interp(const ConstTensor5View yi, const ConstTensor3View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const std::array<Index, 3>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Matrix with selected axis of original data squashed
 */
Matrix interp(const ConstTensor6View yi, const ConstTensor4View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const Lagrange& dim3, const std::array<Index, 4>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] dim4 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Matrix with selected axis of original data squashed
 */
Matrix interp(const ConstTensor7View yi, const ConstTensor5View iw,
              const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
              const Lagrange& dim3, const Lagrange& dim4,
              const std::array<Index, 5>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor3 with selected axis of original data squashed
 */
Tensor3 interp(const ConstTensor4View yi, const ConstVectorView iw,
               const Lagrange& dim0, const std::array<Index, 1>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor3 with selected axis of original data squashed
 */
Tensor3 interp(const ConstTensor5View yi, const ConstMatrixView iw,
               const Lagrange& dim0, const Lagrange& dim1,
               const std::array<Index, 2>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor3 with selected axis of original data squashed
 */
Tensor3 interp(const ConstTensor6View yi, const ConstTensor3View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const std::array<Index, 3>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] dim3 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor3 with selected axis of original data squashed
 */
Tensor3 interp(const ConstTensor7View yi, const ConstTensor4View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const std::array<Index, 4>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor4 with selected axis of original data squashed
 */
Tensor4 interp(const ConstTensor5View yi, const ConstVectorView iw,
               const Lagrange& dim0, const std::array<Index, 1>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor4 with selected axis of original data squashed
 */
Tensor4 interp(const ConstTensor6View yi, const ConstMatrixView iw,
               const Lagrange& dim0, const Lagrange& dim1,
               const std::array<Index, 2>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] dim2 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor4 with selected axis of original data squashed
 */
Tensor4 interp(const ConstTensor7View yi, const ConstTensor3View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const std::array<Index, 3>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor5 with selected axis of original data squashed
 */
Tensor5 interp(const ConstTensor6View yi, const ConstVectorView iw,
               const Lagrange& dim0, const std::array<Index, 1>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] dim1 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor5 with selected axis of original data squashed
 */
Tensor5 interp(const ConstTensor7View yi, const ConstMatrixView iw,
               const Lagrange& dim0, const Lagrange& dim1,
               const std::array<Index, 2>& axis);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weigths along the dimension
 * @param[in] axis - Axis to squash
 * @return Tensor6 with selected axis of original data squashed
 */
Tensor6 interp(const ConstTensor7View yi, const ConstVectorView iw,
               const Lagrange& dim0, const std::array<Index, 1>& axis);
}  // namespace Interpolation

#endif  // lag_interp_h
