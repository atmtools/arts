#ifndef lag_interp_h
#define lag_interp_h

#include <array>
#include <numeric>

#include "matpackVII.h"

namespace Interpolation {

/*! A Lagrange interpolation computer */
struct Lagrange {
  Index pos0;
  Vector lx;
  Vector dlx;
  
  /*! Update the weights for a new x */
  void update(Numeric x, const ConstVectorView xi, const bool do_derivative) noexcept {
    for (Index j=0; j <= lx.nelem(); j++) {
      lx[j] = l(x, xi, j);
      if (do_derivative) dlx[j] = dl(x, xi, j);
    }
  }
  
  /*! Standard and only initializer */
  Lagrange(Index ind0, Index polyorder, Numeric x, const ConstVectorView xi, const bool do_derivative) noexcept : pos0(ind0), lx(polyorder+1), dlx(polyorder+1) {update(x, xi, do_derivative);}
  
private:
  /*! Computes the weights for a given coefficient */
  Numeric l(Numeric x, const ConstVectorView xi, Index j) {
    Numeric val = 1.0;
    for (Index m=pos0; m <= pos0 + lx.nelem(); m++) {
      if (m not_eq j) {
        val *= (x - xi[m]) / (xi[j] - xi[m]);
      }
    }
    return val;
  }
  
  /*! Computes the derivatives of the weights for a given coefficient */
  Numeric dl(Numeric x, const ConstVectorView xi, Index j) {
    Numeric dval = 0.0;
    for (Index i=pos0; i <= pos0 + lx.nelem(); i++) {
      if (i not_eq j) {
        Numeric val = 1.0;
        for (Index m=pos0; m <= pos0 + lx.nelem(); m++) {
          if (m not_eq j and m not_eq i) {
            val *= (x - xi[m]) / (xi[j] - xi[m]);
          }
        }
        dval += val / (xi[j] - xi[i]);
      }
    }
    return dval;
  }
};

/*! Interpolation weights for a 1D reduction
 * 
 * @param[in] dim0 - Interpolation along dimension0
 * @return Vector - interpweights
 */
Vector interpweights(const Lagrange& dim0);

/*! Interpolation weights derivatives for a 1D reduction
 * 
 * @param[in] dim0 - Interpolation along dimension0
 * @return Vector - interpweights derivatives along 0th dimension
 */
Vector dinterpweights(const Lagrange& dim0);

/*! Interpolation weights for a 2D reduction
 * 
 * @param[in] dim0 - Interpolation along dimension0
 * @param[in] dim1 - Interpolation along dimension1
 * @return Matrix - interpweights
 */
Matrix interpweights(const Lagrange& dim0, const Lagrange& dim1);

/*! Interpolation weights derivatives for a 2D reduction
 * 
 * @param[in] dim0 - Interpolation along dimension0
 * @param[in] dim1 - Interpolation along dimension1
 * @return Vector - interpweights derivatives along dim dimension
 */
Matrix dinterpweights(const Lagrange& dim0, const Lagrange& dim1, Index dim);

Tensor3 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2);
Tensor3 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, Index dim);

Tensor4 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3);
Tensor4 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, Index dim);

Tensor5 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4);
Tensor5 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, Index dim);

Tensor6 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5);
Tensor6 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, Index dim);

Tensor7 interpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const Lagrange& dim6);
Tensor7 dinterpweights(const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const Lagrange& dim6, Index dim);

Numeric interp(const ConstVectorView yi, const ConstVectorView iw, const Lagrange& dim0);
Numeric interp(const ConstMatrixView yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1);
Numeric interp(const ConstTensor3View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2);
Numeric interp(const ConstTensor4View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3);
Numeric interp(const ConstTensor5View yi, const ConstTensor5View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4);
Numeric interp(const ConstTensor6View yi, const ConstTensor6View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5);
Numeric interp(const ConstTensor7View yi, const ConstTensor7View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const Lagrange& dim6);

Vector interp(const ConstMatrixView yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis);
Vector interp(const ConstTensor3View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis);
Vector interp(const ConstTensor4View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis);
Vector interp(const ConstTensor5View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const std::array<Index, 4>& axis);
Vector interp(const ConstTensor6View yi, const ConstTensor5View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const std::array<Index, 5>& axis);
Vector interp(const ConstTensor7View yi, const ConstTensor6View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5, const std::array<Index, 6>& axis);

Matrix interp(const ConstTensor3View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis);
Matrix interp(const ConstTensor4View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis);
Matrix interp(const ConstTensor5View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis);
Matrix interp(const ConstTensor6View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const std::array<Index, 4>& axis);
Matrix interp(const ConstTensor7View yi, const ConstTensor5View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const Lagrange& dim4, const std::array<Index, 5>& axis);

Tensor3 interp(const ConstTensor4View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis);
Tensor3 interp(const ConstTensor5View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis);
Tensor3 interp(const ConstTensor6View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis);
Tensor3 interp(const ConstTensor7View yi, const ConstTensor4View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const Lagrange& dim3, const std::array<Index, 4>& axis);

Tensor4 interp(const ConstTensor5View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis);
Tensor4 interp(const ConstTensor6View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis);
Tensor4 interp(const ConstTensor7View yi, const ConstTensor3View iw, const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2, const std::array<Index, 3>& axis);

Tensor5 interp(const ConstTensor6View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis);
Tensor5 interp(const ConstTensor7View yi, const ConstMatrixView iw, const Lagrange& dim0, const Lagrange& dim1, const std::array<Index, 2>& axis);

Tensor6 interp(const ConstTensor7View yi, const ConstVectorView iw, const Lagrange& dim0, const std::array<Index, 1>& axis);
}  // namespace Interpolation

#endif  // lag_interp_h
