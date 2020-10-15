#ifndef interpolation_lagrange_h
#define interpolation_lagrange_h

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#include "matpackVII.h"

namespace Interpolation {

/*! A Lagrange interpolation computer */
struct Lagrange {
  Index pos;
  Vector lx;
  Vector dlx;

  /* Number of weights */
  Index size() const noexcept { return lx.nelem(); }

  /*! Update the weights for a new x
   *
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] extrapol Level of extrapolation
   */
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

      // Adjust back so x is between the two Xs (if possible and not at the end)
      if (pos > 0 and xi[pos] > x) pos--;

      // Set weights
      for (Index j = 0; j < p; j++) {
        lx[j] = l(x, xi, j, p);
        dlx[j] = dl(x, xi, j, p);
      }
    }
  }

  /*! Standard and only initializer, assumes sorted xi and atleast 2 of them
   *
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] polyorder Polynominal degree
   * @param[in] extrapol Level of extrapolation
   */
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
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstVectorView yi, const ConstVectorView iw,
               const Lagrange& dim0);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstMatrixView yi, const ConstMatrixView iw,
               const Lagrange& dim0, const Lagrange& dim1);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
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
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
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
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @param[in] dim4 - Lagrange weights along the dimension
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
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @param[in] dim4 - Lagrange weights along the dimension
 * @param[in] dim5 - Lagrange weights along the dimension
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
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @param[in] dim4 - Lagrange weights along the dimension
 * @param[in] dim5 - Lagrange weights along the dimension
 * @param[in] dim6 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstTensor7View yi, const ConstTensor7View iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5,
               const Lagrange& dim6);

/*! Compute the multiplication of all inds */
template <typename... Inds>
constexpr std::size_t mul(Inds... inds) {
  return (std::size_t(inds) * ...);
}

/*! Compute the multiplication of all inds in arr */
template <std::size_t N>
constexpr std::size_t mul(const std::array<std::size_t, N>& arr) {
  if constexpr (N == 0) {
    return 0;
  } else {
    std::size_t out = 1;
    for (auto i : arr) out *= i;
    return out;
  }
}

/** Row-major grid creation */
template <typename b, std::size_t n>
class Grid {
  std::unique_ptr<b[]> ptr;
  std::array<std::size_t, n> gridsize;

  std::size_t nelem() const { return mul(gridsize); }

 public:
  static constexpr std::size_t N = n;
  using base = b;
  static_assert(N, "Must have size");

  template <typename... Inds>
  Grid(Inds... inds)
      : ptr(std::make_unique<base[]>(mul(inds...))),
        gridsize({std::size_t(inds)...}) {
    static_assert(sizeof...(Inds) == N,
                  "Must have same size for initialization");
  }

  template <typename... Inds>
  base& operator()(Inds... inds) noexcept {
    return ptr[index(std::array<std::size_t, N>{std::size_t(inds)...})];
  }

  template <typename... Inds>
  const base& operator()(Inds... inds) const noexcept {
    return ptr[index(std::array<std::size_t, N>{std::size_t(inds)...})];
  }

  friend std::ostream& operator<<(std::ostream& os, const Grid& g) {
    std::size_t i = 0;
    const std::size_t nel = g.nelem();
    while (i < nel) {
      os << g.ptr[i];
      i++;
      if (i not_eq 0 and i not_eq nel and i % g.gridsize.back() == 0)
        os << '\n';
      else if (i not_eq nel)
        os << ' ';
    }

    return os;
  }

 private:
  base& operator()(std::array<std::size_t, N> inds) noexcept {
    return ptr[index(inds)];
  }

  const base& operator()(std::array<std::size_t, N> inds) const noexcept {
    return ptr[index(inds)];
  }

  std::size_t index(std::array<std::size_t, N> ind) const noexcept {
    std::size_t posmul{gridsize.back()};
    std::size_t pos{ind.back()};
    for (std::size_t i{N - 2}; i < N; i--) {
      pos += posmul * ind[i];
      posmul *= gridsize[i];
      if (ind[i] >= gridsize[i]) {
        std::cerr << "Out of range\n";
        std::terminate();
      }
    }
    return pos;
  }
};  // Grid

/*! Gets a vector of Lagrange interpolation points
 *
 * @param[in] x New grid positions
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] extrapol Level of extrapolation
 * @return vector of Lagrange
 */
std::vector<Lagrange> LagrangeVector(const ConstVectorView x,
                                     const ConstVectorView xi,
                                     const Index polyorder = 1,
                                     const Numeric extrapol = 0.5);

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
Grid<Vector, 1> interpweights(const std::vector<Lagrange>& dim0);

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
Grid<Vector, 1> dinterpweights(const std::vector<Lagrange>& dim0);

/*! Interpolation weights for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @return Matrix - interpweights
 */
Grid<Matrix, 2> interpweights(const std::vector<Lagrange>& dim0,
                              const std::vector<Lagrange>& dim1);

/*! Interpolation weights derivative for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Matrix - interpweights derivative along dim dimension
 */
Grid<Matrix, 2> dinterpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1, Index dim);

/*! Interpolation weights for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @return Tensor3 - interpweights
 */
Grid<Tensor3, 3> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2);

/*! Interpolation weights derivative for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor3 - interpweights derivative along dim dimension
 */
Grid<Tensor3, 3> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2, Index dim);

/*! Interpolation weights for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @return Tensor4 - interpweights
 */
Grid<Tensor4, 4> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3);

/*! Interpolation weights derivative for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor4 - interpweights derivative along dim dimension
 */
Grid<Tensor4, 4> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3, Index dim);

/*! Interpolation weights for a 5D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @return Tensor5 - interpweights
 */
Grid<Tensor5, 5> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3,
                               const std::vector<Lagrange>& dim4);

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
Grid<Tensor5, 5> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3,
                                const std::vector<Lagrange>& dim4, Index dim);

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
Grid<Tensor6, 6> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3,
                               const std::vector<Lagrange>& dim4,
                               const std::vector<Lagrange>& dim5);

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
Grid<Tensor6, 6> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3,
                                const std::vector<Lagrange>& dim4,
                                const std::vector<Lagrange>& dim5, Index dim);

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
Grid<Tensor7, 7> interpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1,
                               const std::vector<Lagrange>& dim2,
                               const std::vector<Lagrange>& dim3,
                               const std::vector<Lagrange>& dim4,
                               const std::vector<Lagrange>& dim5,
                               const std::vector<Lagrange>& dim6);

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
Grid<Tensor7, 7> dinterpweights(const std::vector<Lagrange>& dim0,
                                const std::vector<Lagrange>& dim1,
                                const std::vector<Lagrange>& dim2,
                                const std::vector<Lagrange>& dim3,
                                const std::vector<Lagrange>& dim4,
                                const std::vector<Lagrange>& dim5,
                                const std::vector<Lagrange>& dim6, Index dim);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Vector of interpolated value
 */
Vector reinterp(const ConstVectorView iy, const Grid<Vector, 1>& iw,
                const std::vector<Lagrange>& dim0);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Matrix of interpolated value
 */
Matrix reinterp(const ConstMatrixView iy, const Grid<Matrix, 2>& iw,
                const std::vector<Lagrange>& dim0,
                const std::vector<Lagrange>& dim1);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @return Tensor3 of interpolated value
 */
Tensor3 reinterp(const ConstTensor3View iy, const Grid<Tensor3, 3>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @return Tensor4 of interpolated value
 */
Tensor4 reinterp(const ConstTensor4View iy, const Grid<Tensor4, 4>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @param[in] dim4 - Lagrange weights along the dimension
 * @return Tensor5 of interpolated value
 */
Tensor5 reinterp(const ConstTensor5View iy, const Grid<Tensor5, 5>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @param[in] dim4 - Lagrange weights along the dimension
 * @param[in] dim5 - Lagrange weights along the dimension
 * @return Tensor6 of interpolated value
 */
Tensor6 reinterp(const ConstTensor6View iy, const Grid<Tensor6, 6>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4,
                 const std::vector<Lagrange>& dim5);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @param[in] dim4 - Lagrange weights along the dimension
 * @param[in] dim5 - Lagrange weights along the dimension
 * @param[in] dim6 - Lagrange weights along the dimension
 * @return Tensor7 of interpolated value
 */
Tensor7 reinterp(const ConstTensor7View iy, const Grid<Tensor7, 7>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4,
                 const std::vector<Lagrange>& dim5,
                 const std::vector<Lagrange>& dim6);
}  // namespace Interpolation

#endif  // interpolation_lagrange_h
