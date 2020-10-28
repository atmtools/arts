#ifndef interpolation_lagrange_h
#define interpolation_lagrange_h

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#include "enums.h"
#include "matpackVII.h"

namespace Interpolation {

ENUMCLASS(LagrangeType, unsigned char,
          Linear,
          Log
         )

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

/*! A Lagrange interpolation computer */
struct Lagrange {
  Index pos;
  Vector lx;
  Vector dlx;

  /* Number of weights */
  Index size() const noexcept { return lx.nelem(); }

  /*! Standard and only initializer, assumes sorted xi and atleast 2 of them
   *
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] polyorder Polynominal degree
   * @param[in] extrapol Level of extrapolation
   * @param[in] type Type of Lagrange (Normal or Log)
   */
  Lagrange(Index p0, const Numeric x, const ConstVectorView& xi, const Index polyorder = 1,
           const Numeric extrapol = 0.5, const LagrangeType type=LagrangeType::Linear)
      : lx(polyorder + 1), dlx(polyorder + 1) {
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
      for (; p0 < n - p; p0++)
        if (xi[p0] > x) break;

      // Adjust back so x is between the two Xs (if possible and not at the end)
      if (p0 > 0 and xi[p0] > x) p0--;
      
      // Set the position
      pos = p0;

      // Set weights
      for (Index j = 0; j < p; j++) lx[j] = l(x, xi, j, p, pos, type);
      
      // Set derivatives after the weights
      for (Index j = 0; j < p; j++) dlx[j] = dl(x, xi, lx, j, p, pos, type);
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const Lagrange& l) {
    return os << l.pos << ' ' << l.lx << ' ' << l.dlx;
  }

 private:
  /*! Computes the weights for a given coefficient */
  static Numeric l(const Numeric x, const ConstVectorView& xi, const Index j, const Index n, const Index p0, const LagrangeType type) noexcept {
    Numeric val = 1.0;
    if (type == LagrangeType::Log) {
      for (Index m = 0; m < n; m++) {
        if (m not_eq j) {
          val *= (std::log(x) - std::log(xi[m + p0])) / (std::log(xi[j +  p0]) - std::log(xi[m +  p0]));
        }
      }
    } else if (type == LagrangeType::Linear) {
      for (Index m = 0; m < n; m++) {
        if (m not_eq j) {
          val *= (x - xi[m +  p0]) / (xi[j +  p0] - xi[m +  p0]);
        }
      }
    }
    return val;
  }

  /*! Computes the derivatives of the weights for a given coefficient
   * 
   * Must be called with all lx known
   */
  static Numeric dl(const Numeric x, const ConstVectorView& xi, const ConstVectorView& li, const Index j, const Index n, const Index p0, const LagrangeType type) noexcept {
    Numeric dval = 0.0;
    if (type == LagrangeType::Log) {
      for (Index i = 0; i < n; i++) {
        if (i not_eq j) {
          dval += li[j] / (std::log(x) - std::log(xi[i+p0]));
        }
      }
    } else if (type == LagrangeType::Linear) {
      for (Index i = 0; i < n; i++) {
        if (i not_eq j) {
          dval += li[j] / (x - xi[i+p0]);
        }
      }
    }
    return dval;
  }
};

/*! A Fixed Lagrange interpolation computer */
template <std::size_t PolyOrder>
struct FixedLagrange {
  Index pos;
  std::array<Numeric, PolyOrder + 1> lx;
  std::array<Numeric, PolyOrder + 1> dlx;

  /* Number of weights */
  static constexpr Index size() noexcept { return PolyOrder + 1; }

  /*! Standard initializer from Vector-types
   *
   * @param[in] p0 Guess of first position
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] type Type of grid (Linear or Log)
   */
  template <class SortedVectorType>
  constexpr FixedLagrange(const Index p0, const Numeric x, const SortedVectorType& xi,
                          const LagrangeType type=LagrangeType::Linear) : pos(pos_finder(x, xi, p0)),
                          lx(lx_finder(pos, x, xi, type)), dlx(dlx_finder(pos, x, xi, lx, type)) {
  }
            
 private:
  /*! Finds the position */
  template <class SortedVectorType>
  static constexpr Index pos_finder(const Numeric x,
                                    const SortedVectorType& xi,
                                    const decltype(xi.size()) pos0=0) noexcept {
    decltype(xi.size()) p = pos0;
    for (; p < xi.size() - size(); p++)
      if (xi[p] > x) break;

    // Adjust back so x is between the two Xs (if possible and not at the end)
    if (p > 0 and xi[p] > x) p--;
    return p;
  }

  /*! Finds lx */
  template <class SortedVectorType> static
  constexpr std::array<Numeric, PolyOrder + 1> lx_finder(
    const Index p0,
    const Numeric x, const SortedVectorType& xi,
    const LagrangeType type) noexcept {
    std::array<Numeric, PolyOrder + 1> out{};
    if (type == LagrangeType::Linear) {
      for (Index j = 0; j < size(); j++)
        out[j] = l<LagrangeType::Linear>(p0, x, xi, j);
    } else if (type == LagrangeType::Log) {
      for (Index j = 0; j < size(); j++)
        out[j] = l<LagrangeType::Log>(p0, x, xi, j);
    }
    return out;
  }

  /*! Finds dlx */
  template <class SortedVectorType> static
  constexpr std::array<Numeric, PolyOrder + 1> dlx_finder(
    const Index p0,
    const Numeric x, const SortedVectorType& xi, const std::array<Numeric, PolyOrder + 1>& li,
    const LagrangeType type) noexcept {
    std::array<Numeric, PolyOrder + 1> out{};
    if (type == LagrangeType::Linear) {
      for (Index j = 0; j < size(); j++)
        out[j] = dl<LagrangeType::Linear>(p0, x, xi, li, j);
    } else if (type == LagrangeType::Log) {
      for (Index j = 0; j < size(); j++)
        out[j] = dl<LagrangeType::Log>(p0, x, xi, li, j);
    }
    return out;
  }

  /*! Computes the weights for a given coefficient */
  template <LagrangeType type, typename SortedVectorType>
  static constexpr Numeric l(const Index p0, const Numeric x, const SortedVectorType& xi,
                      const Index j) noexcept {
    Numeric val = 1.0;
    for (Index m = 0; m < size(); m++) {
      if (m not_eq j) {
        if constexpr (type == LagrangeType::Log) {
          val *= (std::log(x) - std::log(xi[m + p0])) / (std::log(xi[j + p0]) - std::log(xi[m + p0]));
        } else if constexpr (type == LagrangeType::Linear) {
          val *= (x - xi[m + p0]) / (xi[j + p0] - xi[m + p0]);
        }
      }
    }
    return val;
  }
  
  /*! Computes the derivatives of the weights for a given coefficient
   * 
   * Must be called with all lx known
   */
  template <LagrangeType type, typename SortedVectorType>
  static constexpr Numeric dl(const Index p0, const Numeric x, const SortedVectorType& xi, const std::array<Numeric, PolyOrder + 1>& li,
                              const Index j) noexcept {
    Numeric dval = 0.0;
    for (Index i = 0; i < size(); i++) {
      if (i not_eq j) {
        if (type == LagrangeType::Log) {
          dval += li[j] / (std::log(x) - std::log(xi[i+p0]));
        } else if (type == LagrangeType::Linear) {
          dval += li[j] / (x - xi[i+p0]);
        }
      }
    }
    return dval;
  }
};

/** Row-major grid creation */
template <typename b, std::size_t n>
class Grid {
  Array<b> ptr;
  std::array<std::size_t, n> gridsize;

  std::size_t nelem() const { return ptr.size(); }

 public:
  static constexpr std::size_t N = n;
  using base = b;
  static_assert(N, "Must have size");

  template <typename... Inds>
  Grid(Inds... inds)
      : ptr(mul(inds...)),
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

/** Row-major fixed grid creation */
template <typename b, std::size_t... Sizes>
class FixedGrid {
  std::array<b, mul(Sizes...)> ptr;
  static constexpr std::array<std::size_t, sizeof...(Sizes)> gridsize{Sizes...};

  constexpr std::size_t nelem() const { return mul(Sizes...); }

 public:
  static constexpr std::size_t N = sizeof...(Sizes);
  using base = b;
  static_assert(N, "Must have size");

  template <typename... Inds>
  base& operator()(Inds... inds) noexcept {
    return ptr[index(std::array<std::size_t, N>{std::size_t(inds)...})];
  }

  template <typename... Inds>
  const base& operator()(Inds... inds) const noexcept {
    return ptr[index(std::array<std::size_t, N>{std::size_t(inds)...})];
  }

  friend std::ostream& operator<<(std::ostream& os, const FixedGrid& g) {
    std::size_t i = 0;
    constexpr std::size_t nel = g.nelem();
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
};  // FixedGrid

////////////////////////////////////////////////////////
////////////////////// For reinterpreting interpolations
////////////////////////////////////////////////////////

/*! Gets a vector of Lagrange interpolation points
 *
 * @param[in] x New grid positions
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] extrapol Level of extrapolation
 * @return vector of Lagrange
 */
std::vector<Lagrange> LagrangeVector(const ConstVectorView& x,
                                     const ConstVectorView& xi,
                                     const Index polyorder,
                                     const Numeric extrapol);

/*! Gets a vector of Lagrange interpolation points
 *
 * @param[in] x New grid positions
 * @param[in] xi Old grid positions
 * @param[in] extrapol Level of extrapolation
 * @return vector of FixedLagrange
 */
template <std::size_t PolyOrder>
std::vector<FixedLagrange<PolyOrder>> FixedLagrangeVector(const ConstVectorView& xs,
                                                          const ConstVectorView& xi,
                                                          const LagrangeType type) {
  std::vector<FixedLagrange<PolyOrder>> out;
  out.reserve(xs.nelem());
  for (Index i = 0; i < xs.nelem(); i++) {
    if (i) {
      out.push_back(FixedLagrange<PolyOrder>(out[i-1].pos, xs[i], xi, type));
    } else {
      out.push_back(FixedLagrange<PolyOrder>(0, xs[i], xi, type));
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

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
const Vector& interpweights(const Lagrange& dim0);

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
Grid<Vector, 1> interpweights(const std::vector<Lagrange>& dim0);

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
template <std::size_t PolyOrder>
constexpr std::array<Numeric, PolyOrder + 1> interpweights(
    const FixedLagrange<PolyOrder>& dim0) {
  return dim0.lx;
}

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
template <std::size_t PolyOrder>
Grid<std::array<Numeric, PolyOrder + 1>, 1> interpweights(
    const std::vector<FixedLagrange<PolyOrder>>& dim0) {
  Grid<std::array<Numeric, PolyOrder + 1>, 1> out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++) out(i) = interpweights(dim0[i]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
const Vector& dinterpweights(const Lagrange& dim0);

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
Grid<Vector, 1> dinterpweights(const std::vector<Lagrange>& dim0);

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
template <std::size_t PolyOrder>
constexpr std::array<Numeric, PolyOrder + 1> dinterpweights(
    const FixedLagrange<PolyOrder>& dim0) {
  return dim0.dlx;
}

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
template <std::size_t PolyOrder>
Grid<std::array<Numeric, PolyOrder + 1>, 1> dinterpweights(
    const std::vector<FixedLagrange<PolyOrder>>& dim0) {
  Grid<std::array<Numeric, PolyOrder + 1>, 1> out(dim0.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    out(i) = dinterpweights(dim0[i]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstVectorView& yi, const ConstVectorView& iw,
               const Lagrange& dim0);

/*! Squashing fixed interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
template <std::size_t PolyOrder, class VectorType>
constexpr Numeric interp(const VectorType& yi,
                         const std::array<Numeric, PolyOrder + 1>& iw,
                         const FixedLagrange<PolyOrder>& dim0) {
  Numeric out(0.0);
  for (Index i = 0; i < dim0.size(); i++) out += iw[i] * yi[i + dim0.pos];
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Vector of interpolated value
 */
Vector reinterp(const ConstVectorView& iy, const Grid<Vector, 1>& iw,
                const std::vector<Lagrange>& dim0);

/*! Reinterpreting fixed interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Vector of interpolated value
 */
template <std::size_t PolyOrder>
Vector reinterp(const ConstVectorView& iy,
                const Grid<std::array<Numeric, PolyOrder + 1>, 1>& iw,
                const std::vector<FixedLagrange<PolyOrder>>& dim0) {
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

/*! Interpolation weights for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @return Matrix - interpweights
 */
Matrix interpweights(const Lagrange& dim0, const Lagrange& dim1);

/*! Interpolation weights for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @return Matrix - interpweights
 */
Grid<Matrix, 2> interpweights(const std::vector<Lagrange>& dim0,
                              const std::vector<Lagrange>& dim1);
/*! Interpolation weights for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @return Matrix - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1> interpweights(
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1) {
  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1> out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++) out(i, j) = dim0.lx[i] * dim1.lx[j];
  return out;
}

/*! Interpolation weights for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @return Matrix - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1>, 2> interpweights(
    const std::vector<FixedLagrange<PolyOrder0>>& dim0,
    const std::vector<FixedLagrange<PolyOrder1>>& dim1) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1>, 2> out(dim0.size(),
                                                                  dim1.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = interpweights(dim0[i], dim1[j]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

/*! Interpolation weights derivative for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Matrix - interpweights derivative along dim dimension
 */
Matrix dinterpweights(const Lagrange& dim0, const Lagrange& dim1, Index dim);

/*! Interpolation weights derivative for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Matrix - interpweights derivative along dim dimension
 */
Grid<Matrix, 2> dinterpweights(const std::vector<Lagrange>& dim0,
                               const std::vector<Lagrange>& dim1, Index dim);

/*! Interpolation weights derivative for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Matrix - interpweights derivative along dim dimension
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t DerivativeDim>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1> dinterpweights(
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1) {
  static_assert(DerivativeDim < 2);

  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1> out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      out(i, j) = (DerivativeDim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                  (DerivativeDim == 1 ? dim1.dlx[j] : dim1.lx[j]);
  return out;
}

/*! Interpolation weights derivative for a 2D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Matrix - interpweights derivative along dim dimension
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t DerivativeDim>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1>, 2> dinterpweights(
    const std::vector<FixedLagrange<PolyOrder0>>& dim0,
    const std::vector<FixedLagrange<PolyOrder1>>& dim1) {
  static_assert(DerivativeDim < 2);

  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1>, 2> out(dim0.size(),
                                                                  dim1.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      out(i, j) = dinterpweights<PolyOrder0, PolyOrder1, DerivativeDim>(
          dim0[i], dim1[j]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
Numeric interp(const ConstMatrixView& yi, const ConstMatrixView& iw,
               const Lagrange& dim0, const Lagrange& dim1);

/*! Squashing interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights or their derivatives from the Lagrange
 * routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Numeric of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1, class MatrixType>
constexpr Numeric interp(
    const MatrixType& yi,
    const FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1>& iw,
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1) {
  Numeric out(0.0);

  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      out += iw(i, j) * yi(i + dim0.pos, j + dim1.pos);
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Matrix of interpolated value
 */
Matrix reinterp(const ConstMatrixView& iy, const Grid<Matrix, 2>& iw,
                const std::vector<Lagrange>& dim0,
                const std::vector<Lagrange>& dim1);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Matrix of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1>
Matrix reinterp(
    const ConstMatrixView& iy,
    const Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1>, 2>& iw,
    const std::vector<FixedLagrange<PolyOrder0>>& dim0,
    const std::vector<FixedLagrange<PolyOrder1>>& dim1) {
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

/*! Interpolation weights for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @return Tensor3 - interpweights
 */
Tensor3 interpweights(const Lagrange& dim0, const Lagrange& dim1,
                      const Lagrange& dim2);

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

/*! Interpolation weights for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @return Tensor3 - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>
interpweights(const FixedLagrange<PolyOrder0>& dim0,
              const FixedLagrange<PolyOrder1>& dim1,
              const FixedLagrange<PolyOrder2>& dim2) {
  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1> out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out(i, j, k) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k];
  return out;
}

/*! Interpolation weights for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @return Tensor3 - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>, 3>
interpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
              const std::vector<FixedLagrange<PolyOrder1>>& dim1,
              const std::vector<FixedLagrange<PolyOrder2>>& dim2) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>, 3>
      out(dim0.size(), dim1.size(), dim2.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        out(i, j, k) = interpweights(dim0[i], dim1[j], dim2[k]);
  return out;
}

////////////////////////////////////////////////
/////////// Derivatives of Interpolation Weights
////////////////////////////////////////////////

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

/*! Interpolation weights derivative for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor3 - interpweights derivative along dim dimension
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t DerivativeDim>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>
dinterpweights(const FixedLagrange<PolyOrder0>& dim0,
               const FixedLagrange<PolyOrder1>& dim1,
               const FixedLagrange<PolyOrder2>& dim2) {
  static_assert(DerivativeDim < 3);

  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1> out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out(i, j, k) = (DerivativeDim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                       (DerivativeDim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                       (DerivativeDim == 2 ? dim2.dlx[k] : dim2.lx[k]);
  return out;
}

/*! Interpolation weights derivative for a 3D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor3 - interpweights derivative along dim dimension
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t DerivativeDim>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>, 3>
dinterpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
               const std::vector<FixedLagrange<PolyOrder1>>& dim1,
               const std::vector<FixedLagrange<PolyOrder2>>& dim2) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>, 3>
      out(dim0.size(), dim0.size(), dim2.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        out(i, j, k) =
            dinterpweights<PolyOrder0, PolyOrder1, PolyOrder2, DerivativeDim>(
                dim0[i], dim1[j], dim2[k]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

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
Numeric interp(const ConstTensor3View& yi, const ConstTensor3View& iw,
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
 * @return Numeric of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, class Tensor3Type>
constexpr Numeric interp(const Tensor3Type& yi,
                         const FixedGrid<Numeric, PolyOrder0 + 1,
                                         PolyOrder1 + 1, PolyOrder2 + 1>& iw,
                         const std::vector<FixedLagrange<PolyOrder0>>& dim0,
                         const std::vector<FixedLagrange<PolyOrder1>>& dim1,
                         const std::vector<FixedLagrange<PolyOrder2>>& dim2) {
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
Tensor3 reinterp(const ConstTensor3View& iy, const Grid<Tensor3, 3>& iw,
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
 * @return Tensor3 of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2>
Tensor3 reinterp(
    const ConstTensor3View& iy,
    const Grid<
        FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1>, 3>&
        iw,
    const std::vector<FixedLagrange<PolyOrder0>>& dim0,
    const std::vector<FixedLagrange<PolyOrder1>>& dim1,
    const std::vector<FixedLagrange<PolyOrder2>>& dim2) {
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

/*! Interpolation weights for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @return Tensor4 - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1>
interpweights(const FixedLagrange<PolyOrder0>& dim0,
              const FixedLagrange<PolyOrder1>& dim1,
              const FixedLagrange<PolyOrder2>& dim2,
              const FixedLagrange<PolyOrder3>& dim3) {
  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1>
      out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          out(i, j, k, l) = dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l];
  return out;
}

/*! Interpolation weights for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @return Tensor4 - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1>,
     4>
interpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
              const std::vector<FixedLagrange<PolyOrder1>>& dim1,
              const std::vector<FixedLagrange<PolyOrder2>>& dim2,
              const std::vector<FixedLagrange<PolyOrder3>>& dim3) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                 PolyOrder3 + 1>,
       4>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
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

/*! Interpolation weights derivative for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor4 - interpweights derivative along dim dimension
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t DerivativeDim>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1>
dinterpweights(const FixedLagrange<PolyOrder0>& dim0,
               const FixedLagrange<PolyOrder1>& dim1,
               const FixedLagrange<PolyOrder2>& dim2,
               const FixedLagrange<PolyOrder3>& dim3) {
  static_assert(DerivativeDim < 4);

  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1>
      out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          out(i, j, k, l) = (DerivativeDim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                            (DerivativeDim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                            (DerivativeDim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                            (DerivativeDim == 3 ? dim3.dlx[l] : dim3.lx[l]);
  return out;
}

/*! Interpolation weights derivative for a 4D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim - Axis along which to compute the derivatives
 * @return Tensor4 - interpweights derivative along dim dimension
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t DerivativeDim>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1>,
     4>
dinterpweights(const std::vector<Lagrange>& dim0,
               const std::vector<Lagrange>& dim1,
               const std::vector<Lagrange>& dim2,
               const std::vector<Lagrange>& dim3) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                 PolyOrder3 + 1>,
       4>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          out(i, j, k, l) =
              dinterpweights<PolyOrder0, PolyOrder1, PolyOrder2, PolyOrder3,
                             DerivativeDim>(dim0[i], dim1[j], dim2[k], dim3[l]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

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
Numeric interp(const ConstTensor4View& yi, const ConstTensor4View& iw,
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
 * @return Numeric of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3, class Tensor4Type>
constexpr Numeric interp(
    const Tensor4Type& yi,
    const FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1>& iw,
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1,
    const FixedLagrange<PolyOrder2>& dim2,
    const FixedLagrange<PolyOrder3>& dim3) {
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
Tensor4 reinterp(const ConstTensor4View& iy, const Grid<Tensor4, 4>& iw,
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
 * @return Tensor4 of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3>
Tensor4 reinterp(const ConstTensor4View& iy,
                 const Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1,
                                      PolyOrder2 + 1, PolyOrder3 + 1>,
                            4>& iw,
                 const std::vector<FixedLagrange<PolyOrder0>>& dim0,
                 const std::vector<FixedLagrange<PolyOrder1>>& dim1,
                 const std::vector<FixedLagrange<PolyOrder2>>& dim2,
                 const std::vector<FixedLagrange<PolyOrder3>>& dim3) {
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

/*! Interpolation weights for a 5D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @return Tensor5 - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1>
interpweights(const FixedLagrange<PolyOrder0>& dim0,
              const FixedLagrange<PolyOrder1>& dim1,
              const FixedLagrange<PolyOrder2>& dim2,
              const FixedLagrange<PolyOrder3>& dim3,
              const FixedLagrange<PolyOrder4>& dim4) {
  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1, PolyOrder4 + 1>
      out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) =
                dim0.lx[i] * dim1.lx[j] * dim2.lx[k] * dim3.lx[l] * dim4.lx[m];
  return out;
}

/*! Interpolation weights for a 5D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @param[in] dim1 - Interpolation along dimension 1
 * @param[in] dim2 - Interpolation along dimension 2
 * @param[in] dim3 - Interpolation along dimension 3
 * @param[in] dim4 - Interpolation along dimension 4
 * @return Tensor5 - interpweights
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1, PolyOrder4 + 1>,
     5>
interpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
              const std::vector<FixedLagrange<PolyOrder1>>& dim1,
              const std::vector<FixedLagrange<PolyOrder2>>& dim2,
              const std::vector<FixedLagrange<PolyOrder3>>& dim3,
              const std::vector<FixedLagrange<PolyOrder4>>& dim4) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                 PolyOrder3 + 1, PolyOrder4 + 1>,
       5>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size());
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t DerivativeDim>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1>
dinterpweights(const FixedLagrange<PolyOrder0>& dim0,
               const FixedLagrange<PolyOrder1>& dim1,
               const FixedLagrange<PolyOrder2>& dim2,
               const FixedLagrange<PolyOrder3>& dim3,
               const FixedLagrange<PolyOrder4>& dim4) {
  static_assert(DerivativeDim < 5);

  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1, PolyOrder4 + 1>
      out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) =
                (DerivativeDim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                (DerivativeDim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                (DerivativeDim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                (DerivativeDim == 3 ? dim3.dlx[l] : dim3.lx[l]) *
                (DerivativeDim == 4 ? dim4.dlx[m] : dim4.lx[m]);
  return out;
}

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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t DerivativeDim>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1, PolyOrder4 + 1>,
     5>
dinterpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
               const std::vector<FixedLagrange<PolyOrder1>>& dim1,
               const std::vector<FixedLagrange<PolyOrder2>>& dim2,
               const std::vector<FixedLagrange<PolyOrder3>>& dim3,
               const std::vector<FixedLagrange<PolyOrder4>>& dim4) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                 PolyOrder3 + 1, PolyOrder4 + 1>,
       5>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            out(i, j, k, l, m) =
                dinterpweights<PolyOrder0, PolyOrder1, PolyOrder2, PolyOrder3,
                               PolyOrder4, DerivativeDim>(
                    dim0[i], dim1[j], dim2[k], dim3[l], dim4[m]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

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
Numeric interp(const ConstTensor5View& yi, const ConstTensor5View& iw,
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
 * @return Numeric of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, class Tensor5Type>
constexpr Numeric interp(
    const Tensor5Type& yi,
    const FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1>& iw,
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1,
    const FixedLagrange<PolyOrder2>& dim2,
    const FixedLagrange<PolyOrder3>& dim3,
    const FixedLagrange<PolyOrder4>& dim4) {
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
Tensor5 reinterp(const ConstTensor5View& iy, const Grid<Tensor5, 5>& iw,
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
 * @return Tensor5 of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4>
Tensor5 reinterp(
    const ConstTensor5View& iy,
    const Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1,
                         PolyOrder2 + 1, PolyOrder3 + 1, PolyOrder4 + 1>,
               5>& iw,
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1,
    const FixedLagrange<PolyOrder2>& dim2,
    const FixedLagrange<PolyOrder3>& dim3,
    const FixedLagrange<PolyOrder4>& dim4) {
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>
interpweights(const FixedLagrange<PolyOrder0>& dim0,
              const FixedLagrange<PolyOrder1>& dim1,
              const FixedLagrange<PolyOrder2>& dim2,
              const FixedLagrange<PolyOrder3>& dim3,
              const FixedLagrange<PolyOrder4>& dim4,
              const FixedLagrange<PolyOrder5>& dim5) {
  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>
      out;
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>,
     6>
interpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
              const std::vector<FixedLagrange<PolyOrder1>>& dim1,
              const std::vector<FixedLagrange<PolyOrder2>>& dim2,
              const std::vector<FixedLagrange<PolyOrder3>>& dim3,
              const std::vector<FixedLagrange<PolyOrder4>>& dim4,
              const std::vector<FixedLagrange<PolyOrder5>>& dim5) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                 PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>,
       6>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
          dim5.size());
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t DerivativeDim>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>
dinterpweights(const FixedLagrange<PolyOrder0>& dim0,
               const FixedLagrange<PolyOrder1>& dim1,
               const FixedLagrange<PolyOrder2>& dim2,
               const FixedLagrange<PolyOrder3>& dim3,
               const FixedLagrange<PolyOrder4>& dim4,
               const FixedLagrange<PolyOrder5>& dim5) {
  static_assert(DerivativeDim < 6);

  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>
      out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) =
                  (DerivativeDim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                  (DerivativeDim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                  (DerivativeDim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                  (DerivativeDim == 3 ? dim3.dlx[l] : dim3.lx[l]) *
                  (DerivativeDim == 4 ? dim4.dlx[m] : dim4.lx[m]) *
                  (DerivativeDim == 5 ? dim5.dlx[n] : dim5.lx[n]);
  return out;
}

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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t DerivativeDim>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>,
     6>
dinterpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
               const std::vector<FixedLagrange<PolyOrder1>>& dim1,
               const std::vector<FixedLagrange<PolyOrder2>>& dim2,
               const std::vector<FixedLagrange<PolyOrder3>>& dim3,
               const std::vector<FixedLagrange<PolyOrder4>>& dim4,
               const std::vector<FixedLagrange<PolyOrder5>>& dim5) {
  Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                 PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>,
       6>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
          dim5.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              out(i, j, k, l, m, n) =
                  dinterpweights<PolyOrder0, PolyOrder1, PolyOrder2, PolyOrder3,
                                 PolyOrder4, PolyOrder5, DerivativeDim>(
                      dim0[i], dim1[j], dim2[k], dim3[l], dim4[m], dim5[n]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

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
Numeric interp(const ConstTensor6View& yi, const ConstTensor6View& iw,
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
 * @return Numeric of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5, class Tensor6Type>
constexpr Numeric interp(
    const Tensor6Type& yi,
    const FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1>& iw,
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1,
    const FixedLagrange<PolyOrder2>& dim2,
    const FixedLagrange<PolyOrder3>& dim3,
    const FixedLagrange<PolyOrder4>& dim4,
    const FixedLagrange<PolyOrder5>& dim5) {
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
Tensor6 reinterp(const ConstTensor6View& iy, const Grid<Tensor6, 6>& iw,
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
 * @return Tensor6 of interpolated value
 */
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5>
Tensor6 reinterp(const ConstTensor6View& iy,
                 const Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1,
                                      PolyOrder2 + 1, PolyOrder3 + 1,
                                      PolyOrder4 + 1, PolyOrder5 + 1>,
                            6>& iw,
                 const std::vector<FixedLagrange<PolyOrder0>>& dim0,
                 const std::vector<FixedLagrange<PolyOrder1>>& dim1,
                 const std::vector<FixedLagrange<PolyOrder2>>& dim2,
                 const std::vector<FixedLagrange<PolyOrder3>>& dim3,
                 const std::vector<FixedLagrange<PolyOrder4>>& dim4,
                 const std::vector<FixedLagrange<PolyOrder5>>& dim5) {
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t PolyOrder6>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1,
                    PolyOrder6 + 1>
interpweights(const FixedLagrange<PolyOrder0>& dim0,
              const FixedLagrange<PolyOrder1>& dim1,
              const FixedLagrange<PolyOrder2>& dim2,
              const FixedLagrange<PolyOrder3>& dim3,
              const FixedLagrange<PolyOrder4>& dim4,
              const FixedLagrange<PolyOrder5>& dim5,
              const FixedLagrange<PolyOrder6>& dim6) {
  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1, PolyOrder6 + 1>
      out;
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t PolyOrder6>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1, PolyOrder6 + 1>,
     7>
interpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
              const std::vector<FixedLagrange<PolyOrder1>>& dim1,
              const std::vector<FixedLagrange<PolyOrder2>>& dim2,
              const std::vector<FixedLagrange<PolyOrder3>>& dim3,
              const std::vector<FixedLagrange<PolyOrder4>>& dim4,
              const std::vector<FixedLagrange<PolyOrder5>>& dim5,
              const std::vector<FixedLagrange<PolyOrder6>>& dim6) {
  Grid<
      FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1, PolyOrder6 + 1>,
      7>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
          dim5.size(), dim6.size());
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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t PolyOrder6, std::size_t DerivativeDim>
constexpr FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1,
                    PolyOrder6 + 1>
dinterpweights(const FixedLagrange<PolyOrder0>& dim0,
               const FixedLagrange<PolyOrder1>& dim1,
               const FixedLagrange<PolyOrder2>& dim2,
               const FixedLagrange<PolyOrder3>& dim3,
               const FixedLagrange<PolyOrder4>& dim4,
               const FixedLagrange<PolyOrder5>& dim5,
               const FixedLagrange<PolyOrder6>& dim6) {
  static_assert(DerivativeDim < 7);

  FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
            PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1, PolyOrder6 + 1>
      out;
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        for (Index l = 0; l < dim3.size(); l++)
          for (Index m = 0; m < dim4.size(); m++)
            for (Index n = 0; n < dim5.size(); n++)
              for (Index o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) =
                    (DerivativeDim == 0 ? dim0.dlx[i] : dim0.lx[i]) *
                    (DerivativeDim == 1 ? dim1.dlx[j] : dim1.lx[j]) *
                    (DerivativeDim == 2 ? dim2.dlx[k] : dim2.lx[k]) *
                    (DerivativeDim == 3 ? dim3.dlx[l] : dim3.lx[l]) *
                    (DerivativeDim == 4 ? dim4.dlx[m] : dim4.lx[m]) *
                    (DerivativeDim == 5 ? dim5.dlx[n] : dim5.lx[n]) *
                    (DerivativeDim == 6 ? dim6.dlx[o] : dim6.lx[o]);
  return out;
}

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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t PolyOrder6, std::size_t DerivativeDim>
Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
               PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1, PolyOrder6 + 1>,
     7>
dinterpweights(const std::vector<FixedLagrange<PolyOrder0>>& dim0,
               const std::vector<FixedLagrange<PolyOrder1>>& dim1,
               const std::vector<FixedLagrange<PolyOrder2>>& dim2,
               const std::vector<FixedLagrange<PolyOrder3>>& dim3,
               const std::vector<FixedLagrange<PolyOrder4>>& dim4,
               const std::vector<FixedLagrange<PolyOrder5>>& dim5,
               const std::vector<FixedLagrange<PolyOrder6>>& dim6) {
  Grid<
      FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1, PolyOrder6 + 1>,
      7>
      out(dim0.size(), dim1.size(), dim2.size(), dim3.size(), dim4.size(),
          dim5.size(), dim6.size());
  for (std::size_t i = 0; i < dim0.size(); i++)
    for (std::size_t j = 0; j < dim1.size(); j++)
      for (std::size_t k = 0; k < dim2.size(); k++)
        for (std::size_t l = 0; l < dim3.size(); l++)
          for (std::size_t m = 0; m < dim4.size(); m++)
            for (std::size_t n = 0; n < dim5.size(); n++)
              for (std::size_t o = 0; o < dim6.size(); o++)
                out(i, j, k, l, m, n, o) =
                    dinterpweights<PolyOrder0, PolyOrder1, PolyOrder2,
                                   PolyOrder3, PolyOrder4, PolyOrder5,
                                   PolyOrder6, DerivativeDim>(
                        dim0[i], dim1[j], dim2[k], dim3[l], dim4[m], dim5[n],
                        dim6[o]);
  return out;
}

////////////////////////////////////////////////
////////////////////////////////// Interpolation
////////////////////////////////////////////////

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
Numeric interp(const ConstTensor7View& yi, const ConstTensor7View& iw,
               const Lagrange& dim0, const Lagrange& dim1, const Lagrange& dim2,
               const Lagrange& dim3, const Lagrange& dim4, const Lagrange& dim5,
               const Lagrange& dim6);

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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t PolyOrder6, class Tensor7Type>
constexpr Numeric interp(
    const Tensor7Type& yi,
    const FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1, PolyOrder2 + 1,
                    PolyOrder3 + 1, PolyOrder4 + 1, PolyOrder5 + 1,
                    PolyOrder6 + 1>& iw,
    const FixedLagrange<PolyOrder0>& dim0,
    const FixedLagrange<PolyOrder1>& dim1,
    const FixedLagrange<PolyOrder2>& dim2,
    const FixedLagrange<PolyOrder3>& dim3,
    const FixedLagrange<PolyOrder4>& dim4,
    const FixedLagrange<PolyOrder5>& dim5,
    const FixedLagrange<PolyOrder6>& dim6) {
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
Tensor7 reinterp(const ConstTensor7View& iy, const Grid<Tensor7, 7>& iw,
                 const std::vector<Lagrange>& dim0,
                 const std::vector<Lagrange>& dim1,
                 const std::vector<Lagrange>& dim2,
                 const std::vector<Lagrange>& dim3,
                 const std::vector<Lagrange>& dim4,
                 const std::vector<Lagrange>& dim5,
                 const std::vector<Lagrange>& dim6);

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
template <std::size_t PolyOrder0, std::size_t PolyOrder1,
          std::size_t PolyOrder2, std::size_t PolyOrder3,
          std::size_t PolyOrder4, std::size_t PolyOrder5,
          std::size_t PolyOrder6>
Tensor7 reinterp(
    const ConstTensor7View& iy,
    const Grid<FixedGrid<Numeric, PolyOrder0 + 1, PolyOrder1 + 1,
                         PolyOrder2 + 1, PolyOrder3 + 1, PolyOrder4 + 1,
                         PolyOrder5 + 1, PolyOrder6 + 1>,
               7>& iw,
    const std::vector<FixedLagrange<PolyOrder0>>& dim0,
    const std::vector<FixedLagrange<PolyOrder1>>& dim1,
    const std::vector<FixedLagrange<PolyOrder2>>& dim2,
    const std::vector<FixedLagrange<PolyOrder3>>& dim3,
    const std::vector<FixedLagrange<PolyOrder4>>& dim4,
    const std::vector<FixedLagrange<PolyOrder5>>& dim5,
    const std::vector<FixedLagrange<PolyOrder6>>& dim6) {
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

#endif  // interpolation_lagrange_h
