#ifndef interpolation_lagrange_h
#define interpolation_lagrange_h

#include "array.h"
#include "arts_conversions.h"
#include "artstime.h"
#include "debug.h"
#include "enums.h"
#include "grids.h"
#include "matpack_data.h"
#include "nonstd.h"
#include <algorithm>
#include <array>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

namespace Interpolation {
/*! Cycle once back through a list
 *
 * @param[in] n Index in a list 0 <= n < 2*N
 * @param[in] N Index size of a list
 * @return n - N if n >= N else n
 */
constexpr Index cycler(const Index n, const Index N) noexcept { return n >= N ? n - N : n; }

/*! Clamp the value within a cycle by recursion
 *
 * Note that your compiler will have a hard limit on
 * recursion.  If this is reached, this function will
 * fail.  If you even encounter this issue, the code
 * below must be updated
 * 
 * @param[in] x Value to clamp
 * @param[in] xlim [Lower, Upper) bound of cycle
 * @return Value of x in the cycle [Lower, Upper)
 */
constexpr Numeric cyclic_clamp(const Numeric x,
                               const std::pair<Numeric, Numeric> xlim) noexcept {
  if (x < xlim.first)
    return cyclic_clamp(x + xlim.second - xlim.first, xlim);
  if (x >= xlim.second)
    return cyclic_clamp(x - xlim.second + xlim.first, xlim);
  return x;
}

/*! Find the absolute minimum in a cycle
 *
 * @param[in] x A position relative to a cycle
 * @param[in] xlim [Lower, Upper) bound of cycle
 * @return x-dx, x, or x+dx, whichever absolute is the smallest, where dx=Upper-Lower
 */
constexpr Numeric min_cyclic(const Numeric x,
                             const std::pair<Numeric, Numeric> xlim) noexcept {
  const Numeric dx = xlim.second - xlim.first;
  const bool lo = nonstd::abs(x) < nonstd::abs(x - dx);
  const bool hi = nonstd::abs(x) < nonstd::abs(x + dx);
  const bool me = nonstd::abs(x + dx) < nonstd::abs(x - dx);
  return (lo and hi) ? x : me ? x + dx : x - dx;
}

/*! Find if the end points represent a full cycle
 * 
 * @param[in] xo Start position of grid
 * @param[in] xn End position of grid
 * @param[in] xlim [Lower, Upper) bound of cycle
 * @return true if xlim = [xo, xn] or xlim = [xn, xo]
 */
constexpr bool full_cycle(const Numeric xo,
                          const Numeric xn,
                          const std::pair<Numeric, Numeric> xlim) noexcept {
  return (xo == xlim.first and xn == xlim.second) or
         (xo == xlim.second and xn == xlim.first);
}

/*! Find an estimation of the start position in a linearly
 * separated grid (useful as a start position esitmated
 * 
 * @param[in] x The position
 * @param[in] xvec The grid
 * @return Estimated position of x [0, xvec.size())
*/
template <class SortedVectorType>
constexpr Index start_pos_finder(const Numeric x, const SortedVectorType& xvec) noexcept {
  if (const Index n = xvec.size(); n > 1) {
    const Numeric x0 = xvec[    0];
    const Numeric x1 = xvec[n - 1];
    const Numeric frac = (x - x0) / (x1 - x0);
    const auto start_pos = Index(frac * (Numeric)(n - 2));
    return start_pos > 0 ? (start_pos < n ? start_pos : n - 1) : 0;
  }
  return 0;
}

//! Return the maximum of two integer numbers.
/*! 
 T his function is based on a macro from Numerical Rec*eipes. The
 original macro:
 
 static Index imaxarg1, imaxarg2;
 #define IMAX(a,b) (imaxarg1=(a), imaxarg2=(b),(imaxarg1) > (imaxarg2) ? \
 (imaxarg1) : (imaxarg2))
 
 The macro can cause trouble if used in parallel regions, so we use this
 function instead.  
 
 \param a Input a.
 \param b Input b.
 
 \return The maximum of a and b.
 */
constexpr Index IMAX(const Index a, const Index b) noexcept { return a > b ? a : b; }

//! Return the minimum of two integer numbers.
/*! 
 T his function is based on a macro from Numerical Rec*eipes. The
 original macro:
 
 static Index iminarg1, iminarg2;
 #define IMIN(a,b) (iminarg1=(a), iminarg2=(b),(iminarg1) < (iminarg2) ? \
 (iminarg1) : (iminarg2))
 
 The macro can cause trouble if used in parallel regions, so we use this
 function instead.  
 
 \param a Input a.
 \param b Input b.
 
 \return The minimum of a and b.
 */
constexpr Index IMIN(const Index a, const Index b) noexcept { return a < b ? a : b; }

/*! Finds the position of interpolation of x in xi
 * 
 * Will first find the first position with one lower value to the
 * side of the point before adjusting so that x is in the center
 * of a multiple interpolation order curve
 *
 * @param[in] pos0 Estimation of the first position, must be [0, xi.size())
 * @param[in] x Coordinate to find a position for
 * @param[in] xi Original sorted grid
 * @param[in] polyorder Polynominal orders
 * @param[in] cyclic The sorting is cyclic (1, 2, 3, 0.5, 1.5...)
 * @param[in] ascending The sorting is ascending (1, 2, 3...)
 * @param[in] cycle The size of a cycle (optional; increasing first->second)
 */
template <class SortedVectorType>
constexpr Index pos_finder(const Index pos0, const Numeric x, const SortedVectorType& xi,
                           const Index polyorder, const bool cyclic,
                           const bool ascending,
                           const std::pair<Numeric, Numeric> cycle = {
                               -180, 180}) noexcept {
  if (cyclic and (x < cycle.first or x > cycle.second))
    // We are below or above the cycle so we must redo calculations after clamping x to the cycle
    return pos_finder(pos0, cyclic_clamp(x, cycle), xi, polyorder, cyclic, ascending, cycle);
  const Index N = xi.size()-1;
  Index p0=pos0;
  
  // Loops to find the first position with a neighbor larger or smaller
  if (ascending) {
    while (p0 < N and xi[p0] < x) ++p0;
    while (p0 > 0 and xi[p0] > x) --p0;
  } else {
    while (p0 < N and xi[p0] > x) ++p0;
    while (p0 > 0 and xi[p0] < x) --p0;
  }
  
  // Adjustment for higher and lower polynominal orders so that x is in the middle
  if (polyorder) {
    if (cyclic)  // Max N since we can overstep bounds
      return IMIN(IMAX(p0 - polyorder / 2, 0), N);
    // Max N-polyorder since we cannot overstep bounds
    return IMIN(IMAX(p0 - polyorder / 2, 0), N-polyorder);
  }
  // In the nearest neighbor case, we mus choose the closest neighbor
  if (p0 < N and (nonstd::abs(xi[p0] - x) >= nonstd::abs(xi[p0 + 1] - x))) 
    return p0 + 1;
  return p0;
}

/*! Type of Lagrange interpolation weights
 * 
 * Note that Cyclic interpolation has reduced
 * poly-order crossing the cycle if both the
 * start and the end of the cycle are in the
 * original x-data and the polyorder is above 1.  
 * 
 * Note that all derivatives are linear regardless
 * of the type
 */
ENUMCLASS(GridType, char,
          Standard,   /* 1-to-1 interpolation grid */
          Cyclic,     /* Cyclic interpolation grid */
          Log,        /* Natural logarithm interpolation grid */
          Log10,      /* 10-base logarithm interpolation grid */
          Log2,       /* 2-base logarithm interpolation grid */
          SinDeg,     /* Cosine in degrees interpolation grid, grid only defined [-90, 90] */
          SinRad,     /* Cosine in radians interpolation grid, grid only defined [-PI/2, PI/2] */
          CosDeg,     /* Cosine in degrees interpolation grid, grid only defined [0, 180] */
          CosRad      /* Cosine in radians interpolation grid, grid only defined [0,  PI] */
         );

/*! Computes the weights for a given coefficient
 *
 * @param[in] p0 The origin position
 * @param[in] n The number of weights
 * @param[in] x The position for the weights
 * @param[in] xi The sorted vector of values
 * @param[in] j The current coefficient
 * @param[in] cycle The size of a cycle (optional)
 */
template <GridType type, class SortedVectorType>
constexpr Numeric l(const Index p0, const Index n, const Numeric x,
                    const SortedVectorType& xi, const Index j,
                    [[maybe_unused]] const std::pair<Numeric, Numeric> cycle = {
                        -180, 180}) noexcept {
  Numeric val = 1.0;
  for (Index m = 0; m < n; m++) {
    if (m not_eq j) {
      if constexpr (type == GridType::Log) {
        // Natural log weights
        val *= (std::log(x) - std::log(xi[m + p0])) /
               (std::log(xi[j + p0]) - std::log(xi[m + p0]));
      } else if constexpr (type == GridType::Log10) {
        // Common log weights
        val *= (std::log10(x) - std::log10(xi[m + p0])) /
               (std::log10(xi[j + p0]) - std::log10(xi[m + p0]));
      } else if constexpr (type == GridType::Log2) {
        // Binary log weights
        val *= (std::log2(x) - std::log2(xi[m + p0])) /
               (std::log2(xi[j + p0]) - std::log2(xi[m + p0]));
      } else if constexpr (type == GridType::SinDeg) {
        // Sine in degrees weights
        using Conversion::sind;
        val *= (sind(x) - sind(xi[m + p0])) / (sind(xi[j + p0]) - sind(xi[m + p0]));
      } else if constexpr (type == GridType::SinRad) {
        // Sine in radians weights
        using std::sin;
        val *= (sin(x) - sin(xi[m + p0])) / (sin(xi[j + p0]) - sin(xi[m + p0]));
      } else if constexpr (type == GridType::CosDeg) {
        // Cosine in degrees weights (nb. order changed)
        using Conversion::cosd;
        val *= (cosd(xi[m + p0]) - cosd(x)) / (cosd(xi[m + p0]) - cosd(xi[j + p0]));
      } else if constexpr (type == GridType::CosRad) {
        // Cosine in radians weights (nb. order changed)
        using std::cos;
        val *= (cos(xi[m + p0]) - cos(x)) / (cos(xi[m + p0]) - cos(xi[j + p0]));
      } else if constexpr (type == GridType::Standard) {
        // Linear weights, simple and straightforward
        val *= (x - xi[m + p0]) / (xi[j + p0] - xi[m + p0]);
      } else if constexpr (type == GridType::Cyclic) {
        // Cyclic weights
        // We have to ensure that all weights are cyclic (e.g., 355 degrees < -6 degrees)
        const Index N = Index(xi.size());
        const Index m_pos = cycler(m + p0, N);
        const Index j_pos = cycler(j + p0, N);
        const Numeric x_val = cyclic_clamp(x, cycle);
        if (full_cycle(xi[0], xi[N-1], cycle)) {
          // We ignore the last point in full cycles
          if (j_pos == N - 1)
            return 0;
          if (m_pos not_eq N - 1)
            val *= min_cyclic(x_val - xi[m_pos], cycle) / min_cyclic(xi[j_pos] - xi[m_pos], cycle);
        } else {
          val *= min_cyclic(x_val - xi[m_pos], cycle) / min_cyclic(xi[j_pos] - xi[m_pos], cycle);
        }
      }
    }
  }
  return val;
}

/*! Computes the derivatives of the weights for a given coefficient for a given
 * weight
 *
 * If x is on the grid, this is more expensive than if x is not on the grid
 *
 * @param[in] p0 The origin position
 * @param[in] n The number of weights
 * @param[in] x The position for the weights
 * @param[in] xi The sorted vector of values
 * @param[in] li The Lagrange weights
 * @param[in] j The current coefficient
 * @param[in] i The current wight Index
 * @param[in] cycle The size of a cycle (optional)
 */
template <GridType type, class SortedVectorType,
          class LagrangeVectorType>
constexpr double dl_dval(
    const Index p0, const Index n, const Numeric x, const SortedVectorType& xi,
    [[maybe_unused]] const LagrangeVectorType& li, const Index j, const Index i,
    [[maybe_unused]] const std::pair<Numeric, Numeric> cycle) noexcept {
  if constexpr (type == GridType::Standard) {
    // Linear weights, simple and straightforward
    if (x not_eq xi[i + p0]) {
      // A simple case when x is not on the grid
      return li[j] / (x - xi[i + p0]);
    }
    // We have to resort to the full recalculations
    Numeric val = 1.0 / (xi[j + p0] - xi[i + p0]);
    for (Index m = 0; m < n; m++) {
      if (m not_eq j and m not_eq i) {
        val *= (x - xi[m + p0]) / (xi[j + p0] - xi[m + p0]);
      }
    }
    return val;
  } else if constexpr (type == GridType::Cyclic) {
    // Cyclic weights
    // We have to ensure that all weights are cyclic (e.g., 355 degrees < -6 degrees)
    const decltype(i + p0) N = xi.size();
    const Index i_pos = cycler(i + p0, N);
    const Index j_pos = cycler(j + p0, N);
    const Numeric x_val = cyclic_clamp(x, cycle);
    if (full_cycle(xi[0], xi[N-1], cycle)) {
      // We ignore the last point in full cycles
      if (i_pos == N - 1 or j_pos == N - 1)
        return 0;
      if (x_val not_eq xi[i_pos])
        // A simple case when x is not on the grid
        return li[j] / min_cyclic(x_val - xi[i_pos], cycle);
      
      // We have to resort to the full recalculations
      Numeric val = 1.0 / min_cyclic(xi[j_pos] - xi[i_pos], cycle);
      for (Index m = 0; m < n; m++) {
        if (m not_eq j and m not_eq i) {
          const Index m_pos = cycler(m + p0, N);
          if (m_pos not_eq N - 1)  {
            val *= min_cyclic(x_val - xi[m_pos], cycle) / min_cyclic(xi[j_pos] - xi[m_pos], cycle);
          }
        }
      }
      return val;
    }
    if (x_val not_eq xi[i_pos])
      // A simple case when x is not on the grid
      return li[j] / min_cyclic(x_val - xi[i_pos], cycle);
    // We have to resort to the full recalculations
    Numeric val = 1.0 / min_cyclic(xi[j_pos] - xi[i_pos], cycle);
    for (Index m = 0; m < n; m++) {
      if (m not_eq j and m not_eq i) {
        const Index m_pos = cycler(m + p0, N);
        val *= min_cyclic(x_val - xi[m_pos], cycle) / min_cyclic(xi[j_pos] - xi[m_pos], cycle);
      }
    }
    return val;
  } else /*if any other case */ {
    // All other cases have to use full calculations because we need a linear derivative
    Numeric val = 1.0 / (xi[j + p0] - xi[i + p0]);
    for (Index m = 0; m < n; m++) {
      if (m not_eq j and m not_eq i) {
        val *= (x - xi[m + p0]) / (xi[j + p0] - xi[m + p0]);
      }
    }
    return val;
  }
}

/*! Computes the derivatives of the weights for a given coefficient
 *
 * @param[in] p0 The origin position
 * @param[in] n The number of weights
 * @param[in] x The position for the weights
 * @param[in] xi The sorted vector of values
 * @param[in] li The Lagrange weights
 * @param[in] j The current coefficient
 * @param[in] cycle The size of a cycle (optional)
 */
template <GridType type, class SortedVectorType,
          class LagrangeVectorType>
constexpr Numeric dl(const Index p0, const Index n, const Numeric x,
                     const SortedVectorType& xi, const LagrangeVectorType& li,
                     const Index j,
                     const std::pair<Numeric, Numeric> cycle = {-180,
                                                                180}) noexcept {
  Numeric dval = 0.0;
  for (Index i = 0; i < n; i++) {
    if (i not_eq j) {
      dval += dl_dval<type>(p0, n, x, xi, li, j, i, cycle);
    }
  }
  return dval;
}

//! Checks whether the Sorted Vector is ascending or not by checking its first two elements
template <class SortedVectorType>
constexpr bool is_ascending(const SortedVectorType& xi) noexcept {
  if (xi.nelem() > 1)
    return xi[0] < xi[1];
  return false;
}

/*! Checks the interpolation grid and throws if it is bad
 * 
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] x {Min new x, Max new x}
 * @param[in] extrapol Level of extrapolation
 * @param[in] type Type of Lagrange
 * @param[in] cycle Size of a cycle if Cyclic type
 */
template <class SortedVectorType>
void check_lagrange_interpolation([[maybe_unused]] const SortedVectorType& xi,
                                  [[maybe_unused]] const Index polyorder = 1,
                                  [[maybe_unused]] const std::pair<Numeric, Numeric> x={std::numeric_limits<Numeric>::infinity(), -std::numeric_limits<Numeric>::infinity()},
                                  [[maybe_unused]] const Numeric extrapol = 0.5,
                                  [[maybe_unused]] const GridType type = GridType::Standard,
                                  [[maybe_unused]] const std::pair<Numeric, Numeric> cycle = {-180, 180}) {
  const Index n = Index(xi.size());
  
  ARTS_USER_ERROR_IF (polyorder >= n,"Interpolation setup has failed!\n"
      "\tRequesting greater interpolation order than possible with given input grid")
  ARTS_USER_ERROR_IF(type == GridType::Cyclic and cycle.first >= cycle.second,
    "Interpolation setup has failed!\n" "\tBad cycle, must be [first, second)")
  if (polyorder and extrapol > 0 and GridType::Cyclic not_eq type) {
    const bool ascending = is_ascending(xi);
    const Numeric xmin = ascending ? xi[0  ] - extrapol * nonstd::abs(xi[1  ] - xi[0  ]) :
                                     xi[n-1] - extrapol * nonstd::abs(xi[n-2] - xi[n-1]) ;
    const Numeric xmax = ascending ? xi[n-1] + extrapol * nonstd::abs(xi[n-2] - xi[n-1]) :
                                     xi[0  ] + extrapol * nonstd::abs(xi[1  ] - xi[0  ]) ;
    ARTS_USER_ERROR_IF(x.first < xmin or x.second > xmax,
      "Interpolation setup has failed!\n"
      "\tThe new grid has limits: ", x.first, ' ', x.second, '\n',
      "\tThe old grid has limits: ", xmin, ' ', xmax)
  }
}


/*! Checks the interpolation grid and throws if it is bad
 * 
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] x New grid position
 * @param[in] extrapol Level of extrapolation
 * @param[in] type Type of Lagrange
 * @param[in] cycle Size of a cycle if Cyclic type
 */
template <class SortedVectorType>
void check_lagrange_interpolation([[maybe_unused]] const SortedVectorType& xi,
                                  [[maybe_unused]] const Index polyorder,
                                  [[maybe_unused]] const Numeric x,
                                  [[maybe_unused]] const Numeric extrapol = 0.5,
                                  [[maybe_unused]] const GridType type = GridType::Standard,
                                  [[maybe_unused]] const std::pair<Numeric, Numeric> cycle = {-180, 180}) {
  check_lagrange_interpolation(xi, polyorder, {x, x}, extrapol, type, cycle);
}

/** A Lagrange interpolation computer */
struct Lagrange {
  /*! The first position of the Lagrange interpolation grid */
  Index pos{0};
  
  /*! The Lagrange interpolation weights at each point */
  std::vector<Numeric> lx;
  
  /*! The Lagrange interpolation weights derivatives at each point */
  std::vector<Numeric> dlx;

  /* Number of weights */
  [[nodiscard]] Index size() const noexcept { return lx.size(); }
  
  //! Ensure that the move constructor exists
  Lagrange(Lagrange&& l) noexcept : pos(l.pos), lx(std::move(l.lx)), dlx(std::move(l.dlx)) {}
  
  //! Ensure that the copy constructor exists
  Lagrange(const Lagrange& l) = default;
  
  //! Ensure that the move operator exists
  Lagrange& operator=(Lagrange&& l) noexcept {
    pos = l.pos;
    lx = std::move(l.lx);
    dlx = std::move(l.dlx);
    return *this;
  }
  
  //! Ensure that the copy operator exists
  Lagrange& operator=(const Lagrange& l) = default;

  /*! Standard initializer, assumes sorted xi
   * 
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] polyorder Polynominal degree
   * @param[in] do_derivs Compute derivatives?
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  Lagrange(const Index pos0, const Numeric x, const SortedVectorType& xi,
           const Index polyorder = 1,
           const bool do_derivs = false,
           const GridType type = GridType::Standard,
           const std::pair<Numeric, Numeric> cycle = {-180, 180}) noexcept {
    const Index n = xi.size();
    if (n > 1) {
      const Index p = polyorder + 1;
      const bool ascending = is_ascending(xi);
      
      pos = pos_finder(pos0, x, xi, polyorder, type == GridType::Cyclic, ascending, cycle);
      lx = lx_finder(pos, p, x, xi, type, cycle);
      if (do_derivs) dlx = dlx_finder(pos, p, x, xi, lx, type, cycle);
    } else {
      pos = 0;
      lx = std::vector<Numeric>(1, 1);
      dlx = std::vector<Numeric>(1, 0);
    }
  }
  
  //! Default constructor for zero-length elements
  Lagrange() noexcept : lx(1, 1), dlx(1, 0) {}
  
  /*! Friendly stream operator */
  friend std::ostream& operator<<(std::ostream& os, const Lagrange& l) {
    os << "pos: " << l.pos << '\n' << "lx: ";
    for (auto x : l.lx) os << ' ' << x;
    os << '\n' << "dlx: ";
    for (auto x : l.dlx) os << ' ' << x;
    return os << '\n';
  }

 private:
  /*! Finds lx 
   * 
   * @param[in] p0 Estimation of original position, must be [0, xi.size())
   * @param[in] n The number of lx to find
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  static std::vector<Numeric> lx_finder(
      const Index p0, const Index n, const Numeric x,
      const SortedVectorType& xi, const GridType type,
      const std::pair<Numeric, Numeric> cycle) noexcept {
    std::vector<Numeric> out(n);
    switch (type) {
      case GridType::Standard:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Standard>(p0, n, x, xi, j);
        break;
      case GridType::Log:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Log>(p0, n, x, xi, j);
        break;
      case GridType::Log10:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Log10>(p0, n, x, xi, j);
        break;
      case GridType::Log2:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Log2>(p0, n, x, xi, j);
        break;
      case GridType::Cyclic:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Cyclic>(p0, n, x, xi, j, cycle);
        break;
      case GridType::SinDeg:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::SinDeg>(p0, n, x, xi, j);
        break;
      case GridType::SinRad:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::SinRad>(p0, n, x, xi, j);
        break;
      case GridType::CosDeg:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::CosDeg>(p0, n, x, xi, j);
        break;
      case GridType::CosRad:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::CosRad>(p0, n, x, xi, j);
        break;
      case GridType::FINAL: { /* pass */
      }
    }
    return out;
  }

  /*! Finds dlx
   * 
   * @param[in] p0 Estimation of original position, must be [0, xi.size())
   * @param[in] n The number of lx to find
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] li The lx produced by lx_finder
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  static std::vector<Numeric> dlx_finder(
      const Index p0, const Index n, const Numeric x,
      const SortedVectorType& xi, const std::vector<Numeric>& li,
      const GridType type,
      const std::pair<Numeric, Numeric> cycle) noexcept {
    std::vector<Numeric> out(n);
    switch (type) {
      case GridType::Standard:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Standard>(p0, n, x, xi, li, j);
        break;
      case GridType::Log:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Log>(p0, n, x, xi, li, j);
        break;
      case GridType::Log10:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Log10>(p0, n, x, xi, li, j);
        break;
      case GridType::Log2:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Log2>(p0, n, x, xi, li, j);
        break;
      case GridType::Cyclic:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Cyclic>(p0, n, x, xi, li, j, cycle);
        break;
      case GridType::SinDeg:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::SinDeg>(p0, n, x, xi, li, j);
        break;
      case GridType::SinRad:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::SinRad>(p0, n, x, xi, li, j);
        break;
      case GridType::CosDeg:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::CosDeg>(p0, n, x, xi, li, j);
        break;
      case GridType::CosRad:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::CosRad>(p0, n, x, xi, li, j);
        break;
      case GridType::FINAL: { /* pass */
      }
    }
    return out;
  }
};  // Lagrange

/*! A Fixed Lagrange interpolation computer */
template <std::size_t PolyOrder>
struct FixedLagrange {
  /*! The first position of the Lagrange interpolation grid */
  Index pos;
  
  /*! The Lagrange interpolation weights at each point */
  std::array<Numeric, PolyOrder + 1> lx;
  
  /*! The Lagrange interpolation weights derivatives at each point */
  std::array<Numeric, PolyOrder + 1> dlx;

  /* Number of weights */
  static constexpr Index size() noexcept { return PolyOrder + 1; }
  
  //! Enusre that the move constructor exists
  constexpr FixedLagrange(FixedLagrange&& l) noexcept : pos(l.pos), lx(std::move(l.lx)), dlx(std::move(l.dlx)) {}
  
  //! Enusre that the move operator exists
  constexpr FixedLagrange& operator=(FixedLagrange&& l) noexcept {
    pos = l.pos;
    lx = std::move(l.lx);
    dlx = std::move(l.dlx);
    return *this;
  }

  /*! Standard initializer from Vector-types
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] do_derivs Compute derivatives?
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  constexpr FixedLagrange(const Index p0, const Numeric x,
                          const SortedVectorType& xi,
                          const bool do_derivs = false,
                          const GridType type = GridType::Standard,
                          const std::pair<Numeric, Numeric> cycle = {-180, 180}) noexcept
      : pos(pos_finder(p0, x, xi, PolyOrder,
                       type == GridType::Cyclic,
                       xi.size() > 1 ? xi[0] < xi[1] : false, cycle)),
        lx(lx_finder(pos, x, xi, type, cycle)),
        dlx(do_derivs ? dlx_finder(pos, x, xi, lx, type, cycle)
                      : std::array<Numeric, PolyOrder + 1>{}) {}

  /*! Standard initializer from Vector-types
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] do_derivs Compute derivatives?
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  FixedLagrange(const Index p0,
                const Time x,
                const SortedVectorType& xi,
                const bool do_derivs = false,
                const GridType type = GridType::Standard,
                const std::pair<Numeric, Numeric> cycle = {-180, 180}) noexcept {
    Vector xs(static_cast<Index>(xi.size()));
    std::transform(
        xi.begin(), xi.end(), xs.begin(), [](auto& v) { return v.Seconds(); });
    *this = FixedLagrange(p0, x.Seconds(), xs, do_derivs, type, cycle);
  }

  /*! Friendly stream operator */
  friend std::ostream& operator<<(std::ostream& os, const FixedLagrange& l) {
    os << "pos: " << l.pos << '\n' << "lx: ";
    for (auto x : l.lx) os << ' ' << x;
    os << '\n' << "dlx: ";
    for (auto x : l.dlx) os << ' ' << x;
    return os << '\n';
  }

 private:
  /*! Finds lx 
   * 
   * @param[in] p0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  static constexpr std::array<Numeric, PolyOrder + 1> lx_finder(
      const Index p0, const Numeric x, const SortedVectorType& xi,
      const GridType type,
      const std::pair<Numeric, Numeric> cycle) noexcept {
    std::array<Numeric, PolyOrder + 1> out{};
    constexpr Index n = PolyOrder + 1;
    switch (type) {
      case GridType::Standard:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Standard>(p0, n, x, xi, j);
        break;
      case GridType::Log:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Log>(p0, n, x, xi, j);
        break;
      case GridType::Log10:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Log10>(p0, n, x, xi, j);
        break;
      case GridType::Log2:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Log2>(p0, n, x, xi, j);
        break;
      case GridType::Cyclic:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::Cyclic>(p0, n, x, xi, j, cycle);
        break;
      case GridType::SinDeg:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::SinDeg>(p0, n, x, xi, j);
        break;
      case GridType::SinRad:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::SinRad>(p0, n, x, xi, j);
        break;
      case GridType::CosDeg:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::CosDeg>(p0, n, x, xi, j);
        break;
      case GridType::CosRad:
        for (Index j = 0; j < n; j++)
          out[j] = l<GridType::CosRad>(p0, n, x, xi, j);
        break;
      case GridType::FINAL: { /* pass */
      }
    }
    return out;
  }
  
  /*! Finds dlx
   * 
   * @param[in] p0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] li The lx produced by lx_finder
   * @param[in] type Type of Lagrange
   * @param[in] cycle Size of a cycle if Cyclic type
   */
  template <class SortedVectorType>
  static constexpr std::array<Numeric, PolyOrder + 1> dlx_finder(
      const Index p0, const Numeric x, const SortedVectorType& xi,
      const std::array<Numeric, PolyOrder + 1>& li, const GridType type,
      const std::pair<Numeric, Numeric> cycle) noexcept {
    std::array<Numeric, PolyOrder + 1> out{};
    constexpr Index n = PolyOrder + 1;
    switch (type) {
      case GridType::Standard:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Standard>(p0, n, x, xi, li, j);
        break;
      case GridType::Log:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Log>(p0, n, x, xi, li, j);
        break;
      case GridType::Log10:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Log10>(p0, n, x, xi, li, j);
        break;
      case GridType::Log2:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Log2>(p0, n, x, xi, li, j);
        break;
      case GridType::SinDeg:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::SinDeg>(p0, n, x, xi, li, j);
        break;
      case GridType::SinRad:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::SinRad>(p0, n, x, xi, li, j);
        break;
      case GridType::CosDeg:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::CosDeg>(p0, n, x, xi, li, j);
        break;
      case GridType::CosRad:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::CosRad>(p0, n, x, xi, li, j);
        break;
      case GridType::Cyclic:
        for (Index j = 0; j < n; j++)
          out[j] = dl<GridType::Cyclic>(p0, n, x, xi, li, j, cycle);
        break;
      case GridType::FINAL: { /* pass */
      }
    }
    return out;
  }
};  // FixedLagrange

////////////////////////////////////////////////////////
////////////////////// For reinterpreting interpolations
////////////////////////////////////////////////////////

/*! Gets a vector of Lagrange interpolation points
 * 
 * Estimates start position using start_pos_finder
 *
 * @param[in] x New grid positions
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] extrapol Level of extrapolation
 * @param[in] do_derivs Flag to activate derivative calculations
 * @param[in] type Type of Lagrange
 * @param[in] cycle Size of a cycle if Cyclic type
 * @return vector of Lagrange
 */
Array<Lagrange> LagrangeVector(
    const ConstVectorView& x, const ConstVectorView& xi, const Index polyorder,
    const Numeric extrapol=0.5, const bool do_derivs=false, const GridType type=GridType::Standard, const std::pair<Numeric, Numeric> cycle={-180, 180});

/*! Gets a vector of Lagrange interpolation points
 * 
 * Estimates start position using start_pos_finder
 *
 * @param[in] x New grid positions
 * @param[in] xi Old grid positions
 * @param[in] do_derivs Flag to activate derivative calculations
 * @param[in] type Type of Lagrange
 * @param[in] cycle Size of a cycle if Cyclic type
 * @return vector of FixedLagrange
 */
template <std::size_t PolyOrder, class UnsortedVectorType,
          class SortedVectorType>
Array<FixedLagrange<PolyOrder>> FixedLagrangeVector(
    const UnsortedVectorType& xs, const SortedVectorType& xi,
    const bool do_derivs=false, const GridType type=GridType::Standard, const std::pair<Numeric, Numeric> cycle={-180, 180}) {
  Array<FixedLagrange<PolyOrder>> out;
  out.reserve(xs.size());
  bool has_one = false;
  for (auto x : xs) {
    if (has_one) {
      out.emplace_back(out.back().pos, x, xi, do_derivs, type, cycle);
    } else {
      out.emplace_back(start_pos_finder(x, xi), x, xi, do_derivs, type, cycle);
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

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
Vector interpweights(const Lagrange& dim0);

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
Grid<Vector, 1> interpweights(const Array<Lagrange>& dim0);

/*! Interpolation weights for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights
 */
template <std::size_t PolyOrder>
constexpr const std::array<Numeric, PolyOrder + 1>& interpweights(
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
    const Array<FixedLagrange<PolyOrder>>& dim0) {
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
Vector dinterpweights(const Lagrange& dim0);

/*! Interpolation weights derivative for a 1D reduction
 *
 * @param[in] dim0 - Interpolation along dimension 0
 * @return Vector - interpweights derivative along 0th dimension
 */
Grid<Vector, 1> dinterpweights(const Array<Lagrange>& dim0);

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
    const Array<FixedLagrange<PolyOrder>>& dim0) {
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
  const Index I = yi.size();
  for (Index i = 0; i < dim0.size(); i++)
    out += iw[i] * yi[cycler(i + dim0.pos, I)];
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Vector of interpolated value
 */
void reinterp(VectorView out, const ConstVectorView& iy,
              const Grid<Vector, 1>& iw, const Array<Lagrange>& dim0);

/*! Reinterpreting interpolation routine
 *
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @return Vector of interpolated value
 */
Vector reinterp(const ConstVectorView& iy, const Grid<Vector, 1>& iw,
                const Array<Lagrange>& dim0);

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
                const Array<FixedLagrange<PolyOrder>>& dim0) {
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
Grid<Matrix, 2> interpweights(const Array<Lagrange>& dim0,
                              const Array<Lagrange>& dim1);
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
    const Array<FixedLagrange<PolyOrder0>>& dim0,
    const Array<FixedLagrange<PolyOrder1>>& dim1) {
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
Grid<Matrix, 2> dinterpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1, Index dim);

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
    const Array<FixedLagrange<PolyOrder0>>& dim0,
    const Array<FixedLagrange<PolyOrder1>>& dim1) {
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

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @return Matrix of interpolated value
 */
void reinterp(MatrixView out, const ConstMatrixView& iy,
              const Grid<Matrix, 2>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1);

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
                const Array<Lagrange>& dim0,
                const Array<Lagrange>& dim1);

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
    const Array<FixedLagrange<PolyOrder0>>& dim0,
    const Array<FixedLagrange<PolyOrder1>>& dim1) {
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
Grid<Tensor3, 3> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2);

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
interpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
              const Array<FixedLagrange<PolyOrder1>>& dim1,
              const Array<FixedLagrange<PolyOrder2>>& dim2) {
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
Grid<Tensor3, 3> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2, Index dim);

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
dinterpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
               const Array<FixedLagrange<PolyOrder1>>& dim1,
               const Array<FixedLagrange<PolyOrder2>>& dim2) {
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
                         const FixedLagrange<PolyOrder0>& dim0,
                         const FixedLagrange<PolyOrder1>& dim1,
                         const FixedLagrange<PolyOrder2>& dim2) {
  Numeric out(0.0);
  const Index I = yi.npages();
  const Index J = yi.nrows();
  const Index K = yi.ncols();
  for (Index i = 0; i < dim0.size(); i++)
    for (Index j = 0; j < dim1.size(); j++)
      for (Index k = 0; k < dim2.size(); k++)
        out +=
            iw(i, j, k) * yi(cycler(i + dim0.pos, I), cycler(j + dim1.pos, J),
                             cycler(k + dim2.pos, K));
  return out;
}

////////////////////////////////////////////////
/////////////////////////////// Re-Interpolation
////////////////////////////////////////////////

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @return Tensor3 of interpolated value
 */
void reinterp(Tensor3View out, const ConstTensor3View& iy,
              const Grid<Tensor3, 3>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2);

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
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2);

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
    const Array<FixedLagrange<PolyOrder0>>& dim0,
    const Array<FixedLagrange<PolyOrder1>>& dim1,
    const Array<FixedLagrange<PolyOrder2>>& dim2) {
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
Grid<Tensor4, 4> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3);

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
interpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
              const Array<FixedLagrange<PolyOrder1>>& dim1,
              const Array<FixedLagrange<PolyOrder2>>& dim2,
              const Array<FixedLagrange<PolyOrder3>>& dim3) {
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
Grid<Tensor4, 4> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3, Index dim);

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
dinterpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
               const Array<FixedLagrange<PolyOrder1>>& dim1,
               const Array<FixedLagrange<PolyOrder2>>& dim2,
               const Array<FixedLagrange<PolyOrder3>>& dim3) {
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

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
 * @param[in] yi - Original data to squash
 * @param[in] iw - Interpolation weights grid or their derivatives from the
 * Lagrange routines
 * @param[in] dim0 - Lagrange weights along the dimension
 * @param[in] dim1 - Lagrange weights along the dimension
 * @param[in] dim2 - Lagrange weights along the dimension
 * @param[in] dim3 - Lagrange weights along the dimension
 * @return Tensor4 of interpolated value
 */
void reinterp(Tensor4View out, const ConstTensor4View& iy,
              const Grid<Tensor4, 4>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3);

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
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3);

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
                 const Array<FixedLagrange<PolyOrder0>>& dim0,
                 const Array<FixedLagrange<PolyOrder1>>& dim1,
                 const Array<FixedLagrange<PolyOrder2>>& dim2,
                 const Array<FixedLagrange<PolyOrder3>>& dim3) {
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
Grid<Tensor5, 5> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3,
                               const Array<Lagrange>& dim4);

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
interpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
              const Array<FixedLagrange<PolyOrder1>>& dim1,
              const Array<FixedLagrange<PolyOrder2>>& dim2,
              const Array<FixedLagrange<PolyOrder3>>& dim3,
              const Array<FixedLagrange<PolyOrder4>>& dim4) {
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
Grid<Tensor5, 5> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3,
                                const Array<Lagrange>& dim4, Index dim);

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
dinterpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
               const Array<FixedLagrange<PolyOrder1>>& dim1,
               const Array<FixedLagrange<PolyOrder2>>& dim2,
               const Array<FixedLagrange<PolyOrder3>>& dim3,
               const Array<FixedLagrange<PolyOrder4>>& dim4) {
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

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
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
void reinterp(Tensor5View out, const ConstTensor5View& iy,
              const Grid<Tensor5, 5>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3,
              const Array<Lagrange>& dim4);

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
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3,
                 const Array<Lagrange>& dim4);

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
Grid<Tensor6, 6> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3,
                               const Array<Lagrange>& dim4,
                               const Array<Lagrange>& dim5);

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
interpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
              const Array<FixedLagrange<PolyOrder1>>& dim1,
              const Array<FixedLagrange<PolyOrder2>>& dim2,
              const Array<FixedLagrange<PolyOrder3>>& dim3,
              const Array<FixedLagrange<PolyOrder4>>& dim4,
              const Array<FixedLagrange<PolyOrder5>>& dim5) {
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
Grid<Tensor6, 6> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3,
                                const Array<Lagrange>& dim4,
                                const Array<Lagrange>& dim5, Index dim);

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
dinterpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
               const Array<FixedLagrange<PolyOrder1>>& dim1,
               const Array<FixedLagrange<PolyOrder2>>& dim2,
               const Array<FixedLagrange<PolyOrder3>>& dim3,
               const Array<FixedLagrange<PolyOrder4>>& dim4,
               const Array<FixedLagrange<PolyOrder5>>& dim5) {
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

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
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
void reinterp(Tensor6View out, const ConstTensor6View& iy,
              const Grid<Tensor6, 6>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3,
              const Array<Lagrange>& dim4,
              const Array<Lagrange>& dim5);

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
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3,
                 const Array<Lagrange>& dim4,
                 const Array<Lagrange>& dim5);

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
                 const Array<FixedLagrange<PolyOrder0>>& dim0,
                 const Array<FixedLagrange<PolyOrder1>>& dim1,
                 const Array<FixedLagrange<PolyOrder2>>& dim2,
                 const Array<FixedLagrange<PolyOrder3>>& dim3,
                 const Array<FixedLagrange<PolyOrder4>>& dim4,
                 const Array<FixedLagrange<PolyOrder5>>& dim5) {
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
Grid<Tensor7, 7> interpweights(const Array<Lagrange>& dim0,
                               const Array<Lagrange>& dim1,
                               const Array<Lagrange>& dim2,
                               const Array<Lagrange>& dim3,
                               const Array<Lagrange>& dim4,
                               const Array<Lagrange>& dim5,
                               const Array<Lagrange>& dim6);

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
interpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
              const Array<FixedLagrange<PolyOrder1>>& dim1,
              const Array<FixedLagrange<PolyOrder2>>& dim2,
              const Array<FixedLagrange<PolyOrder3>>& dim3,
              const Array<FixedLagrange<PolyOrder4>>& dim4,
              const Array<FixedLagrange<PolyOrder5>>& dim5,
              const Array<FixedLagrange<PolyOrder6>>& dim6) {
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
Grid<Tensor7, 7> dinterpweights(const Array<Lagrange>& dim0,
                                const Array<Lagrange>& dim1,
                                const Array<Lagrange>& dim2,
                                const Array<Lagrange>& dim3,
                                const Array<Lagrange>& dim4,
                                const Array<Lagrange>& dim5,
                                const Array<Lagrange>& dim6, Index dim);

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
dinterpweights(const Array<FixedLagrange<PolyOrder0>>& dim0,
               const Array<FixedLagrange<PolyOrder1>>& dim1,
               const Array<FixedLagrange<PolyOrder2>>& dim2,
               const Array<FixedLagrange<PolyOrder3>>& dim3,
               const Array<FixedLagrange<PolyOrder4>>& dim4,
               const Array<FixedLagrange<PolyOrder5>>& dim5,
               const Array<FixedLagrange<PolyOrder6>>& dim6) {
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

/*! Reinterpreting interpolation routine
 *
 * @param[in,out] out - Reinterpreted field
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
void reinterp(Tensor7View out, const ConstTensor7View& iy,
              const Grid<Tensor7, 7>& iw, const Array<Lagrange>& dim0,
              const Array<Lagrange>& dim1,
              const Array<Lagrange>& dim2,
              const Array<Lagrange>& dim3,
              const Array<Lagrange>& dim4,
              const Array<Lagrange>& dim5,
              const Array<Lagrange>& dim6);

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
                 const Array<Lagrange>& dim0,
                 const Array<Lagrange>& dim1,
                 const Array<Lagrange>& dim2,
                 const Array<Lagrange>& dim3,
                 const Array<Lagrange>& dim4,
                 const Array<Lagrange>& dim5,
                 const Array<Lagrange>& dim6);

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
    const Array<FixedLagrange<PolyOrder0>>& dim0,
    const Array<FixedLagrange<PolyOrder1>>& dim1,
    const Array<FixedLagrange<PolyOrder2>>& dim2,
    const Array<FixedLagrange<PolyOrder3>>& dim3,
    const Array<FixedLagrange<PolyOrder4>>& dim4,
    const Array<FixedLagrange<PolyOrder5>>& dim5,
    const Array<FixedLagrange<PolyOrder6>>& dim6) {
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

/*! The Lagrange interpolation base point */
using LagrangeInterpolation = Interpolation::Lagrange;

/*! A list of Lagrange interpolation points */
using ArrayOfLagrangeInterpolation = Array<LagrangeInterpolation>;

/*! A Vector of Vector for Lagrange Interpolation */
using VectorOfVector = Grid<Vector, 1>;

/*! A Matrix of Matrix for Lagrange Interpolation */
using MatrixOfMatrix = Grid<Matrix, 2>;

/*! A Tensor3 of Tensor3 for Lagrange Interpolation */
using Tensor3OfTensor3 = Grid<Tensor3, 3>;

/*! The fixed Lagrange interpolation base point */
template <std::size_t N>
using FixedLagrangeInterpolation = Interpolation::FixedLagrange<N>;

#endif  // interpolation_lagrange_h
