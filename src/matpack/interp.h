#pragma once

#include "array.h"
#include "arts_conversions.h"
#include "artstime.h"
#include "debug.h"
#include "enums.h"
#include "grids.h"
#include "interpolation.h"
#include "matpack_concepts.h"
#include "matpack_constexpr.h"
#include "matpack_data.h"
#include "matpack_iter.h"
#include "nonstd.h"

#include <__concepts/same_as.h>
#include <__iterator/concepts.h>
#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <vector>

namespace my_interp {
//! A concept for a type that CAN be sorted and can thus be a grid
template <typename T>
concept sortable_grid_t = matpack::rank<T>() == 1 and requires (T a) {
  { a[0] } -> matpack::arithmetic;  // Has access operator returning an arithmetic
  { a.size() } -> matpack::integral;  // Has a size
  { a.begin() } -> std::random_access_iterator;  // Can iterate
  { a.end() } -> std::random_access_iterator;  // Can iterate
};

/*! Cycle once back through a list
 *
 * @param[in] n Index in a list 0 <= n < 2*N
 * @param[in] N Index size of a list
 * @return n - N if n >= N else n
 */
constexpr Index cycler(const Index n, const Index N) noexcept { return n >= N ? n - N : n; }

//! A helper class to select bounds
enum class cycle_limit {lower, upper};

//! A helper class to denote a longitude cycle
template <cycle_limit lim> requires (lim == cycle_limit::lower or lim == cycle_limit::upper)
struct longitude_cycle {
  static constexpr Numeric bound = lim == cycle_limit::upper ? 180 : -180;
};

//! A helper class to denote no cyclic limit
template <cycle_limit lim> requires (lim == cycle_limit::lower or lim == cycle_limit::upper)
struct no_cycle {
  static constexpr Numeric bound = lim == cycle_limit::upper ? std::numeric_limits<Numeric>::infinity() : -std::numeric_limits<Numeric>::infinity();
};

//! A helper concept to denote a cyclic limit
template <typename T>
concept cyclic_limit =
requires(T a) {
  { a.bound } -> matpack::arithmetic;
};

template <template <cycle_limit lim> class Limit>
constexpr bool test_cyclic_limit()
  requires(cyclic_limit<Limit<cycle_limit::lower>> and
           cyclic_limit<Limit<cycle_limit::upper>>)
{
  return Limit<cycle_limit::lower>::bound < Limit<cycle_limit::upper>::bound;
}

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
template <template <cycle_limit lim> class Limit>
constexpr Numeric cyclic_clamp(Numeric x) noexcept
  requires(test_cyclic_limit<Limit>())
{
  constexpr auto lb = Limit<cycle_limit::lower>::bound;
  constexpr auto ub = Limit<cycle_limit::upper>::bound;
  constexpr auto db = ub - lb;

  while (x < lb) x += db;
  while (x >= ub) x -= db;
  return x;
}

/*! Find the absolute minimum in a cycle
 *
 * @param[in] x A position relative to a cycle
 * @param[in] xlim [Lower, Upper) bound of cycle
 * @return x-dx, x, or x+dx, whichever absolute is the smallest, where dx=Upper-Lower
 */
template <template <cycle_limit lim> class Limit>
constexpr Numeric min_cyclic(const Numeric x) noexcept
  requires(test_cyclic_limit<Limit>()) {
  constexpr auto lb = Limit<cycle_limit::lower>::bound;
  constexpr auto ub = Limit<cycle_limit::upper>::bound;
  constexpr auto dx = ub - lb;

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
template <template <cycle_limit lim> class Limit>
constexpr bool full_cycle(const Numeric xo,
                          const Numeric xn) noexcept
  requires(test_cyclic_limit<Limit>()) {
  constexpr auto lb = Limit<cycle_limit::lower>::bound;
  constexpr auto ub = Limit<cycle_limit::upper>::bound;

  return (xo == lb and xn == ub) or
         (xo == ub and xn == lb);
}

/*! Find an estimation of the start position in a linearly
 * separated grid (useful as a start position esitmated
 * 
 * @param[in] x The position
 * @param[in] xvec The grid
 * @return Estimated position of x [0, xvec.size())
*/
template <sortable_grid_t Vec>
constexpr Index start_pos_finder(const Numeric x, const Vec& xvec) noexcept {
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
 * @param[in] ascending The sorting is ascending (1, 2, 3...)
 * @param[in] cycle The size of a cycle (optional; increasing first->second)
 */
template <template <cycle_limit lim> class Limit=no_cycle, sortable_grid_t Vec=Vector>
constexpr Index pos_finder(const Index pos0, const Numeric x, const Vec& xi,
                           const Index polyorder,
                           const bool ascending) noexcept requires(test_cyclic_limit<Limit>()) {
  constexpr Numeric lb = Limit<cycle_limit::lower>::bound;
  constexpr Numeric ub = Limit<cycle_limit::upper>::bound;
  constexpr bool cyclic = (ub - lb) < std::numeric_limits<Numeric>::infinity();

  if constexpr (cyclic) if (x < lb or x > ub) {
    // We are below or above the cycle so we must redo calculations after clamping x to the cycle
    return pos_finder<no_cycle>(pos0, cyclic_clamp<Limit>(x), xi, polyorder, ascending);
  }
  
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
    if constexpr (cyclic)  // Max N since we can overstep bounds
      return IMIN(IMAX(p0 - polyorder / 2, 0), N);
    else // Max N-polyorder since we cannot overstep bounds
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
template <GridType type, template <cycle_limit lim> class Limit, sortable_grid_t Vec>
constexpr Numeric l(const Index p0, const Index n, const Numeric x,
                    const Vec& xi, const Index j) noexcept requires(test_cyclic_limit<Limit>()) {
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
        const Numeric x_val = cyclic_clamp<Limit>(x);
        if (full_cycle<Limit>(xi[0], xi[N-1])) {
          // We ignore the last point in full cycles
          if (j_pos == N - 1)
            return 0;
          if (m_pos not_eq N - 1)
            val *= min_cyclic<Limit>(x_val - xi[m_pos]) / min_cyclic(xi[j_pos] - xi[m_pos]);
        } else {
          val *= min_cyclic<Limit>(x_val - xi[m_pos]) / min_cyclic(xi[j_pos] - xi[m_pos]);
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
template <GridType type, template <cycle_limit lim> class Limit, sortable_grid_t Vec,
          class LagrangeVectorType>
constexpr double dl_dval(
    const Index p0, const Index n, const Numeric x, const Vec& xi,
    [[maybe_unused]] const LagrangeVectorType& li, const Index j, const Index i) noexcept requires(test_cyclic_limit<Limit>()) {
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
    const Numeric x_val = cyclic_clamp<Limit>(x);
    if (full_cycle<Limit>(xi[0], xi[N-1])) {
      // We ignore the last point in full cycles
      if (i_pos == N - 1 or j_pos == N - 1)
        return 0;
      if (x_val not_eq xi[i_pos])
        // A simple case when x is not on the grid
        return li[j] / min_cyclic<Limit>(x_val - xi[i_pos]);
      
      // We have to resort to the full recalculations
      Numeric val = 1.0 / min_cyclic<Limit>(xi[j_pos] - xi[i_pos]);
      for (Index m = 0; m < n; m++) {
        if (m not_eq j and m not_eq i) {
          const Index m_pos = cycler(m + p0, N);
          if (m_pos not_eq N - 1)  {
            val *= min_cyclic<Limit>(x_val - xi[m_pos]) / min_cyclic<Limit>(xi[j_pos] - xi[m_pos]);
          }
        }
      }
      return val;
    }
    if (x_val not_eq xi[i_pos])
      // A simple case when x is not on the grid
      return li[j] / min_cyclic<Limit>(x_val - xi[i_pos]);
    // We have to resort to the full recalculations
    Numeric val = 1.0 / min_cyclic<Limit>(xi[j_pos] - xi[i_pos]);
    for (Index m = 0; m < n; m++) {
      if (m not_eq j and m not_eq i) {
        const Index m_pos = cycler(m + p0, N);
        val *= min_cyclic<Limit>(x_val - xi[m_pos]) / min_cyclic<Limit>(xi[j_pos] - xi[m_pos]);
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
template <GridType type, template <cycle_limit lim> class Limit, sortable_grid_t Vec,
          class LagrangeVectorType>
constexpr Numeric dl(const Index p0, const Index n, const Numeric x,
                     const Vec& xi, const LagrangeVectorType& li,
                     const Index j) noexcept requires(test_cyclic_limit<Limit>()) {
  Numeric dval = 0.0;
  for (Index i = 0; i < n; i++) {
    if (i not_eq j) {
      dval += dl_dval<type, Limit>(p0, n, x, xi, li, j, i);
    }
  }
  return dval;
}

//! Checks whether the Sorted Vector is ascending or not by checking its first two elements
template <sortable_grid_t Vec>
constexpr bool is_ascending(const Vec& xi) ARTS_NOEXCEPT {
  if (xi.size() > 1)
    return xi[0] < xi[1];
  return false;
}

/*! Checks the interpolation grid and throws if it is bad
 * 
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] x {Min new x, Max new x}
 * @param[in] extrapol Level of extrapolation
 */
template <GridType type, template <cycle_limit lim> class Limit,
          sortable_grid_t Vec>
void check_lagrange_interpolation(
    [[maybe_unused]] const Vec &xi,
    [[maybe_unused]] const Index polyorder = 1,
    [[maybe_unused]] const std::pair<Numeric, Numeric> x =
        {std::numeric_limits<Numeric>::infinity(),
         -std::numeric_limits<Numeric>::infinity()},
    [[maybe_unused]] const Numeric extrapol = 0.5)
  requires(test_cyclic_limit<Limit>())
{
  constexpr Numeric lb = Limit<cycle_limit::lower>::bound;
  constexpr Numeric ub = Limit<cycle_limit::upper>::bound;

  const Index n = Index(xi.size());

  ARTS_USER_ERROR_IF(polyorder >= n, "Interpolation setup has failed!\n"
                                     "\tRequesting greater interpolation order "
                                     "than possible with given input grid")
  ARTS_USER_ERROR_IF(type == GridType::Cyclic and lb >= ub,
                     "Interpolation setup has failed!\n"
                     "\tBad cycle, must be [first, second)")
  if constexpr (GridType::Cyclic not_eq type) {
    if (polyorder and extrapol > 0) {
      const bool ascending = is_ascending(xi);
      const Numeric xmin =
          ascending ? xi[0] - extrapol * nonstd::abs(xi[1] - xi[0])
                    : xi[n - 1] - extrapol * nonstd::abs(xi[n - 2] - xi[n - 1]);
      const Numeric xmax =
          ascending ? xi[n - 1] + extrapol * nonstd::abs(xi[n - 2] - xi[n - 1])
                    : xi[0] + extrapol * nonstd::abs(xi[1] - xi[0]);
      ARTS_USER_ERROR_IF(x.first < xmin or x.second > xmax,
                         "Interpolation setup has failed!\n"
                         "\tThe new grid has limits: ",
                         x.first, ' ', x.second, '\n',
                         "\tThe old grid has limits: ", xmin, ' ', xmax)
    }
  }
}

/*! Checks the interpolation grid and throws if it is bad
 * 
 * @param[in] xi Old grid positions
 * @param[in] polyorder Polynominal degree
 * @param[in] x New grid position
 * @param[in] extrapol Level of extrapolation
 */
template <GridType type, template <cycle_limit lim> class Limit, sortable_grid_t Vec> constexpr
void check_lagrange_interpolation([[maybe_unused]] const Vec& xi,
                                  [[maybe_unused]] const Index polyorder,
                                  [[maybe_unused]] const Numeric x,
                                  [[maybe_unused]] const Numeric extrapol = 0.5)
  requires(test_cyclic_limit<Limit>()) {
  check_lagrange_interpolation<type, Limit>(xi, polyorder, {x, x}, extrapol);
}

//! Completely empty struct that may store as 0 bytes when used with [[no_unique_address]] on most compilers
struct Empty {
  [[nodiscard]] static constexpr Index size() noexcept {return 0;}
  [[nodiscard]] static constexpr std::nullptr_t begin() noexcept {return nullptr;}
  [[nodiscard]] static constexpr std::nullptr_t end() noexcept {return nullptr;}
  [[nodiscard]] constexpr Numeric operator[](Index) noexcept {return 0;}
};

/*! A Lagrange interpolation computer */
template <Index PolyOrder=-1,
bool do_derivs=false, GridType type=GridType::Standard,
template <cycle_limit lim> class Limit=no_cycle> requires(test_cyclic_limit<Limit>())
struct Lagrange {
  static constexpr bool runtime_polyorder() noexcept {return PolyOrder < 0;}
  static constexpr bool has_derivatives() noexcept {return do_derivs;}

  //! std::vector if runtime_polyorder, otherwise std::array
  using lx_type = std::conditional_t<runtime_polyorder(), std::vector<Numeric>, std::array<Numeric, PolyOrder + 1>>;

  //! std::vector if runtime_polyorder and has_derivatives, std::array if has_derivatives, otherwise Empty
  using dlx_type = std::conditional_t<has_derivatives(), std::conditional_t<runtime_polyorder(), std::vector<Numeric>, std::array<Numeric, PolyOrder + 1>>, Empty>;

  /*! The first position of the Lagrange interpolation grid */
  Index pos{0};
  
  /*! The Lagrange interpolation weights at each point */
  lx_type lx{};
  
  /*! The Lagrange interpolation weights derivatives at each point */
  [[no_unique_address]] dlx_type dlx{};

  /* Number of weights */
  static constexpr Index size() noexcept requires(not runtime_polyorder()) { return PolyOrder + 1; }

  /* Number of weights */
  [[nodiscard]] constexpr Index size() const noexcept requires(runtime_polyorder()) { return lx.size(); }

  /*! Get the index position in the original grid when applying an offset
   * 
   *  Note that this cycles the position iff the grid-type is cyclic
   *
   * @param offset The offset of the new index
   * @param maxsize The maximum size of the original grid; only used by the cyclic path
   * @return constexpr Index The position of a point in an original grid
   */
  [[nodiscard]] constexpr Index index_pos(Index offset, Index maxsize [[maybe_unused]]) const noexcept {
    if constexpr (type == GridType::Cyclic) return cycler(pos+offset, maxsize);
    else return pos+offset;
  }
  
  //! Enusre that the move constructor exists
  constexpr Lagrange() = default;
  constexpr Lagrange(const Lagrange& l) = default;
  constexpr Lagrange(Lagrange&& l) noexcept = default;
  constexpr Lagrange& operator=(const Lagrange& l) = default;
  constexpr Lagrange& operator=(Lagrange&& l) noexcept = default;

  /*! Standard initializer from Vector-types for runtime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] polyorder Polynominal degree
   */
  template <sortable_grid_t Vec>
  constexpr Lagrange(const Index p0, const Numeric x,
                     const Vec &xi, Index polyorder=1) noexcept
    requires(runtime_polyorder())
      : pos(pos_finder<Limit>(p0, x, xi, polyorder,
                              xi.size() > 1 ? xi[0] < xi[1] : false)),
        lx(polyorder + 1), dlx(do_derivs ? polyorder + 1 : 0) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Standard initializer from Vector-types for compiletime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   */
  template <sortable_grid_t Vec>
  constexpr Lagrange(const Index p0, const Numeric x,
                     const Vec &xi) noexcept
    requires(not runtime_polyorder())
      : pos(pos_finder<Limit>(p0, x, xi, PolyOrder,
                              xi.size() > 1 ? xi[0] < xi[1] : false)) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Friendly stream operator */
  friend std::ostream &operator<<(std::ostream &os, const Lagrange &l) {
    os << "Lagrange interpolation ";
    if constexpr (not runtime_polyorder())
      os << "of constant ";
    else
      os << "of runtime ";
    os << "polynominal order: " << (l.size() - 1) << '\n';

    os << "Grid type is: " << std::quoted(var_string(type));
    if constexpr (type == GridType::Cyclic)
      os << " in range [" << Limit<cycle_limit::lower>::bound << ", "
         << Limit<cycle_limit::upper>::bound << ')';
    os << '\n';

    os << "pos: " << l.pos << '\n';

    os << "weights lx: ";
    for (auto x : l.lx)
      os << ' ' << x;

    if constexpr (do_derivs) {
      os << "\nweights dlx:";
      for (auto x : l.dlx)
        os << ' ' << x;
    }

    return os;
  }

 private:
   /*! Finds lx
    *
    * @param[in] x New grid position
    * @param[in] xi Old grid positions
    */
   template <sortable_grid_t Vec>
   constexpr void lx_finder(const Numeric x,
                            const Vec &xi) noexcept {
    for (Index j = 0; j < size(); j++)
      lx[j] = l<type, Limit>(pos, size(), x, xi, j);
   }

   /*! Finds dlx
    *
    * @param[in] x New grid position
    * @param[in] xi Old grid positions
    */
   template <sortable_grid_t Vec>
   constexpr void dlx_finder(const Numeric x, const Vec &xi) noexcept {
    if constexpr (do_derivs) {
      for (Index j = 0; j < size(); j++)
        dlx[j] = dl<type, Limit>(pos, size(), x, xi, lx, j);
    }
   }
};  // Lagrange

namespace internal {
/**  Checks if the Lagrange type has runtime or compile time polynominal order
 * 
 * @tparam T A type that has a static bool-convertible runtime_polyorder() method
 * @return true If runtime polynominal order
 * @return false Otherwise
 */
template <typename T>
constexpr bool runtime_polyorder() {
  return std::remove_cvref_t<T>::runtime_polyorder();
}

/**  Get the compile-time size of the object
 * 
 * @tparam T A type that has a static size() method returning an Index
 * @return constexpr Index The size of the object at compile time
 */
template <typename T>
constexpr Index compile_time_size() {
  return std::remove_cvref_t<T>::size();
}

/**  Checks if the Lagrange type has runtime or compile time polynominal order
 * 
 * @tparam T A type that has a static bool-convertible has_derivatives() method
 * @return true If derivatives are computed
 * @return false Otherwise
 */
template <typename T>
constexpr bool has_derivatives() {
  return std::remove_cvref_t<T>::has_derivatives();
}
}  // namespace internal

//! Test to make sure that the type can be used as a lagrange key type
template <typename T>
concept lagrange_type = requires(const T a) {
  /* Normal lagrange weight */
  { a.pos } -> matpack::integral;
  { a.size() } -> matpack::integral;
  { a.index_pos(0, 0) } -> matpack::integral;
} and std::same_as<decltype(internal::has_derivatives<T>()), bool> and
      std::same_as<decltype(internal::runtime_polyorder<T>()), bool> and
requires (T a) {
   /* The lagrange weights can be accessed */
  { a.lx.size() } -> matpack::integral;  // Size is a size-type
  { a.lx[0] } -> matpack::arithmetic;  // May contain some values
} and
requires (T a) {
  { a.dlx.size() } -> matpack::integral;  // Size is a size-type
  { a.dlx[0] } -> matpack::arithmetic;  // May contain some values
};

/** Get a list of Lagrange types for use with reinterp
 *
 * This is the runtime polyorder version
 * 
 * @tparam T A Lagrage type
 * @tparam NewVec A grid-type that is sortable
 * @tparam Vec A grid-type that is sortable
 * @param[in] xs The new grid, does not have to be sorted
 * @param[in] xi The original grid, should be sorted
 * @param[in] order The runtime polynominal order of the interpolation
 * @return An array of Lagrange types
 */
template <lagrange_type T, sortable_grid_t NewVec, sortable_grid_t Vec>
constexpr Array<T> lagrange_interpolation_list(NewVec &&xs,
                                               Vec &&xi, Index order=1) requires(std::remove_cvref_t<T>::runtime_polyorder()) {
   Array<T> out;
   out.reserve(xs.size());

   for (auto &x : xs) {
    if (out.size())
      out.emplace_back(out.back().pos, x, xi, order);
    else
      out.emplace_back(start_pos_finder(x, xi), x, xi, order);
   }
   
   return out;
}

template <lagrange_type T, sortable_grid_t NewVec, sortable_grid_t Vec>
constexpr Array<T> lagrange_interpolation_list(NewVec &&xs, Vec &&xi)
  requires(not std::remove_cvref_t<T>::runtime_polyorder()) {
   Array<T> out;
   out.reserve(xs.size());

   for (auto &x : xs) {
    if (out.size())
      out.emplace_back(out.back().pos, x, xi);
    else
      out.emplace_back(start_pos_finder(x, xi), x, xi);
   }
   
   return out;
}

namespace internal {
/** Helper struct that selects dlx or lx from the inputs
 * 
 * @tparam dlx The index at which the derivative is computed
 * @tparam T A list of types that are lagrange_type compatible
 */
template <Index dlx, lagrange_type... T> struct select_derivative {
private:
  /** Selection mechanism
   * 
   * @tparam selection_flag If true, return the derivative, otherwise the pure weights 
   * @tparam lag_t A Lagrange type
   * @param lag The Lagrange value
   * @return constexpr const auto& The derivative or pure weights 
   */
  template <bool selection_flag, lagrange_type lag_t>
  static constexpr const auto &one_by_one(const lag_t &lag) {
    if constexpr (selection_flag) {
      static_assert(has_derivatives<lag_t>(), "Your type lacks derivatives");
      return lag.dlx;
    } else {
      return lag.lx;
    }
  }

public:
  /** Only interface to the helper struct
   * 
   * @tparam inds A compile-time list of Index generated by std::make_integet_sequence<sizeof...(T)>()
   * @param all A list of Lagrange values
   * @return constexpr auto An iterable list
   */
  template <Index... inds>
  static constexpr auto as_elemwise(std::integer_sequence<Index, inds...>,
                                    T &&...all) requires (sizeof...(inds) == sizeof...(T)) {
    return matpack::elemwise{one_by_one<inds == dlx>(all)...};
  }
};
} // namespace internal

template<lagrange_type... lags, Index N = sizeof...(lags)>
constexpr auto interpweights2(lags&&... lag) requires (N > 0) {
  if constexpr (N > 1) {
    if constexpr ((internal::runtime_polyorder<lags>() or ...)) {
      const auto in = matpack::elemwise{lag.lx...};
      matpack::matpack_data<Numeric, N> out(lag.size()...);

      std::transform(in.begin(), in.end(), out.elem_begin(), [](auto&& v){
        return std::apply([](auto&&... x){return (x * ...);}, v);
      });

      return out;
    } else {
      const auto in = matpack::elemwise{lag.lx...};
      matpack::matpack_constant_data<Numeric, internal::compile_time_size<lags>()...> out{};

      std::transform(in.begin(), in.end(), out.elem_begin(), [](auto&& v){
        return std::apply([](auto&&... x){return (x * ...);}, v);
      });
      
      return out;
    }
  } else {
    return std::get<0>(std::tuple{lag...}).lx;
  }
}

template<Index dlx, lagrange_type... lags, Index N = sizeof...(lags)>
constexpr auto dinterpweights2(lags&&... lag) requires (N > 0 and dlx >= 0 and dlx < N) {
  if constexpr (N > 1) {
    if constexpr ((internal::runtime_polyorder<lags>() or ...)) {
      const auto in = internal::select_derivative<dlx, lags...>::as_elemwise(std::make_integer_sequence<Index, N>{}, lag...);
      matpack::matpack_data<Numeric, N> out(lag.size()...);

      std::transform(in.begin(), in.end(), out.elem_begin(), [](auto&& v){
        return std::apply([](auto&&... x){return (x * ...);}, v);
      });

      return out;
    } else {
      const auto in = internal::select_derivative<dlx, lags...>::as_elemwise(std::make_integer_sequence<Index, N>{}, lag...);
      matpack::matpack_constant_data<Numeric, internal::compile_time_size<lags>()...> out{};

      std::transform(in.begin(), in.end(), out.elem_begin(), [](auto&& v){
        return std::apply([](auto&&... x){return (x * ...);}, v);
      });

      return out;
    }
  } else {
    return std::get<0>(std::tuple{lag...}).dlx;
  }
}

//! Test to make sure that the type is some list of a type that can be used as a lagrange key type
template <typename T>
concept list_of_lagrange_type = requires(T a) {
  { a.size() } -> matpack::integral;
  { a[0] } -> lagrange_type;
};

template <list_of_lagrange_type... list_lags, Index N = sizeof...(list_lags)>
constexpr auto interpweights2(list_lags &&...lags)
  requires(N > 0)
{
  using internal_type =
      std::remove_cvref_t<decltype(interpweights2(lags[0]...))>;

  const auto in = matpack::elemwise{lags...};
  matpack::matpack_data<internal_type, N> out{
      static_cast<Index>(lags.size())...};

  std::transform(in.begin(), in.end(), out.elem_begin(), [](auto &&v) {
    return std::apply([](auto &&...lag) { return interpweights2(lag...); }, v);
  });

  return out;
}

template <Index dlx, list_of_lagrange_type... list_lags,
          Index N = sizeof...(list_lags)>
constexpr auto dinterpweights2(list_lags &&...lags)
  requires(N > 0 and dlx >= 0 and dlx < N)
{
  using internal_type =
      std::remove_cvref_t<decltype(dinterpweights2<dlx>(lags[0]...))>;

  const auto in = matpack::elemwise{lags...};
  matpack::matpack_data<internal_type, N> out{static_cast<Index>(lags.size())...};

  std::transform(in.begin(), in.end(), out.elem_begin(), [](auto &&v) {
    return std::apply([](auto &&...lag) { return dinterpweights2<dlx>(lag...); }, v);
  });

  return out;
}

template <typename FieldType, typename InterpWeights, lagrange_type... lags>
constexpr auto interp2(FieldType &&field, InterpWeights &&iw, lags &&...lag) {
  typename std::remove_cvref_t<FieldType>::value_type out{0};

  //! Deal with the cyclicity by a applying wrapping lambda
  auto val_fn = std::apply(
      [&](auto &&...maxsize) {
        return [&](auto &&...offset) {
          return iw(offset...) * field(lag.index_pos(offset, maxsize)...);
        };
      },
      matpack::mdshape(field));

  //! Now the cyclicity is wrapped, we can just add the values up
  for (matpack::flat_shape_pos<sizeof...(lags)> pos{matpack::mdshape(iw)};
       pos.pos.front() < pos.shp.front(); ++pos) {
    out += std::apply(val_fn, pos.pos);
  }

  return out;
}

template <typename FieldType, typename ListOfInterpWeights,
          list_of_lagrange_type... lags>
constexpr auto reinterp2(FieldType &&field, ListOfInterpWeights &&iw,
                         lags &&...list_lag) {
  const auto in = matpack::elemwise{list_lag...};
  matpack::matpack_data<typename std::remove_cvref_t<FieldType>::value_type,
                        sizeof...(lags)>
  out(list_lag.size()...);

  std::transform(
      iw.elem_begin(), iw.elem_end(), in.begin(), out.elem_begin(),
      [&](auto &&internal_iw, auto &&lag_t) {
        return std::apply(
            [&](auto &&...lag) { return interp2(field, internal_iw, lag...); },
            lag_t);
      });

  return out;
}
}  // namespace my_interp
