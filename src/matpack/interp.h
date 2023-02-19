#pragma once

#include "array.h"
#include "debug.h"
#include "enums.h"
#include "nonstd.h"

#include "arts_constants.h"
#include "arts_conversions.h"

#include "matpack_concepts.h"
#include "matpack_constexpr.h"
#include "matpack_data.h"
#include "matpack_iter.h"
#include "matpack_math.h"
#include "matpack_math_extra.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <exception>
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
/*! Cycle once through a list
 *
 * @param[in] n Index in a list -N <= n < 2*N
 * @param[in] N Index size of a list
 * @return n - N if n >= N else n + N if n < 0 else n
 */
constexpr Index cycler(const Index n, const Index N) noexcept { return n >= N ? n - N : n < 0 ? n + N : n; }

//! A helper class to select bounds
enum class cycle_limit {lower, upper};

//! A helper class to denote a longitude cycle
template <cycle_limit lim> requires (lim == cycle_limit::lower or lim == cycle_limit::upper)
struct cycle_m180_p180 {
  static constexpr Numeric bound = lim == cycle_limit::upper ? 180 : -180;
};

//! A helper class to denote a longitude cycle
template <cycle_limit lim> requires (lim == cycle_limit::lower or lim == cycle_limit::upper)
struct cycle_0_p360 {
  static constexpr Numeric bound = lim == cycle_limit::upper ? 360 : 0;
};

//! A helper class to denote a longitude cycle
template <cycle_limit lim> requires (lim == cycle_limit::lower or lim == cycle_limit::upper)
struct cycle_0_p2pi {
  static constexpr Numeric bound = lim == cycle_limit::upper ? Constant::two_pi : 0;
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
template <matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr Index start_pos_finder(const Numeric x, const Vec& xvec) {
  if (const Index n = xvec.size(); n > 1) {
    const Numeric x0 = xvec[    0];
    const Numeric x1 = xvec[n - 1];
    const Numeric frac = (x - x0) / (x1 - x0);
    const auto start_pos = static_cast<Index>(frac * (Numeric)(n - 2));
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
 * @tparam ascending The sorting is ascending (1, 2, 3...)
 * @param[in] pos0 Estimation of the first position, must be [0, xi.size())
 * @param[in] x Coordinate to find a position for
 * @param[in] xi Original sorted grid
 * @param[in] polyorder Polynominal orders
 */
template <bool ascending, template <cycle_limit lim> class Limit=no_cycle, matpack::ranked_matpack_type<Numeric, 1> Vec=Vector>
constexpr Index pos_finder(Index p0, const Numeric x, const Vec& xi,
                           const Index polyorder) requires(test_cyclic_limit<Limit>()) {
  constexpr Numeric lb = Limit<cycle_limit::lower>::bound;
  constexpr Numeric ub = Limit<cycle_limit::upper>::bound;
  constexpr bool cyclic = (ub - lb) < std::numeric_limits<Numeric>::infinity();

  if constexpr (cyclic) if (x < lb or x > ub) {
    // We are below or above the cycle so we must redo calculations after clamping x to the cycle
    return pos_finder<ascending, no_cycle>(p0, cyclic_clamp<Limit>(x), xi, polyorder);
  }
  
  const Index N = xi.size()-1;
  
  // Loops to find the first position with a neighbor larger or smaller
  if constexpr (ascending) {
    while (p0 < N and xi[p0] < x) ++p0;
    while (p0 > 0 and xi[p0] > x) --p0;
  } else {
    while (p0 < N and xi[p0] > x) ++p0;
    while (p0 > 0 and xi[p0] < x) --p0;
  }
  
  // Adjustment for higher and lower polynominal orders so that x is in the middle
  if (polyorder > 0) {
    if constexpr (cyclic)  // Max N since we can overstep bounds
      return std::clamp<Index>(p0 - polyorder / 2, 0, N);
    else // Max N-polyorder since we cannot overstep bounds
      return std::clamp<Index>(p0 - polyorder / 2, 0, N-polyorder);
  }

  // In the nearest neighbor case, we must choose the closest neighbor
  return p0 + (p0 < N and (nonstd::abs(xi[p0] - x) >= nonstd::abs(xi[p0 + 1] - x)));
}

/*! Finds the position of interpolation of x in xi
 * 
 * Will first find the first position with one lower value to the
 * side of the point before adjusting so that x is in the center
 * of a multiple interpolation order curve
 *
 * @tparam ascending The sorting is ascending (1, 2, 3...)
 * @tparam polyorder Polynominal orders
 * @param[in] p0 Estimation of the first position, must be [0, xi.size())
 * @param[in] x Coordinate to find a position for
 * @param[in] xi Original sorted grid
 */
template <bool ascending, Index polyorder, template <cycle_limit lim> class Limit=no_cycle, matpack::ranked_matpack_type<Numeric, 1> Vec=Vector>
constexpr Index pos_finder(Index p0, const Numeric x, const Vec& xi) requires(test_cyclic_limit<Limit>()) {
  constexpr Numeric lb = Limit<cycle_limit::lower>::bound;
  constexpr Numeric ub = Limit<cycle_limit::upper>::bound;
  constexpr bool cyclic = (ub - lb) < std::numeric_limits<Numeric>::infinity();

  // Ensure we are within the limits for cyclic orders
  if constexpr (cyclic) {
   if (x < lb or x > ub) {
      return pos_finder<ascending, polyorder, Limit>(p0, cyclic_clamp<Limit>(x), xi);
    }
  }
  
  /*
  Polynominal Layout
  polyorder   indices-offsets-p0
  0           0
  1           0 1
  2           -1 0 1
  3           -1 0 1 2
  4           -2 -1 0 1 2
  ---
  The output of this function is the left-most index position
  The input of this function is also the left-most index position
  */
  constexpr Index p0_offset = polyorder / 2;

  //! Add the mean offset so the tracing loops put this close to the center
  if constexpr (p0_offset > 0) p0 += p0_offset;

  //! Cycle the position if necessary
  if constexpr (cyclic) p0 = cycler(p0, xi.size());
 
  /*
  Minimum values of the output index
  polyorder linear cyclic
  0         0      0
  1         0      0
  2         1      0
  3         1      0
  4         2      0
  ---
  */
  constexpr Index min_p0_offset = (not cyclic) * p0_offset;

  /*
  Maximum values
  polyorder max-linear max-cyclic
  0         N-1        N-1
  1         N-2        N-1
  2         N-3        N-1
  3         N-4        N-1
  4         N-5        N-1
  ---
  */
  constexpr Index max_p0_offset = 1 + (not cyclic) * polyorder;

  /* The maximum offset also gives the maximum p0 */
  const Index N = xi.size() - max_p0_offset;

  /* Clamp higher order indices to ensure we are in a good state */
  if constexpr (not cyclic and p0_offset > 0) p0 = std::clamp<Index>(p0, 0, N);

  //! The upper index is 1 or 0 if polyorder exist
  constexpr Index upper_offset = polyorder not_eq 0;

  /** Loop and find a position */
  if constexpr (ascending) {
    if (x < xi[p0] and p0 > min_p0_offset) {
      do --p0; while (x < xi[p0] and p0 > min_p0_offset);
    } else if (x >= xi[p0 + upper_offset] and p0 < N) {
      do ++p0; while (x >= xi[p0 + upper_offset] and p0 < N);
    }
  } else /* if constexpr descending */ {
    if (x > xi[p0] and p0 > min_p0_offset) {
      do --p0; while (x > xi[p0] and p0 > min_p0_offset);
    } else if (x <= xi[p0 + upper_offset] and p0 < N) {
      do ++p0; while (x <= xi[p0 + upper_offset] and p0 < N);
    }
  }

  if constexpr (cyclic) return cycler(p0 - p0_offset, N + 1);
  else if constexpr (p0_offset > 0) return p0 - p0_offset;
  else return p0;
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

template <GridType type, template <cycle_limit lim> class Limit, matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr Numeric l_factor(const Numeric x, const Vec& xi, const Index j, const Index m) {
  if constexpr (type == GridType::Log) {
    return (std::log(x) - std::log(xi[m])) /
           (std::log(xi[j]) - std::log(xi[m]));
  } else if constexpr (type == GridType::Log10) {
    return (std::log10(x) - std::log10(xi[m])) /
           (std::log10(xi[j]) - std::log10(xi[m]));
  } else if constexpr (type == GridType::Log2) {
    // Binary log weights
    return (std::log2(x) - std::log2(xi[m])) /
            (std::log2(xi[j]) - std::log2(xi[m]));
  } else if constexpr (type == GridType::SinDeg) {
    // Sine in degrees weights
    using Conversion::sind;
    return (sind(x) - sind(xi[m])) / (sind(xi[j]) - sind(xi[m]));
  } else if constexpr (type == GridType::SinRad) {
    // Sine in radians weights
    using std::sin;
    return (sin(x) - sin(xi[m])) / (sin(xi[j]) - sin(xi[m]));
  } else if constexpr (type == GridType::CosDeg) {
    // Cosine in degrees weights (nb. order changed)
    using Conversion::cosd;
    return (cosd(xi[m]) - cosd(x)) / (cosd(xi[m]) - cosd(xi[j]));
  } else if constexpr (type == GridType::CosRad) {
    // Cosine in radians weights (nb. order changed)
    using std::cos;
    return (cos(xi[m]) - cos(x)) / (cos(xi[m]) - cos(xi[j]));
  } else if constexpr (type == GridType::Standard) {
    // Linear weights, simple and straightforward
    return (x - xi[m]) / (xi[j] - xi[m]);
  } else if constexpr (type == GridType::Cyclic) {
    // Cyclic weights
    // We have to ensure that all weights are cyclic (e.g., 355 degrees < -6 degrees)
    const auto N = static_cast<Index>(xi.size());
    const Index m_pos = cycler(m, N);
    const Index j_pos = cycler(j, N);
    const Numeric x_val = cyclic_clamp<Limit>(x);
    
    // We ignore the last point in full cycles
    if (full_cycle<Limit>(xi[0], xi[N-1])) {
      if (j_pos == N - 1) return 0.0;
      if (m_pos == N - 1) return 1.0;
    }
    return min_cyclic<Limit>(x_val - xi[m_pos]) / min_cyclic<Limit>(xi[j_pos] - xi[m_pos]);
  } else {
    std::terminate();
  }
}

/*! Computes the weights for a given coefficient
 *
 * @param[in] p0 The origin position
 * @param[in] n The number of weights
 * @param[in] x The position for the weights
 * @param[in] xi The sorted vector of values
 * @param[in] j The current coefficient
 * @param[in] cycle The size of a cycle (optional)
 */
template <GridType type, template <cycle_limit lim> class Limit, matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr Numeric l(const Index p0, const Index order, const Numeric x,
                    const Vec& xi, const Index j) requires(test_cyclic_limit<Limit>()) {
  Numeric val = 1.0;
  for (Index m = 0; m < order; m++) {
    val *= l_factor<type, Limit, Vec>(x, xi, j+p0, m+p0+(m>=j));
  }
  return val;
}

/*! Computes the weights for a given coefficient
 *
 * @param[in] p0 The origin position
 * @param[in] n The number of weights
 * @param[in] x The position for the weights
 * @param[in] xi The sorted vector of values
 * @param[in] j The current coefficient
 * @param[in] cycle The size of a cycle (optional)
 */
template <Index order, GridType type, template <cycle_limit lim> class Limit, matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr Numeric l(const Index p0, const Numeric x,
                    const Vec& xi, const Index j) requires(test_cyclic_limit<Limit>()) {
  if constexpr (order == 0) {
    return 1.0;
  } else if constexpr (order == 1) {
    return l_factor<type, Limit, Vec>(x, xi, j+p0, p0+(0>=j));
  } else if constexpr (order == 2) {
    return l_factor<type, Limit, Vec>(x, xi, j+p0, p0+(0>=j)) *
           l_factor<type, Limit, Vec>(x, xi, j+p0, p0+(1>=j));
  } else {
    Numeric val = 1.0;
    for (Index m = 0; m < order; m++) {
      val *= l_factor<type, Limit, Vec>(x, xi, j+p0, m+p0+(m>=j));
    }
    return val;
  }
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
template <GridType type, template <cycle_limit lim> class Limit, matpack::ranked_matpack_type<Numeric, 1> Vec,
          class LagrangeVectorType>
constexpr double dl_dval(
    const Index p0, const Index n, const Numeric x, const Vec& xi,
    [[maybe_unused]] const LagrangeVectorType& li, const Index j, const Index i) requires(test_cyclic_limit<Limit>()) {
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
template <GridType type, template <cycle_limit lim> class Limit, matpack::ranked_matpack_type<Numeric, 1> Vec,
          class LagrangeVectorType>
constexpr Numeric dl(const Index p0, const Index n, const Numeric x,
                     const Vec& xi, const LagrangeVectorType& li,
                     const Index j) requires(test_cyclic_limit<Limit>()) {
  Numeric dval = 0.0;
  for (Index i = 0; i < n; i++) {
    if (i not_eq j) {
      dval += dl_dval<type, Limit>(p0, n, x, xi, li, j, i);
    }
  }
  return dval;
}

//! Checks whether the Sorted Vector is ascending or not by checking its first two elements
template <matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr bool is_ascending(const Vec& xi) {
  return (xi.size() > 1) and xi[0] < xi[1];
}

//! Completely empty struct that may store as 0 bytes when used with [[no_unique_address]] on most compilers
struct Empty {
  constexpr Empty(auto&&...) {}
  [[nodiscard]] static constexpr Index size() noexcept {return 0;}
  [[nodiscard]] static constexpr std::nullptr_t begin() noexcept {return nullptr;}
  [[nodiscard]] static constexpr std::nullptr_t end() noexcept {return nullptr;}
  [[nodiscard]] constexpr Numeric operator[](Index) noexcept {return 0;}
};

//! Ascending flag struct
struct AscendingOrder {constexpr AscendingOrder(auto&&...) noexcept {}};

//! Descending flag struct
struct DescendingOrder {constexpr DescendingOrder(auto&&...) noexcept {}};

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
  [[nodiscard]] constexpr Index size() const requires(runtime_polyorder()) { return lx.size(); }

  /*! Get the index position in the original grid when applying an offset
   * 
   *  Note that this cycles the position iff the grid-type is cyclic
   *
   * @param[in] offset The offset of the new index
   * @param[in] maxsize The maximum size of the original grid; only used by the cyclic path
   * @return constexpr Index The position of a point in an original grid
   */
  [[nodiscard]] constexpr Index index_pos(Index offset, Index maxsize [[maybe_unused]]) const noexcept {
    if constexpr (type == GridType::Cyclic) return cycler(pos+offset, maxsize);
    else return pos+offset;
  }
  
  //! Enusre that the standard constructors and set-operators exist
  constexpr Lagrange(std::size_t order=0) requires(runtime_polyorder()) : lx(order+1, 0), dlx(do_derivs * (order+1), 0) {lx.front() = 1;}
  constexpr Lagrange() requires(not runtime_polyorder()) : lx({1}), dlx({0}) {}
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
  template <matpack::ranked_matpack_type<Numeric, 1> Vec>
  constexpr Lagrange(const Index p0, const Numeric x, const Vec &xi,
                     Index polyorder = 1)
    requires(runtime_polyorder())
      : pos(is_ascending(xi)
                ? pos_finder<true, Limit>(p0, x, xi, polyorder)
                : pos_finder<false, Limit>(p0, x, xi, polyorder)),
        lx(polyorder + 1), dlx(polyorder + 1) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Standard initializer from Vector-types for runtime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] polyorder Polynominal degree
   * @param[in] flag For the order of xi
   */
  template <matpack::ranked_matpack_type<Numeric, 1> Vec>
  constexpr Lagrange(const Index p0, const Numeric x, const Vec &xi,
                     Index polyorder, AscendingOrder)
    requires(runtime_polyorder())
      : pos(pos_finder<true, Limit>(p0, x, xi, polyorder)),
        lx(polyorder + 1), dlx(polyorder + 1) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Standard initializer from Vector-types for runtime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] polyorder Polynominal degree
   * @param[in] flag For the order of xi
   */
  template <matpack::ranked_matpack_type<Numeric, 1> Vec>
  constexpr Lagrange(const Index p0, const Numeric x, const Vec &xi,
                     Index polyorder, DescendingOrder)
    requires(runtime_polyorder())
      : pos(pos_finder<false, Limit>(p0, x, xi, polyorder)),
        lx(polyorder + 1), dlx(polyorder + 1) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Standard initializer from Vector-types for compiletime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   */
  template <matpack::ranked_matpack_type<Numeric, 1> Vec>
  constexpr Lagrange(const Index p0, const Numeric x, const Vec &xi)
    requires(not runtime_polyorder())
      : pos(is_ascending(xi)
                ? pos_finder<true, PolyOrder, Limit>(p0, x, xi)
                : pos_finder<false, PolyOrder, Limit>(p0, x, xi)) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Standard initializer from Vector-types for compiletime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] flag For the order of xi
   */
  template <matpack::ranked_matpack_type<Numeric, 1> Vec>
  constexpr Lagrange(const Index p0, const Numeric x, const Vec &xi, AscendingOrder)
    requires(not runtime_polyorder())
      : pos(pos_finder<true, PolyOrder, Limit>(p0, x, xi)) {
    lx_finder(x, xi);
    dlx_finder(x, xi);
  }

  /*! Standard initializer from Vector-types for compiletime polyorder
   *
   * @param[in] pos0 Estimation of original position, must be [0, xi.size())
   * @param[in] x New grid position
   * @param[in] xi Old grid positions
   * @param[in] flag For the order of xi
   */
  template <matpack::ranked_matpack_type<Numeric, 1> Vec>
  constexpr Lagrange(const Index p0, const Numeric x, const Vec &xi, DescendingOrder)
    requires(not runtime_polyorder())
      : pos(pos_finder<false, PolyOrder, Limit>(p0, x, xi)) {
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
    * Note that the sum(lx) == 1, and this is guaranteed by
    * a reduction of the first N-1 elements of lx in this
    * method
    *
    * @param[in] x New grid position
    * @param[in] xi Old grid positions
    */
   template <matpack::ranked_matpack_type<Numeric, 1> Vec>
   constexpr void lx_finder(const Numeric x, const Vec &xi) {
     if constexpr (runtime_polyorder()) {
       for (Index j = 0; j < size() - 1; j++)
         lx[j] = l<type, Limit>(pos, size() - 1, x, xi, j);
     } else {
       if constexpr (PolyOrder == 0) {
       } else if constexpr (PolyOrder == 1) {
         lx[0] = l<PolyOrder, type, Limit>(pos, x, xi, 0);
       } else {
         for (Index j = 0; j < PolyOrder; j++)
           lx[j] = l<PolyOrder, type, Limit>(pos, x, xi, j);
       }
     }
     lx.back() = 1.0 - std::reduce(lx.begin(), lx.end() - 1);
   }

   /*! Finds dlx
    *
    * FIXME: Remove this comment if sum(dlx) == 0 or not.
    * I think it is, and if it is the loop should be fixed
    * to set dlx.back = - std::reduce(dlx.begin(), dlx.end() - 1);
    *
    * @param[in] x New grid position
    * @param[in] xi Old grid positions
    */
   template <matpack::ranked_matpack_type<Numeric, 1> Vec>
   constexpr void dlx_finder(const Numeric x, const Vec &xi) {
     if constexpr (do_derivs) {
       for (Index j = 0; j < size(); j++)
         dlx[j] = dl<type, Limit>(pos, size(), x, xi, lx, j);
     }
   }

 public:
   /*! Get the polynominal order of a compile time type */
   static constexpr Index polyorder()
     requires(not runtime_polyorder())
   {
    return PolyOrder;
   }

   /*! Checks the interpolation grid and throws if it is bad
    *
    * @param[in] xi Old grid positions
    * @param[in] polyorder Polynominal degree
    * @param[in] x {Min new x, Max new x} if non-cyclic
    * @param[in] extrapol Level of extrapolation allowed if non-cyclic
    */
   template <matpack::ranked_matpack_type<Numeric, 1> Vec>
   static void check(const Vec &xi, const Index polyorder,
                     [[maybe_unused]] const std::pair<Numeric, Numeric> x =
                         {std::numeric_limits<Numeric>::infinity(),
                          -std::numeric_limits<Numeric>::infinity()},
                     [[maybe_unused]] const Numeric extrapol = 0.5)
     requires(test_cyclic_limit<Limit>())
   {
    constexpr Numeric lb = Limit<cycle_limit::lower>::bound;
    constexpr Numeric ub = Limit<cycle_limit::upper>::bound;

    const Index n = Index(xi.size());

    ARTS_USER_ERROR_IF(polyorder >= n,
                       "Interpolation setup has failed!\n"
                       "\tRequesting greater interpolation order "
                       "than possible with given input grid")
    ARTS_USER_ERROR_IF(type == GridType::Cyclic and lb >= ub,
                       "Interpolation setup has failed!\n"
                       "\tBad cycle, must be [first, second)")
    if constexpr (GridType::Cyclic not_eq type) {
      if (polyorder and extrapol > 0) {
        const bool ascending = is_ascending(xi);
        const Numeric xmin =
            ascending
                ? xi[0] - extrapol * nonstd::abs(xi[1] - xi[0])
                : xi[n - 1] - extrapol * nonstd::abs(xi[n - 2] - xi[n - 1]);
        const Numeric xmax =
            ascending
                ? xi[n - 1] + extrapol * nonstd::abs(xi[n - 2] - xi[n - 1])
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
    * @param[in] x New grid position if non-cyclic
    * @param[in] extrapol Level of extrapolation allowed if non-cyclic
    */
   template <matpack::ranked_matpack_type<Numeric, 1> Vec>
   static constexpr void check(const Vec &xi, const Index polyorder,
                               const Numeric x = 0,
                               const Numeric extrapol = 0.5)
     requires(test_cyclic_limit<Limit>())
   {
    check(xi, polyorder, std::pair{x, x}, extrapol);
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
concept lagrange_type = requires(std::remove_cvref_t<T> a) {
  /* Normal lagrange weight */
  { a.pos } -> matpack::integral;
  { a.size() } -> matpack::integral;
  { a.index_pos(0, 0) } -> matpack::integral;
} and std::same_as<decltype(internal::has_derivatives<T>()), bool> and
      std::same_as<decltype(internal::runtime_polyorder<T>()), bool> and
requires (std::remove_cvref_t<T> a) {
   /* The lagrange weights can be accessed */
  { a.lx.size() } -> matpack::integral;  // Size is a size-type
  { a.lx[0] } -> matpack::arithmetic;  // May contain some values
} and
requires (std::remove_cvref_t<T> a) {
  { a.dlx.size() } -> matpack::integral;  // Size is a size-type
  { a.dlx[0] } -> matpack::arithmetic;  // May contain some values
};

/** Get a list of Lagrange types for use with reinterp
 *
 * This is the runtime polyorder version
 *
 * @param[in] xs The new grid, does not have to be sorted
 * @param[in] xi The original grid, should be sorted
 * @param[in] order The runtime polynominal order of the interpolation
 * @param[in] extrapol The allowed extrapolation (inf means all extrapolations
 * are allowed)
 * @return An array of Lagrange types
 */
template <lagrange_type T, matpack::ranked_matpack_type<Numeric, 1> NewVec, matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr Array<T> lagrange_interpolation_list(
    NewVec &&xs, Vec &&xi, Index order = 1,
    Numeric extrapol = 0.5)
  requires(std::remove_cvref_t<T>::runtime_polyorder())
{
  if (extrapol < std::numeric_limits<Numeric>::infinity()) {
    T::check(xi, order, minmax(xs), extrapol);
  }

  Array<T> out;
  out.reserve(xs.size());

  if (is_ascending(xi)) {
    for (auto x : xs) {
      out.emplace_back(out.size() ? out.back().pos : start_pos_finder(x, xi), x, xi, order, AscendingOrder{});
    }
  } else {
    for (auto x : xs) {
      out.emplace_back(out.size() ? out.back().pos : start_pos_finder(x, xi), x, xi, order, DescendingOrder{});
    }
  }

  return out;
}

/** Get a list of Lagrange types for use with reinterp
 *
 * This is the runtime polyorder version
 *
 * @param[in] xs The new grid, does not have to be sorted
 * @param[in] xi The original grid, should be sorted
 * @param[in] extrapol The allowed extrapolation (inf means all extrapolations
 * are allowed)
 * @return An array of Lagrange types
 */
template <lagrange_type T, matpack::ranked_matpack_type<Numeric, 1> NewVec, matpack::ranked_matpack_type<Numeric, 1> Vec>
constexpr Array<T> lagrange_interpolation_list(
    NewVec &&xs, Vec &&xi,
    Numeric extrapol = 0.5)
  requires(not std::remove_cvref_t<T>::runtime_polyorder())
{
  if (extrapol < std::numeric_limits<Numeric>::infinity()) {
    T::check(xi, T::polyorder(), minmax(xs), extrapol);
  }

  Array<T> out;
  out.reserve(xs.size());

  if (is_ascending(xi)) {
    for (auto x : xs) {
      out.emplace_back(out.size() ? out.back().pos : start_pos_finder(x, xi), x, xi, AscendingOrder{});
    }
  } else {
    for (auto x : xs) {
      out.emplace_back(out.size() ? out.back().pos : start_pos_finder(x, xi), x, xi, DescendingOrder{});
    }
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
   * @param[in] lag The Lagrange value
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
   * @param[in] all A list of Lagrange values
   * @return constexpr auto An iterable list
   */
  template <Index... inds>
  static constexpr auto as_elemwise(std::integer_sequence<Index, inds...>,
                                    const T &...all) requires (sizeof...(inds) == sizeof...(T)) {
    return matpack::elemwise{one_by_one<inds == dlx>(all)...};
  }
};
} // namespace internal

//! Conceptually, the type is assignable
template <typename T, typename U>
concept matpack_assignable = matpack::any_matpack_type<T> and requires(T& a) {
  std::apply(a, std::array<Index, matpack::rank<T>()>{}) = U{};
};

/** Precompute the interpolation weights
 *
 * For use with the interp method when re-usability is required
 *
 * The size of out is already set before calling this method
 *
 * @param[inout] out A writable matpack type 
 * @param[in] lag... Several (at least 1) Lagrange value
 */
template <lagrange_type... lags, Index N = sizeof...(lags)>
void interpweights(matpack::ranked_matpack_type<Numeric, N> auto&& out, const lags &...lag) 
requires(N > 0 and matpack_assignable<decltype(out), Numeric>) {
  ARTS_ASSERT(std::array{lag.size()...} == out.shape(),
              matpack::shape_help{std::array{lag.size()...}}, " vs ",
              matpack::shape_help{out.shape()})

  const auto in = matpack::elemwise{lag.lx...};
  std::transform(in.begin(), in.end(), out.elem_begin(), [](auto&& v){
    return std::apply([](auto&&... x){return (x * ...);}, v);
  });
}

/** Precompute the interpolation weights
 *
 * For use with the interp method when re-usability is required
 *
 * The size of out is already set before calling this method
 *
 * @param[inout] out A writable matpack type 
 * @param[in] lag... Several (at least 1) Lagrange value
 */
template <Index dlx, lagrange_type... lags, Index N = sizeof...(lags)>
void dinterpweights(matpack::ranked_matpack_type<Numeric, N> auto&& out, const lags &...lag) 
requires(N > 0 and matpack_assignable<decltype(out), Numeric> and dlx >= 0 and dlx < N) {
  ARTS_ASSERT(std::array{lag.size()...} == out.shape(),
              matpack::shape_help{std::array{lag.size()...}}, " vs ",
              matpack::shape_help{out.shape()})

  const auto in = internal::select_derivative<dlx, lags...>::as_elemwise(std::make_integer_sequence<Index, N>{}, lag...);
  std::transform(in.begin(), in.end(), out.elem_begin(), [](auto&& v){
    return std::apply([](auto&&... x){return (x * ...);}, v);
  });
}

/** Precompute the interpolation derivative weights
 *
 * For use with the interp method when re-usability is required
 * 
 * @param[inout] out A writable matpack type 
 * @param[in] lag... Several (at least 1) Lagrange value
 */
template<lagrange_type... lags, Index N = sizeof...(lags)>
constexpr auto interpweights(const lags&... lag) requires (N > 0) {
  if constexpr (N > 1) {
    if constexpr ((internal::runtime_polyorder<lags>() or ...)) {
      matpack::matpack_data<Numeric, N> out(lag.size()...);
      interpweights(out, lag...);
      return out;
    } else {
      matpack::matpack_constant_data<Numeric, internal::compile_time_size<lags>()...> out{};
      interpweights(out, lag...);
      return out;
    }
  } else {
    return std::get<0>(std::tuple{lag...}).lx;
  }
}

/** Precompute the interpolation derivative weights
 *
 * For use with the interp method when re-usability is required
 * 
 * @tparam dlx And Index indicating which dimension's derivative is wanted
 * @tparam lags... Several Lagrange types
 * @param[in] lag... Several Lagrange values, where the dlx:th one has a derivative
 */
template<Index dlx, lagrange_type... lags, Index N = sizeof...(lags)>
constexpr auto dinterpweights(const lags&... lag) requires (N > 0 and dlx >= 0 and dlx < N) {
  if constexpr (N > 1) {
    if constexpr ((internal::runtime_polyorder<lags>() or ...)) {
      matpack::matpack_data<Numeric, N> out(lag.size()...);
      dinterpweights<dlx>(out, lag...);
      return out;
    } else {
      matpack::matpack_constant_data<Numeric, internal::compile_time_size<lags>()...> out{};
      dinterpweights<dlx>(out, lag...);
      return out;
    }
  } else {
    return std::get<0>(std::tuple{lag...}).dlx;
  }
}

//! Test to make sure that the type is some list of a type that can be used as a lagrange key type
template <typename T>
concept list_of_lagrange_type = requires(T a) {
  { a.size() } -> std::convertible_to<Index>;
  { a[0] } -> lagrange_type;
};

/** Precompute the interpolation weights
 *
 * For use with the reinterp method when re-usability is required
 * 
 * @param[in] lags... Several lists of Lagrange values
 */
template <list_of_lagrange_type... list_lags, Index N = sizeof...(list_lags)>
constexpr auto interpweights(const list_lags &...lags)
  requires(N > 0)
{
  matpack::matpack_data<matpack::matpack_data<Numeric, N>, N> out(
      static_cast<Index>(lags.size())...,
      matpack::matpack_data<Numeric, N>(
          std::array{
              static_cast<Index>(lags.size() == 0 ? 0 : lags[0].size())...},
          0.0));

  ARTS_ASSERT((std::all_of(lags.begin(), lags.end(),
                           [SZ = lags.size() == 0 ? 0 : lags[0].size()](
                               auto &l) { return SZ == l.size(); }) and
               ...),
              "All input lags must have the same interpolation order")

  if (out.empty()) return out;

  for (matpack::flat_shape_pos<N> pos{
           std::array{static_cast<Index>(lags.size())...}};
       pos.pos.front() < pos.shp.front(); ++pos) {
    std::apply([&](auto&&... ind) { interpweights(out(ind...), lags[ind]...); },
               pos.pos);
  }

  return out;
}

/** Precompute the interpolation derivative weights
 *
 * For use with the reinterp method when re-usability is required
 * 
 * @param[in] lags... Several lists of Lagrange values, where the dlx:th one has a derivatives
 */
template <Index dlx, list_of_lagrange_type... list_lags,
          Index N = sizeof...(list_lags)>
constexpr auto dinterpweights(const list_lags &...lags)
  requires(N > 0 and dlx >= 0 and dlx < N)
{
  matpack::matpack_data<matpack::matpack_data<Numeric, N>, N> out(
      static_cast<Index>(lags.size())...,
      matpack::matpack_data<Numeric, N>(
          std::array{
              static_cast<Index>(lags.size() == 0 ? 0 : lags[0].size())...},
          0.0));

  ARTS_ASSERT((std::all_of(lags.begin(), lags.end(),
                           [SZ = lags.size() == 0 ? 0 : lags[0].size()](
                               auto &l) { return SZ == l.size(); }) and
               ...),
              "All input lags must have the same interpolation order")

  if (out.empty()) return out;

  for (matpack::flat_shape_pos<N> pos{
           std::array{static_cast<Index>(lags.size())...}};
       pos.pos.front() < pos.shp.front(); ++pos) {
    std::apply(
        [&](auto &&...ind) { dinterpweights<dlx>(out(ind...), lags[ind]...); },
        pos.pos);
    }

    return out;
  }

/** Interpolate a single output value with re-usable weights
 * 
 * @tparam lags... Several Lagrange types
 * @param[in] field A field value of the same rank as the number of lags
 * @param[in] iw An interpolation weight of the same rank as the number of lags
 * @param[in] lag... Several Lagrange values
 * @return constexpr auto The value of the interpolation
 */
template <lagrange_type... lags, std::size_t N = sizeof...(lags)>
constexpr auto interp(const matpack::strict_rank_matpack_type<N> auto &field,
                      const auto &iw,  // FIXME: MUST BE MATPACK TYPE IN THE END
                      const lags &...lag) {
  matpack::matpack_value_type<decltype(field)> out{0};

  //! Now the cyclicity is wrapped, we can just add the values up
  for (matpack::flat_shape_pos<sizeof...(lags)> pos{matpack::mdshape(iw)};
       pos.pos.front() < pos.shp.front(); ++pos) {
    out += std::apply(std::apply(
                          [&](auto &&...maxsize) {
                            return [&](auto &&...offset) {
                              return matpack::mdvalue(iw, {offset...}) *
                                     matpack::mdvalue(
                                         field,
                                         {lag.index_pos(offset, maxsize)...});
                            };
                          },
                          matpack::mdshape(field)),
                      pos.pos);
  }

  return out;
}

/** Reinterpolates a field as another field with re-usable weights
 *
 *  This is done by calling interp for all the combinations of the inputs
 * 
 * @tparam lags... Several lists of Lagrange types
 * @param[inout] out A writable output-field of the size of the list_lag-variables
 * @param[in] field A field value of the same rank as the number of lags
 * @param[in] iw_field A field of interpolation weights of the same rank as the number of lags
 * @param[in] list_lag... Several lists of Lagrange values
 * @return constexpr auto A new field
 */
template <list_of_lagrange_type... lags, std::size_t N = sizeof...(lags)>
constexpr void reinterp(auto &&out,
                        const matpack::strict_rank_matpack_type<N> auto &field,
                        const matpack::strict_rank_matpack_type<N> auto &iw_field, const lags &...list_lag) {
  const auto in = matpack::elemwise{list_lag...};
  std::transform(
      iw_field.elem_begin(), iw_field.elem_end(), in.begin(), out.elem_begin(),
      [&](auto &&internal_iw, auto &&lag_t) {
        return std::apply(
            [&](auto &&...lag) { return interp(field, internal_iw, lag...); },
            lag_t);
      });
}

/** Reinterpolates a field as another field with re-usable weights
 *
 *  This is done by calling interp for all the combinations of the inputs
 * 
 * @tparam lags... Several lists of Lagrange types
 * @param[in] field A field value of the same rank as the number of lags
 * @param[in] iw_field A field of interpolation weights of the same rank as the number of lags
 * @param[in] lag... Several lists of Lagrange values
 * @return constexpr auto A new field
 */
template <list_of_lagrange_type... lags, std::size_t N = sizeof...(lags)>
constexpr auto reinterp(const matpack::strict_rank_matpack_type<N> auto &field,
                        const matpack::strict_rank_matpack_type<N> auto &iw_field, const lags &...list_lag) {
  matpack::matpack_data<matpack::matpack_value_type<decltype(field)>,
                        sizeof...(lags)>
  out(list_lag.size()...);
  reinterp(out, field, iw_field, list_lag...);
  return out;
}

/** Interpolate a single output value
 * 
 * @tparam lags... Several Lagrange types
 * @param[in] field A field value of the same rank as the number of lags
 * @param[in] lag... Several Lagrange values
 * @return constexpr auto The value of the interpolation
 */
template <lagrange_type... lags, std::size_t N = sizeof...(lags)>
constexpr auto interp(const matpack::strict_rank_matpack_type<N> auto &field,
                      const lags &...lag) {
  matpack::matpack_value_type<decltype(field)> out{0};

  for (matpack::flat_shape_pos<sizeof...(lags)> pos{std::array{lag.size()...}};
       pos.pos.front() < pos.shp.front(); ++pos) {
    out += std::apply(std::apply(
                          [&](auto &&...maxsize) {
                            return [&](auto &&...offset) {
                              return (lag.lx[offset] * ...) *
                                     matpack::mdvalue(
                                         field,
                                         {lag.index_pos(offset, maxsize)...});
                            };
                          },
                          matpack::mdshape(field)),
                      pos.pos);
  }

  return out;
}

/** Reinterpolates a field as another field
 *
 *  This is done by calling interp for all the combinations of the inputs
 * 
 * @tparam lags... Several lists of Lagrange types
 * @param[in] field A field value of the same rank as the number of lags
 * @param[in] lag... Several lists of Lagrange values
 * @return constexpr auto A new field
 */
template <list_of_lagrange_type... lags, std::size_t N = sizeof...(lags)>
constexpr auto reinterp(const matpack::strict_rank_matpack_type<N> auto &field,
                        const lags &...list_lag) {
  const auto in = matpack::elemwise{list_lag...};
  matpack::matpack_data<matpack::matpack_value_type<decltype(field)>,
                        sizeof...(lags)>
  out(list_lag.size()...);

  std::transform(
      in.begin(), in.end(), out.elem_begin(),
      [&](auto &&lag_t) {
        return std::apply(
            [&](auto &&...lag) { return interp(field, lag...); },
            lag_t);
      });

  return out;
}
}  // namespace my_interp

// For linear interpolation
using LagrangeInterpolation = my_interp::Lagrange<>;
using ArrayOfLagrangeInterpolation = Array<LagrangeInterpolation>;

// For logarithmic interpolation
using LagrangeLogInterpolation = my_interp::Lagrange<-1, false, my_interp::GridType::Log>;
using ArrayOfLagrangeLogInterpolation = Array<LagrangeLogInterpolation>;

// For cyclic interpolation between 0 and 360
using LagrangeCyclic0to360Interpolation = my_interp::Lagrange<-1, false, my_interp::GridType::Cyclic, my_interp::cycle_0_p360>;
using ArrayOfLagrangeCyclic0to360Interpolation = Array<LagrangeCyclic0to360Interpolation>;

template <std::size_t sz, bool deriv=false>
using FixedLagrangeInterpolation = my_interp::Lagrange<sz, deriv>;
