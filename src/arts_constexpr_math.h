/*!
 * \file   arts_constexpr_math.h
 * \brief  Simple constexpr math that we can make use of in Arts
 *
 * The main advantage of maning things like here is to avoid repetition
 * later x*x*x*x is more difficult than pow4(x).
 * 
 * The main point of these functions is to avoid repetition as above in
 * context where we also want to support compile time computations
 * 
 * \author Richard Larsson
 * \date   2019-04-01
 */

#ifndef CONSTANT_MATHS_IN_ARTS_H
#define CONSTANT_MATHS_IN_ARTS_H

/** Namespace containing several constants, physical and mathematical **/
namespace Math {
/** power of two */
constexpr auto pow2(auto x) noexcept { return x * x; }

/** power of three */
constexpr auto pow3(auto x) noexcept { return pow2(x) * x; }

/** power of four */
constexpr auto pow4(auto x) noexcept { return pow2(pow2(x)); }
}  // namespace Math

#endif 
