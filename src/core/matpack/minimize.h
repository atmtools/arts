#ifndef minimize_h
#define minimize_h

#include <matpack.h>

#include <optional>

namespace Minimize {
/*! Fit a curve to a polynomial.
 * 
 * Wraps curve_fit(Polynom(X, Y, order)). See it for more details.
 *
 * 
 * @param[in] X The X-Grid
 * @param[in] Y The computed/measured/simulated values at X
 * @param[in] order The order of the polynomial to fit, e.g., 3 means x**3 is largest factor
 * @return A Vector if successful, otherwise an empty optional
 */
std::optional<Vector> polyfit(const StridedConstVectorView& X,
                              const StridedConstVectorView& Y,
                              const Index& order);
}  // namespace Minimize

#endif  // minimize_wrap_h
