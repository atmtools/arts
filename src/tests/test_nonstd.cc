/* nonstd::pow against std::pow

   exp(v log(x)) rounds v log(x), which the exponential amplifies, besides the
   rounding of exp and log themselves: a relative error of about
   epsilon (|v log(x)| + 2), bounded here by 4 epsilon (|v log(x)| + 2).  Its
   domain is x >= 0: zero gives 0 for v > 0 and +inf for v < 0, as std::pow
   does, and a negative x asserts in debug builds and is NaN otherwise.
   Integral arguments do not compile.
*/
#include <nonstd.h>

#include <cmath>
#include <iostream>
#include <limits>

namespace {
template <class T, class U>
concept powable = requires(T x, U v) { nonstd::pow(x, v); };
static_assert(powable<double, double> and powable<float, double>);
static_assert(not powable<int, double> and not powable<double, int> and not powable<long, long>);
}  // namespace

int main() {
  constexpr double eps   = std::numeric_limits<double>::epsilon();
  double           worst = 0.0;
  for (int i = -300; i <= 300; i += 7) {
    const double x = 1.37 * std::pow(10.0, i / 10.0);
    for (const double v : {-7.3, -2.0, -0.5, 0.0, 0.25, 1.0, 3.7, 11.0}) {
      const double expected = std::pow(x, v);
      if (not std::isfinite(expected) or expected == 0.0) continue;
      const double bound = 4 * eps * (std::abs(v * std::log(x)) + 2);
      worst              = std::max(worst, std::abs(nonstd::pow(x, v) / expected - 1) / bound);
    }
  }
  if (not(worst <= 1.0)) {
    std::cerr << "nonstd::pow is " << worst << " times its rounding bound from std::pow\n";
    return 1;
  }
  if (nonstd::pow(0.0, 2.5) != 0.0 or nonstd::pow(0.0, -2.5) != std::numeric_limits<double>::infinity()) {
    std::cerr << "nonstd::pow of zero must be 0 for v > 0 and +inf for v < 0\n";
    return 1;
  }
#ifdef NDEBUG
  if (not std::isnan(nonstd::pow(-2.0, 2.0)) or not std::isnan(nonstd::pow(-2.0, 0.5))) {
    std::cerr << "nonstd::pow of a negative x must be NaN\n";
    return 1;
  }
#endif
  std::cout << "nonstd::pow at " << worst << " of its rounding bound\n";
}
