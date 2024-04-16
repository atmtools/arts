#include <iomanip>
#include <iostream>
#include <numeric>
#include <ranges>

#include "debug.h"
#include "legendre.h"

void sumlegtest() {
  const Vector c = {1.,          3.,          5.,          7.,
                    9.,          11.,         13.,         15.,
                    17.,         19.,         21.,         23.,
                    25.,         27.,         29.,         31.,
                    33.,         32.8125,     29.91796875, 25.96508789,
                    21.84068298, 17.98672199, 14.59263057, 11.71004179,
                    9.31995848,  7.37111478,  5.80115833,  4.54775559,
                    3.55391083,  2.77005711,  2.15442342,  1.67254542};

  const Vector nu{0.01985507,
                  0.10166676,
                  0.2372338,
                  0.40828268,
                  0.59171732,
                  0.7627662,
                  0.89833324,
                  0.98014493};

  // From numpy.polynomial.legendre.Legendre
  const Vector res{0.048454,
                   0.38029738,
                   -0.52891716,
                   0.44314246,
                   -1.33650515,
                   1.33803551,
                   8.08787746,
                   -29.07110553};

  Vector sum(8);
  for (Index i = 0; i < 8; i++) {
    sum[i] = Legendre::legendre_sum(c, nu[i]);
  }

  const Numeric sum_squared_diff =
      std::transform_reduce(res.begin(),
                            res.end(),
                            sum.begin(),
                            0.,
                            std::plus<>(),
                            [](auto a, auto b) { return (a - b) * (a - b); });

  if (sum_squared_diff > 1e-10) {
    throw std::runtime_error("Legendre sum test failed");
  }
}

int main() { sumlegtest(); }
