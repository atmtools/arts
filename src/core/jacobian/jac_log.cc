#include "jac_log.h"

#include <algorithm>
#include <cmath>

#include "jacobian.h"

Vector logfwd::operator()(ConstVectorView xx) const {
  if (N != xx.size())
    arts::user_error(R"(Mismatch size.

Cannot produce a logarithmic profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size());

  constexpr auto log = [](Numeric v) { return std::log(v); };

  Vector x(N);

  std::transform(xx.begin(), xx.end(), x.begin(), log);

  return x;
}

Vector loginv::operator()(ConstVectorView xx) const {
  if (N != xx.size())
    arts::user_error(R"(Mismatch size.

Cannot produce a logarithmic profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size());

  constexpr auto exp = [](Numeric v) { return std::exp(v); };

  Vector x(N);

  std::transform(xx.begin(), xx.end(), x.begin(), exp);

  return x;
}

Matrix loginv::operator()(ConstMatrixView dyy, ConstVectorView x) const {
  if (N != static_cast<Size>(dyy.ncols()))
    arts::user_error(R"(Mismatch size.

Cannot produce a logarithmic profile.

The original data     has {} elements
The incoming jacobian has {} elements
)",
                     N,
                     dyy.ncols());
  Matrix dy(dyy);

  for (Size i = 0; i < N; i++) {
    const Numeric ex  = std::exp(x[i]);
    dy[joker, i]     *= ex;
  }

  return dy;
}
