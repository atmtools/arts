#include "jac_logrel.h"

#include "jacobian.h"

Vector logrelfwd::operator()(ConstVectorView xx) const {
  if (not sorig) arts::user_error("No original-vector provided for logrelfwd.");
  const Vector& orig = *sorig;

  const Size N = orig.size();
  if (N != xx.size())
    arts::user_error(R"(Mismatch size.

Cannot produce a log-relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size());

  Vector x(N);
  for (Size i = 0; i < N; i++) x[i] = std::log(xx[i] / orig[i]);
  return x;
}

Vector logrelinv::operator()(ConstVectorView xx) const {
  if (not sorig) arts::user_error("No original-vector provided for logrelinv.");
  const Vector& orig = *sorig;

  const Size N = orig.size();
  if (N != xx.size())
    arts::user_error(R"(Mismatch size.

Cannot produce a log-relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size());

  Vector x(N);
  for (Size i = 0; i < N; i++) x[i] = std::exp(xx[i]) * orig[i];
  return x;
}

Matrix logrelinv::operator()(ConstMatrixView dyy, ConstVectorView x) const {
  if (not sorig) arts::user_error("No original-vector provided for logrelinv.");
  const Vector& orig = *sorig;

  const Size N = orig.size();

  if (N != static_cast<Size>(dyy.ncols()) or x.size() != N)
    arts::user_error(R"(Mismatch size.

Cannot produce a log-relative profile.

The original data     has {} elements
The incoming jacobian has {} elements
The incoming x        has {} elements
)",
                     N,
                     dyy.ncols(),
                     x.size());

  Matrix dy(dyy);
  for (Size j = 0; j < N; j++) dy[joker, j] *= std::exp(x[j]) * orig[j];
  return dy;
}
