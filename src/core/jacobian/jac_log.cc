#include "jac_log.h"

#include <algorithm>
#include <cmath>

#include "jacobian.h"

Vector logfwd::operator()(ConstVectorView xx) const {
  ARTS_USER_ERROR_IF(N != xx.size(),
                     R"(Mismatch size.

Cannot produce a logarithmic profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size())

  constexpr auto log = [](Numeric v) { return std::log(v); };

  Vector x(N);

  std::transform(xx.begin(), xx.end(), x.begin(), log);

  return x;
}

Vector logfwd::operator()(ConstVectorView x, const AtmField&) const {
  return this->operator()(x);
}

Vector logfwd::operator()(ConstVectorView x, const SurfaceField&) const {
  return this->operator()(x);
}

Vector logfwd::operator()(ConstVectorView x, const SubsurfaceField&) const {
  return this->operator()(x);
}

Vector loginv::operator()(ConstVectorView xx) const {
  ARTS_USER_ERROR_IF(N != xx.size(),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size())

  constexpr auto exp = [](Numeric v) { return std::exp(v); };

  Vector x(N);

  std::transform(xx.begin(), xx.end(), x.begin(), exp);

  return x;
}

Vector loginv::operator()(ConstVectorView x, const AtmField&) const {
  return this->operator()(x);
}

Vector loginv::operator()(ConstVectorView x, const SurfaceField&) const {
  return this->operator()(x);
}

Vector loginv::operator()(ConstVectorView x, const SubsurfaceField&) const {
  return this->operator()(x);
}

Matrix loginv::operator()(ConstMatrixView dyy, ConstVectorView x) const {
  ARTS_USER_ERROR_IF(N != static_cast<Size>(dyy.ncols()),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data     has {} elements
The incoming jacobian has {} elements
)",
                     N,
                     dyy.ncols())
  Matrix dy(dyy);

  for (Size i = 0; i < N; i++) {
    const Numeric ex  = std::exp(x[i]);
    dy[joker, i]     *= ex;
  }

  return dy;
}

Matrix loginv::operator()(ConstMatrixView dy,
                          ConstVectorView x,
                          const AtmField&) const {
  return this->operator()(dy, x);
}

Matrix loginv::operator()(ConstMatrixView dy,
                          ConstVectorView x,
                          const SurfaceField&) const {
  return this->operator()(dy, x);
}

Matrix loginv::operator()(ConstMatrixView dy,
                          ConstVectorView x,
                          const SubsurfaceField&) const {
  return this->operator()(dy, x);
}

void make_logfit(Jacobian::AtmTarget& x, const AtmField& atm) {
  const logfwd rfwd{.N = atm[x.type].flat_view().size()};
  const loginv rinv{.N = atm[x.type].flat_view().size()};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}

void make_logfit(Jacobian::SurfaceTarget& x, const SurfaceField& surf) {
  const logfwd rfwd{.N = surf[x.type].flat_view().size()};
  const loginv rinv{.N = surf[x.type].flat_view().size()};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}

void make_logfit(Jacobian::SubsurfaceTarget& x,
                 const SubsurfaceField& subsurf) {
  const logfwd rfwd{.N = subsurf[x.type].flat_view().size()};
  const loginv rinv{.N = subsurf[x.type].flat_view().size()};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}
