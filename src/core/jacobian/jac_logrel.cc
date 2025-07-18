#include "jac_logrel.h"

#include "jacobian.h"

Vector logrelfwd::operator()(ConstVectorView xx) const {
  ARTS_USER_ERROR_IF(not sorig, "No original-vector provided for polyinv.")
  const Vector& orig = *sorig;

  const Size N = orig.size();
  ARTS_USER_ERROR_IF(N != xx.size(),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size())

  Vector x(N);
  for (Size i = 0; i < N; i++) x[i] = std::log(xx[i] / orig[i]);
  return x;
}

Vector logrelfwd::operator()(ConstVectorView x, const AtmField&) const {
  return this->operator()(x);
}

Vector logrelfwd::operator()(ConstVectorView x, const SurfaceField&) const {
  return this->operator()(x);
}

Vector logrelfwd::operator()(ConstVectorView x, const SubsurfaceField&) const {
  return this->operator()(x);
}

Vector logrelinv::operator()(ConstVectorView xx) const {
  ARTS_USER_ERROR_IF(not sorig, "No original-vector provided for polyinv.")
  const Vector& orig = *sorig;

  const Size N = orig.size();
  ARTS_USER_ERROR_IF(N != xx.size(),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     N,
                     xx.size())

  Vector x(N);
  for (Size i = 0; i < N; i++) x[i] = std::exp(xx[i]) * orig[i];
  return x;
}

Vector logrelinv::operator()(ConstVectorView x, const AtmField&) const {
  return this->operator()(x);
}

Vector logrelinv::operator()(ConstVectorView x, const SurfaceField&) const {
  return this->operator()(x);
}

Vector logrelinv::operator()(ConstVectorView x, const SubsurfaceField&) const {
  return this->operator()(x);
}

Matrix logrelinv::operator()(ConstMatrixView dyy, ConstVectorView x) const {
  ARTS_USER_ERROR_IF(not sorig, "No original-vector provided for polyinv.")
  const Vector& orig = *sorig;

  const Size N = orig.size();

  ARTS_USER_ERROR_IF(N != static_cast<Size>(dyy.ncols()) or x.size() != N,
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data     has {} elements
The incoming jacobian has {} elements
The incoming x        has {} elements
)",
                     N,
                     dyy.ncols(),
                     x.size())

  Matrix dy(dyy);
  for (Size j = 0; j < N; j++) dy[joker, j] *= std::exp(x[j]) * orig[j];
  return dy;
}

Matrix logrelinv::operator()(ConstMatrixView dy,
                             ConstVectorView x,
                             const AtmField&) const {
  return this->operator()(dy, x);
}

Matrix logrelinv::operator()(ConstMatrixView dy,
                             ConstVectorView x,
                             const SurfaceField&) const {
  return this->operator()(dy, x);
}

Matrix logrelinv::operator()(ConstMatrixView dy,
                             ConstVectorView x,
                             const SubsurfaceField&) const {
  return this->operator()(dy, x);
}

void make_logrelfit(Jacobian::AtmTarget& x, const AtmField& atm) {
  auto sorig = std::make_shared<Vector>(atm[x.type].flat_view());
  const logrelfwd rfwd{.sorig = sorig};
  const logrelinv rinv{.sorig = sorig};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}

void make_logrelfit(Jacobian::SurfaceTarget& x, const SurfaceField& surf) {
  auto sorig = std::make_shared<Vector>(surf[x.type].flat_view());
  const logrelfwd rfwd{.sorig = sorig};
  const logrelinv rinv{.sorig = sorig};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}

void make_logrelfit(Jacobian::SubsurfaceTarget& x,
                    const SubsurfaceField& subsurf) {
  auto sorig = std::make_shared<const Vector>(subsurf[x.type].flat_view());
  const logrelfwd rfwd{.sorig = sorig};
  const logrelinv rinv{.sorig = sorig};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}
