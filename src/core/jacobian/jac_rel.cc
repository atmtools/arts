#include "jac_rel.h"

#include "matpack_einsum.h"

Vector relfwd::operator()(ConstVectorView xx) const {
  ARTS_USER_ERROR_IF(not sorig, "No original-vector provided for polyinv.")
  const Vector& orig = *sorig;

  ARTS_USER_ERROR_IF(orig.size() != xx.size(),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     orig.size(),
                     xx.size())

  Vector x{xx};
  x /= orig;
  return x;
}

Vector relfwd::operator()(ConstVectorView x, const AtmField&) const {
  return this->operator()(x);
}

Vector relfwd::operator()(ConstVectorView x, const SurfaceField&) const {
  return this->operator()(x);
}

Vector relfwd::operator()(ConstVectorView x, const SubsurfaceField&) const {
  return this->operator()(x);
}

Vector relinv::operator()(ConstVectorView xx) const {
  ARTS_USER_ERROR_IF(not sorig, "No original-vector provided for polyinv.")
  const Vector& orig = *sorig;

  ARTS_USER_ERROR_IF(orig.size() != xx.size(),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     orig.size(),
                     xx.size())
  return matpack::einsum<Vector, "i", "i", "i">(xx.shape(), xx, orig);
}

Vector relinv::operator()(ConstVectorView x, const AtmField&) const {
  return this->operator()(x);
}

Vector relinv::operator()(ConstVectorView x, const SurfaceField&) const {
  return this->operator()(x);
}

Vector relinv::operator()(ConstVectorView x, const SubsurfaceField&) const {
  return this->operator()(x);
}

Matrix relinv::operator()(ConstMatrixView dyy) const {
  ARTS_USER_ERROR_IF(not sorig, "No original-vector provided for polyinv.")
  const Vector& orig = *sorig;
  const Size N       = orig.size();

  ARTS_USER_ERROR_IF(N != static_cast<Size>(dyy.ncols()),
                     R"(Mismatch size.

Cannot produce a relative profile.

The original data     has {} elements
The incoming jacobian has {} elements
)",
                     N,
                     dyy.ncols())

  // Matrix dy(dyy);

  // for (Size i = 0; i < N; i++) {
  //   dy[joker, i] *= y[i] / x[i];  // ??? y-> e.g., real vmr, x-> the ratio
  // }

  // return dy;

  // This is what I suspect is the truth, but ARTS2 does it as above
  return matpack::einsum<Matrix, "ij", "ij", "j">(dyy.shape(), dyy, orig);
}

Matrix relinv::operator()(ConstMatrixView dy,
                          ConstVectorView,
                          const AtmField&) const {
  return this->operator()(dy);
}

Matrix relinv::operator()(ConstMatrixView dy,
                          ConstVectorView,
                          const SurfaceField&) const {
  return this->operator()(dy);
}

Matrix relinv::operator()(ConstMatrixView dy,
                          ConstVectorView,
                          const SubsurfaceField&) const {
  return this->operator()(dy);
}

void make_relfit(Jacobian::AtmTarget& x, const AtmField& atm) {
  auto sorig = std::make_shared<Vector>(atm[x.type].flat_view());
  const relfwd rfwd{.sorig = sorig};
  const relinv rinv{.sorig = sorig};  // , .atm = x.type};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}

void make_relfit(Jacobian::SurfaceTarget& x, const SurfaceField& surf) {
  auto sorig = std::make_shared<Vector>(surf[x.type].flat_view());
  const relfwd rfwd{.sorig = sorig};
  const relinv rinv{.sorig = sorig};  // , .surf = x.type};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}

void make_relfit(Jacobian::SubsurfaceTarget& x,
                 const SubsurfaceField& subsurf) {
  auto sorig = std::make_shared<const Vector>(subsurf[x.type].flat_view());
  const relfwd rfwd{.sorig = sorig};
  const relinv rinv{.sorig = sorig};  // , .subsurf = x.type};
  x.inverse_state    = rinv;
  x.inverse_jacobian = rinv;
  x.transform_state  = rfwd;
}
