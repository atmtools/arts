#include "jac_rel.h"

#include "jacobian.h"

Vector relfwd::operator()(ConstVectorView xx) const {
  if (not sorig) arts::user_error("No original-vector provided for relfwd.");
  const Vector& orig = *sorig;

  if (orig.size() != xx.size())
    arts::user_error(R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     orig.size(),
                     xx.size());

  Vector x{xx};
  x /= orig;
  return x;
}

Vector relinv::operator()(ConstVectorView xx) const {
  if (not sorig) arts::user_error("No original-vector provided for relinv.");
  const Vector& orig = *sorig;

  if (orig.size() != xx.size())
    arts::user_error(R"(Mismatch size.

Cannot produce a relative profile.

The original data has {} elements
The incoming data has {} elements
)",
                     orig.size(),
                     xx.size());
  return matpack::einsum<Vector, "i", "i", "i">(xx, orig);
}

Matrix relinv::operator()(ConstMatrixView dyy) const {
  if (not sorig) arts::user_error("No original-vector provided for relinv.");
  const Vector& orig = *sorig;
  const Size    N    = orig.size();

  if (N != static_cast<Size>(dyy.ncols()))
    arts::user_error(R"(Mismatch size.

Cannot produce a relative profile.

The original data     has {} elements
The incoming jacobian has {} elements
)",
                     N,
                     dyy.ncols());

  // y = x * orig, so dy/dx = orig (ARTS2 used y / x, which is the same)
  return matpack::einsum<Matrix, "ij", "ij", "j">(dyy, orig);
}
