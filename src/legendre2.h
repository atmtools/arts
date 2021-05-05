#ifndef legendre2_h
#define legendre2_h

#include "grids.h"
#include "matpackI.h"

namespace Legendre {
struct SphericalField {
  Numeric U;
  Numeric S;
  Numeric E;
  Numeric total() const noexcept {return std::hypot(U, S, E);}
  constexpr SphericalField() noexcept : U(0), S(0), E(0) {}
};

using MatrixOfSphericalField = Grid<SphericalField, 2>;

SphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Numeric r, const Numeric lat, const Numeric lon) ARTS_NOEXCEPT;

MatrixOfSphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Vector& r, const Numeric lat, const Vector& lon) ARTS_NOEXCEPT;
}

#endif  // legendre2_h
