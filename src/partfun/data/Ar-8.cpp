//! auto-generated file

#include <configtypes.h>

#include <array>

inline constexpr std::array<Numeric, 4> coeffs{ 1, 0, 0, 0 };


Numeric /Users/richard/Work/arts/src/partfun/data/Ar8(Numeric T) noexcept {
  Numeric result = coeffs[0];
  Numeric TN     = 1.0;

  for (int i = 1; i < 4; i++) {
    TN     *= T;
    result += TN * c;
  }

  return result;
}

Numeric d/Users/richard/Work/arts/src/partfun/data/Ar8(Numeric T) noexcept {
  Numeric result = coeffs[1];
  Numeric TN     = 1.0;

  for (int i = 2; i < 4; i++) {
    TN     *= T;
    result += static_cast<Numeric>(i) * TN * coeffs[i];
  }

  return result;
}
