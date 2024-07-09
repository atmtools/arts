#include "lbl_temperature_model.h"

#include <limits>

#include "debug.h"
#include "double_imanip.h"

namespace lbl::temperature {
namespace model {
using std::log;
using std::pow;

Numeric T1(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1);
}

Numeric dT1_dX0(Numeric, Numeric X1, Numeric T0, Numeric T) {
  return pow(T0 / T, X1);
}

Numeric dT1_dX1(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * log(T0 / T);
}

Numeric dT1_dT0(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * X1 * pow(T0 / T, X1) / T0;
}

Numeric dT1_dT(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return -X0 * X1 * pow(T0 / T, X1) / T;
}

Numeric T2(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * (1 + X2 * log(T / T0));
}

Numeric dT2_dX0(
    Numeric, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X1) * (1 + X2 * log(T / T0));
}

Numeric dT2_dX1(
    Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) * log(T0 / T);
}

Numeric dT2_dX2(
    Numeric X0, Numeric X1, Numeric, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * log(T / T0);
}

Numeric dT2_dT0(
    Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X0 * X1 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) / T0 -
         X0 * X2 * pow(T0 / T, X1) / T0;
}

Numeric dT2_dT(
    Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return -X0 * X1 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) / T +
         X0 * X2 * pow(T0 / T, X1) / T;
}

Numeric T4(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return (X0 + X1 * (T0 / T - 1)) * pow(T0 / T, X2);
}

Numeric dT4_dX0(Numeric, Numeric, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X2);
}

Numeric dT4_dX1(Numeric, Numeric, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X2) * (T0 / T - 1.);
}

Numeric dT4_dX2(
    Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1)) * log(T0 / T);
}

Numeric dT4_dT0(
    Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return X2 * pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1.)) / T0 +
         X1 * pow(T0 / T, X2) / T;
}

Numeric dT4_dT(
    Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T) {
  return -X2 * pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1.)) / T -
         T0 * X1 * pow(T0 / T, X2) / (T * T);
}

Numeric T5(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, 0.25 + 1.5 * X1);
}

Numeric dT5_dX0(Numeric, Numeric X1, Numeric T0, Numeric T) {
  return pow(T0 / T, 1.5 * X1 + 0.25);
}

Numeric dT5_dX1(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return 1.5 * X0 * pow(T0 / T, 1.5 * X1 + 0.25) * log(T0 / T);
}

Numeric dT5_dT0(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, 1.5 * X1 + 0.25) * (1.5 * X1 + 0.25) / T0;
}

Numeric dT5_dT(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return -X0 * pow(T0 / T, 1.5 * X1 + 0.25) * (1.5 * X1 + 0.25) / T;
}

Numeric DPL(Numeric X0,
            Numeric X1,
            Numeric X2,
            Numeric X3,
            Numeric T0,
            Numeric T) {
  return X0 * pow(T0 / T, X1) + X2 * pow(T0 / T, X3);
}

Numeric dDPL_dX0(
    Numeric, Numeric X1, Numeric, Numeric, Numeric T0, Numeric T) {
  return pow(T0 / T, X1);
}

Numeric dDPL_dX1(
    Numeric X0, Numeric X1, Numeric, Numeric, Numeric T0, Numeric T) {
  return X0 * pow(T0 / T, X1) * log(T0 / T);
}

Numeric dDPL_dX2(
    Numeric, Numeric, Numeric, Numeric X3, Numeric T0, Numeric T) {
  return pow(T0 / T, X3);
}

Numeric dDPL_dX3(
    Numeric, Numeric, Numeric X2, Numeric X3, Numeric T0, Numeric T) {
  return X2 * pow(T0 / T, X3) * log(T0 / T);
}

Numeric dDPL_dT0(Numeric X0,
                 Numeric X1,
                 Numeric X2,
                 Numeric X3,
                 Numeric T0,
                 Numeric T) {
  return X0 * X1 * pow(T0 / T, X1) / T0 + X2 * X3 * pow(T0 / T, X3) / T0;
}

Numeric dDPL_dT(Numeric X0,
                Numeric X1,
                Numeric X2,
                Numeric X3,
                Numeric T0,
                Numeric T) {
  return -X0 * X1 * pow(T0 / T, X1) / T + -X2 * X3 * pow(T0 / T, X3) / T;
}

Numeric POLY(const ExhaustiveConstVectorView& x, Numeric T) {
  Numeric poly_fac = 1.0;
  Numeric poly_sum = 0.0;
  for (auto X : x) {
    poly_sum += X * poly_fac;
    poly_fac *= T;
  }
  return poly_sum;
}

Numeric dPOLY_dT(const ExhaustiveConstVectorView& x, Numeric T) {
  Numeric poly_fac = 1.0;
  Numeric poly_sum = 0.0;
  for (Index i = 1; i < x.size(); ++i) {
    poly_sum += static_cast<Numeric>(i) * x[i] * poly_fac;
    poly_fac *= T;
  }
  return poly_sum;
}
}  // namespace model

Numeric data::operator()(Numeric T0, Numeric T) const {
#define SWITCHCASE(mod) \
  case LineShapeModelType::mod: \
    return operator()<LineShapeModelType::mod>(T0, T)

  switch (t) {
    SWITCHCASE(T0);
    SWITCHCASE(T1);
    SWITCHCASE(T2);
    SWITCHCASE(T3);
    SWITCHCASE(T4);
    SWITCHCASE(T5);
    SWITCHCASE(AER);
    SWITCHCASE(DPL);
    SWITCHCASE(POLY);
  }

  return NAN;

#undef SWITCHCASE
}

#define SWITCHCASE(name, mod) \
  case LineShapeModelType::mod:       \
    return d##name<LineShapeModelType::mod>(T0, T)

#define DERIVATIVES(name)                                            \
  Numeric data::d##name(Numeric T0, Numeric T) const { \
    switch (t) {                                                     \
      SWITCHCASE(name, T0);                                          \
      SWITCHCASE(name, T1);                                          \
      SWITCHCASE(name, T2);                                          \
      SWITCHCASE(name, T3);                                          \
      SWITCHCASE(name, T4);                                          \
      SWITCHCASE(name, T5);                                          \
      SWITCHCASE(name, AER);                                         \
      SWITCHCASE(name, DPL);                                         \
      SWITCHCASE(name, POLY);                                        \
    }                                                                \
                                                                     \
    return std::numeric_limits<Numeric>::quiet_NaN();                \
  }

DERIVATIVES(X0)
DERIVATIVES(X1)
DERIVATIVES(X2)
DERIVATIVES(X3)
DERIVATIVES(T0)
DERIVATIVES(T)

#undef DERIVATIVES
#undef SWITCHCASE

std::ostream& operator<<(std::ostream& os, const temperature::data& x) {
  os << x.t;

  if (model_size[static_cast<Size>(x.t)] == std::numeric_limits<Size>::max()) {
    os << ' ' << x.x.size();
  }

  return os << ' ' << x.x;
}

std::istream& operator>>(std::istream& is, temperature::data& x) {
  String name;
  is >> name;
  x.t = to<LineShapeModelType>(name);

  Size n = model_size[static_cast<Size>(x.t)];
  if (n == std::numeric_limits<Size>::max()) {
    is >> n;
  }
  x.x.resize(n);

  for (auto& v : x.x) is >> double_imanip() >> v;
  return is;
}

LineShapeModelType data::Type() const { return t; }

const Vector& data::X() const { return x; }

bool data::is_zero() const {
  switch (t) {
    using enum LineShapeModelType;
    case T0:
      [[fallthrough]];
    case T1:
      [[fallthrough]];
    case T2:
      [[fallthrough]];
    case T5:
      return x[0] == 0;
    case T3:
      [[fallthrough]];
    case T4:
      return x[0] == 0 and x[1] == 0;
    case AER:
      [[fallthrough]];
    case POLY:
      for (auto y : x)
        if (y != 0) return false;
      return true;
    case DPL:
      return x[0] == 0 and x[2] == 0;
  }
  return true;
}
}  // namespace lbl::temperature
