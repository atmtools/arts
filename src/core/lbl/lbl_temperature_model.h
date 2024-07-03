#pragma once

#include <matpack.h>

#include <limits>

#include "debug.h"
#include "enums.h"
#include "matpack_view.h"

namespace lbl::temperature {
inline constexpr std::size_t LineShapeModelTypeSize = 9;
static_assert(LineShapeModelTypeSize == enumsize::LineShapeModelTypeSize,
              "Invalid number of models");

//! The number of parameters per model type (max() means and arbitrary number)
inline constexpr std::array<Size, LineShapeModelTypeSize> model_size{
    1,                                // T0
    2,                                // T1
    3,                                // T2
    2,                                // T3
    3,                                // T4
    2,                                // T5
    4,                                // AER
    4,                                // DPL
    std::numeric_limits<Size>::max()  // POLY
};

namespace model {
#define EMPTY(name, deriv)                                          \
  template <typename... T>                                          \
  [[nodiscard]] constexpr Numeric d##name##_d##deriv(const T&...) { \
    return 0;                                                       \
  }

constexpr Numeric T0(Numeric X0) { return X0; }
constexpr Numeric dT0_dX0(Numeric) { return 1; }
EMPTY(T0, X1)
EMPTY(T0, X2)
EMPTY(T0, X3)
EMPTY(T0, T0)
EMPTY(T0, T)

Numeric T1(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT1_dX0(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT1_dX1(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT1_dT0(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT1_dT(Numeric X0, Numeric X1, Numeric T0, Numeric T);
EMPTY(T1, X2)
EMPTY(T1, X3)

Numeric T2(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT2_dX0(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT2_dX1(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT2_dX2(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT2_dT0(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT2_dT(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
EMPTY(T2, X3)

constexpr Numeric T3(Numeric X0, Numeric X1, Numeric T0, Numeric T) {
  return X0 + X1 * (T - T0);
}
constexpr Numeric dT3_dX0(Numeric, Numeric, Numeric, Numeric) { return 1; }
constexpr Numeric dT3_dX1(Numeric, Numeric, Numeric T0, Numeric T) {
  return T - T0;
}
constexpr Numeric dT3_dT0(Numeric, Numeric X1, Numeric, Numeric) { return -X1; }
constexpr Numeric dT3_dT(Numeric, Numeric X1, Numeric, Numeric) { return X1; }
EMPTY(T3, X2)
EMPTY(T3, X3)

Numeric T4(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT4_dX0(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT4_dX1(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT4_dX2(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT4_dT0(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
Numeric dT4_dT(Numeric X0, Numeric X1, Numeric X2, Numeric T0, Numeric T);
EMPTY(T4, X2)
EMPTY(T4, X3)

Numeric T5(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT5_dX0(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT5_dX1(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT5_dT0(Numeric X0, Numeric X1, Numeric T0, Numeric T);
Numeric dT5_dT(Numeric X0, Numeric X1, Numeric T0, Numeric T);
EMPTY(T5, X2)
EMPTY(T5, X3)

constexpr Numeric AER(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T) {
  if (T < 250.0) return X0 + (T - 200.0) * (X1 - X0) / (250.0 - 200.0);
  if (T > 296.0) return X2 + (T - 296.0) * (X3 - X2) / (340.0 - 296.0);
  return X1 + (T - 250.0) * (X2 - X1) / (296.0 - 250.0);
}
constexpr Numeric dAER_dX0(Numeric, Numeric, Numeric, Numeric, Numeric T) {
  if (T < 250.0) return 1 - (T - 200.0) / (250.0 - 200.0);
  return 0;
}
constexpr Numeric dAER_dX1(Numeric, Numeric, Numeric, Numeric, Numeric T) {
  if (T < 250.0) return (T - 200.0) / (250.0 - 200.0);
  if (T > 296.0) return 0;
  return 1 - (T - 250.0) / (296.0 - 250.0);
}
constexpr Numeric dAER_dX2(Numeric, Numeric, Numeric, Numeric, Numeric T) {
  if (T < 250.0) return 0;
  if (T > 296.0) return 1 - (T - 296.0) / (340.0 - 296.0);
  return (T - 250.0) / (296.0 - 250.0);
}
constexpr Numeric dAER_dX3(Numeric, Numeric, Numeric, Numeric, Numeric T) {
  if (T > 296.0) return (T - 296.0) / (340.0 - 296.0);
  return 0;
}
constexpr Numeric dAER_dT(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T) {
  if (T < 250.0) return (X1 - X0) / (250.0 - 200.0);
  if (T > 296.0) return (X3 - X2) / (340.0 - 296.0);
  return (X2 - X1) / (296.0 - 250.0);
}
EMPTY(AER, T0)

Numeric DPL(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);
Numeric dDPL_dX0(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);
Numeric dDPL_dX1(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);
Numeric dDPL_dX2(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);
Numeric dDPL_dX3(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);
Numeric dDPL_dT0(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);
Numeric dDPL_dT(
    Numeric X0, Numeric X1, Numeric X2, Numeric X3, Numeric T0, Numeric T);

Numeric POLY(const ExhaustiveConstVectorView& x, Numeric T);
constexpr Numeric dPOLY_dX0(const ExhaustiveConstVectorView&, Numeric) {
  return 1.0;
}
constexpr Numeric dPOLY_dX1(const ExhaustiveConstVectorView&, Numeric T) {
  return T;
}
constexpr Numeric dPOLY_dX2(const ExhaustiveConstVectorView&, Numeric T) {
  return T * T;
}
constexpr Numeric dPOLY_dX3(const ExhaustiveConstVectorView&, Numeric T) {
  return T * T * T;
}
Numeric dPOLY_dT(const ExhaustiveConstVectorView& x, Numeric T);
EMPTY(POLY, T0)

#undef EMPTY
}  // namespace model

class data {
  LineShapeModelType t{LineShapeModelType::T0};
  Vector x{};

 public:
  [[nodiscard]] LineShapeModelType Type() const;
  [[nodiscard]] const Vector& X() const;

  friend std::ostream& operator<<(std::ostream& os, const temperature::data& x);
  friend std::istream& operator>>(std::istream& is, temperature::data& x);

  constexpr data(LineShapeModelType type = LineShapeModelType::T0,
                 Vector X                = {0.0})
      : t(type), x(std::move(X)) {
    ARTS_USER_ERROR_IF(not good_enum(t), "Invalid model type")
    ARTS_USER_ERROR_IF(auto n = model_size[static_cast<Size>(t)];
                       n != std::numeric_limits<Size>::max() and
                           n != static_cast<Size>(x.size()),
                       "Invalid number of parameters for model ",
                       t,
                       "\nExpected ",
                       n,
                       " parameters but got ",
                       x.size(),
                       " parameters.")
  }

  template <LineShapeModelType mod>
  [[nodiscard]] constexpr Numeric operator()(Numeric T0 [[maybe_unused]],
                                             Numeric T [[maybe_unused]]) const {
    if constexpr (mod == LineShapeModelType::T0)
      return model::T0(x[0]);
    else if constexpr (mod == LineShapeModelType::T1)
      return model::T1(x[0], x[1], T0, T);
    else if constexpr (mod == LineShapeModelType::T2)
      return model::T2(x[0], x[1], x[2], T0, T);
    else if constexpr (mod == LineShapeModelType::T3)
      return model::T3(x[0], x[1], T0, T);
    else if constexpr (mod == LineShapeModelType::T4)
      return model::T4(x[0], x[1], x[2], T0, T);
    else if constexpr (mod == LineShapeModelType::T5)
      return model::T5(x[0], x[1], T0, T);
    else if constexpr (mod == LineShapeModelType::AER)
      return model::AER(x[0], x[1], x[2], x[3], T);
    else if constexpr (mod == LineShapeModelType::DPL)
      return model::DPL(x[0], x[1], x[2], x[3], T0, T);
    else if constexpr (mod == LineShapeModelType::POLY)
      return model::POLY(x, T);
  }

#define DERIVATIVE(name)                                                      \
  template <LineShapeModelType mod>                                           \
  [[nodiscard]] constexpr Numeric d##name(Numeric T0 [[maybe_unused]],        \
                                          Numeric T [[maybe_unused]]) const { \
    if constexpr (mod == LineShapeModelType::T0)                              \
      return model::dT0_d##name(x[0]);                                        \
    else if constexpr (mod == LineShapeModelType::T1)                         \
      return model::dT1_d##name(x[0], x[1], T0, T);                           \
    else if constexpr (mod == LineShapeModelType::T2)                         \
      return model::dT2_d##name(x[0], x[1], x[2], T0, T);                     \
    else if constexpr (mod == LineShapeModelType::T3)                         \
      return model::dT3_d##name(x[0], x[1], T0, T);                           \
    else if constexpr (mod == LineShapeModelType::T4)                         \
      return model::dT4_d##name(x[0], x[1], x[2], T0, T);                     \
    else if constexpr (mod == LineShapeModelType::T5)                         \
      return model::dT5_d##name(x[0], x[1], T0, T);                           \
    else if constexpr (mod == LineShapeModelType::AER)                        \
      return model::dAER_d##name(x[0], x[1], x[2], x[3], T);                  \
    else if constexpr (mod == LineShapeModelType::DPL)                        \
      return model::dDPL_d##name(x[0], x[1], x[2], x[3], T0, T);              \
    else if constexpr (mod == LineShapeModelType::POLY)                       \
      return model::dPOLY_d##name(x, T);                                      \
  }

  DERIVATIVE(X0)
  DERIVATIVE(X1)
  DERIVATIVE(X2)
  DERIVATIVE(X3)
  DERIVATIVE(T0)
  DERIVATIVE(T)

#undef DERIVATIVE

  [[nodiscard]] Numeric operator()(Numeric T0, Numeric T) const;
  [[nodiscard]] Numeric dX0(Numeric T0, Numeric T) const;
  [[nodiscard]] Numeric dX1(Numeric T0, Numeric T) const;
  [[nodiscard]] Numeric dX2(Numeric T0, Numeric T) const;
  [[nodiscard]] Numeric dX3(Numeric T0, Numeric T) const;
  [[nodiscard]] Numeric dT0(Numeric T0, Numeric T) const;
  [[nodiscard]] Numeric dT(Numeric T0, Numeric T) const;

  [[nodiscard]] bool is_zero() const;
};
}  // namespace lbl::temperature

template <>
struct std::formatter<lbl::temperature::data> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::temperature::data& v,
                              FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');

    std::formatter<LineShapeModelType> t{};
    std::formatter<Vector> x{};
    make_compat(t, x);

    t.format(v.Type(), ctx);
    std::format_to(ctx.out(), "{}", tags.sep());
    x.format(v.X(), ctx);

    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};
