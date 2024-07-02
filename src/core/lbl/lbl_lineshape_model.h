#pragma once

#include <matpack.h>

#include <format>
#include <ostream>
#include <vector>

#include "array.h"
#include "atm.h"
#include "enums.h"
#include "lbl_temperature_model.h"

namespace lbl::line_shape {
struct species_model {
  SpeciesEnum species{};

  std::vector<std::pair<LineShapeModelVariable, temperature::data>> data{};

#define VARIABLE(name) \
  [[nodiscard]] Numeric name(Numeric T0, Numeric T, Numeric P) const

  VARIABLE(G0);
  VARIABLE(D0);
  VARIABLE(G2);
  VARIABLE(D2);
  VARIABLE(ETA);
  VARIABLE(FVC);
  VARIABLE(Y);
  VARIABLE(G);
  VARIABLE(DV);

#undef VARIABLE

#define DERIVATIVE(name)                                                      \
  [[nodiscard]] Numeric dG0_d##name(Numeric T0, Numeric T, Numeric P) const;  \
  [[nodiscard]] Numeric dD0_d##name(Numeric T0, Numeric T, Numeric P) const;  \
  [[nodiscard]] Numeric dG2_d##name(Numeric T0, Numeric T, Numeric P) const;  \
  [[nodiscard]] Numeric dD2_d##name(Numeric T0, Numeric T, Numeric P) const;  \
  [[nodiscard]] Numeric dETA_d##name(Numeric T0, Numeric T, Numeric P) const; \
  [[nodiscard]] Numeric dFVC_d##name(Numeric T0, Numeric T, Numeric P) const; \
  [[nodiscard]] Numeric dY_d##name(Numeric T0, Numeric T, Numeric P) const;   \
  [[nodiscard]] Numeric dG_d##name(Numeric T0, Numeric T, Numeric P) const;   \
  [[nodiscard]] Numeric dDV_d##name(Numeric T0, Numeric T, Numeric P) const

  DERIVATIVE(P);
  DERIVATIVE(T);
  DERIVATIVE(T0);
  DERIVATIVE(X0);
  DERIVATIVE(X1);
  DERIVATIVE(X2);
  DERIVATIVE(X3);

#undef DERIVATIVE

  [[nodiscard]] Numeric dG0_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dD0_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dG2_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dD2_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dETA_dX(Numeric T0,
                                Numeric T,
                                Numeric P,
                                LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dFVC_dX(Numeric T0,
                                Numeric T,
                                Numeric P,
                                LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dY_dX(Numeric T0,
                              Numeric T,
                              Numeric P,
                              LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dG_dX(Numeric T0,
                              Numeric T,
                              Numeric P,
                              LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dDV_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const;

  friend std::ostream& operator<<(std::ostream& os, const species_model& x);
  friend std::istream& operator>>(std::istream& os, species_model& x);
};

struct model {
  bool one_by_one{false};

  Numeric T0{0};

  std::vector<species_model> single_models{};

  friend std::ostream& operator<<(std::ostream& os, const model& x);
  friend std::istream& operator>>(std::istream& is, model& x);

#define VARIABLE(name)                                      \
  [[nodiscard]] Numeric name(const AtmPoint& atm) const;    \
                                                            \
  [[nodiscard]] Numeric d##name##_dVMR(const AtmPoint& atm, \
                                       SpeciesEnum species) const

  VARIABLE(G0);
  VARIABLE(D0);
  VARIABLE(G2);
  VARIABLE(D2);
  VARIABLE(ETA);
  VARIABLE(FVC);
  VARIABLE(Y);
  VARIABLE(G);
  VARIABLE(DV);

#undef VARIABLE

#define DERIVATIVE(name)                                         \
  [[nodiscard]] Numeric dG0_d##name(const AtmPoint& atm) const;  \
  [[nodiscard]] Numeric dD0_d##name(const AtmPoint& atm) const;  \
  [[nodiscard]] Numeric dG2_d##name(const AtmPoint& atm) const;  \
  [[nodiscard]] Numeric dD2_d##name(const AtmPoint& atm) const;  \
  [[nodiscard]] Numeric dETA_d##name(const AtmPoint& atm) const; \
  [[nodiscard]] Numeric dFVC_d##name(const AtmPoint& atm) const; \
  [[nodiscard]] Numeric dY_d##name(const AtmPoint& atm) const;   \
  [[nodiscard]] Numeric dG_d##name(const AtmPoint& atm) const;   \
  [[nodiscard]] Numeric dDV_d##name(const AtmPoint& atm) const

  DERIVATIVE(T);
  DERIVATIVE(T0);
  DERIVATIVE(P);

#undef DERIVATIVE

#define DERIVATIVE(name)                                                   \
  [[nodiscard]] Numeric dG0_d##name(const AtmPoint& atm, const Size spec)  \
      const;                                                               \
  [[nodiscard]] Numeric dD0_d##name(const AtmPoint& atm, const Size spec)  \
      const;                                                               \
  [[nodiscard]] Numeric dG2_d##name(const AtmPoint& atm, const Size spec)  \
      const;                                                               \
  [[nodiscard]] Numeric dD2_d##name(const AtmPoint& atm, const Size spec)  \
      const;                                                               \
  [[nodiscard]] Numeric dETA_d##name(const AtmPoint& atm, const Size spec) \
      const;                                                               \
  [[nodiscard]] Numeric dFVC_d##name(const AtmPoint& atm, const Size spec) \
      const;                                                               \
  [[nodiscard]] Numeric dY_d##name(const AtmPoint& atm, const Size spec)   \
      const;                                                               \
  [[nodiscard]] Numeric dG_d##name(const AtmPoint& atm, const Size spec)   \
      const;                                                               \
  [[nodiscard]] Numeric dDV_d##name(const AtmPoint& atm, const Size spec) const

  DERIVATIVE(X0);
  DERIVATIVE(X1);
  DERIVATIVE(X2);
  DERIVATIVE(X3);

#undef DERIVATIVE

#undef DERIVATIVE

  [[nodiscard]] Numeric dG0_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dD0_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dG2_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dD2_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dETA_dX(const AtmPoint& atm,
                                const Size spec,
                                LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dFVC_dX(const AtmPoint& atm,
                                const Size spec,
                                LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dY_dX(const AtmPoint& atm,
                              const Size spec,
                              LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dG_dX(const AtmPoint& atm,
                              const Size spec,
                              LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dDV_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const;

  //! Remove all line shape variables that evaluate unconditionally to 0
  void clear_zeroes();
};

std::ostream& operator<<(std::ostream& os, const std::vector<species_model>& x);
std::istream& operator>>(std::istream& is, std::vector<species_model>& x);
}  // namespace lbl::line_shape

std::ostream& operator<<(
    std::ostream& os,
    const std::vector<
        std::pair<LineShapeModelVariable, lbl::temperature::data>>& x);

std::istream& operator>>(
    std::istream& is,
    std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>& x);

template <>
struct std::formatter<
    std::pair<LineShapeModelVariable, lbl::temperature::data>> {
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
  FmtContext::iterator format(
      const std::pair<LineShapeModelVariable, lbl::temperature::data>& v,
      FmtContext& ctx) const {
    if (tags.bracket) std::ranges::copy("["sv, ctx.out());

    std::formatter<LineShapeModelVariable> p1{};
    std::formatter<lbl::temperature::data> p2{};
    make_compat(p1, p2);

    p1.format(v.first, ctx);
    if (tags.comma) std::ranges::copy(","sv, ctx.out());
    std::ranges::copy(" "sv, ctx.out());
    p2.format(v.second, ctx);

    if (tags.bracket) std::ranges::copy("]"sv, ctx.out());
    return ctx.out();
  }
};

template <>
struct std::formatter<lbl::line_shape::species_model> {
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
  FmtContext::iterator format(const lbl::line_shape::species_model& v,
                              FmtContext& ctx) const {
    if (tags.bracket) std::ranges::copy("["sv, ctx.out());

    std::formatter<SpeciesEnum> species{};
    std::formatter<
        std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>>
        data{};
    make_compat(species, data);

    species.format(v.species, ctx);
    if (tags.comma) std::ranges::copy(","sv, ctx.out());
    std::ranges::copy(" "sv, ctx.out());
    data.format(v.data, ctx);

    if (tags.bracket) std::ranges::copy("]"sv, ctx.out());
    return ctx.out();
  }
};

template <>
struct std::formatter<lbl::line_shape::model> {
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
  FmtContext::iterator format(const lbl::line_shape::model& v,
                              FmtContext& ctx) const {
    if (tags.bracket) std::ranges::copy("["sv, ctx.out());

    if (tags.comma) {
      std::format_to(ctx.out(), "{}, {}, ", v.one_by_one, v.T0);
    } else {
      std::format_to(ctx.out(), "{} {} ", v.one_by_one, v.T0);
    }

    std::formatter<std::vector<lbl::line_shape::species_model>> s;
    make_compat(s);
    s.format(v.single_models, ctx);

    if (tags.bracket) std::ranges::copy("]"sv, ctx.out());
    return ctx.out();
  }
};
