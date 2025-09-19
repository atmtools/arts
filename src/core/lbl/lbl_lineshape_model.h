#pragma once

#include <array.h>
#include <atm.h>
#include <enumsLineShapeModelCoefficient.h>
#include <enumsLineShapeModelVariable.h>
#include <enumsSpeciesEnum.h>
#include <matpack.h>

#include <format>
#include <iosfwd>
#include <unordered_map>

#include "lbl_temperature_model.h"

namespace lbl::line_shape {
struct species_model {
  using map_t = std::unordered_map<LineShapeModelVariable, temperature::data>;
  map_t data{};

  //! Removes the variables from the model.
  template <LineShapeModelVariable... V>
  map_t::size_type remove_variables() {
    return (... + std::erase_if(
                      data, [v = V](const auto& x) { return x.first == v; }));
  }

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

  friend std::istream& operator>>(std::istream& os, species_model& x);
};

struct model {
  bool one_by_one{false};

  Numeric T0{0};

  using map_t = std::unordered_map<SpeciesEnum, species_model>;
  map_t single_models{};

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

#define DERIVATIVE(name)                                            \
  [[nodiscard]] Numeric dG0_d##name(const AtmPoint& atm,            \
                                    const SpeciesEnum spec) const;  \
  [[nodiscard]] Numeric dD0_d##name(const AtmPoint& atm,            \
                                    const SpeciesEnum spec) const;  \
  [[nodiscard]] Numeric dG2_d##name(const AtmPoint& atm,            \
                                    const SpeciesEnum spec) const;  \
  [[nodiscard]] Numeric dD2_d##name(const AtmPoint& atm,            \
                                    const SpeciesEnum spec) const;  \
  [[nodiscard]] Numeric dETA_d##name(const AtmPoint& atm,           \
                                     const SpeciesEnum spec) const; \
  [[nodiscard]] Numeric dFVC_d##name(const AtmPoint& atm,           \
                                     const SpeciesEnum spec) const; \
  [[nodiscard]] Numeric dY_d##name(const AtmPoint& atm,             \
                                   const SpeciesEnum spec) const;   \
  [[nodiscard]] Numeric dG_d##name(const AtmPoint& atm,             \
                                   const SpeciesEnum spec) const;   \
  [[nodiscard]] Numeric dDV_d##name(const AtmPoint& atm,            \
                                    const SpeciesEnum spec) const

  DERIVATIVE(X0);
  DERIVATIVE(X1);
  DERIVATIVE(X2);
  DERIVATIVE(X3);

#undef DERIVATIVE

#undef DERIVATIVE

  [[nodiscard]] Numeric dG0_dX(const AtmPoint& atm,
                               const SpeciesEnum spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dD0_dX(const AtmPoint& atm,
                               const SpeciesEnum spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dG2_dX(const AtmPoint& atm,
                               const SpeciesEnum spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dD2_dX(const AtmPoint& atm,
                               const SpeciesEnum spec,
                               LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dETA_dX(const AtmPoint& atm,
                                const SpeciesEnum spec,
                                LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dFVC_dX(const AtmPoint& atm,
                                const SpeciesEnum spec,
                                LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dY_dX(const AtmPoint& atm,
                              const SpeciesEnum spec,
                              LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dG_dX(const AtmPoint& atm,
                              const SpeciesEnum spec,
                              LineShapeModelCoefficient coeff) const;

  [[nodiscard]] Numeric dDV_dX(const AtmPoint& atm,
                               const SpeciesEnum spec,
                               LineShapeModelCoefficient coeff) const;

  //! Remove all line shape variables that evaluate unconditionally to 0
  void clear_zeroes();
};

std::istream& operator>>(std::istream& is,
                         std::unordered_map<SpeciesEnum, species_model>& x);
}  // namespace lbl::line_shape

template <>
struct std::formatter<lbl::line_shape::species_model> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line_shape::species_model& v,
                              FmtContext& ctx) const {
    if (tags.help) {
      tags.format(ctx, "Data: "sv, v.data);
    } else if (tags.io) {
      tags.format(ctx, v.data.size(), ' ', v.data);
    } else {
      tags.add_if_bracket(ctx, '[');
      tags.format(ctx, v.data);
      tags.add_if_bracket(ctx, ']');
    }

    return ctx.out();
  }
};

template <>
struct std::formatter<lbl::line_shape::model> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::line_shape::model& v,
                              FmtContext& ctx) const {
    if (tags.help) {
      return tags.format(ctx,
                         "Reference temperature: "sv,
                         v.T0,
                         " K; One-by-one: "sv,
                         v.one_by_one ? "<on>"sv : "<off>"sv,
                         "; Single models: "sv,
                         v.single_models);
    }

    if (tags.io) {
      return tags.format(ctx,
                         v.T0,
                         ' ',
                         Index{v.one_by_one},
                         ' ',
                         v.single_models.size(),
                         ' ',
                         v.single_models);
    }

    const auto sep = tags.sep();
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.one_by_one, sep, v.T0, sep, v.single_models);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
std::optional<std::string>
to_helper_string<lbl::line_shape::species_model::map_t>(
    const lbl::line_shape::species_model::map_t&);
