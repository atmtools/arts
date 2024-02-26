#pragma once

#include <matpack.h>

#include <ostream>
#include <vector>

#include "atm.h"
#include "enums.h"
#include "lbl_temperature_model.h"
#include "species.h"

namespace lbl::line_shape {
struct species_model {
  SpeciesEnum species{};

  std::vector<std::pair<LineShapeModelVariable, temperature::data>> data{};

#define VARIABLE(name) \
  [[nodiscard]] Numeric name(Numeric T0, Numeric T, Numeric P) const noexcept

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

#define DERIVATIVE(name)                                               \
  [[nodiscard]] Numeric dG0_d##name(Numeric T0, Numeric T, Numeric P)  \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dD0_d##name(Numeric T0, Numeric T, Numeric P)  \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dG2_d##name(Numeric T0, Numeric T, Numeric P)  \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dD2_d##name(Numeric T0, Numeric T, Numeric P)  \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dETA_d##name(Numeric T0, Numeric T, Numeric P) \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dFVC_d##name(Numeric T0, Numeric T, Numeric P) \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dY_d##name(Numeric T0, Numeric T, Numeric P)   \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dG_d##name(Numeric T0, Numeric T, Numeric P)   \
      const noexcept;                                                  \
  [[nodiscard]] Numeric dDV_d##name(Numeric T0, Numeric T, Numeric P)  \
      const noexcept

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
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dD0_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dG2_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dD2_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dETA_dX(Numeric T0,
                                Numeric T,
                                Numeric P,
                                LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dFVC_dX(Numeric T0,
                                Numeric T,
                                Numeric P,
                                LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dY_dX(Numeric T0,
                              Numeric T,
                              Numeric P,
                              LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dG_dX(Numeric T0,
                              Numeric T,
                              Numeric P,
                              LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dDV_dX(Numeric T0,
                               Numeric T,
                               Numeric P,
                               LineShapeModelCoefficient coeff) const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const species_model& x);
  friend std::istream& operator>>(std::istream& os, species_model& x);
};

struct model {
  bool one_by_one{false};

  Numeric T0{0};

  std::vector<species_model> single_models{};

  friend std::ostream& operator<<(std::ostream& os, const model& x);
  friend std::istream& operator>>(std::istream& is, model& x);

#define VARIABLE(name)                                            \
  [[nodiscard]] Numeric name(const AtmPoint& atm) const noexcept; \
                                                                  \
  [[nodiscard]] Numeric d##name##_dVMR(const AtmPoint& atm,       \
                                       SpeciesEnum species) const noexcept

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

#define DERIVATIVE(name)                                                  \
  [[nodiscard]] Numeric dG0_d##name(const AtmPoint& atm) const noexcept;  \
  [[nodiscard]] Numeric dD0_d##name(const AtmPoint& atm) const noexcept;  \
  [[nodiscard]] Numeric dG2_d##name(const AtmPoint& atm) const noexcept;  \
  [[nodiscard]] Numeric dD2_d##name(const AtmPoint& atm) const noexcept;  \
  [[nodiscard]] Numeric dETA_d##name(const AtmPoint& atm) const noexcept; \
  [[nodiscard]] Numeric dFVC_d##name(const AtmPoint& atm) const noexcept; \
  [[nodiscard]] Numeric dY_d##name(const AtmPoint& atm) const noexcept;   \
  [[nodiscard]] Numeric dG_d##name(const AtmPoint& atm) const noexcept;   \
  [[nodiscard]] Numeric dDV_d##name(const AtmPoint& atm) const noexcept

  DERIVATIVE(T);
  DERIVATIVE(T0);
  DERIVATIVE(P);

#undef DERIVATIVE

#define DERIVATIVE(name)                                                   \
  [[nodiscard]] Numeric dG0_d##name(const AtmPoint& atm, const Size spec)  \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dD0_d##name(const AtmPoint& atm, const Size spec)  \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dG2_d##name(const AtmPoint& atm, const Size spec)  \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dD2_d##name(const AtmPoint& atm, const Size spec)  \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dETA_d##name(const AtmPoint& atm, const Size spec) \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dFVC_d##name(const AtmPoint& atm, const Size spec) \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dY_d##name(const AtmPoint& atm, const Size spec)   \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dG_d##name(const AtmPoint& atm, const Size spec)   \
      const noexcept;                                                      \
  [[nodiscard]] Numeric dDV_d##name(const AtmPoint& atm, const Size spec)  \
      const noexcept

  DERIVATIVE(X0);
  DERIVATIVE(X1);
  DERIVATIVE(X2);
  DERIVATIVE(X3);

#undef DERIVATIVE

#undef DERIVATIVE

  [[nodiscard]] Numeric dG0_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dD0_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dG2_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dD2_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dETA_dX(const AtmPoint& atm,
                                const Size spec,
                                LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dFVC_dX(const AtmPoint& atm,
                                const Size spec,
                                LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dY_dX(const AtmPoint& atm,
                              const Size spec,
                              LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dG_dX(const AtmPoint& atm,
                              const Size spec,
                              LineShapeModelCoefficient coeff) const noexcept;

  [[nodiscard]] Numeric dDV_dX(const AtmPoint& atm,
                               const Size spec,
                               LineShapeModelCoefficient coeff) const noexcept;

  //! Remove all line shape variables that evaluate unconditionally to 0
  void clear_zeroes();
};

std::ostream& operator<<(std::ostream& os, const std::vector<species_model>& x);
std::istream& operator>>(std::istream& is, std::vector<species_model>& x);
}  // namespace lbl::line_shape

std::ostream& operator<<(
    std::ostream& os,
    const std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>& x);

std::istream& operator>>(
    std::istream& is,
    std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>& x);
