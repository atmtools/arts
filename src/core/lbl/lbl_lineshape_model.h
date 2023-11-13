#pragma once

#include <matpack.h>

#include <ostream>
#include <vector>

#include "atm.h"
#include "enums.h"
#include "lbl_temperature_model.h"
#include "species.h"

namespace lbl::line_shape {
ENUMCLASS(variable,
          char,
          G0,   // Pressure broadening speed-independent
          D0,   // Pressure f-shifting speed-dependent
          G2,   // Pressure broadening speed-dependent
          D2,   // Pressure f-shifting speed-independent
          FVC,  // Frequency of velocity-changing collisions
          ETA,  // Correlation
          Y,    // First order line mixing coefficient
          G,    // Second order line mixing coefficient
          DV    // Second order line mixing f-shifting
)

struct species_model {
  Species::Species species{};

  std::vector<std::pair<variable, temperature::data>> data{};

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
};

struct model {
  bool one_by_one{false};

  Numeric T0{0};

  std::vector<species_model> single_models{};

  friend std::ostream& operator<<(std::ostream& os, const model& x);

#define VARIABLE(name)                                            \
  [[nodiscard]] Numeric name(const AtmPoint& atm) const noexcept; \
                                                                  \
  [[nodiscard]] Numeric d##name##_dVMR(                           \
      const AtmPoint& atm, Species::Species species) const noexcept

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
  DERIVATIVE(X0);
  DERIVATIVE(X1);
  DERIVATIVE(X2);
  DERIVATIVE(X3);

#undef DERIVATIVE
};

  std::ostream& operator<<(std::ostream& os, const std::vector<species_model>& x);
}  // namespace lbl::line_shape
