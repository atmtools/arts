#include "lbl_lineshape_model.h"

namespace lbl::line_shape {
#define VARIABLE(name, PVAR, DPVAR)                                       \
  Numeric species_model::name(                                            \
      Numeric T0, Numeric T, Numeric P [[maybe_unused]]) const noexcept { \
    auto ptr = std::ranges::find_if(                                      \
        data, [](auto& x) { return x.first == variable::name; });         \
    return ptr != data.end() ? PVAR * ptr->second(T0, T) : 0.0;           \
  }                                                                       \
                                                                          \
  Numeric species_model::d##name##_dP(                                    \
      Numeric T0, Numeric T, Numeric P [[maybe_unused]]) const noexcept { \
    auto ptr = std::ranges::find_if(                                      \
        data, [](auto& x) { return x.first == variable::name; });         \
    return ptr != data.end() ? DPVAR * ptr->second(T0, T) : 0.0;          \
  }

VARIABLE(G0, P, 1);
VARIABLE(D0, P, 1);
VARIABLE(G2, P, 1);
VARIABLE(D2, P, 1);
VARIABLE(FVC, P, 1);
VARIABLE(ETA, 1, 0);
VARIABLE(Y, P, 1);
VARIABLE(G, P* P, 2 * P);
VARIABLE(DV, P* P, 2 * P);

#undef VARIABLE

#define DERIVATIVE(name, deriv, PVAR)                                     \
  Numeric species_model::d##name##_d##deriv(                              \
      Numeric T0, Numeric T, Numeric P [[maybe_unused]]) const noexcept { \
    auto ptr = std::ranges::find_if(                                      \
        data, [](auto& x) { return x.first == variable::name; });         \
    return ptr != data.end() ? PVAR * ptr->second.d##deriv(T0, T) : 0.0;  \
  }

#define VARIABLE(deriv)        \
  DERIVATIVE(G0, deriv, P);    \
  DERIVATIVE(D0, deriv, P);    \
  DERIVATIVE(G2, deriv, P);    \
  DERIVATIVE(D2, deriv, P);    \
  DERIVATIVE(FVC, deriv, P);   \
  DERIVATIVE(Y, deriv, P);     \
  DERIVATIVE(G, deriv, P* P);  \
  DERIVATIVE(DV, deriv, P* P); \
  DERIVATIVE(ETA, deriv, 1);

VARIABLE(T);
VARIABLE(T0);
VARIABLE(X0);
VARIABLE(X1);
VARIABLE(X2);
VARIABLE(X3);

#undef VARIABLE

#undef DERIVATIVE

/////////////////////////////////////////////////////////////////////////////

#define VARIABLE(mod)                                                         \
  Numeric model::mod(const AtmPoint& atm) const noexcept {                    \
    Numeric vmr = 0.0;                                                        \
    Numeric out = 0.0;                                                        \
                                                                              \
    const auto compute = [&](const species_model& m) {                        \
      const Numeric this_vmr = atm[m.species];                                \
      vmr += this_vmr;                                                        \
      out += this_vmr * m.mod(T0, atm.temperature, atm.pressure);             \
    };                                                                        \
                                                                              \
    for (auto& m :                                                            \
         std::ranges::take_view(single_models, single_models.size() - 1)) {   \
      compute(m);                                                             \
    }                                                                         \
                                                                              \
    if (const auto& m = single_models.back();                                 \
        m.species == Species::Species::Bath) {                                \
      out += (1.0 - vmr) * m.mod(T0, atm.temperature, atm.pressure);          \
    } else {                                                                  \
      compute(m);                                                             \
      out /= vmr;                                                             \
    }                                                                         \
                                                                              \
    return out;                                                               \
  }                                                                           \
                                                                              \
  [[nodiscard]] Numeric model::d##mod##_dVMR(                                 \
      const AtmPoint& atm, Species ::Species species) const noexcept {        \
    auto ptr = std::ranges::find_if(                                          \
        single_models, [&](auto& m) { return m.species == species; });        \
                                                                              \
    if (ptr == single_models.end()) {                                         \
      return 0.0;                                                             \
    }                                                                         \
                                                                              \
    const Numeric x = ptr->mod(T0, atm.temperature, atm.pressure);            \
                                                                              \
    if (species == Species::Species::Bath) return -x;                         \
    if (single_models.back().species == Species::Species::Bath) return x;     \
    return x /                                                                \
           std::transform_reduce(single_models.begin(),                       \
                                 single_models.end(),                         \
                                 0.0,                                         \
                                 std::plus<>{},                               \
                                 [&atm](auto& m) { return atm[m.species]; }); \
  }

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

#define DERIVATIVE(mod, deriv)                                               \
  Numeric model::d##mod##_d##deriv(const AtmPoint& atm) const noexcept {     \
    Numeric vmr = 0.0;                                                       \
    Numeric out = 0.0;                                                       \
                                                                             \
    const auto compute = [&](const species_model& m) {                       \
      const Numeric this_vmr = atm[m.species];                               \
      vmr += this_vmr;                                                       \
      out +=                                                                 \
          this_vmr * m.d##mod##_d##deriv(T0, atm.temperature, atm.pressure); \
    };                                                                       \
                                                                             \
    for (auto& m :                                                           \
         std::ranges::take_view(single_models, single_models.size() - 1)) {  \
      compute(m);                                                            \
    }                                                                        \
                                                                             \
    if (const auto& m = single_models.back();                                \
        m.species == Species::Species::Bath) {                               \
      out += (1.0 - vmr) * m.mod(T0, atm.temperature, atm.pressure);         \
    } else {                                                                 \
      compute(m);                                                            \
      out /= vmr;                                                            \
    }                                                                        \
                                                                             \
    return out;                                                              \
  }

#define VARIABLE(deriv)   \
  DERIVATIVE(G0, deriv);  \
  DERIVATIVE(D0, deriv);  \
  DERIVATIVE(G2, deriv);  \
  DERIVATIVE(D2, deriv);  \
  DERIVATIVE(FVC, deriv); \
  DERIVATIVE(Y, deriv);   \
  DERIVATIVE(G, deriv);   \
  DERIVATIVE(DV, deriv);  \
  DERIVATIVE(ETA, deriv);

VARIABLE(P);
VARIABLE(T);
VARIABLE(T0);
VARIABLE(X0);
VARIABLE(X1);
VARIABLE(X2);
VARIABLE(X3);

#undef VARIABLE

#undef DERIVATIVE
}  // namespace lbl::line_shape
