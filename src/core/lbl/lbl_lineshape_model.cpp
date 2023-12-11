#include "lbl_lineshape_model.h"

#include "debug.h"
#include "double_imanip.h"
#include "enums.h"
#include "species.h"

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
    if (single_models.back().species == Species::Species::Bath)               \
      return x - single_models.back().mod(T0, atm.temperature, atm.pressure); \
    const Numeric t =                                                         \
        std::transform_reduce(single_models.begin(),                          \
                              single_models.end(),                            \
                              0.0,                                            \
                              std::plus<>{},                                  \
                              [&atm](auto& m) { return atm[m.species]; });    \
    return (t - x) / t * t;                                                   \
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
      out += (1.0 - vmr) *                                                   \
             m.d##mod##_d##deriv(T0, atm.temperature, atm.pressure);         \
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

#undef VARIABLE

#undef DERIVATIVE

#define DERIVATIVE(mod, deriv)                                                 \
  Numeric model::d##mod##_d##deriv(const AtmPoint& atm, const Size spec)       \
      const noexcept {                                                         \
    if (spec >= single_models.size()) return 0.0;                              \
    const auto& m = single_models[spec];                                       \
                                                                               \
    const Numeric x = m.d##mod##_d##deriv(T0, atm.temperature, atm.pressure);  \
                                                                               \
    if (spec == single_models.size() - 1 and                                   \
        single_models.back().species == Species::Species::Bath) {              \
      const Numeric vmr =                                                      \
          std::transform_reduce(single_models.begin(),                         \
                                single_models.end() - 1,                       \
                                0.0,                                           \
                                std::plus<>{},                                 \
                                [&atm](auto& ml) { return atm[ml.species]; }); \
      return (1 - vmr) * x;                                                    \
    }                                                                          \
                                                                               \
    if (single_models.back().species == Species::Species::Bath) {              \
      return atm[m.species] * x;                                               \
    }                                                                          \
                                                                               \
    const Numeric vmr =                                                        \
        std::transform_reduce(single_models.begin(),                           \
                              single_models.end(),                             \
                              0.0,                                             \
                              std::plus<>{},                                   \
                              [&atm](auto& ml) { return atm[ml.species]; });   \
    return x * atm[m.species] / vmr;                                           \
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

VARIABLE(X0);
VARIABLE(X1);
VARIABLE(X2);
VARIABLE(X3);

#undef VARIABLE

#undef DERIVATIVE

#define VARIABLE(name)                                                      \
  Numeric model::d##name##_dX(                                              \
      const AtmPoint& atm, const Size spec, temperature::coefficient coeff) \
      const noexcept {                                                      \
    switch (coeff) {                                                        \
      case temperature::coefficient::X0:                                    \
        return d##name##_dX0(atm, spec);                                    \
      case temperature::coefficient::X1:                                    \
        return d##name##_dX1(atm, spec);                                    \
      case temperature::coefficient::X2:                                    \
        return d##name##_dX2(atm, spec);                                    \
      case temperature::coefficient::X3:                                    \
        return d##name##_dX3(atm, spec);                                    \
      case temperature::coefficient::FINAL:;                                \
    }                                                                       \
    return 0.0;                                                             \
  }                                                                         \
  Numeric species_model::d##name##_dX(                                      \
      Numeric T0, Numeric T, Numeric P, temperature::coefficient coeff)     \
      const noexcept {                                                      \
    switch (coeff) {                                                        \
      case temperature::coefficient::X0:                                    \
        return d##name##_dX0(T0, T, P);                                     \
      case temperature::coefficient::X1:                                    \
        return d##name##_dX1(T0, T, P);                                     \
      case temperature::coefficient::X2:                                    \
        return d##name##_dX2(T0, T, P);                                     \
      case temperature::coefficient::X3:                                    \
        return d##name##_dX3(T0, T, P);                                     \
      case temperature::coefficient::FINAL:;                                \
    }                                                                       \
    return 0.0;                                                             \
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

std::ostream& operator<<(std::ostream& os,
                         const std::pair<variable, temperature::data>& x) {
  return os << x.first << ' ' << x.second;
}

std::istream& operator>>(std::istream& is,
                         std::pair<variable, temperature::data>& x) {
  String name;
  is >> name >> x.second;

  x.first = tovariableOrThrow(name);

  return is;
}

std::ostream& operator<<(
    std::ostream& os,
    const std::vector<std::pair<variable, temperature::data>>& x) {
  os << x.size();
  for (auto& y : x) os << ' ' << y;
  return os;
}

std::istream& operator>>(
    std::istream& is, std::vector<std::pair<variable, temperature::data>>& x) {
  Size n;
  is >> n;
  x.resize(n);
  for (auto& y : x) is >> y;
  return is;
}

std::ostream& operator<<(std::ostream& os, const species_model& x) {
  return os << toShortName(x.species) << ' ' << x.data;
}

std::istream& operator>>(std::istream& is, species_model& x) {
  String name;
  is >> name >> x.data;
  x.species = Species::toSpeciesEnumOrThrow(name);
  return is;
}

std::ostream& operator<<(std::ostream& os,
                         const std::vector<species_model>& x) {
  os << x.size();
  for (auto& y : x) os << ' ' << y;
  return os;
}

std::istream& operator>>(std::istream& is, std::vector<species_model>& x) {
  Size n;
  is >> n;
  x.resize(n);
  for (auto& y : x) is >> y;
  return is;
}

std::ostream& operator<<(std::ostream& os, const model& x) {
  return os << x.T0 << ' ' << x.one_by_one << ' ' << x.single_models;
}

std::istream& operator>>(std::istream& is, model& x) {
  Index i;
  is >> double_imanip() >> x.T0;
  is >> i;
  x.one_by_one = static_cast<bool>(i);
  return is >> x.single_models;
}

void model::clear_zeroes() {
  auto is_zero = [](auto& x) { return x.second.is_zero(); };
  for (auto& m : single_models) {
    auto ptr = std::ranges::find_if(m.data, is_zero);
    while (ptr != m.data.end()) {
      m.data.erase(ptr);
      ptr = std::ranges::find_if(m.data, is_zero);
    }
  }
}
}  // namespace lbl::line_shape
