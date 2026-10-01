#include <arts_constants.h>
#include <arts_conversions.h>
#include <atm.h>
#include <jacobian.h>
#include <lbl_lineshape_voigt_ecs.h>
#include <wigner_functions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using lbl::voigt::ecs::ComputeData;
using Samples = std::vector<Complex>;

void require(bool condition, const std::string& message) {
  if (not condition) throw std::runtime_error(message);
}

void compare(
    const Samples& actual, const Samples& expected, Numeric relative, Numeric absolute, const std::string& context) {
  require(actual.size() == expected.size(), context + ": size mismatch");
  Numeric scale = 0;
  for (const auto value : expected) scale = std::max(scale, std::abs(value));
  for (Size i = 0; i < actual.size(); ++i) {
    const Numeric tolerance = absolute + relative * (std::abs(expected[i]) + 1e-3 * scale);
    require(std::isfinite(actual[i].real()) and std::isfinite(actual[i].imag()) and
                std::isfinite(expected[i].real()) and std::isfinite(expected[i].imag()) and
                std::abs(actual[i] - expected[i]) <= tolerance,
            std::format("{}[{}]: analytic={}, reference={}, error={}, tolerance={}",
                        context,
                        i,
                        actual[i],
                        expected[i],
                        std::abs(actual[i] - expected[i]),
                        tolerance));
  }
}

struct Fixture {
  std::string              name;
  QuantumIdentifier        id;
  lbl::band_data           band;
  LinemixingSpeciesEcsData rates;
};

Fixture oxygen() {
  Fixture out{"O2", QuantumIdentifier{"O2-66"_isot}, {}, {}};
  out.id.state[QuantumNumberType::S] = {.upper = Rational{1}, .lower = Rational{1}};
  out.band.lineshape                 = LineByLineLineshape::VP_ECS_MAKAROV;
  // Six microwave lines copied from arts-cat-data/lines/O2-66.xml, ground
  // electronic/vibrational band, 2026-09-24. Units: Hz, s^-1, J, Hz/Pa.
  // AIR widths are used as an explicit N2 test surrogate; this truncated
  // fixture tests derivatives, not the accuracy of atmospheric O2 absorption.
  // Catalogue Rosenkranz Y parameters are omitted because ECS supplies mixing.
  struct Entry {
    Numeric f, a, e, air, self;
    Index   n, jl;
  };
  constexpr std::array entries{Entry{58446595363.984505,
                                     7.680324949634841e-10,
                                     3.255308012861378e-22,
                                     15207.82861208981,
                                     15385.351903281518,
                                     3,
                                     4},
                               Entry{59164208870.49157,
                                     8.208551850869097e-10,
                                     1.5813499535005478e-21,
                                     13610.118991364421,
                                     13610.118991364421,
                                     7,
                                     6},
                               Entry{59590989816.19087,
                                     8.337540976469208e-10,
                                     8.387568987225636e-22,
                                     14231.450510535406,
                                     14201.863295336789,
                                     5,
                                     6},
                               Entry{60306062281.82392,
                                     8.588706898275419e-10,
                                     8.3828213816270505e-22,
                                     14231.450510535406,
                                     14201.863295336789,
                                     5,
                                     4},
                               Entry{60434784469.30783,
                                     8.771434358812755e-10,
                                     1.5805077004571166e-21,
                                     13610.118991364421,
                                     13610.118991364421,
                                     7,
                                     8},
                               Entry{62486259462.72251,
                                     9.145437479925664e-10,
                                     3.2285505871655825e-22,
                                     15207.82861208981,
                                     15385.351903281518,
                                     3,
                                     2}};
  for (const auto& entry : entries) {
    lbl::line line;
    line.f0                       = entry.f;
    line.a                        = entry.a;
    line.e0                       = entry.e;
    line.gu                       = Numeric(2 * entry.n + 1);
    line.gl                       = Numeric(2 * entry.jl + 1);
    line.ls.T0                    = 296;
    line.qn[QuantumNumberType::J] = {.upper = Rational{entry.n}, .lower = Rational{entry.jl}};
    line.qn[QuantumNumberType::N] = {.upper = Rational{entry.n}, .lower = Rational{entry.n}};
    line.ls.single_models[SpeciesEnum::Nitrogen].data[LineShapeModelVariable::G0] = {LineShapeModelType::T1,
                                                                                     {entry.air, 0.72}};
    line.ls.single_models[SpeciesEnum::Oxygen].data[LineShapeModelVariable::G0]   = {LineShapeModelType::T1,
                                                                                     {entry.self, 0.72}};
    out.band.lines.push_back(line);
  }
  // Same collision parameters as abs_ecs_dataAddMakarov2020.
  for (const auto partner : {SpeciesEnum::Nitrogen, SpeciesEnum::Oxygen}) {
    auto& rate                = out.rates[partner];
    rate.scaling              = {LineShapeModelType::T0, {1}};
    rate.beta                 = {LineShapeModelType::T0, {0.567}};
    rate.lambda               = {LineShapeModelType::T0, {0.39}};
    rate.collisional_distance = {LineShapeModelType::T0, {0.61e-10}};
  }
  return out;
}

Fixture carbon_dioxide() {
  Fixture out{"CO2", QuantumIdentifier{"CO2-626"_isot}, {}, {}};
  out.id.state[QuantumNumberType::l2] = {.upper = Rational{0}, .lower = Rational{0}};
  out.band.lineshape                  = LineByLineLineshape::VP_ECS_HARTMANN;
  // A deliberately synthetic coupled R-branch fixture. Nonzero, temperature
  // dependent shifts exercise real eigenvalue and Doppler-center derivatives.
  for (Index i = 0; i < 4; ++i) {
    const Index j = 2 * i;
    lbl::line   line;
    line.f0                       = 2e13 + Numeric(i) * 8e7;
    line.a                        = 1e-3 * Numeric(i + 1);
    line.e0                       = Conversion::kaycm2joule(0.39021 * Numeric(j * (j + 1)));
    line.gu                       = Numeric(2 * j + 3);
    line.gl                       = Numeric(2 * j + 1);
    line.ls.T0                    = 296;
    line.qn[QuantumNumberType::J] = {.upper = Rational{j + 1}, .lower = Rational{j}};
    for (const auto partner : {SpeciesEnum::Nitrogen, SpeciesEnum::Oxygen}) {
      const Numeric factor                   = partner == SpeciesEnum::Nitrogen ? 1 : 0.8;
      auto&         model                    = line.ls.single_models[partner];
      model.data[LineShapeModelVariable::G0] = {LineShapeModelType::T1,
                                                {factor * (1.5e4 + 500 * Numeric(i)), 0.7 + 0.03 * Numeric(i)}};
      model.data[LineShapeModelVariable::D0] = {LineShapeModelType::T3,
                                                {factor * (50 + 10 * Numeric(i)), factor * (0.2 + 0.03 * Numeric(i))}};
    }
    out.band.lines.push_back(line);
  }
  for (const auto partner : {SpeciesEnum::Nitrogen, SpeciesEnum::Oxygen}) {
    const bool nitrogen = partner == SpeciesEnum::Nitrogen;
    auto&      rate     = out.rates[partner];
    // Rodrigues1997 scaling/lambda, with synthetic temperature dependence of
    // beta and collision distance to cover every ECS coefficient derivative.
    rate.scaling = {LineShapeModelType::T1,
                    {Conversion::kaycm_per_atm2hz_per_pa(nitrogen ? 0.0180 : 0.0168), nitrogen ? 0.85 : 0.5}};
    rate.lambda  = {LineShapeModelType::T1, {nitrogen ? 0.81 : 0.82, nitrogen ? 0.0152 : -0.091}};
    rate.beta    = {LineShapeModelType::T1, {nitrogen ? 0.008 : 0.007, 0.17}};
    rate.collisional_distance = {LineShapeModelType::T1, {nitrogen ? 2.2e-10 : 2.4e-10, -0.12}};
  }
  return out;
}

AtmPoint atmosphere(const Fixture& fixture, Numeric temperature, Numeric pressure) {
  AtmPoint atm;
  atm.temperature                 = temperature;
  atm.pressure                    = pressure;
  atm[SpeciesEnum::Nitrogen]      = 0.79;
  atm[SpeciesEnum::Oxygen]        = 0.21;
  atm[SpeciesEnum::CarbonDioxide] = 4e-4;
  atm[fixture.id.isot]            = 1;
  return atm;
}

Vector frequencies(const Fixture& fixture, const AtmPoint& atm) {
  std::vector<Numeric> values;
  for (const auto& line : fixture.band) {
    const Numeric doppler =
        line.f0 * std::sqrt(Constant::doppler_broadening_const_squared * atm.temperature / fixture.id.isot.mass);
    const Numeric width = line.ls.G0(atm);
    for (Numeric offset : {-5 * width, -width, -0.7 * doppler, 0.0, 0.7 * doppler, width, 5 * width})
      values.push_back(line.f0 + offset);
  }
  for (Size i = 1; i < fixture.band.size(); ++i)
    values.push_back(0.5 * (fixture.band.lines[i - 1].f0 + fixture.band.lines[i].f0));
  std::sort(values.begin(), values.end());
  Vector grid(values.size());
  std::copy(values.begin(), values.end(), grid.begin());
  return grid;
}

struct Parameter {
  std::string                                       name;
  Numeric                                           step;
  std::function<void(Fixture&, AtmPoint&, Numeric)> perturb;
};
struct Request {
  Jacobian::Targets      targets;
  std::vector<Parameter> parameters;
};

Request requested_derivatives(const Fixture& fixture, const AtmPoint& atm) {
  Request request;
  auto    atmosphere_target = [&](const AtmKeyVal& key, std::string name, Numeric step, auto perturb) {
    request.targets.emplace_back(key);
    request.parameters.push_back({std::move(name), step, perturb});
  };
  auto line_target = [&](LineByLineVariable variable, std::string name, Numeric step, auto perturb) {
    LblLineKey key;
    key.band = fixture.id;
    key.line = 1;
    key.var  = variable;
    request.targets.emplace_back(key);
    request.parameters.push_back({std::move(name), step, perturb});
  };
  atmosphere_target(
      AtmKey::t, "temperature", 0.01, [](Fixture&, AtmPoint& state, Numeric delta) { state.temperature += delta; });
  line_target(LineByLineVariable::f0, "line f0", 100, [](Fixture& data, AtmPoint&, Numeric delta) {
    data.band.lines[1].f0 += delta;
  });
  atmosphere_target(AtmKey::p, "pressure", 1e-4 * atm.pressure, [](Fixture&, AtmPoint& state, Numeric delta) {
    state.pressure += delta;
  });
  line_target(LineByLineVariable::e0,
              "line e0",
              1e-4 * Constant::k * atm.temperature,
              [](Fixture& data, AtmPoint&, Numeric delta) { data.band.lines[1].e0 += delta; });
  for (const auto species : {SpeciesEnum::Nitrogen, SpeciesEnum::Oxygen, SpeciesEnum::CarbonDioxide}) {
    atmosphere_target(species,
                      std::format("VMR {}", species),
                      species == SpeciesEnum::CarbonDioxide ? 1e-5 : (atm.pressure < 1e3 ? 7.5e-4 : 1e-4),
                      [species](Fixture&, AtmPoint& state, Numeric delta) { state[species] += delta; });
  }
  line_target(
      LineByLineVariable::a, "Einstein A", fixture.band.lines[1].a * 1e-4, [](Fixture& data, AtmPoint&, Numeric delta) {
        data.band.lines[1].a += delta;
      });
  atmosphere_target(fixture.id.isot,
                    "isotopologue ratio",
                    atm.pressure < 1e3 ? 1e-4 : 2e-4,
                    [isot = fixture.id.isot](Fixture&, AtmPoint& state, Numeric delta) { state[isot] += delta; });
  for (const auto species : {SpeciesEnum::Nitrogen, SpeciesEnum::Oxygen}) {
    for (const auto variable : {LineShapeModelVariable::G0, LineShapeModelVariable::D0}) {
      if (not fixture.band.lines[1].ls.single_models.contains(species)) continue;
      const auto& model = fixture.band.lines[1].ls.single_models.at(species).data;
      if (not model.contains(variable)) continue;
      for (const auto coefficient : {LineShapeModelCoefficient::X0, LineShapeModelCoefficient::X1}) {
        LblLineKey key;
        key.band     = fixture.id;
        key.line     = 1;
        key.spec     = species;
        key.ls_var   = variable;
        key.ls_coeff = coefficient;
        request.targets.emplace_back(key);
        // Keep the relative width perturbation resolved even near T0, where
        // differentiating its exponent otherwise changes almost no digits.
        const Numeric width_step = atm.pressure < 1e3 ? 3e-4 : 3e-5;
        Numeric       step       = coefficient == LineShapeModelCoefficient::X0
                                       ? width_step * model.at(variable).X(coefficient)
                                       : width_step / std::abs(std::log(296 / atm.temperature));
        // Pressure shifts are added to a large optical carrier. Resolve their
        // finite differences without materially changing the spectral shape.
        if (variable == LineShapeModelVariable::D0) {
          const Numeric shift_step = atm.pressure < 1e3 ? 20 : 1;
          step = coefficient == LineShapeModelCoefficient::X0 ? shift_step
                                                              : shift_step / (1 + std::abs(atm.temperature - 296));
        }
        request.parameters.push_back({std::format("{} {} {}", species, variable, coefficient),
                                      step,
                                      [species, variable, coefficient](Fixture& data, AtmPoint&, Numeric delta) {
                                        data.band.lines[1].ls.single_models.at(species).data.at(variable).X(
                                            coefficient) += delta;
                                      }});
      }
    }
  }
  return request;
}

ComputeData evaluate(const Fixture&           fixture,
                     const AtmPoint&          atm,
                     const Vector&            grid,
                     const Jacobian::Targets* targets = nullptr,
                     const ArrayOfIndex*      order   = nullptr) {
  ComputeData data(grid, atm);
  if (order) data.sort = *order;
  if (targets)
    data.adapt_single(fixture.id, fixture.band, fixture.rates, atm, *targets, order != nullptr);
  else
    data.adapt_single(fixture.id, fixture.band, fixture.rates, atm, order != nullptr);
  data.core_calc(grid);
  return data;
}

std::vector<Size> match_modes(const ComputeData& reference, const ComputeData& side) {
  const Size        n = reference.pop.size();
  std::vector<Size> matched(n);
  std::vector<bool> used(n, false);
  for (Size i = 0; i < n; ++i) {
    Numeric distance = std::numeric_limits<Numeric>::infinity();
    Size    closest  = n;
    for (Size j = 0; j < n; ++j) {
      const Numeric candidate = std::abs(reference.eqv_vals[0, i] - side.eqv_vals[0, j]);
      if (not used[j] and candidate < distance) {
        distance = candidate;
        closest  = j;
      }
    }
    require(closest != n, "Could not match perturbed ECS modes");
    matched[i]    = closest;
    used[closest] = true;
  }
  return matched;
}

Numeric max_abs(const Samples& values) {
  Numeric result = 0;
  for (const auto value : values) result = std::max(result, std::abs(value));
  return result;
}

Numeric     largest_profile_error = 0;
std::string largest_profile_context;

void compare_profile(const Samples&     actual,
                     const Samples&     expected,
                     const std::string& context,
                     Numeric            reference_roundoff = 0) {
  // ECS uses the existing LTE dF facility, which is a complex forward
  // difference. Its approximation error also contributes near derivative zeros.
  // The O2 T derivative's peak-normalized error is 3.19e-4 for both
  // 0.01 K and 0.05 K central steps, distinguishing this bias from FD error.
  Numeric error = 0;
  for (Size i = 0; i < actual.size(); ++i) error = std::max(error, std::abs(actual[i] - expected[i]));
  const Numeric normalized = max_abs(expected) == 0 ? error : error / max_abs(expected);
  if (normalized > largest_profile_error) {
    largest_profile_error   = normalized;
    largest_profile_context = context;
  }
  if (context.starts_with("O2 T=220.25 P=10 temperature"))
    std::cout << context << " peak-normalized error=" << normalized << '\n';
  compare(actual, expected, 1.5e-3, 1.2e-4 * max_abs(expected) + reference_roundoff + 1e-60, context);
}

struct Absorption {
  Samples              value;
  std::vector<Samples> derivative;
};

Absorption calculate(const Fixture&           fixture,
                     const AtmPoint&          atm,
                     const Vector&            grid,
                     const Jacobian::Targets* targets = nullptr) {
  ComputeData             data(grid, atm);
  const Jacobian::Targets none;
  const Size              count = targets ? targets->target_count() : 0;
  PropmatVector           pm(grid.size());
  PropmatMatrix           dpm(count, grid.size());
  pm  = 0;
  dpm = 0;
  lbl::voigt::ecs::calculate(pm,
                             dpm,
                             data,
                             grid,
                             Range(0, grid.size()),
                             targets ? *targets : none,
                             fixture.id,
                             fixture.band,
                             fixture.rates,
                             atm,
                             ZeemanPolarization::no,
                             false);
  Absorption result;
  result.derivative.resize(count);
  for (Size i = 0; i < grid.size(); ++i) {
    result.value.push_back(pm[i].A());
    for (Size target = 0; target < count; ++target) result.derivative[target].push_back(dpm[target, i].A());
  }
  return result;
}

void derivative_case(const Fixture& fixture, Numeric temperature, Numeric pressure) {
  const auto atm        = atmosphere(fixture, temperature, pressure);
  const auto grid       = frequencies(fixture, atm);
  const auto request    = requested_derivatives(fixture, atm);
  const auto base       = evaluate(fixture, atm, grid, &request.targets);
  const auto plain      = evaluate(fixture, atm, grid, nullptr, &base.sort);
  const auto absorption = calculate(fixture, atm, grid, &request.targets);
  const auto context    = std::format("{} T={} P={}", fixture.name, temperature, pressure);
  const Size n = base.pop.size(), count = request.parameters.size();
  require(base.dW.npages() == static_cast<Index>(count) and base.dW.nrows() == static_cast<Index>(n) and
              base.dW.ncols() == static_cast<Index>(n),
          context + " incorrect batched matrix dimensions");
  Samples active_shape, plain_shape;
  for (Size i = 0; i < grid.size(); ++i) {
    active_shape.push_back(base.shape[i]);
    plain_shape.push_back(plain.shape[i]);
  }
  compare(active_shape, plain_shape, 2e-13, 0, context + " primal spectrum");
  Numeric offdiagonal = 0, damping_scale = 0;
  for (Size i = 0; i < n; ++i) {
    for (Size j = 0; j < n; ++j) {
      require(base.Ws[0, i, j] == plain.Ws[0, i, j], context + " derivative request changed matrix");
      damping_scale = std::max(damping_scale, std::abs(base.Ws[0, i, j].imag()));
      if (i != j) offdiagonal = std::max(offdiagonal, std::abs(base.Ws[0, i, j]));
    }
  }
  require(offdiagonal > 1, context + " fixture must exercise nonzero line coupling");
  // Targets deliberately interleave atmospheric and line parameters, exercising
  // global target_pos indexing rather than the storage order of each category.
  for (Size target = 0; target < count; ++target) {
    const auto& parameter = request.parameters[target];
    for (Numeric multiplier : {1.0, parameter.name == "temperature" ? 5.0 : 2.0}) {
      const Numeric step     = multiplier * parameter.step;
      auto          high_atm = atm, low_atm = atm;
      auto          high_fixture = fixture, low_fixture = fixture;
      parameter.perturb(high_fixture, high_atm, step);
      parameter.perturb(low_fixture, low_atm, -step);
      const auto high  = evaluate(high_fixture, high_atm, grid, nullptr, &base.sort);
      const auto low   = evaluate(low_fixture, low_atm, grid, nullptr, &base.sort);
      const auto label = context + std::format(" {} h={}", parameter.name, step);
      Samples    pop, pop_fd, dip, dip_fd;
      for (Size i = 0; i < n; ++i) {
        pop.push_back(base.dpop[target, i]);
        pop_fd.push_back((high.pop[i] - low.pop[i]) / (2 * step));
        dip.push_back(base.ddip[target, i]);
        dip_fd.push_back((high.dip[i] - low.dip[i]) / (2 * step));
      }
      compare(pop, pop_fd, 2e-6, 1e-60, label + " population");
      Numeric dipole_scale = 0;
      for (Size i = 0; i < n; ++i) dipole_scale = std::max(dipole_scale, std::abs(base.dip[i]));
      compare(dip, dip_fd, 2e-6, 8 * std::numeric_limits<Numeric>::epsilon() * dipole_scale / step, label + " dipole");
      Samples matrix_real, matrix_real_fd, matrix_imag, matrix_imag_fd;
      for (Size i = 0; i < n; ++i) {
        for (Size j = 0; j < n; ++j) {
          const Complex difference = (high.Ws[0, i, j] - low.Ws[0, i, j]) / (2 * step);
          matrix_real.push_back(base.dW[target, i, j].real());
          matrix_real_fd.push_back(difference.real());
          matrix_imag.push_back(base.dW[target, i, j].imag());
          matrix_imag_fd.push_back(difference.imag());
        }
      }
      const Numeric eps              = std::numeric_limits<Numeric>::epsilon();
      const Numeric carrier_roundoff = 8 * eps * fixture.band.front().f0 / step;
      compare(matrix_real, matrix_real_fd, 3e-6, carrier_roundoff, label + " matrix real");
      compare(matrix_imag, matrix_imag_fd, 3e-6, 8 * eps * damping_scale / step, label + " matrix imaginary");
      const auto hi = match_modes(base, high), lo = match_modes(base, low);
      Samples    values, values_fd, strengths, strengths_fd, shape, shape_fd;
      for (Size i = 0; i < n; ++i) {
        values.push_back(base.deqv_vals[target, i]);
        values_fd.push_back((high.eqv_vals[0, hi[i]] - low.eqv_vals[0, lo[i]]) / (2 * step));
        strengths.push_back(base.deqv_strs[target, i]);
        strengths_fd.push_back((high.eqv_strs[0, hi[i]] - low.eqv_strs[0, lo[i]]) / (2 * step));
      }
      compare(values, values_fd, 3e-5, carrier_roundoff, label + " equivalent positions");
      Numeric strength_scale = 0;
      for (Size i = 0; i < n; ++i) strength_scale = std::max(strength_scale, std::abs(base.eqv_strs[0, i]));
      // A central difference of very small residue changes subtracts primal
      // strengths. Allow their matrix-size-scaled floating-point roundoff,
      // while keeping the relative derivative tolerance unchanged.
      const Numeric strength_roundoff = 8 * Numeric(n) * eps * strength_scale / step;
      compare(strengths, strengths_fd, 3e-5, strength_roundoff, label + " equivalent strengths");
      for (Size i = 0; i < grid.size(); ++i) {
        shape.push_back(base.dshape[target, i]);
        shape_fd.push_back((high.shape[i] - low.shape[i]) / (2 * step));
      }
      // At 10 Pa the isotope-dependent O2 collision shifts are below a
      // 60 GHz carrier ULP: h=1e-4 and h=1e-3 yield nonconvergent profile FDs.
      // Check this weak operator direction at higher pressure; its matrix,
      // residues and full absorption-amplitude derivative remain tested here.
      if (not(fixture.id.isot.spec == SpeciesEnum::Oxygen and pressure < 1e3 and
              parameter.name == "isotopologue ratio"))
        compare_profile(shape, shape_fd, label + " complex Voigt spectrum");
      compare(
          {base.dgd_fac[target]}, {(high.gd_fac - low.gd_fac) / (2 * step)}, 2e-7, 1e-30, label + " Doppler factor");
      const auto high_absorption = calculate(high_fixture, high_atm, grid);
      const auto low_absorption  = calculate(low_fixture, low_atm, grid);
      Samples    absorption_fd;
      for (Size i = 0; i < grid.size(); ++i)
        absorption_fd.push_back((high_absorption.value[i] - low_absorption.value[i]) / (2 * step));
      // At low pressure the CO2 VMR-induced center perturbations approach the
      // 20 THz carrier ULP. Their central differences oscillate with h. Bound
      // that reference noise by one carrier ULP / h times peak absorption /
      // Doppler width; the carrier cancels against the Doppler factor below.
      const Numeric absorption_roundoff = fixture.id.isot.spec == SpeciesEnum::CarbonDioxide and pressure < 1e3 and
                                                  (parameter.name == "VMR N2" or parameter.name == "VMR O2")
                                              ? eps * max_abs(absorption.value) / (step * base.gd_fac)
                                              : 0;
      compare_profile(absorption.derivative[target], absorption_fd, label + " calculate Jacobian", absorption_roundoff);
    }
  }
}
void frequency_and_zero_targets(const Fixture& fixture) {
  const auto        atm  = atmosphere(fixture, 280.25, 1e4);
  const auto        grid = frequencies(fixture, atm);
  Jacobian::Targets targets;
  for (const auto key : {AtmKey::wind_u, AtmKey::wind_v, AtmKey::wind_w, AtmKey::mag_u}) targets.emplace_back(key);
  const auto base       = evaluate(fixture, atm, grid, &targets);
  const auto absorption = calculate(fixture, atm, grid, &targets);
  // Other LBL profiles also return d/df for each wind target; propagation
  // subsequently applies its line-of-sight projection and Doppler conversion.
  for (Numeric step : {100.0, 500.0}) {
    Vector high_grid = grid, low_grid = grid;
    high_grid                  += step;
    low_grid                   -= step;
    const auto high             = evaluate(fixture, atm, high_grid);
    const auto low              = evaluate(fixture, atm, low_grid);
    const auto high_absorption  = calculate(fixture, atm, high_grid);
    const auto low_absorption   = calculate(fixture, atm, low_grid);
    Samples    expected_shape, expected_absorption;
    for (Size i = 0; i < grid.size(); ++i) {
      expected_shape.push_back((high.shape[i] - low.shape[i]) / (2 * step));
      expected_absorption.push_back((high_absorption.value[i] - low_absorption.value[i]) / (2 * step));
    }
    for (Size target = 0; target < 3; ++target) {
      Samples analytic;
      for (Size i = 0; i < grid.size(); ++i) analytic.push_back(base.dshape[target, i]);
      compare_profile(analytic, expected_shape, fixture.name + " wind frequency spectrum");
      compare_profile(absorption.derivative[target], expected_absorption, fixture.name + " wind frequency absorption");
    }
  }
  // With no Zeeman polarization, the magnetic target is exactly inactive.
  for (Size i = 0; i < grid.size(); ++i)
    require(base.dshape[3, i] == Complex{} and absorption.derivative[3][i] == Complex{},
            fixture.name + " inactive magnetic target must be zero");
}

void degenerate_fixed_operator() {
  AtmPoint atm;
  atm.temperature = 296.25;
  atm.pressure    = 1e5;
  const Vector grid{0.999e9, 1e9, 1.001e9};
  ComputeData  data(grid, atm);
  data.pop    = Vector{2.0, 0.5};
  data.dip    = Vector{1.0, -0.8};
  data.vmrs   = Vector{1.0};
  data.gd_fac = 1e-4;
  data.Ws.resize(1, 2, 2);
  data.Ws          = 0;
  data.Ws[0, 0, 0] = data.Ws[0, 1, 1] = Complex(1e9, 1e5);
  data.dW.resize(3, 2, 2);
  data.dpop.resize(3, 2);
  data.ddip.resize(3, 2);
  data.dgd_fac.resize(3);
  data.df.resize(3);
  data.dW         = 0;
  data.dpop       = 0;
  data.ddip       = 0;
  data.dgd_fac    = 0;
  data.df         = 0;
  data.dpop[1, 0] = 0.2;
  data.dpop[1, 1] = -0.1;
  data.ddip[2, 0] = 0.03;
  data.ddip[2, 1] = 0.07;
  // These directions leave the degenerate operator fixed. Their derivatives
  // are well-defined and must not demand derivatives of its eigenvectors.
  data.core_calc(grid);
  constexpr Numeric                strength = 2.0 + 0.5 * 0.8 * 0.8;
  constexpr std::array<Numeric, 3> derivatives{
      0.0, 0.2 - 0.1 * 0.8 * 0.8, 2 * 2.0 * 1.0 * 0.03 + 2 * 0.5 * (-0.8) * 0.07};
  for (Size target = 0; target < derivatives.size(); ++target) {
    Complex total = 0;
    for (Size i = 0; i < 2; ++i) {
      require(data.deqv_vals[target, i] == Complex{}, "Fixed degenerate eigenvalue derivative must be zero");
      total += data.deqv_strs[target, i];
    }
    compare({total}, {derivatives[target]}, 2e-13, 0, "Degenerate fixed-operator strength derivative");
    for (Size i = 0; i < grid.size(); ++i)
      compare({data.dshape[target, i]},
              {data.shape[i] * derivatives[target] / strength},
              2e-13,
              0,
              "Degenerate fixed-operator spectrum derivative");
  }
}

void bath_case() {
  auto fixture = oxygen();
  fixture.name = "O2 Bath";
  // The AIR-width/N2-rate surrogate is explicitly assigned to Bath here.
  // N2/CO2 VMRs then affect Bath's mean mass without changing its fraction;
  // O2 VMR changes both its explicit fraction and the complementary Bath.
  for (auto& line : fixture.band.lines) {
    line.ls.single_models[SpeciesEnum::Bath] = line.ls.single_models.at(SpeciesEnum::Nitrogen);
    line.ls.single_models.erase(SpeciesEnum::Nitrogen);
  }
  fixture.rates[SpeciesEnum::Bath] = fixture.rates.at(SpeciesEnum::Nitrogen);
  fixture.rates.erase(SpeciesEnum::Nitrogen);
  derivative_case(fixture, 280.25, 3e4);
}

}  // namespace

int main() try {
  WignerInformation wigner(100, 0, true, true);
  for (const auto& fixture : {oxygen(), carbon_dioxide()}) {
    // Keep central differences inside one partition-function interpolation cell.
    for (Numeric temperature : {220.25, 296.25, 330.25}) {
      for (Numeric pressure : {10.0, 1e4, 1e5}) derivative_case(fixture, temperature, pressure);
    }
    frequency_and_zero_targets(fixture);
  }
  bath_case();
  degenerate_fixed_operator();
  wigner.unload();
  std::cout << "Largest peak-normalized profile error=" << largest_profile_error << " (" << largest_profile_context
            << ")\n";
  std::cout << "ECS batched derivative tests passed\n";
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
