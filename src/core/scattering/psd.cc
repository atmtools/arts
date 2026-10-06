#include "psd.h"

#include <arts_constants.h>

#include <algorithm>
#include <cmath>

namespace scattering {

MonodispersePSD::MonodispersePSD(ScatteringSpeciesProperty number_density_, Numeric t_min_, Numeric t_max_)
    : number_density(std::move(number_density_)), t_min(t_min_), t_max(t_max_) {}

Vector MonodispersePSD::evaluate(const AtmPoint& point, const Vector& particle_sizes, Numeric, Numeric) const {
  ARTS_USER_ERROR_IF(particle_sizes.size() != 1, "MonodispersePSD requires exactly one particle.")
  ARTS_USER_ERROR_IF(not point.has(number_density), "The PSD requires {}.", number_density)
  return Vector{point.temperature < t_min || point.temperature > t_max ? 0.0 : point[number_density]};
}

PSDData MonodispersePSD::evaluate_with_derivatives(const AtmPoint& point,
                                                   const Vector&   sizes,
                                                   Numeric         a,
                                                   Numeric         b) const {
  const Numeric derivative = point.temperature < t_min || point.temperature > t_max ? 0.0 : 1.0;
  return {.values = evaluate(point, sizes, a, b), .derivatives = {{number_density, Vector{derivative}}}};
}
namespace {
bool temperature_allowed(Numeric temperature, Numeric t_min, Numeric t_max, bool picky) {
  if (temperature >= t_min and temperature <= t_max) return true;
  if (picky) { ARTS_USER_ERROR("PSD temperature {} K is outside [{}, {}] K.", temperature, t_min, t_max) }
  return false;
}

template <typename Evaluator> PSDData numerical_derivatives(const AtmPoint&                               point,
                                                            const Vector&                                 sizes,
                                                            Numeric                                       a,
                                                            Numeric                                       b,
                                                            const std::vector<ScatteringSpeciesProperty>& properties,
                                                            const Evaluator&                              evaluate) {
  PSDData result{.values = evaluate(point, sizes, a, b), .derivatives = {}};
  for (const auto& property : properties) {
    if (property.species_name.empty() or not point.has(property)) continue;

    const Numeric value  = point[property];
    const Numeric step   = 1e-5 * std::max(std::abs(value), 1e-8);
    AtmPoint      plus   = point;
    plus[property]      += step;
    Vector plus_value    = evaluate(plus, sizes, a, b);
    Vector derivative(sizes.size());

    if (value > step) {
      AtmPoint minus            = point;
      minus[property]          -= step;
      const Vector minus_value  = evaluate(minus, sizes, a, b);
      for (Size i = 0; i < sizes.size(); ++i) derivative[i] = (plus_value[i] - minus_value[i]) / (2.0 * step);
    } else {
      for (Size i = 0; i < sizes.size(); ++i) derivative[i] = (plus_value[i] - result.values[i]) / step;
    }
    result.derivatives[property] = std::move(derivative);
  }
  return result;
}

Vector modified_gamma(const Vector& sizes, Numeric n0, Numeric mu, Numeric lambda, Numeric gamma) {
  Vector result(sizes.size());
  Matrix unused(4, sizes.size());
  mgd_with_derivatives(result, unused, sizes, n0, mu, lambda, gamma, false, false, false, false);
  return result;
}

Numeric delanoe_n0_from_temperature(Numeric t) { return std::exp(-0.076586 * (t - 273.15) + 17.948); }

Numeric delanoe_dm_from_mass_n0(Numeric mass, Numeric n0, Numeric rho) {
  return mass == 0.0 ? 1e-9 : std::pow(256.0 * mass / (Constant::pi * rho * n0), 0.25);
}

Numeric delanoe_n0_from_mass_dm(Numeric mass, Numeric dm, Numeric rho) {
  return dm > 1e-9 ? 256.0 * mass / (Constant::pi * rho * std::pow(dm, 4)) : 0.0;
}
}  // namespace

MGDSingleMoment::MGDSingleMoment(ScatteringSpeciesProperty moment_,
                                 Numeric                   n_alpha_,
                                 Numeric                   n_b_,
                                 Numeric                   mu_,
                                 Numeric                   gamma_,
                                 Numeric                   t_min_,
                                 Numeric                   t_max_,
                                 bool                      picky_)
    : moment(std::move(moment_)),
      n_alpha(n_alpha_),
      n_b(n_b_),
      mu(mu_),
      gamma(gamma_),
      t_min(t_min_),
      t_max(t_max_),
      picky(picky_) {}

Vector MGDSingleMoment::evaluate(const AtmPoint& point,
                                 const Vector&   particle_sizes,
                                 const Numeric&  scat_species_a,
                                 const Numeric&  scat_species_b) const {
  if (!point.has(moment)) {
    std::ostringstream os;
    os << "The PSD requires water content from '" << moment << "' but it is not part of the AtmPoint.";
    throw std::runtime_error(os.str());
  }

  Numeric water_content = point[moment];
  Numeric t             = point[AtmKey::t];

  // Outside of [t_min,tmax]?
  if (not temperature_allowed(t, t_min, t_max, picky)) return Vector(particle_sizes.size(), 0.0);

  // Negative wc?
  const Numeric sign = water_content < 0 ? -1.0 : 1.0;
  water_content      = std::abs(water_content);

  auto   nsi = particle_sizes.size();
  Vector psd(nsi);
  psd = 0.0;
  if (water_content == 0.0) return psd;

  // Calculate PSD
  // Calculate lambda for modified gamma distribution from mass density
  Numeric k     = (scat_species_b + mu + 1 - gamma) / gamma;
  Numeric expo  = 1.0 / (n_b - k - 1);
  Numeric denom = scat_species_a * n_alpha * tgamma(k + 1);
  Numeric lam   = pow(water_content * gamma / denom, expo);
  Numeric n_0   = n_alpha * pow(lam, n_b);
  Matrix  jac_data(4, nsi);

  mgd_with_derivatives(psd,
                       jac_data,
                       particle_sizes,
                       n_0,
                       mu,
                       lam,
                       gamma,
                       false,   // n_0 jacobian
                       false,   // mu jacobian
                       false,   // lambda jacobian
                       false);  // gamma jacobian
  psd *= sign;
  return psd;
}

PSDData MGDSingleMoment::evaluate_with_derivatives(const AtmPoint& point,
                                                   const Vector&   sizes,
                                                   Numeric         a,
                                                   Numeric         b) const {
  return numerical_derivatives(point, sizes, a, b, {moment}, [this](const auto& p, const auto& s, auto aa, auto bb) {
    return evaluate(p, s, aa, bb);
  });
}

MGDSingleMoment::MGDSingleMoment(
    ScatteringSpeciesProperty moment_, std::string name, Numeric t_min_, Numeric t_max_, bool picky_)
    : moment(std::move(moment_)), t_min(t_min_), t_max(t_max_), picky(picky_) {
  if (name == "Abel12") {
    n_alpha = 0.22;
    n_b     = 2.2;
    mu      = 0.0;
    gamma   = 1.0;
  } else if (name == "Wang16") {
    // Wang 16 parameters converted to SI units
    n_alpha = 14.764;
    n_b     = 1.49;
    mu      = 0.0;
    gamma   = 1.0;
  } else if (name == "Field19") {
    n_alpha = 7.9e9;
    n_b     = -2.58;
    mu      = 0.0;
    gamma   = 1.0;
  } else {
    std::ostringstream os;
    os << "The PSD configuration '" << name << "' is currently not supported."
       << " Supported config names are 'Abel12', 'Wang16', 'Field19'.";
    throw std::runtime_error(os.str());
  }
}

MGDMass::MGDMass(ScatteringSpeciesProperty mass_,
                 Numeric                   n0_,
                 Numeric                   mu_,
                 Numeric                   lambda_,
                 Numeric                   gamma_,
                 Numeric                   t_min_,
                 Numeric                   t_max_,
                 bool                      picky_)
    : mass(std::move(mass_)),
      n0(n0_),
      mu(mu_),
      lambda(lambda_),
      gamma(gamma_),
      t_min(t_min_),
      t_max(t_max_),
      picky(picky_) {
  if (std::isnan(n0) == std::isnan(lambda)) { ARTS_USER_ERROR("Exactly one of n0 and lambda must be NaN.") }
}

Vector MGDMass::evaluate(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  ARTS_USER_ERROR_IF(not point.has(mass), "The PSD requires {}.", mass)
  ARTS_USER_ERROR_IF(a <= 0.0 or b <= 0.0 or gamma <= 0.0, "Invalid modified-gamma or mass-size parameters.")
  if (not temperature_allowed(point.temperature, t_min, t_max, picky)) return Vector(sizes.size(), 0.0);

  const Numeric content = std::abs(point[mass]);
  if (content == 0.0) return Vector(sizes.size(), 0.0);
  const Numeric sign  = point[mass] < 0.0 ? -1.0 : 1.0;
  const Numeric order = (mu + b + 1.0) / gamma;
  ARTS_USER_ERROR_IF(order <= 0.0, "Modified-gamma mass moment requires mu + b + 1 > 0.")

  Numeric local_n0     = n0;
  Numeric local_lambda = lambda;
  if (std::isnan(local_n0)) {
    local_n0 = gamma * std::pow(local_lambda, order) * content / (a * std::tgamma(order));
  } else {
    local_lambda =
        std::pow(gamma / (a * local_n0 * std::tgamma(order)), -1.0 / order) * std::pow(content, -1.0 / order);
  }
  Vector result  = modified_gamma(sizes, local_n0, mu, local_lambda, gamma);
  result        *= sign;
  return result;
}

PSDData MGDMass::evaluate_with_derivatives(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  return numerical_derivatives(point, sizes, a, b, {mass}, [this](const auto& p, const auto& s, auto aa, auto bb) {
    return evaluate(p, s, aa, bb);
  });
}

MGDTwoMoment::MGDTwoMoment(ScatteringSpeciesProperty mass_,
                           ScatteringSpeciesProperty second_moment_,
                           MGDTwoMomentType          type_,
                           Numeric                   mu_,
                           Numeric                   gamma_,
                           Numeric                   t_min_,
                           Numeric                   t_max_,
                           bool                      picky_)
    : mass(std::move(mass_)),
      second_moment(std::move(second_moment_)),
      type(type_),
      mu(mu_),
      gamma(gamma_),
      t_min(t_min_),
      t_max(t_max_),
      picky(picky_) {}

Vector MGDTwoMoment::evaluate(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  ARTS_USER_ERROR_IF(
      not point.has(mass) or not point.has(second_moment), "The PSD requires both {} and {}.", mass, second_moment)
  ARTS_USER_ERROR_IF(a <= 0.0 or b <= 0.0 or gamma <= 0.0, "Invalid modified-gamma or mass-size parameters.")
  if (not temperature_allowed(point.temperature, t_min, t_max, picky)) return Vector(sizes.size(), 0.0);

  const Numeric content = std::abs(point[mass]);
  const Numeric second  = point[second_moment];
  ARTS_USER_ERROR_IF(second <= 0.0, "The second PSD moment {} must be positive.", second_moment)
  if (content == 0.0) return Vector(sizes.size(), 0.0);

  const Numeric order = (mu + b + 1.0) / gamma;
  ARTS_USER_ERROR_IF(order <= 0.0, "Modified-gamma mass moment requires mu + b + 1 > 0.")
  const Numeric g_order = std::tgamma(order);
  Numeric       lambda;
  if (type == MGDTwoMomentType::MassMeanSize) {
    lambda = std::pow(std::tgamma(order + 1.0 / gamma) / g_order, gamma) * std::pow(second, -gamma);
  } else {
    ARTS_USER_ERROR_IF(mu + 1.0 <= 0.0, "Number-constrained PSD requires mu + 1 > 0.")
    const Numeric exponent = gamma / b;
    lambda = std::pow(a * g_order / std::tgamma((mu + 1.0) / gamma), exponent) * std::pow(second / content, exponent);
  }
  const Numeric local_n0 = gamma * std::pow(lambda, order) * content / (a * g_order);
  Vector        result   = modified_gamma(sizes, local_n0, mu, lambda, gamma);
  if (point[mass] < 0.0) result *= -1.0;
  return result;
}

PSDData MGDTwoMoment::evaluate_with_derivatives(const AtmPoint& point,
                                                const Vector&   sizes,
                                                Numeric         a,
                                                Numeric         b) const {
  return numerical_derivatives(
      point, sizes, a, b, {mass, second_moment}, [this](const auto& p, const auto& s, auto aa, auto bb) {
        return evaluate(p, s, aa, bb);
      });
}

DelanoeEtAl14::DelanoeEtAl14(ScatteringSpeciesProperty mass_,
                             Numeric                   rho_,
                             Numeric                   alpha_,
                             Numeric                   beta_,
                             Numeric                   t_min_,
                             Numeric                   t_max_,
                             Numeric                   dm_min_,
                             bool                      picky_)
    : DelanoeEtAl14(std::move(mass_), {}, {}, rho_, alpha_, beta_, t_min_, t_max_, dm_min_, picky_) {}

DelanoeEtAl14::DelanoeEtAl14(ScatteringSpeciesProperty mass_,
                             ScatteringSpeciesProperty intercept_parameter_,
                             ScatteringSpeciesProperty mean_size_,
                             Numeric                   rho_,
                             Numeric                   alpha_,
                             Numeric                   beta_,
                             Numeric                   t_min_,
                             Numeric                   t_max_,
                             Numeric                   dm_min_,
                             bool                      picky_)
    : mass(std::move(mass_)),
      intercept_parameter(std::move(intercept_parameter_)),
      mean_size(std::move(mean_size_)),
      rho(rho_),
      alpha(alpha_),
      beta(beta_),
      t_min(t_min_),
      t_max(t_max_),
      dm_min(dm_min_),
      picky(picky_) {}

Vector DelanoeEtAl14::evaluate(const AtmPoint& point, const Vector& sizes, Numeric, Numeric) const {
  ARTS_USER_ERROR_IF(not point.has(mass), "The PSD requires {}.", mass)
  if (not temperature_allowed(point.temperature, t_min, t_max, picky)) return Vector(sizes.size(), 0.0);
  const Numeric iwc = point[mass];
  if (iwc == 0.0) return Vector(sizes.size(), 0.0);
  ARTS_USER_ERROR_IF(iwc < 0.0, "Delanoe14 requires non-negative ice mass density.")

  const bool has_n0 = not intercept_parameter.species_name.empty() and point.has(intercept_parameter);
  const bool has_dm = not mean_size.species_name.empty() and point.has(mean_size);
  Numeric    n0     = has_n0 ? point[intercept_parameter] : delanoe_n0_from_temperature(point.temperature);
  Numeric    dm     = has_dm ? point[mean_size] : delanoe_dm_from_mass_n0(iwc, n0, rho);
  if (has_dm and not has_n0) n0 = delanoe_n0_from_mass_dm(iwc, dm, rho);
  ARTS_USER_ERROR_IF(dm <= 0.0 or dm < dm_min, "Delanoe14 Dm={} m is below the allowed minimum {} m.", dm, dm_min)

  Vector scaled  = sizes;
  scaled        *= 1.0 / dm;
  if (not scaled.empty() and scaled[0] <= std::numeric_limits<Numeric>::epsilon()) {
    ARTS_USER_ERROR_IF(scaled.size() < 2, "Delanoe14 cannot evaluate a single zero-sized particle.")
    scaled[0] = 0.1 * scaled[1];
  }
  Vector result(sizes.size());
  Matrix derivative(1, sizes.size());
  delanoe_shape_with_derivative(result, derivative, scaled, alpha, beta);
  result *= n0;
  return result;
}

PSDData DelanoeEtAl14::evaluate_with_derivatives(const AtmPoint& point,
                                                 const Vector&   sizes,
                                                 Numeric         a,
                                                 Numeric         b) const {
  return numerical_derivatives(
      point,
      sizes,
      a,
      b,
      {mass, intercept_parameter, mean_size},
      [this](const auto& p, const auto& s, auto aa, auto bb) { return evaluate(p, s, aa, bb); });
}

FieldEtAl07::FieldEtAl07(ScatteringSpeciesProperty mass_,
                         std::string               regime,
                         Numeric                   t_min_,
                         Numeric                   t_max_,
                         Numeric                   t_min_psd_,
                         Numeric                   t_max_psd_,
                         bool                      picky_)
    : mass(std::move(mass_)),
      tropical(regime == "TR"),
      t_min(t_min_),
      t_max(t_max_),
      t_min_psd(t_min_psd_),
      t_max_psd(t_max_psd_),
      picky(picky_) {
  if (regime != "TR" and regime != "ML") { ARTS_USER_ERROR("Field07 regime must be 'TR' or 'ML'.") }
}

Vector FieldEtAl07::evaluate(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  ARTS_USER_ERROR_IF(not point.has(mass), "The PSD requires {}.", mass)
  ARTS_USER_ERROR_IF(a <= 0.0 or b < 1.01 or b > 4.0, "Field07 received an invalid mass-size relationship.")
  if (not temperature_allowed(point.temperature, t_min, t_max, picky)) return Vector(sizes.size(), 0.0);
  Numeric content = point[mass];
  if (std::abs(content) < 1e-15) return Vector(sizes.size(), 0.0);
  const Numeric sign = content < 0.0 ? -1.0 : 1.0;
  content            = std::abs(content);
  const Numeric tc   = std::clamp(point.temperature, t_min_psd, t_max_psd) - 273.15;

  const std::array<Numeric, 5>     q = tropical ? std::array<Numeric, 5>{152., -12.4, 3.28, -0.78, -1.94}
                                                : std::array<Numeric, 5>{141., -16.8, 102., 2.07, -4.82};
  constexpr std::array<Numeric, 3> aq{13.6, -7.76, 0.479};
  constexpr std::array<Numeric, 3> bq{-0.0361, 0.0151, 0.00149};
  constexpr std::array<Numeric, 3> cq{0.807, 0.00581, 0.0457};
  auto                             moment_factors = [&](Numeric n) {
    return std::array{std::exp(aq[0] + aq[1] * n + aq[2] * n * n),
                      bq[0] + bq[1] * n + bq[2] * n * n,
                      cq[0] + cq[1] * n + cq[2] * n * n};
  };

  Numeric m2 = content / a;
  if (b != 2.0) {
    const auto [an, bn, cn] = moment_factors(b);
    m2                      = std::pow(m2 * std::exp(-bn * tc) / an, 1.0 / cn);
  }
  const auto [an, bn, cn] = moment_factors(3.0);
  const Numeric m3        = an * std::exp(bn * tc) * std::pow(m2, cn);
  const Numeric scale     = std::pow(m2, 4) / std::pow(m3, 3);
  Vector        result(sizes.size());
  for (Size i = 0; i < sizes.size(); ++i) {
    const Numeric x   = sizes[i] * m2 / m3;
    const Numeric phi = q[0] * std::exp(q[1] * x) + q[2] * std::pow(x, q[3]) * std::exp(q[4] * x);
    result[i]         = sign * phi * scale;
  }
  return result;
}

PSDData FieldEtAl07::evaluate_with_derivatives(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  return numerical_derivatives(point, sizes, a, b, {mass}, [this](const auto& p, const auto& s, auto aa, auto bb) {
    return evaluate(p, s, aa, bb);
  });
}

McFarquharHeymsfield97::McFarquharHeymsfield97(ScatteringSpeciesProperty mass_,
                                               Numeric                   t_min_,
                                               Numeric                   t_max_,
                                               Numeric                   t_min_psd_,
                                               Numeric                   t_max_psd_,
                                               bool                      picky_)
    : mass(std::move(mass_)),
      t_min(t_min_),
      t_max(t_max_),
      t_min_psd(t_min_psd_),
      t_max_psd(t_max_psd_),
      picky(picky_) {}

Vector McFarquharHeymsfield97::evaluate(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  ARTS_USER_ERROR_IF(not point.has(mass), "The PSD requires {}.", mass)
  ARTS_USER_ERROR_IF(b < 2.9 or b > 3.1 or a < 460.0 or a > 500.0,
                     "MH97 requires an ice-sphere mass-size relationship near a=480, b=3.")
  if (not temperature_allowed(point.temperature, t_min, t_max, picky)) return Vector(sizes.size(), 0.0);
  Numeric iwc = point[mass];
  if (iwc == 0.0) return Vector(sizes.size(), 0.0);
  const Numeric sign = iwc < 0.0 ? -1.0 : 1.0;
  iwc                = std::abs(iwc) * 1e3;  // kg/m3 to g/m3
  const Numeric tc   = std::clamp(point.temperature, t_min_psd, t_max_psd) - 273.15;
  const Numeric rho  = Constant::density_of_ice_at_0c * 1e3;

  const Numeric small_iwc = std::min(iwc, 0.252 * std::pow(iwc, 0.837));
  const Numeric large_iwc = iwc - small_iwc;
  const Numeric alpha     = -4.99e-3 - 0.0494 * std::log10(small_iwc);
  Vector        result(sizes.size(), 0.0);
  if (alpha > 0.0) {
    const Numeric ns = 6.0 * small_iwc * std::pow(alpha, 5) / (Constant::pi * rho * std::tgamma(5.0));
    for (Size i = 0; i < sizes.size(); ++i) {
      const Numeric d_um = sizes[i] * 1e6;
      result[i]          = 1e24 * ns * d_um * std::exp(-alpha * d_um);
    }
  }
  if (large_iwc > 0.0) {
    const Numeric mu    = 5.20 + 0.0013 * tc + (0.026 - 1.2e-3 * tc) * std::log10(large_iwc);
    const Numeric sigma = 0.47 + 2.1e-3 * tc + (0.018 - 2.1e-4 * tc) * std::log10(large_iwc);
    if (mu > 0.0 and sigma > 0.0) {
      const Numeric denom =
          std::pow(Constant::pi, 1.5) * rho * std::sqrt(2.0) * std::exp(3.0 * mu + 4.5 * sigma * sigma) * sigma;
      for (Size i = 0; i < sizes.size(); ++i) {
        const Numeric d_um = sizes[i] * 1e6;
        if (d_um > 0.0) {
          result[i] +=
              1e24 * 6.0 * large_iwc / (denom * d_um) * std::exp(-0.5 * std::pow((std::log(d_um) - mu) / sigma, 2));
        }
      }
    }
  }
  result *= sign;
  return result;
}

PSDData McFarquharHeymsfield97::evaluate_with_derivatives(const AtmPoint& point,
                                                          const Vector&   sizes,
                                                          Numeric         a,
                                                          Numeric         b) const {
  return numerical_derivatives(point, sizes, a, b, {mass}, [this](const auto& p, const auto& s, auto aa, auto bb) {
    return evaluate(p, s, aa, bb);
  });
}

BinnedPSD::BinnedPSD(SizeParameter size_parameter_, Vector bins_, Vector counts_, Numeric t_min_, Numeric t_max_)
    : size_parameter(size_parameter_), bins(bins_), counts(counts_), t_min(t_min_), t_max(t_max_) {
  if (bins_.size() != (counts_.size() + 1)) {
    ARTS_USER_ERROR("The bin vector must have exactly one element more than the counts vector.");
  }
  if (!std::is_sorted(bins.begin(), bins.end())) { ARTS_USER_ERROR("The bins vector must be strictly increasing."); }
}

Vector BinnedPSD::evaluate(const AtmPoint& point,
                           const Vector&   particle_sizes,
                           const Numeric& /*scat_species_a*/,
                           const Numeric& /*scat_species_b*/) const {
  Index  n_parts = particle_sizes.size();
  Vector pnd     = Vector(n_parts);
  for (Index ind = 0; ind < n_parts; ++ind) {
    if ((point.temperature < t_min) || (t_max < point.temperature)) {
      pnd[ind] = 0.0;
    } else {
      Index bin_ind = digitize(bins, particle_sizes[ind]);
      Index n_bins  = bins.size();
      if (bin_ind < 0) {
        pnd[ind] = 0.0;
      } else if (bin_ind >= n_bins) {
        pnd[ind] = 0.0;
      } else {
        pnd[ind] = counts[bin_ind];
      }
    }
  }
  return pnd;
}

PSDData BinnedPSD::evaluate_with_derivatives(const AtmPoint& point, const Vector& sizes, Numeric a, Numeric b) const {
  return {.values = evaluate(point, sizes, a, b), .derivatives = {}};
}

Vector size_bin_widths(const Vector& sizes) {
  const Size n = sizes.size();
  ARTS_USER_ERROR_IF(n < 2,
                     "A particle size distribution given as a density per unit size needs at least two particle sizes "
                     "to weight them, got {}",
                     n)
  for (Size i = 0; i + 1 < n; ++i)
    ARTS_USER_ERROR_IF(not(sizes[i] < sizes[i + 1]),
                       "The particle sizes must ascend strictly to be weighted, but size {} is {} and size {} is {}",
                       i,
                       sizes[i],
                       i + 1,
                       sizes[i + 1])
  Vector w(n);
  w[0]     = sizes[1] - sizes[0];
  w[n - 1] = sizes[n - 1] - sizes[n - 2];
  for (Size i = 1; i + 1 < n; ++i) w[i] = 0.5 * (sizes[i + 1] - sizes[i - 1]);
  return w;
}
}  // namespace scattering
