#include "disort.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <compare.h>
#include <debug.h>
#include <legendre.h>
#include <matpack.h>
#include <matpack_mdspan_helpers_eigen.h>
#include <time_report.h>
#include <xml.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <ranges>
#include <vector>

namespace dc = disort_common;

namespace disort {
delta_m_scaling delta_m_plus(const ConstMatrixView phase_moments, const Index nleg) {
  ARTS_USER_ERROR_IF(nleg < 1, "nleg must be positive, got {}", nleg);
  ARTS_USER_ERROR_IF(phase_moments.ncols() < nleg + 2,
                     "delta-M-plus requires phase moments through degree {}, got {} columns",
                     nleg + 1,
                     phase_moments.ncols());

  delta_m_scaling out{Vector(phase_moments.nrows()), Matrix(phase_moments.nrows(), nleg, 1.0)};
  const bool      use_classical_delta_m = stdr::any_of(
      phase_moments, [nleg](const auto layer) { return layer[nleg] < 1e-4 || layer[nleg + 1] < 0.7 * layer[nleg]; });

  for (Index layer = 0; layer < phase_moments.nrows(); ++layer) {
    const Numeric pn  = phase_moments[layer, nleg];
    const Numeric pn1 = phase_moments[layer, nleg + 1];

    // DISORT 4.0.99 applies these guards globally across all layers.
    if (use_classical_delta_m) {
      out.fraction[layer] = pn;
      continue;
    }

    ARTS_USER_ERROR_IF(pn <= 0.0 || pn1 <= 0.0 || pn == pn1,
                       "Cannot derive delta-M-plus width from P_{}={} and P_{}={} in layer {}",
                       nleg,
                       pn,
                       nleg + 1,
                       pn1,
                       layer);
    const Numeric n = static_cast<Numeric>(nleg);
    const Numeric sigma_sq =
        (Math::pow2(n + 1.0) - Math::pow2(n)) / (std::log(Math::pow2(pn)) - std::log(Math::pow2(pn1)));
    ARTS_USER_ERROR_IF(!(sigma_sq > 0.0) || !std::isfinite(sigma_sq),
                       "Invalid delta-M-plus Gaussian width {} in layer {}",
                       sigma_sq,
                       layer);

    out.fraction[layer] = pn * std::exp(Math::pow2(n) / (2.0 * sigma_sq));
    for (Index degree = 0; degree < nleg; ++degree) {
      const Numeric l            = static_cast<Numeric>(degree);
      out.moments[layer, degree] = std::exp(-Math::pow2(l) / (2.0 * sigma_sq));
    }
  }
  return out;
}

void radiances::resize(AscendingGrid f_grid, DescendingGrid alt_grid_, AziGrid azi_grid_, ZenGrid zen_grid_) {
  freq_grid = std::move(f_grid);
  alt_grid  = std::move(alt_grid_);
  azi_grid  = std::move(azi_grid_);
  zen_grid  = std::move(zen_grid_);

  data.resize(freq_grid.size(), alt_grid.size() - 1, azi_grid.size(), zen_grid.size());
}

namespace {
enum class CoupledType : char { freq_unique, freq_coupled, alt_coupled };

/** Construct the joint storage and overlap classification for two radiance grids.
 *
 * Grids are extended from the first object into the second object. A shared
 * endpoint is represented once.
 *
 * @param freq_grid The first frequency grid
 * @param freq_grid_second The second frequency grid
 * @param alt_grid The first altitude grid
 * @param alt_grid_second The second altitude grid
 * @param azi_grid The first azimuth grid
 * @param azi_grid_second The second azimuth grid
 * @param zen_grid The first zenith grid
 * @param zen_grid_second The second zenith grid
 *
 * @return A pair of the coupled type and the new radiances object with the coupled size.
 */
std::pair<CoupledType, radiances> coupled_size(AscendingGrid         freq_grid,
                                               const AscendingGrid&  freq_grid_second,
                                               DescendingGrid        alt_grid,
                                               const DescendingGrid& alt_grid_second,
                                               AziGrid               azi_grid,
                                               const AziGrid&        azi_grid_second,
                                               ZenGrid               zen_grid,
                                               const ZenGrid&        zen_grid_second) {
  const bool freq_match = freq_grid == freq_grid_second;
  const bool azi_match  = azi_grid == azi_grid_second;
  const bool zen_match  = zen_grid == zen_grid_second;
  const bool alt_match  = alt_grid == alt_grid_second;

  ARTS_USER_ERROR_IF((azi_match + zen_match + freq_match + alt_match) != 3, "Cannot couple disort radiances");

  std::pair<CoupledType, radiances> out;

  if (not freq_match) {
    const bool correct = freq_grid.back() <= freq_grid_second.front();

    if (not correct) {
      throw std::runtime_error(std::format(R"(Incorrectly overlapping frequency grids.
      
original:  {:Bs,}
extension: {:Bs,}
)",
                                           freq_grid,
                                           freq_grid_second));
    }

    const bool coupled = freq_grid.back() == freq_grid_second.front();
    out.first          = coupled ? CoupledType::freq_coupled : CoupledType::freq_unique;

    freq_grid.extend(freq_grid_second | stdv::drop(coupled));
  } else if (not alt_match) {
    const bool correct = alt_grid.back() >= alt_grid_second.front();

    if (not correct) {
      throw std::runtime_error(std::format(R"(Incorrectly overlapping altitude grids.

original:  {:Bs,}
extension: {:Bs,}
)",
                                           alt_grid,
                                           alt_grid_second));
    }

    if (alt_grid.back() != alt_grid_second.front()) {
      throw std::runtime_error(
          "Unique altitude grids not supported for disort radiances as the "
          "altitude grid represent layer boundaries and the two grids would "
          "not be compatible for coupling");
    }

    out.first = CoupledType::alt_coupled;
    alt_grid.extend(alt_grid_second | stdv::drop(1));
  } else {
    throw std::runtime_error(
        "Coupling of disort radiances with different angular grids is "
        "not implemented nor is it intended to be implemented");
  }

  out.second.resize(std::move(freq_grid), std::move(alt_grid), std::move(azi_grid), std::move(zen_grid));

  return out;
}

/** Construct the joint storage and overlap classification for two flux grids.
 *
 * Grids are extended from the first object into the second object. A shared
 * endpoint is represented once.
 *
 * @param freq_grid The first frequency grid
 * @param freq_grid_second The second frequency grid
 * @param alt_grid The first altitude grid
 * @param alt_grid_second The second altitude grid
 * @return A pair of the coupled type and the new fluxes object with the coupled size.
 */
std::pair<CoupledType, fluxes> coupled_size(AscendingGrid         freq_grid,
                                            const AscendingGrid&  freq_grid_second,
                                            DescendingGrid        alt_grid,
                                            const DescendingGrid& alt_grid_second) {
  const bool freq_match = freq_grid == freq_grid_second;
  const bool alt_match  = alt_grid == alt_grid_second;

  ARTS_USER_ERROR_IF((freq_match + alt_match) != 1, "Cannot couple disort fluxes");

  std::pair<CoupledType, fluxes> out;

  if (not freq_match) {
    const bool correct = freq_grid.back() <= freq_grid_second.front();

    if (not correct) {
      throw std::runtime_error(std::format(R"(Incorrectly overlapping frequency grids.
      
original:  {:Bs,}
extension: {:Bs,}
)",
                                           freq_grid,
                                           freq_grid_second));
    }

    const bool coupled = freq_grid.back() == freq_grid_second.front();
    out.first          = coupled ? CoupledType::freq_coupled : CoupledType::freq_unique;

    freq_grid.extend(freq_grid_second | stdv::drop(coupled));
  } else if (not alt_match) {
    const bool correct = alt_grid.back() >= alt_grid_second.front();

    if (not correct) {
      throw std::runtime_error(std::format(R"(Incorrectly overlapping altitude grids.
      
original:  {:Bs,}
extension: {:Bs,}
)",
                                           alt_grid,
                                           alt_grid_second));
    }

    if (alt_grid.back() != alt_grid_second.front()) {
      throw std::runtime_error(
          "Unique altitude grids not supported for disort radiances as the "
          "altitude grid represent layer boundaries and the two grids would "
          "not be compatible for coupling");
    }

    out.first = CoupledType::alt_coupled;

    alt_grid.extend(alt_grid_second | stdv::drop(1));
  } else {
    throw std::runtime_error(
        "Coupling of disort fluxes with different angular grids is "
        "not implemented nor is it intended to be implemented");
  }

  out.second.resize(std::move(freq_grid), std::move(alt_grid));

  return out;
}
}  // namespace

radiances radiances::combine(const radiances& other) const {
  auto [type, out] = coupled_size(
      freq_grid, other.freq_grid, alt_grid, other.alt_grid, azi_grid, other.azi_grid, zen_grid, other.zen_grid);

  switch (type) {
    using enum CoupledType;
    case freq_unique:
      out.data[Range(0, data.nbooks())]                   = data;
      out.data[Range(data.nbooks(), other.data.nbooks())] = other.data;
      break;
    case freq_coupled:
      out.data[Range(0, data.nbooks())]                       = data;
      out.data[Range(data.nbooks() - 1, other.data.nbooks())] = other.data;
      break;
    case alt_coupled:
      out.data[joker, Range(0, data.npages())]                   = data;
      out.data[joker, Range(data.npages(), other.data.npages())] = other.data;
      break;
  }

  return out;
}

void radiances::sort(const Vector& solver_mu) {
  using Conversion::acosd;

  ArrayOfIndex idx(solver_mu.size());
  stdr::iota(idx, 0);
  stdr::sort(idx, [&solver_mu](Index i, Index j) { return solver_mu[i] > solver_mu[j]; });

  Tensor4 data_copy = data;
  for (Size i = 0; i < solver_mu.size(); i++) { data[joker, joker, joker, i] = data_copy[joker, joker, joker, idx[i]]; }
}

void fluxes::resize(AscendingGrid f, DescendingGrid a) {
  freq_grid = std::move(f);
  alt_grid  = std::move(a);
  up.resize(freq_grid.size(), alt_grid.size() - 1);
  down_diffuse.resize(freq_grid.size(), alt_grid.size() - 1);
  down_direct.resize(freq_grid.size(), alt_grid.size() - 1);
  dfdt.resize(freq_grid.size(), alt_grid.size() - 1);
}

[[nodiscard]] fluxes fluxes::combine(const fluxes& other) const {
  auto [type, out] = coupled_size(freq_grid, other.freq_grid, alt_grid, other.alt_grid);

  switch (type) {
    using enum CoupledType;
    case freq_unique: {
      const Range f{0, up.nrows()};
      const Range g{up.nrows(), other.up.nrows()};
      out.up[f]           = up;
      out.up[g]           = other.up;
      out.down_diffuse[f] = down_diffuse;
      out.down_diffuse[g] = other.down_diffuse;
      out.down_direct[f]  = down_direct;
      out.down_direct[g]  = other.down_direct;
      out.dfdt[f]         = dfdt;
      out.dfdt[g]         = other.dfdt;
    } break;
    case freq_coupled: {
      const Range f{0, up.nrows()};
      const Range g{up.nrows() - 1, other.up.nrows()};
      out.up[f]           = up;
      out.up[g]           = other.up;
      out.down_diffuse[f] = down_diffuse;
      out.down_diffuse[g] = other.down_diffuse;
      out.down_direct[f]  = down_direct;
      out.down_direct[g]  = other.down_direct;
      out.dfdt[f]         = dfdt;
      out.dfdt[g]         = other.dfdt;
    } break;
    case alt_coupled: {
      const Range f{0, up.ncols()};
      const Range g{up.ncols(), other.up.ncols()};
      out.up[joker, f]           = up;
      out.up[joker, g]           = other.up;
      out.down_diffuse[joker, f] = down_diffuse;
      out.down_diffuse[joker, g] = other.down_diffuse;
      out.down_direct[joker, f]  = down_direct;
      out.down_direct[joker, g]  = other.down_direct;
      out.dfdt[joker, f]         = dfdt;
      out.dfdt[joker, g]         = other.dfdt;
    } break;
  }

  return out;
}

void BDRF::operator()(MatrixView x, const ConstVectorView& a, const ConstVectorView& b) const {
  ARTS_TIME_REPORT

  assert(static_cast<Size>(x.nrows()) == a.size());
  assert(static_cast<Size>(x.ncols()) == b.size());
  f(x, a, b);
}

Matrix BDRF::operator()(const Vector& a, const Vector& b) const {
  ARTS_TIME_REPORT

  Matrix x(a.size(), b.size());
  f(x, a, b);
  return x;
}

namespace {
/** Return the positive-stream half-range [0, N). */
constexpr Range rf(Size N) { return {0, N}; }
/** Return the negative-stream half-range [N, 2N). */
constexpr Range rb(Size N) { return {N, N}; }

/** Solve the stream-coupling system used by the ordinary polynomial source. */
void ordinary_source_set_k1(mathscr_v_data& data, const ConstMatrixView& G, const ConstVectorView& inv_mu) {
  stdr::copy(inv_mu, data.k1.elem_begin());
  stdr::copy(G | by_elem, data.G.elem_begin());
  solve_inplace(data.k1, data.G, data.solve_work);
}

/** Form eigenmode coefficients for an ordinary polynomial source at one depth. */
void ordinary_source_set_k2(mathscr_v_data&        data,
                            const Numeric          tau,
                            const ConstVectorView& source_poly_coeffs,
                            const ConstVectorView& K) {
  const Index nk = K.size();
  const Index nc = source_poly_coeffs.size();
  const Index n  = nc - 1;

  Numeric power = 1.0;
  for (Index i = n; i >= 0; power *= tau, --i) data.cvec[i] = power;

  for (Index eigen = 0; eigen < nk; ++eigen) {
    data.k2[eigen] = 0.0;
    for (Index i = 0; i < nc; ++i) {
      Numeric factor = 1.0 / (std::pow(K[eigen], i + 1) * Legendre::factorial(n - i));
      for (Index j = 0; j < i; factor *= K[eigen], ++j) {
        data.k2[eigen] += data.cvec[i] * Legendre::factorial(n - j) * source_poly_coeffs[n - j] * factor;
      }
      data.k2[eigen] += data.cvec[i] * Legendre::factorial(n - i) * source_poly_coeffs[n - i] * factor;
    }
  }
}

/** Replace or accumulate an ordinary-source particular field. */
void ordinary_source_update(VectorView             field,
                            mathscr_v_data&        data,
                            const Numeric          tau,
                            const ConstVectorView& source_poly_coeffs,
                            const ConstMatrixView& G,
                            const ConstVectorView& K,
                            const ConstVectorView& inv_mu,
                            const Index            row0 = 0,
                            const Numeric          add  = 0.0) {
  data.resize(K.size(), source_poly_coeffs.size());
  ordinary_source_set_k1(data, G, inv_mu);
  ordinary_source_set_k2(data, tau, source_poly_coeffs, K);
  data.k2 *= data.k1;
  mult(field, G[Range{row0, static_cast<Index>(field.size())}], data.k2, 1.0, add);
}

/** Add an ordinary polynomial-source particular solution to an existing field. */
void ordinary_source_add(VectorView             field,
                         mathscr_v_data&        data,
                         const Numeric          tau,
                         const ConstVectorView& source_poly_coeffs,
                         const ConstMatrixView& G,
                         const ConstVectorView& K,
                         const ConstVectorView& inv_mu) {
  ordinary_source_update(field, data, tau, source_poly_coeffs, G, K, inv_mu, 0, 1.0);
}

/** Assemble ordinary-source boundary and interface terms for the coefficient solve. */
void ordinary_source_terms(MatrixView              SRC0,
                           MatrixView              SRC1,
                           VectorView              SRCB,
                           mathscr_v_data&         data,
                           const ConstVectorView&  tau,
                           const ConstMatrixView&  source_poly_coeffs,
                           const ConstTensor3View& G,
                           const ConstMatrixView&  K,
                           const ConstVectorView&  inv_mu,
                           const Index             N,
                           const Index             nlayers) {
  data.resize(K.ncols(), source_poly_coeffs.ncols());
  ordinary_source_set_k1(data, G[0], inv_mu);
  ordinary_source_set_k2(data, 0.0, source_poly_coeffs[0], K[0]);
  data.k2 *= data.k1;
  mult(SRCB[rf(N)], G[0][Range{N, N}], data.k2);

  for (Index layer = 0; layer < nlayers; ++layer) {
    ordinary_source_set_k2(data, tau[layer], source_poly_coeffs[layer], K[layer]);
    data.k2 *= data.k1;
    mult(SRC0[layer], G[layer], data.k2);

    if (layer < nlayers - 1) {
      ordinary_source_set_k1(data, G[layer + 1], inv_mu);
      ordinary_source_set_k2(data, tau[layer], source_poly_coeffs[layer + 1], K[layer + 1]);
      data.k2 *= data.k1;
      mult(SRC1[layer], G[layer + 1], data.k2);
    }
  }

  const Index last = nlayers - 1;
  ordinary_source_set_k2(data, tau.back(), source_poly_coeffs[last], K[last]);
  data.k2 *= data.k1;
  mult(SRCB[rb(N)], G[last][rf(N)], data.k2);
}
}  // namespace

Numeric main_data::homogeneous(
    const Index m, const Index layer, const Index state, const Index eigen, const Numeric tau) const {
  const Index pair = conservative_pair_index[static_cast<std::size_t>(layer)];
  if (m == 0 and pair >= 0 and (eigen == pair or eigen == pair + N)) {
    const Numeric center = 0.5 * (scaled_tau_arr_with_0[layer] + scaled_tau_arr_with_0[layer + 1]);
    const Numeric s      = tau - center;
    const Numeric kappa  = conservative_pair_kappa[layer];
    const auto    basis  = dc::centered_pair(kappa, s);
    const auto    columns =
        dc::centered_pair_columns(G_collect[m, layer, state, pair], G_collect[m, layer, state, pair + N], basis);
    return columns[eigen == pair ? 0 : 1];
  }

  const Numeric anchor = eigen < N ? scaled_tau_arr_with_0[layer] : scaled_tau_arr_with_0[layer + 1];
  return G_collect[m, layer, state, eigen] * std::exp(K_collect[m, layer, eigen] * (tau - anchor));
}

void main_data::homogeneous_field(VectorView out, const Index m, const Index layer, const Numeric tau) const {
  out = 0.0;

  const Index pair = m == 0 ? conservative_pair_index[static_cast<std::size_t>(layer)] : Index{-1};
  for (Index eigen = 0; eigen < NQuad; ++eigen) {
    if (pair >= 0 and (eigen == pair or eigen == pair + N)) continue;

    const Numeric anchor    = eigen < N ? scaled_tau_arr_with_0[layer] : scaled_tau_arr_with_0[layer + 1];
    const Numeric amplitude = C_collect[m, layer, eigen] * std::exp(K_collect[m, layer, eigen] * (tau - anchor));
    for (Index state = 0; state < NQuad; ++state)
      out[state] = std::fma(amplitude, G_collect[m, layer, state, eigen], out[state]);
  }

  if (pair >= 0) {
    const Numeric center = 0.5 * (scaled_tau_arr_with_0[layer] + scaled_tau_arr_with_0[layer + 1]);
    const Numeric s      = tau - center;
    const Numeric kappa  = conservative_pair_kappa[layer];
    const auto    basis  = dc::centered_pair(kappa, s);
    const auto    amplitudes =
        dc::centered_pair_amplitudes(C_collect[m, layer, pair], C_collect[m, layer, pair + N], basis);

    for (Index state = 0; state < NQuad; ++state) {
      out[state] = std::fma(amplitudes[0],
                            G_collect[m, layer, state, pair],
                            std::fma(amplitudes[1], G_collect[m, layer, state, pair + N], out[state]));
    }
  }
}

Numeric main_data::source_particular(const Index layer, const Index state, const Numeric tau) const {
  return dc::horner_polynomial(Nscoeffs, tau, [&](const Index p) { return source_collect[layer, state, p]; });
}

Numeric main_data::particular(const Index m, const Index layer, const Index state, const Numeric tau) const {
  Numeric value = m == 0 and has_source_poly ? source_particular(layer, state, tau) : 0.0;
  if (has_beam_source) value += B_collect[m, layer, state] * std::exp(-tau / mu0);
  return value;
}

void main_data::solve_for_coefs() {
  ARTS_TIME_REPORT

  const Index ln = NLayers - 1;

  for (Index m = 0; m < NFourier; m++) {
    const bool m_equals_0 = m == 0;
    const bool has_brdf   = m < NBDRF;

    R             = 0.0;
    mathscr_X_pos = 0.0;
    if (has_brdf) {
      brdf_fourier_modes[m](mathscr_D_neg, mu_arr[rf(N)], mu_arr[rb(N)]),
          einsum<"ij", "", "ij", "j", "j">(R, 1 + m_equals_0, mathscr_D_neg, mu_arr[rf(N)], W);
      if (has_beam_source) {
        brdf_fourier_modes[m](mathscr_X_pos.view_as(N, 1), mu_arr[rf(N)], ConstVectorView{-mu0});
        mathscr_X_pos *= mu0 * I0 / Constant::pi;
      }
    }

    const bool has_stable_pair = m == 0 and stdr::any_of(conservative_pair_index, [](const Index i) { return i >= 0; });
    if (not has_stable_pair) {
      // Ordinary-mode storage is an anchored-exponential representation.
      auto RHS_middle = RHS[Range{N, n - NQuad}].view_as(NLayers - 1, NQuad);
      RHS             = 0.0;

      if (has_source_poly and m_equals_0) {
        for (Index i = 0; i < N; ++i) {
          RHS[i]         = -SRCB[i];
          RHS[n - N + i] = -SRCB[i + N];
        }
        std::transform(
            SRC1.elem_begin(), SRC1.elem_end() - NQuad, SRC0.elem_begin(), RHS_middle.elem_begin(), std::minus{});
        if (has_brdf) {
          ordinary_source_update(jvec[rf(N)],
                                 comp_data,
                                 scaled_tau_arr_with_0.back(),
                                 scaled_source_poly_coeffs[ln],
                                 G_collect[0, ln],
                                 K_collect[0, ln],
                                 inv_mu_arr,
                                 N);
          mult(RHS[Range{n - N, N}], R, jvec[rf(N)], 1.0, 1.0);
        }
      }

      if (has_beam_source) {
        if (has_brdf) {
          BDRF_RHS_contribution = mathscr_X_pos;
          mult(BDRF_RHS_contribution, R, B_collect[m, ln, rb(N)], 1.0, 1.0);
        } else {
          BDRF_RHS_contribution = 0.0;
        }
        for (Index l = 0; l < ln; ++l) {
          const Numeric attenuation = std::exp(-scaled_tau_arr_with_0[l + 1] / mu0);
          for (Index state = 0; state < NQuad; ++state)
            RHS_middle[l, state] += (B_collect[m, l + 1, state] - B_collect[m, l, state]) * attenuation;
        }
        for (Index i = 0; i < N; ++i) {
          RHS[i]         += boundary_down[m, i] - B_collect[m, 0, N + i];
          RHS[n - N + i] += boundary_up[m, i] + (BDRF_RHS_contribution[i] - B_collect[m, ln, i]) *
                                                    std::exp(-scaled_tau_arr_with_0.back() / mu0);
        }
      } else {
        RHS[rf(N)]           += boundary_down[m];
        RHS[Range{n - N, N}] += boundary_up[m];
      }

      if (has_brdf)
        mult(BDRF_LHS, R, G_collect[m, ln, rb(N)]);
      else
        BDRF_LHS = 0.0;

      LHSB.zero();
      for (Index j = 0; j < N; ++j) {
        for (Index i = 0; i < N; ++i) {
          LHSB[i, j]                     = G_collect[m, 0, i + N, j];
          LHSB[i, N + j]                 = G_collect[m, 0, i + N, j + N] * expK_collect[m, 0, j];
          LHSB[n - N + i, n - 2 * N + j] = (G_collect[m, ln, i, j] - BDRF_LHS[i, j]) * expK_collect[m, ln, j];
          LHSB[n - N + i, n - N + j]     = G_collect[m, ln, i, j + N] - BDRF_LHS[i, j + N];
        }
      }
      for (Index l = 0; l < ln; ++l) {
        for (Index j = 0; j < N; ++j) {
          const Numeric e1 = 1.0 / expK_collect[m, l, j + N];
          const Numeric e2 = 1.0 / expK_collect[m, l + 1, j + N];
          for (Index i = 0; i < N; ++i) {
            LHSB[N + l * NQuad + i, l * NQuad + j]                     = G_collect[m, l, i, j] * e1;
            LHSB[2 * N + l * NQuad + i, l * NQuad + j]                 = G_collect[m, l, N + i, j] * e1;
            LHSB[N + l * NQuad + i, l * NQuad + 2 * NQuad - N + j]     = -G_collect[m, l + 1, i, N + j] * e2;
            LHSB[2 * N + l * NQuad + i, l * NQuad + 2 * NQuad - N + j] = -G_collect[m, l + 1, N + i, N + j] * e2;
          }
        }
        for (Index state = 0; state < NQuad; ++state) {
          for (Index j = 0; j < N; ++j) {
            LHSB[N + l * NQuad + state, l * NQuad + N + j]     = G_collect[m, l, state, N + j];
            LHSB[N + l * NQuad + state, l * NQuad + 2 * N + j] = -G_collect[m, l + 1, state, j];
          }
        }
      }

      if (LHSB.solve(RHS)) throw std::runtime_error(std::format("Disort failed to converge for Fourier mode {}.", m));
      for (Index l = 0; l < NLayers; ++l) {
        for (Index eigen = 0; eigen < NQuad; ++eigen) {
          const Numeric coefficient = RHS[l * NQuad + eigen];
          C_collect[m, l, eigen]    = coefficient;
          for (Index state = 0; state < NQuad; ++state)
            GC_collect[m, l, state, eigen] = coefficient * G_collect[m, l, state, eigen];
        }
      }
      continue;
    }

    LHSB.zero();
    RHS = 0.0;

    // Prescribed diffuse field at the top, for the downward ordinates.
    for (Index i = 0; i < N; ++i) {
      const Index state = N + i;
      RHS[i]            = boundary_down[m, i] - particular(m, 0, state, 0.0);
      for (Index eigen = 0; eigen < NQuad; ++eigen) LHSB[i, eigen] = homogeneous(m, 0, state, eigen, 0.0);
    }

    // Continuity of the total field at every layer interface.
    for (Index l = 0; l < ln; ++l) {
      const Numeric tau  = scaled_tau_arr_with_0[l + 1];
      const Index   row0 = N + l * NQuad;
      for (Index state = 0; state < NQuad; ++state) {
        const Index row = row0 + state;
        RHS[row]        = particular(m, l + 1, state, tau) - particular(m, l, state, tau);
        for (Index eigen = 0; eigen < NQuad; ++eigen) {
          LHSB[row, l * NQuad + eigen]       = homogeneous(m, l, state, eigen, tau);
          LHSB[row, (l + 1) * NQuad + eigen] = -homogeneous(m, l + 1, state, eigen, tau);
        }
      }
    }

    // Upward field at the lower boundary, including diffuse and direct reflection.
    const Numeric bottom = scaled_tau_arr_with_0.back();
    for (Index i = 0; i < N; ++i) {
      const Index row = n - N + i;
      Numeric     rhs = boundary_up[m, i] - particular(m, ln, i, bottom);
      if (has_beam_source) rhs += mathscr_X_pos[i] * std::exp(-bottom / mu0);
      for (Index j = 0; j < N; ++j) rhs += R[i, j] * particular(m, ln, N + j, bottom);
      RHS[row] = rhs;

      for (Index eigen = 0; eigen < NQuad; ++eigen) {
        Numeric value = homogeneous(m, ln, i, eigen, bottom);
        for (Index j = 0; j < N; ++j) value -= R[i, j] * homogeneous(m, ln, N + j, eigen, bottom);
        LHSB[row, ln * NQuad + eigen] = value;
      }
    }

    if (LHSB.solve(RHS)) { throw std::runtime_error(std::format("Disort failed to converge for Fourier mode {}.", m)); }

    for (Index l = 0; l < NLayers; ++l) {
      for (Index eigen = 0; eigen < NQuad; ++eigen) {
        const Numeric coefficient = RHS[l * NQuad + eigen];
        C_collect[m, l, eigen]    = coefficient;
        for (Index state = 0; state < NQuad; ++state)
          GC_collect[m, l, state, eigen] = coefficient * G_collect[m, l, state, eigen];
      }
    }
  }
}

namespace {
/** Evaluate the rising factorial (x)_n through the Legendre gamma-ratio helper. */
Numeric poch(Index x, Index n) { return Legendre::tgamma_ratio(static_cast<Numeric>(x + n), static_cast<Numeric>(x)); }
}  // namespace

void main_data::diagonalize() {
  ARTS_TIME_REPORT

  stdr::fill(conservative_pair_index, Index{-1});
  conservative_pair_kappa = 0.0;
  B_collect               = 0.0;

  for (Index m = 0; m < NFourier; m++) {
    auto Km = K_collect[m];
    auto Gm = G_collect[m];
    auto Bm = B_collect[m];

    D_temp.resize(N, NLeg - m);
    X_temp.resize(NLeg - m);
    asso_leg_term_pos.resize(NLeg - m, N);
    asso_leg_term_neg.resize(NLeg - m, N);
    asso_leg_term_mu0.resize(NLeg - m);
    weighted_asso_Leg_coeffs_l.resize(NLeg - m);

    auto xtemp = X_temp.view_as(1, NLeg - m);

    const bool m_equals_0_bool = (m == 0);

    fac.resize(NLeg - m);
    for (Index i = m; i < NLeg; i++) fac[i - m] = poch(i + m + 1, -2 * m);

    for (Index i = m; i < NLeg; i++) {
      for (Index j = 0; j < N; j++) {
        asso_leg_term_pos[i - m, j] = Legendre::assoc_legendre(i, m, mu_arr[j]);
        asso_leg_term_neg[i - m, j] = asso_leg_term_pos[i - m, j] * ((i - m) % 2 ? -1.0 : 1.0);
      }
      asso_leg_term_mu0[i - m] = Legendre::assoc_legendre(i, m, -mu0);
    }

    const bool all_asso_leg_term_pos_finite =
        stdr::all_of(asso_leg_term_pos | by_elem, [](auto& x) { return std::isfinite(x); });

    for (Index l = 0; l < NLayers; l++) {
      VectorView K = Km[l];
      MatrixView G = Gm[l];

      for (Index i = 0; i < NLeg - m; i++) {
        weighted_asso_Leg_coeffs_l[i] = weighted_scaled_Leg_coeffs[l, i + m] * fac[i];
      }

      const Numeric scaled_omega_l = scaled_omega_arr[l];

      if (scaled_omega_l != 0.0 or
          (all_asso_leg_term_pos_finite and stdr::any_of(weighted_asso_Leg_coeffs_l, Cmp::gt(0)))) {
        einsum<"ij", "j", "ji">(D_temp, weighted_asso_Leg_coeffs_l, asso_leg_term_pos);
        mult(D_pos, D_temp, asso_leg_term_pos, 0.5 * scaled_omega_l);
        mult(D_neg, D_temp, asso_leg_term_neg, 0.5 * scaled_omega_l);

        einsum<"ij", "i", "ij", "j">(sqr, inv_mu_arr[rf(N)], D_neg, W);
        einsum<"ij", "i", "ij", "j">(apb, inv_mu_arr[rf(N)], D_pos, W);
        diagonal(apb) -= inv_mu_arr[rf(N)];

        // Only the zeroth mode contains a conservative pair or isotropic source polynomial.
        if (m == 0) {
          for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
              const Numeric alpha               = apb[i, j];
              const Numeric beta                = sqr[i, j];
              transport_matrix[l, i, j]         = -alpha;
              transport_matrix[l, i, N + j]     = -beta;
              transport_matrix[l, N + i, j]     = beta;
              transport_matrix[l, N + i, N + j] = alpha;
            }
          }
        }

        amb  = apb;  // still just alpha
        apb += sqr;  // sqr is beta
        amb -= sqr;

        VectorView eval = K[rf(N)];
        MatrixView evec = amb;
        MatrixView AB   = apb;

        mult(sqr, evec, AB);

        diagonalize_inplace(evec, eval, sqr, diag_work);

        Index pair = -1;
        if (m == 0) {
          pair = 0;
          for (Index i = 1; i < N; ++i)
            if (std::abs(eval[i]) < std::abs(eval[pair])) pair = i;

          const Numeric kappa = std::sqrt(std::abs(eval[pair]));
          // Reduced eigenvectors coalesce as 1/kappa near conservation.
          if (dc::use_centered_pair(omega_arr[l], kappa)) {
            conservative_pair_index[static_cast<std::size_t>(l)] = pair;
            conservative_pair_kappa[l]                           = omega_arr[l] == 1.0 ? 0.0 : kappa;
          } else {
            pair = -1;
          }
        }

        for (Index i = 0; i < N; i++) {
          const Numeric sqrt_x = std::sqrt(std::abs(eval[i]));
          K[i]                 = -sqrt_x;
          K[i + N]             = sqrt_x;
        }

        mult(sqr, AB, evec);

        for (Index i = 0; i < N; i++) {
          for (Index j = 0; j < N; j++) {
            if (j == pair) continue;
            const Numeric a = evec[i, j];
            const Numeric b = sqr[i, j] / K[j];
            G[i, j]         = 0.5 * (a - b);
            G[i, j + N]     = 0.5 * (a + b);
            G[i + N, j]     = G[i, j + N];
            G[i + N, j + N] = G[i, j];
          }
        }

        if (pair >= 0) {
          if (has_beam_source) {
            // A stabilized Jordan pair has no invertible eigenvector matrix.
            einsum<"i", "i", "i", "">(X_temp,
                                      weighted_asso_Leg_coeffs_l,
                                      asso_leg_term_mu0,
                                      (scaled_omega_l * I0 * (2 - m_equals_0_bool) / (4 * Constant::pi)));
            mult(jvec[rf(N)].view_as(1, N), xtemp, asso_leg_term_pos, -1);
            jvec[rf(N)] *= inv_mu_arr[rf(N)];
            mult(jvec[rb(N)].view_as(1, N), xtemp, asso_leg_term_neg);
            jvec[rb(N)] *= inv_mu_arr[rf(N)];
            jvec        *= -1.0;

            Gml            = transport_matrix[l];
            diagonal(Gml) += 1.0 / mu0;
            solve_inplace(jvec, Gml, solve_work);
            Bm[l] = jvec;
          }

          // At exact conservation the null vector of alpha+beta is the constant angular field.
          for (Index i = 0; i < N; ++i) jvec[i] = omega_arr[l] == 1.0 ? 1.0 : evec[i, pair];

          // The pair columns store X=[x;x] and -R=[-r;r], where
          // (alpha-beta) r = x.
          for (Index i = 0; i < N; ++i)
            for (Index j = 0; j < N; ++j) D_pos[i, j] = -transport_matrix[l, i, j] + transport_matrix[l, i, N + j];
          solve_inplace(jvec[rf(N)], D_pos);
          for (Index i = 0; i < N; ++i) {
            const Numeric x    = omega_arr[l] == 1.0 ? 1.0 : evec[i, pair];
            const Numeric r    = jvec[i];
            G[i, pair]         = x;
            G[N + i, pair]     = x;
            G[i, N + pair]     = -r;
            G[N + i, N + pair] = r;
          }
          K[pair] = K[N + pair] = 0.0;
          if (omega_arr[l] != 1.0) {
            K[pair]     = -conservative_pair_kappa[l];
            K[N + pair] = conservative_pair_kappa[l];
          }
        } else if (has_beam_source) {
          // Ordinary layers have an invertible eigenbasis.
          einsum<"i", "i", "i", "">(X_temp,
                                    weighted_asso_Leg_coeffs_l,
                                    asso_leg_term_mu0,
                                    (scaled_omega_l * I0 * (2 - m_equals_0_bool) / (4 * Constant::pi)));
          mult(jvec[rf(N)].view_as(1, N), xtemp, asso_leg_term_pos, -1);
          jvec[rf(N)] *= inv_mu_arr[rf(N)];
          mult(jvec[rb(N)].view_as(1, N), xtemp, asso_leg_term_neg);
          jvec[rb(N)] *= inv_mu_arr[rf(N)];
          Gml          = G;
          solve_inplace(jvec, Gml, solve_work);
          for (Index j = 0; j < NQuad; ++j) jvec[j] *= mu0 / (1.0 + K[j] * mu0);
          mult(Bm[l], G, jvec, -1.0);
        }
      } else {
        if (m == 0) {
          transport_matrix[l]           = 0.0;
          diagonal(transport_matrix[l]) = inv_mu_arr;
        }
        G[rf(N), rf(N)]           = 0.0;
        G[rb(N), rb(N)]           = 0.0;
        G[rb(N), rf(N)]           = 0.0;
        G[rf(N), rb(N)]           = 0.0;
        diagonal(G[rb(N), rf(N)]) = 1.0;
        diagonal(G[rf(N), rb(N)]) = 1.0;
        for (Index i = 0; i < N; i++) {
          K[i]     = -inv_mu_arr[i];
          K[i + N] = inv_mu_arr[i];
        }
      }
    }
  }
}

/** Computes the IMS factors
 * 
 * Dependent on:
 * - omega_arr
 * - tau_arr
 * - f_arr
 * - mu0
 * - Leg_coeffs_all
 *
 * Modifies:
 * - scaled_mu0
 * - Leg_coeffs_residue_avg
 * - IMS_scalar
 */
void main_data::set_ims_factors() {
  ARTS_TIME_REPORT

  Numeric cumulative_peak_scattering = 0.0;
  Vector  cumulative_residue(NLeg_all, 0.0);

  for (Index l = 0; l < NLayers; ++l) {
    const Numeric tau_top   = l == 0 ? 0.0 : tau_arr[l - 1];
    const Numeric thickness = tau_arr[l] - tau_top;
    const Numeric weight    = omega_arr[l] * thickness;

    cumulative_peak_scattering += f_arr[l] * weight;

    const Numeric omega_f_avg = cumulative_peak_scattering / tau_arr[l];
    scaled_mu0[l + 1]         = mu0 / (1.0 - omega_f_avg);
    IMS_scalar[l + 1]         = I0 / (4.0 * Constant::pi) * Math::pow2(omega_f_avg) / (1.0 - omega_f_avg);

    for (Index i = 0; i < NLeg_all; ++i) {
      const Numeric residue  = i < NLeg ? f_arr[l] : Leg_coeffs_all[l, i];
      cumulative_residue[i] += residue * weight;

      const Numeric residue_avg =
          cumulative_peak_scattering > 0.0 ? cumulative_residue[i] / cumulative_peak_scattering : 0.0;
      Leg_coeffs_residue_avg[l + 1, i] =
          static_cast<Numeric>(2 * i + 1) * (2.0 * residue_avg - Math::pow2(residue_avg));
    }
  }

  // The top-boundary limiting set samples only the first layer.
  scaled_mu0[0]             = scaled_mu0[1];
  IMS_scalar[0]             = IMS_scalar[1];
  Leg_coeffs_residue_avg[0] = Leg_coeffs_residue_avg[1];
}

void main_data::set_scales() {
  ARTS_TIME_REPORT

  eintra<"i", "i", "i">([](auto om, auto fr) { return 1 - om * fr; }, scale_tau, omega_arr, f_arr);

  scaled_tau_arr_with_0[0] = 0;
  for (Index l = 0; l < NLayers; ++l) {
    const Numeric tau_top        = l == 0 ? 0.0 : tau_arr[l - 1];
    const Numeric thickness      = tau_arr[l] - tau_top;
    scaled_tau_arr_with_0[l + 1] = scaled_tau_arr_with_0[l] + scale_tau[l] * thickness;
  }

  eintra<"i", "i", "i", "i">(
      [](auto om, auto fr, auto st) { return om * (1.0 - fr) / st; }, scaled_omega_arr, omega_arr, f_arr, scale_tau);

  // The transport emission is (1 - omega'_l) B_l(tau(tau')), with physical
  // tau affine in the cumulative delta-M-scaled coordinate tau'.
  scaled_source_poly_coeffs = 0.0;
  if (has_source_poly) {
    for (Index l = 0; l < NLayers; ++l) {
      const Numeric physical_top  = l == 0 ? 0.0 : tau_arr[l - 1];
      const Numeric scaled_top    = scaled_tau_arr_with_0[l];
      const Numeric affine_scale  = 1.0 / scale_tau[l];
      const Numeric affine_offset = physical_top - affine_scale * scaled_top;
      auto          coefficients  = scaled_source_poly_coeffs[l];
      coefficients[0]             = source_poly_coeffs[l, Nscoeffs - 1];
      Index degree                = 0;
      for (Index p = Nscoeffs - 2; p >= 0; --p) {
        coefficients[degree + 1] = affine_scale * coefficients[degree];
        for (Index j = degree; j >= 1; --j)
          coefficients[j] = std::fma(affine_offset, coefficients[j], affine_scale * coefficients[j - 1]);
        coefficients[0] = std::fma(affine_offset, coefficients[0], source_poly_coeffs[l, p]);
        ++degree;
      }
      coefficients *= 1.0 - scaled_omega_arr[l];
    }
  }

  eintra<"ij", "j", "ij", "i", "ij">(
      [](auto j, auto lca, auto f, auto peak) {
        return static_cast<Numeric>(2 * j + 1) * (lca - f * peak) / (1.0 - f);
      },
      weighted_scaled_Leg_coeffs,
      stdv::iota(Index{0}, NLeg),
      Leg_coeffs_all[joker, Range{0, NLeg}],
      f_arr,
      delta_m_peak);
}

void main_data::set_weighted_Leg_coeffs_all() {
  ARTS_TIME_REPORT

  eintra<"ij", "j", "ij">(
      std::multiplies<>{},
      weighted_Leg_coeffs_all,
      stdv::iota(Index{0}, NLeg_all) | stdv::transform([](Index x) { return static_cast<Numeric>(2 * x + 1); }),
      Leg_coeffs_all);
}

void main_data::set_beam_source(const Numeric I0_) {
  ARTS_TIME_REPORT

  has_beam_source = I0_ > 0;

  if (not has_source_poly and has_beam_source and stdr::all_of(boundary_up | by_elem, Cmp::eq<0>()) and
      stdr::all_of(boundary_down | by_elem, Cmp::eq<0>())) {
    I0_orig = I0_;
    I0      = 1;
  } else {
    I0_orig = 1;
    I0      = I0_;
  }
}

void main_data::check_input_size() const {
  ARTS_TIME_REPORT

  dc::check_layer_input_sizes(NLayers, tau_arr, omega_arr);

  ARTS_USER_ERROR_IF((source_poly_coeffs.shape() != std::array{NLayers, Nscoeffs}),
                     "{:B,} vs [{}, {}]",
                     source_poly_coeffs.shape(),
                     NLayers,
                     Nscoeffs);

  ARTS_USER_ERROR_IF(static_cast<Index>(f_arr.size()) != NLayers, "{} vs {}", f_arr.size(), NLayers);

  ARTS_USER_ERROR_IF(
      (delta_m_peak.shape() != std::array{NLayers, NLeg}), "{:B,} vs [{}, {}]", delta_m_peak.shape(), NLayers, NLeg);

  ARTS_USER_ERROR_IF((Leg_coeffs_all.shape() != std::array{NLayers, NLeg_all}),
                     "{:B,} vs [{}, {}]",
                     Leg_coeffs_all.shape(),
                     NLayers,
                     NLeg_all);

  ARTS_USER_ERROR_IF(
      (boundary_up.shape() != std::array{NFourier, N}), "{:B,} vs [{}, {}]", boundary_up.shape(), NFourier, N)

  ARTS_USER_ERROR_IF(
      (boundary_down.shape() != std::array{NFourier, N}), "{:B,} vs [{}, {}]", boundary_down.shape(), NFourier, N)

  ARTS_USER_ERROR_IF(
      brdf_fourier_modes.size() != static_cast<std::size_t>(NBDRF), "{} vs {}", brdf_fourier_modes.size(), NBDRF);
}

void main_data::check_input_value() const {
  ARTS_TIME_REPORT

  dc::check_layer_input_values(tau_arr, omega_arr, mu_arr, mu0, phi0, has_beam_source);

  ARTS_USER_ERROR_IF(stdr::any_of(Leg_coeffs_all,
                                  [](auto&& x) {
                                    return x.size() == 0 or x[0] != 1 or
                                           stdr::any_of(x, [](auto&& u) { return std::abs(u) > 1; });
                                  }),
                     "Leg_coeffs_all must have 1 in the first column and be [-1, 1] elsewhere, got {:B,}",
                     Leg_coeffs_all);

  ARTS_USER_ERROR_IF(I0 < 0, "I0 must be non-negative, got {}", I0);

  ARTS_USER_ERROR_IF(
      stdr::any_of(f_arr, [](auto&& x) { return x >= 1 or x < 0; }), "f_arr must be [0, 1), got {:B,}", f_arr);

  ARTS_USER_ERROR_IF(
      stdr::any_of(delta_m_peak | by_elem, [](const Numeric x) { return !std::isfinite(x) || x < 0.0 || x > 1.0; }),
      "delta_m_peak moments must be finite and [0, 1], got {:B,}",
      delta_m_peak);
}

void main_data::transmission() {
  ARTS_TIME_REPORT

  eintra<"mli", "mli", "l", "l">([](auto k, auto dtau1, auto dtau0) { return std::exp(k * (dtau1 - dtau0)); },
                                 expK_collect,
                                 K_collect,
                                 scaled_tau_arr_with_0[Range(1, NLayers)],
                                 scaled_tau_arr_with_0[Range(0, NLayers)]);

  for (Index j = 0; j < NFourier; j++) { exponent[joker, j, rf(N)] = expK_collect[j, joker, rf(N)]; }
}

void main_data::rad_field() {
  ARTS_TIME_REPORT

  // Ordinary-mode storage contains coefficient-weighted eigenvectors and layer transmissions.
  static_assert(matpack::einsum_optpath<"mi", "mij", "mj">(),
                "On Failure, the einsum has been changed to not use optimal path");
  for (Index l = 0; l < NLayers; ++l) {
    einsum<"mi", "mij", "mj">(um[l], GC_collect[joker, l], exponent[l]);
    if (conservative_pair_index[static_cast<std::size_t>(l)] >= 0)
      homogeneous_field(um[l, 0], 0, l, scaled_tau_arr_with_0[l + 1]);
  }

  if (has_beam_source) {
    for (Index l = 0; l < NLayers; ++l) {
      const Numeric attenuation = std::exp(-scaled_tau_arr_with_0[l + 1] / mu0);
      for (Index m = 0; m < NFourier; ++m)
        for (Index state = 0; state < NQuad; ++state)
          um[l, m, state] = std::fma(attenuation, B_collect[m, l, state], um[l, m, state]);
    }
  }

  if (has_source_poly) {
    for (Index l = 0; l < NLayers; ++l)
      for (Index state = 0; state < NQuad; ++state) um[l, 0, state] += SRC0[l, state];
  }
}

void main_data::source_function() {
  ARTS_TIME_REPORT

  source_collect = 0.0;
  SRC0           = 0.0;
  SRC1           = 0.0;
  SRCB           = 0.0;
  if (not has_source_poly) return;

  const bool all_ordinary = stdr::none_of(conservative_pair_index, [](const Index pair) { return pair >= 0; });
  if (all_ordinary) {
    ordinary_source_terms(SRC0,
                          SRC1,
                          SRCB,
                          comp_data,
                          scaled_tau_arr_with_0[Range{1, NLayers}],
                          scaled_source_poly_coeffs,
                          G_collect[0],
                          K_collect[0],
                          inv_mu_arr,
                          N,
                          NLayers);
  }

  for (Index l = 0; l < NLayers; ++l) {
    // Exact conservative scattering has zero absorption emission and a singular transport matrix.
    if (stdr::all_of(scaled_source_poly_coeffs[l], Cmp::eq<0>())) continue;

    const Index pair = conservative_pair_index[static_cast<std::size_t>(l)];
    if (pair < 0) {
      // Away from the conservative pair, the eigenbasis recurrence and stream solve are equivalent.
      jvec = inv_mu_arr;
      Gml  = G_collect[0, l];
      solve_inplace(jvec, Gml, solve_work);
      comp_data.k2 = 0.0;
      for (Index p = Nscoeffs - 1; p >= 0; --p) {
        for (Index eigen = 0; eigen < NQuad; ++eigen)
          comp_data.k1[eigen] =
              (jvec[eigen] * scaled_source_poly_coeffs[l, p] + static_cast<Numeric>(p + 1) * comp_data.k2[eigen]) /
              K_collect[0, l, eigen];
        mult(source_collect[l, joker, p], G_collect[0, l], comp_data.k1);
        comp_data.k2 = comp_data.k1;
      }
    } else {
      for (Index p = Nscoeffs - 1; p >= 0; --p) {
        for (Index state = 0; state < NQuad; ++state) {
          const Numeric next = p + 1 < Nscoeffs ? source_collect[l, state, p + 1] : 0.0;
          jvec[state]        = inv_mu_arr[state] * scaled_source_poly_coeffs[l, p] + static_cast<Numeric>(p + 1) * next;
        }
        Gml = transport_matrix[l];
        solve_inplace(jvec, Gml, solve_work);
        for (Index state = 0; state < NQuad; ++state) source_collect[l, state, p] = jvec[state];
      }
    }
  }

  if (all_ordinary) return;

  for (Index l = 0; l < NLayers; ++l) {
    const Numeric top    = scaled_tau_arr_with_0[l];
    const Numeric bottom = scaled_tau_arr_with_0[l + 1];
    for (Index state = 0; state < NQuad; ++state) {
      SRC0[l, state] = source_particular(l, state, bottom);
      if (l < NLayers - 1) SRC1[l, state] = source_particular(l + 1, state, bottom);
    }
    if (l == 0)
      for (Index i = 0; i < N; ++i) SRCB[i] = source_particular(l, i + N, top);
    if (l == NLayers - 1)
      for (Index i = 0; i < N; ++i) SRCB[i + N] = source_particular(l, i, bottom);
  }
}

void main_data::update_all(const Numeric I0_) {
  ARTS_TIME_REPORT

  check_input_value();

  set_weighted_Leg_coeffs_all();
  if (I0_ >= 0 or has_beam_source) { set_beam_source(I0_ >= 0 ? I0_ : I0 * I0_orig); }
  set_scales();
  set_ims_factors();
  diagonalize();
  transmission();
  source_function();
  solve_for_coefs();
  rad_field();
}

main_data::main_data(const Index NLayers_,
                     const Index NQuad_,
                     const Index NLeg_,
                     const Index NFourier_,
                     const Index Nscoeffs_,
                     const Index NLeg_all_,
                     const Index NBDRF_)
    : NLayers(NLayers_),
      NQuad(NQuad_),
      NLeg(NLeg_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(Nscoeffs_),
      NLeg_all(NLeg_all_),
      NBDRF(NBDRF_),
      has_source_poly(Nscoeffs > 0),
      // User data
      tau_arr(NLayers),
      omega_arr(NLayers),
      f_arr(NLayers),
      delta_m_peak(NLayers, NLeg, 1.0),
      source_poly_coeffs(NLayers, Nscoeffs),
      Leg_coeffs_all(NLayers, NLeg_all),
      boundary_up(NFourier, N),
      boundary_down(NFourier, N),
      brdf_fourier_modes(NBDRF),
      // Derived data
      scale_tau(NLayers),
      scaled_omega_arr(NLayers),
      scaled_tau_arr_with_0(NLayers + 1),
      scaled_source_poly_coeffs(NLayers, Nscoeffs),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      Leg_coeffs_residue_avg(NLayers + 1, NLeg_all),
      IMS_scalar(NLayers + 1),
      weighted_scaled_Leg_coeffs(NLayers, NLeg),
      weighted_Leg_coeffs_all(NLayers, NLeg_all),
      GC_collect(NFourier, NLayers, NQuad, NQuad),
      C_collect(NFourier, NLayers, NQuad),
      G_collect(NFourier, NLayers, NQuad, NQuad),
      K_collect(NFourier, NLayers, NQuad),
      expK_collect(NFourier, NLayers, NQuad),
      exponent(NLayers, NFourier, NQuad, 1.0),
      um(NLayers, NFourier, NQuad),
      B_collect(NFourier, NLayers, NQuad),
      source_collect(NLayers, NQuad, Nscoeffs),
      transport_matrix(NLayers, NQuad, NQuad),
      conservative_pair_index(static_cast<std::size_t>(NLayers), Index{-1}),
      conservative_pair_kappa(NLayers),
      scaled_mu0(NLayers + 1),
      // Pure compute allocations
      n(NQuad * NLayers),
      RHS(n),
      jvec(NQuad),
      fac(NLeg),
      weighted_asso_Leg_coeffs_l(NLeg),
      asso_leg_term_mu0(NLeg),
      X_temp(NLeg),
      mathscr_X_pos(N),
      E_Lm1L(N),
      E_lm1l(N),
      E_llp1(N),
      BDRF_RHS_contribution(N),
      SRCB(NQuad),
      SRC0(NLayers, NQuad),
      SRC1(NLayers, NQuad),
      Gml(NQuad, NQuad),
      BDRF_LHS(N, NQuad),
      R(N, N),
      mathscr_D_neg(N, N),
      D_pos(N, N),
      D_neg(N, N),
      apb(N, N),
      amb(N, N),
      sqr(N, N),
      asso_leg_term_pos(N, NLeg),
      asso_leg_term_neg(N, NLeg),
      D_temp(N, NLeg),
      solve_work(NQuad),
      diag_work(N),
      LHSB(3 * N - 1, 3 * N - 1, n, n),
      comp_data(NQuad, Nscoeffs) {
  ARTS_TIME_REPORT

  dc::initialize_streams(mu_arr, inv_mu_arr, W);
}

main_data::main_data(const Index       NQuad_,
                     const Index       NLeg_,
                     const Index       NFourier_,
                     AscendingGrid     tau_arr_,
                     Vector            omega_arr_,
                     Matrix            Leg_coeffs_all_,
                     Matrix            boundary_up_,
                     Matrix            boundary_down_,
                     Vector            f_arr_,
                     Matrix            source_poly_coeffs_,
                     std::vector<BDRF> brdf_fourier_modes_,
                     Numeric           mu0_,
                     Numeric           I0_,
                     Numeric           phi0_,
                     Matrix            delta_m_peak_)
    : NLayers(tau_arr_.size()),
      NQuad(NQuad_),
      NLeg(NLeg_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(source_poly_coeffs_.ncols()),
      NLeg_all(Leg_coeffs_all_.ncols()),
      NBDRF(brdf_fourier_modes_.size()),
      has_source_poly(Nscoeffs > 0),
      has_beam_source(I0_ > 0),
      // User data
      tau_arr(std::move(tau_arr_)),
      omega_arr(std::move(omega_arr_)),
      f_arr(std::move(f_arr_)),
      delta_m_peak(std::move(delta_m_peak_)),
      source_poly_coeffs(std::move(source_poly_coeffs_)),
      Leg_coeffs_all(std::move(Leg_coeffs_all_)),
      boundary_up(std::move(boundary_up_)),
      boundary_down(std::move(boundary_down_)),
      brdf_fourier_modes(std::move(brdf_fourier_modes_)),
      mu0(mu0_),
      I0(I0_),
      phi0(phi0_),
      // Derived data
      scale_tau(NLayers),
      scaled_omega_arr(NLayers),
      scaled_tau_arr_with_0(NLayers + 1),
      scaled_source_poly_coeffs(NLayers, Nscoeffs),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      Leg_coeffs_residue_avg(NLayers + 1, NLeg_all),
      IMS_scalar(NLayers + 1),
      weighted_scaled_Leg_coeffs(NLayers, NLeg),
      weighted_Leg_coeffs_all(NLayers, NLeg_all),
      GC_collect(NFourier, NLayers, NQuad, NQuad),
      C_collect(NFourier, NLayers, NQuad),
      G_collect(NFourier, NLayers, NQuad, NQuad),
      K_collect(NFourier, NLayers, NQuad),
      expK_collect(NFourier, NLayers, NQuad),
      exponent(NLayers, NFourier, NQuad, 1.0),
      um(NLayers, NFourier, NQuad),
      B_collect(NFourier, NLayers, NQuad),
      source_collect(NLayers, NQuad, Nscoeffs),
      transport_matrix(NLayers, NQuad, NQuad),
      conservative_pair_index(static_cast<std::size_t>(NLayers), Index{-1}),
      conservative_pair_kappa(NLayers),
      scaled_mu0(NLayers + 1),
      // Pure compute allocations
      n(NQuad * NLayers),
      RHS(n),
      jvec(NQuad),
      fac(NLeg),
      weighted_asso_Leg_coeffs_l(NLeg),
      asso_leg_term_mu0(NLeg),
      X_temp(NLeg),
      mathscr_X_pos(N),
      E_Lm1L(N),
      E_lm1l(N),
      E_llp1(N),
      BDRF_RHS_contribution(N),
      SRCB(NQuad),
      SRC0(NLayers, NQuad),
      SRC1(NLayers, NQuad),
      Gml(NQuad, NQuad),
      BDRF_LHS(N, NQuad),
      R(N, N),
      mathscr_D_neg(N, N),
      D_pos(N, N),
      D_neg(N, N),
      apb(N, N),
      amb(N, N),
      sqr(N, N),
      asso_leg_term_pos(N, NLeg),
      asso_leg_term_neg(N, NLeg),
      D_temp(N, NLeg),
      solve_work(NQuad),
      diag_work(N),
      LHSB(3 * N - 1, 3 * N - 1, n, n),
      comp_data(NQuad, Nscoeffs) {
  ARTS_TIME_REPORT

  dc::initialize_streams(mu_arr, inv_mu_arr, W);

  if (delta_m_peak.empty()) delta_m_peak = Matrix(NLayers, NLeg, 1.0);

  check_input_size();
  update_all(I0_);
}

[[nodiscard]] Index main_data::tau_index(const Numeric tau) const {
  ARTS_TIME_REPORT
  return dc::layer_index(tau_arr, tau);
}

void main_data::u(u_data& data, const Numeric tau, const Numeric phi) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau       = scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.um.resize(NFourier, NQuad);
  for (Index m = 0; m < NFourier; ++m) homogeneous_field(data.um[m], m, l, scaled_tau);

  if (has_beam_source) {
    const Numeric attenuation = std::exp(-scaled_tau / mu0);
    for (Index m = 0; m < NFourier; ++m)
      for (Index state = 0; state < NQuad; ++state) data.um[m, state] += B_collect[m, l, state] * attenuation;
  }
  if (has_source_poly) {
    if (conservative_pair_index[static_cast<std::size_t>(l)] < 0) {
      ordinary_source_add(
          data.um[0], data.src, scaled_tau, scaled_source_poly_coeffs[l], G_collect[0, l], K_collect[0, l], inv_mu_arr);
    } else {
      for (Index state = 0; state < NQuad; ++state) data.um[0, state] += source_particular(l, state, scaled_tau);
    }
  }

  data.intensities.resize(NQuad);
  data.intensities = 0.0;
  for (Index m = 0; m < NFourier; m++) {
    const Numeric cp  = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const auto    umm = data.um[m];
    for (Index i = 0; i < NQuad; i++) { data.intensities[i] += umm[i] * cp; }
  }

  data.intensities *= I0_orig;
}

void main_data::u_user(user_u_data& data, const Numeric tau, const Numeric phi, const ConstVectorView& user_mu) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      tau < 0.0 || tau > tau_arr.back(), "Optical depth must be in [0, {}], got {}", tau_arr.back(), tau);
  ARTS_USER_ERROR_IF(
      stdr::any_of(user_mu, [](const Numeric mu) { return !std::isfinite(mu) || mu == 0.0 || std::abs(mu) > 1.0; }),
      "User polar-angle cosines must be finite, nonzero, and in [-1, 1], got {:B,}",
      user_mu);

  const Index   output_layer = tau_index(tau);
  const Numeric scaled_output =
      scaled_tau_arr_with_0[output_layer + 1] - (tau_arr[output_layer] - tau) * scale_tau[output_layer];

  const auto scattering_source = [&](const Index m, const Index layer, const Numeric mu, const auto& ordinate_values) {
    Numeric result = 0.0;
    for (Index degree = m; degree < NLeg; ++degree) {
      Numeric moment = 0.0;
      for (Index j = 0; j < NQuad; ++j) {
        moment += W[j % N] * Legendre::assoc_legendre(degree, m, mu_arr[j]) * ordinate_values[j];
      }
      result += 0.5 * scaled_omega_arr[layer] * weighted_scaled_Leg_coeffs[layer, degree] *
                poch(degree + m + 1, -2 * m) * Legendre::assoc_legendre(degree, m, mu) * moment;
    }
    return result;
  };

  const auto external_boundary = [&](const Index m, const Numeric mu) {
    const auto nodes  = mu > 0.0 ? mu_arr[rf(N)] : mu_arr[rb(N)];
    const auto values = mu > 0.0 ? boundary_up[m] : boundary_down[m];
    if (data.barycentric_weights.size() != static_cast<Size>(N)) {
      dc::barycentric_weights(data.barycentric_weights, mu_arr[rf(N)]);
    }
    return dc::barycentric_interpolate(nodes, data.barycentric_weights, values, mu);
  };

  static const auto source_quadrature = [] {
    std::pair<Vector, Vector> out{Vector(32), Vector(32)};
    Legendre::GaussLegendre(out.first, out.second);
    return out;
  }();

  data.source.resize(NQuad, Nscoeffs);
  data.particular.resize(NQuad);
  data.intensities.resize(user_mu.size());
  data.intensities = 0.0;

  for (Size iu = 0; iu < user_mu.size(); ++iu) {
    const Numeric mu       = user_mu[iu];
    const Numeric abs_mu   = std::abs(mu);
    const bool    downward = mu < 0.0;

    for (Index m = 0; m < NFourier; ++m) {
      Numeric mode = external_boundary(m, mu);
      if (downward) {
        mode *= std::exp(-scaled_output / abs_mu);
      } else {
        const Index bottom_layer = NLayers - 1;
        Numeric     bottom       = mode;
        if (m < NBDRF) {
          const Vector outgoing{mu};
          Matrix       reflectivity(1, N);
          brdf_fourier_modes[m](reflectivity, outgoing, mu_arr[rb(N)]);
          for (Index j = 0; j < N; ++j) {
            bottom += (1 + (m == 0)) * reflectivity[0, j] * mu_arr[j] * W[j] * um[bottom_layer, m, N + j];
          }
          if (has_beam_source) {
            Matrix       direct_reflectivity(1, 1);
            const Vector beam_direction{-mu0};
            brdf_fourier_modes[m](direct_reflectivity, outgoing, beam_direction);
            bottom +=
                direct_reflectivity[0, 0] * mu0 * I0 / Constant::pi * std::exp(-scaled_tau_arr_with_0.back() / mu0);
          }
        }
        mode = bottom * std::exp(-(scaled_tau_arr_with_0.back() - scaled_output) / abs_mu);
      }

      const Index first_layer = downward ? 0 : output_layer;
      const Index last_layer  = downward ? output_layer : NLayers - 1;
      for (Index layer = first_layer; layer <= last_layer; ++layer) {
        const Numeric layer_top    = scaled_tau_arr_with_0[layer];
        const Numeric layer_bottom = scaled_tau_arr_with_0[layer + 1];
        const Numeric lower        = downward ? layer_top : std::max(layer_top, scaled_output);
        const Numeric upper        = downward ? std::min(layer_bottom, scaled_output) : layer_bottom;
        if (upper <= lower) continue;

        const Index pair = m == 0 ? conservative_pair_index[static_cast<std::size_t>(layer)] : Index{-1};
        for (Index q = 0; q < NQuad; ++q) {
          if (pair >= 0 and (q == pair or q == pair + N)) continue;
          const Numeric source     = scattering_source(m, layer, mu, GC_collect[m, layer, joker, q]);
          const Numeric reference  = q < N ? layer_top : layer_bottom;
          mode                    += source * dc::user_angle_exponential_integral(
                                                  K_collect[m, layer, q], reference, lower, upper, scaled_output, abs_mu, downward);
        }

        if (pair >= 0) {
          const Numeric center = 0.5 * (layer_top + layer_bottom);
          const Numeric kappa  = conservative_pair_kappa[layer];
          const Numeric l0     = scattering_source(m, layer, mu, G_collect[m, layer, joker, pair]);
          const Numeric l1     = scattering_source(m, layer, mu, G_collect[m, layer, joker, pair + N]);
          const Numeric c0     = C_collect[m, layer, pair];
          const Numeric c1     = C_collect[m, layer, pair + N];
          const Numeric ac     = l0 * c0 + l1 * c1;
          const Numeric as     = l0 * c1 + kappa * kappa * l1 * c0;

          const auto [integral_c, integral_s] =
              dc::centered_pair_integrals(kappa, center, lower, upper, scaled_output, abs_mu, downward);
          mode += ac * integral_c + as * integral_s;
        }

        if (has_beam_source) {
          Numeric source = scattering_source(m, layer, mu, B_collect[m, layer]);
          for (Index degree = m; degree < NLeg; ++degree) {
            source += scaled_omega_arr[layer] * I0 * (2 - (m == 0)) / (4 * Constant::pi) *
                      weighted_scaled_Leg_coeffs[layer, degree] * poch(degree + m + 1, -2 * m) *
                      Legendre::assoc_legendre(degree, m, mu) * Legendre::assoc_legendre(degree, m, -mu0);
          }
          mode += source *
                  dc::user_angle_exponential_integral(-1.0 / mu0, 0.0, lower, upper, scaled_output, abs_mu, downward);
        }

        if (has_source_poly && m == 0) {
          const Numeric midpoint  = 0.5 * (lower + upper);
          const Numeric halfwidth = 0.5 * (upper - lower);
          Numeric       integral  = 0.0;
          for (Index k = 0; k < static_cast<Index>(source_quadrature.first.size()); ++k) {
            const Numeric scaled_point = midpoint + halfwidth * source_quadrature.first[k];
            for (Index state = 0; state < NQuad; ++state)
              data.particular[state] = source_particular(layer, state, scaled_point);
            Numeric polynomial = 0.0;
            for (Index coefficient = Nscoeffs - 1; coefficient >= 0; --coefficient) {
              polynomial = std::fma(polynomial, scaled_point, scaled_source_poly_coeffs[layer, coefficient]);
            }
            const Numeric source    = scattering_source(0, layer, mu, data.particular) + polynomial;
            const Numeric distance  = downward ? scaled_output - scaled_point : scaled_point - scaled_output;
            integral               += source_quadrature.second[k] * source * std::exp(-distance / abs_mu) / abs_mu;
          }
          mode += halfwidth * integral;
        }
      }

      data.intensities[iu] += mode * std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    }
    data.intensities[iu] *= I0_orig;
  }
}

void main_data::u0(u0_data& data, const Numeric tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau       = scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.u0.resize(NQuad);
  homogeneous_field(data.u0, 0, l, scaled_tau);

  if (has_beam_source) {
    const Numeric attenuation = std::exp(-scaled_tau / mu0);
    for (Index state = 0; state < NQuad; ++state) data.u0[state] += B_collect[0, l, state] * attenuation;
  }
  if (has_source_poly) {
    if (conservative_pair_index[static_cast<std::size_t>(l)] < 0) {
      ordinary_source_add(
          data.u0, data.src, scaled_tau, scaled_source_poly_coeffs[l], G_collect[0, l], K_collect[0, l], inv_mu_arr);
    } else {
      for (Index state = 0; state < NQuad; ++state) data.u0[state] += source_particular(l, state, scaled_tau);
    }
  }

  data.u0 *= I0_orig;
}

namespace {
/** Return the cosine of the scattering angle between two directions. */
Numeric calculate_nu(const Numeric mu, const Numeric phi, const Numeric mu_p, const Numeric phi_p) {
  const Numeric scl = std::sqrt(1.0 - mu_p * mu_p) * std::cos(phi_p - phi);
  return mu * mu_p + scl * std::sqrt(1.0 - mu * mu);
}

/** Vectorized scattering-angle cosine for several outgoing polar directions. */
void calculate_nu(Vector& nu, const ConstVectorView& mu, const Numeric phi, const Numeric mu_p, const Numeric phi_p) {
  nu.resize(mu.size());

  std::transform(
      mu.begin(), mu.end(), nu.begin(), [mu_p, scl = std::sqrt(1.0 - mu_p * mu_p) * std::cos(phi_p - phi)](auto&& x) {
        return x * mu_p + scl * std::sqrt(1.0 - x * x);
      });
}

}  // namespace

void main_data::TMS(tms_data& data, const Numeric tau, const Numeric phi) const { TMS(data, tau, phi, mu_arr); }

void main_data::TMS(tms_data& data, const Numeric tau, const Numeric phi, const ConstVectorView& mu) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);
  prepare_TMS(data, phi, mu);
  evaluate_TMS(data, tau, mu);
}

void main_data::prepare_TMS(tms_data& data, const Numeric phi, const ConstVectorView& mu) const {
  check_classical_delta_m_correction();

  for (const Numeric x : mu)
    ARTS_USER_ERROR_IF(x == 0.0 || std::abs(x) > 1.0, "Polar-angle cosine ({}) must be nonzero and in [-1, 1]", x);

  // The upward beam factor is regular and remains in mathscr_B.  The
  // downward factor has a removable singularity at mu == mu0 and is folded
  // into the exponential kernel below.
  calculate_nu(data.nu, mu, phi, -mu0, phi0);

  const Index nmu = mu.size();
  data.mathscr_B.resize(NLayers, nmu);
  for (Index j = 0; j < NLayers; j++) {
    for (Index i = 0; i < nmu; i++) {
      const Numeric p_true      = Legendre::legendre_sum(weighted_Leg_coeffs_all[j], data.nu[i]);
      const Numeric p_trun      = Legendre::legendre_sum(weighted_scaled_Leg_coeffs[j], data.nu[i]);
      const Numeric beam_factor = mu[i] > 0.0 ? mu0 / (mu0 + mu[i]) : 1.0;
      data.mathscr_B[j, i] =
          (scaled_omega_arr[j] * I0) / (4 * Constant::pi) * beam_factor * (p_true / (1.0 - f_arr[j]) - p_trun);
    }
  }
}

void main_data::check_classical_delta_m_correction() const {
  for (Index layer = 0; layer < NLayers; ++layer) {
    if (f_arr[layer] == 0.0) continue;
    ARTS_USER_ERROR_IF(stdr::any_of(delta_m_peak[layer], [](const Numeric x) { return x != 1.0; }),
                       "IMS/TMS corrections are unavailable for a non-classical delta-M removed peak; "
                       "use u() or u_user() for delta-M-plus, matching DISORT 4.0.99");
  }
}

void main_data::evaluate_TMS(tms_data& data, const Numeric tau, const ConstVectorView& mu) const {
  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau         = scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  const Index nmu = mu.size();
  data.TMS.resize(nmu);
  data.TMS             = 0.0;
  const Numeric exptau = std::exp(-scaled_tau / mu0);
  for (Index j = 0; j < nmu; ++j) {
    const Numeric abs_mu = std::abs(mu[j]);
    if (mu[j] > 0.0) {
      const Numeric at_bottom = std::exp((scaled_tau - scaled_tau_arr_l) / abs_mu - scaled_tau_arr_l / mu0);
      data.TMS[j]             = data.mathscr_B[l, j] * (exptau - at_bottom);
      for (Index i = l + 1; i < NLayers; ++i) {
        const Numeric top          = scaled_tau_arr_with_0[i];
        const Numeric bottom       = scaled_tau_arr_with_0[i + 1];
        const Numeric at_top       = std::exp((scaled_tau - top) / abs_mu - top / mu0);
        const Numeric at_bottom_i  = std::exp((scaled_tau - bottom) / abs_mu - bottom / mu0);
        data.TMS[j]               += data.mathscr_B[i, j] * (at_top - at_bottom_i);
      }
    } else {
      const Numeric at_top = std::exp((scaled_tau_arr_lm1 - scaled_tau) / abs_mu - scaled_tau_arr_lm1 / mu0);
      data.TMS[j] =
          data.mathscr_B[l, j] * dc::downward_tms_kernel(abs_mu, mu0, scaled_tau - scaled_tau_arr_lm1, exptau, at_top);
      for (Index i = 0; i < l; ++i) {
        const Numeric top       = scaled_tau_arr_with_0[i];
        const Numeric bottom    = scaled_tau_arr_with_0[i + 1];
        const Numeric at_bottom = std::exp((bottom - scaled_tau) / abs_mu - bottom / mu0);
        const Numeric at_top_i  = std::exp((top - scaled_tau) / abs_mu - top / mu0);
        data.TMS[j] += data.mathscr_B[i, j] * dc::downward_tms_kernel(abs_mu, mu0, bottom - top, at_bottom, at_top_i);
      }
    }
  }
}

void main_data::IMS(Vector& ims, const Numeric tau, const Numeric phi, const ims_convention convention) const {
  IMS(ims, tau, phi, mu_arr[rb(N)], convention);
}

void main_data::IMS(Vector&                ims,
                    const Numeric          tau,
                    const Numeric          phi,
                    const ConstVectorView& mu,
                    const ims_convention   convention) const {
  ARTS_TIME_REPORT

  check_classical_delta_m_correction();

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  const Index   l         = tau_index(tau);
  const Numeric tau_top   = l == 0 ? 0.0 : tau_arr[l - 1];
  const Numeric thickness = tau_arr[l] - tau_top;

  ims.resize(mu.size());
  ims             = 0.0;
  const Index nmu = mu.size();
  for (Index i = 0; i < nmu; i++) {
    ARTS_USER_ERROR_IF(
        mu[i] == 0.0 || std::abs(mu[i]) > 1.0, "Polar-angle cosine ({}) must be nonzero and in [-1, 1]", mu[i]);
    if (mu[i] > 0.0) continue;
    const Numeric abs_mu = -mu[i];
    const Numeric nu     = calculate_nu(mu[i], phi, -mu0, phi0);

    if (convention == ims_convention::disort) {
      // The DISORT convention defines IMS only inside the 10-degree incident-beam aureole.
      const Numeric beam_theta = std::acos(-mu0);
      const Numeric ray_theta  = std::acos(mu[i]);
      if (std::abs(beam_theta - ray_theta) > Constant::pi / 18.0) continue;
    }

    const auto correction_for_boundary = [&](const Index k) {
      const Numeric sign = convention == ims_convention::disort ? -1.0 : 1.0;
      return sign * IMS_scalar[k] * Legendre::legendre_sum(Leg_coeffs_residue_avg[k], nu) *
             dc::ims_chi(tau, abs_mu, scaled_mu0[k]);
    };

    if (NLayers == 1) {
      ims[i] = correction_for_boundary(1);
    } else {
      const Numeric lower  = correction_for_boundary(l);
      const Numeric upper  = correction_for_boundary(l + 1);
      const Numeric weight = (tau - tau_top) / thickness;
      ims[i]               = std::lerp(lower, upper, weight);
    }
  }
}

void main_data::u_user_corr(user_u_data&           data,
                            Vector&                ims,
                            tms_data&              tms,
                            const Numeric          tau,
                            const Numeric          phi,
                            const ConstVectorView& user_mu) const {
  ARTS_TIME_REPORT

  u_user(data, tau, phi, user_mu);
  if (!has_beam_source) return;

  TMS(tms, tau, phi, user_mu);
  IMS(ims, tau, phi, user_mu);
  const Index nmu = user_mu.size();
  for (Index i = 0; i < nmu; ++i) data.intensities[i] += I0_orig * (tms.TMS[i] + ims[i]);
}

void main_data::u_corr(u_data&              u_data,
                       Vector&              ims,
                       tms_data&            tms_data,
                       const Numeric        tau,
                       const Numeric        phi,
                       const ims_convention convention) const {
  ARTS_TIME_REPORT

  u(u_data, tau, phi);
  if (!has_beam_source) return;

  TMS(tms_data, tau, phi);
  IMS(ims, tau, phi, convention);

  for (Index i = 0; i < N; i++) { u_data.intensities[i] += I0_orig * tms_data.TMS[i]; }

  for (Index i = N; i < NQuad; i++) { u_data.intensities[i] += I0_orig * (tms_data.TMS[i] + ims[i - N]); }
}

void main_data::gridded_TMS(Tensor3View tms, const Vector& phi) const {
  ARTS_TIME_REPORT

  const Index M = phi.size();
  tms_data    t{};

  Matrix ray_transport(N, NLayers);
  Vector beam_at_boundary(NLayers + 1);
  for (Index l = 0; l <= NLayers; ++l) beam_at_boundary[l] = std::exp(-scaled_tau_arr_with_0[l] / mu0);
  for (Index angle = 0; angle < N; ++angle)
    for (Index l = 0; l < NLayers; ++l)
      ray_transport[angle, l] = std::exp(-(scaled_tau_arr_with_0[l + 1] - scaled_tau_arr_with_0[l]) / mu_arr[angle]);

  // Boundary-to-boundary transport is a downward/upward recurrence with fixed
  // phase coefficients for each azimuth.
  for (Index j = 0; j < M; j++) {
    prepare_TMS(t, phi[j], mu_arr);
    for (Index angle = 0; angle < N; ++angle) {
      const Numeric mu = mu_arr[angle];

      Numeric downward = 0.0;
      for (Index l = 0; l < NLayers; ++l) {
        const Numeric top       = scaled_tau_arr_with_0[l];
        const Numeric bottom    = scaled_tau_arr_with_0[l + 1];
        const Numeric thickness = bottom - top;
        const Numeric transport = ray_transport[angle, l];
        const Numeric at_bottom = beam_at_boundary[l + 1];
        const Numeric at_top    = transport * beam_at_boundary[l];
        downward = transport * downward +
                   t.mathscr_B[l, angle + N] * dc::downward_tms_kernel(mu, mu0, thickness, at_bottom, at_top);
        tms[l, j, angle + N] = downward;
      }

      Numeric upward             = 0.0;
      tms[NLayers - 1, j, angle] = 0.0;
      for (Index l = NLayers - 1; l > 0; --l) {
        const Numeric transport = ray_transport[angle, l];
        const Numeric at_top    = beam_at_boundary[l];
        const Numeric at_bottom = transport * beam_at_boundary[l + 1];
        upward                  = transport * upward + t.mathscr_B[l, angle] * (at_top - at_bottom);
        tms[l - 1, j, angle]    = upward;
      }
    }
  }
}

void main_data::gridded_IMS(Tensor3View ims, const Vector& phi, const ims_convention convention) const {
  ARTS_TIME_REPORT

  check_classical_delta_m_correction();

  const Index M = phi.size();

  Vector        nu;
  const auto    downward_mu = mu_arr[rb(N)];
  const Numeric beam_theta  = std::acos(-mu0);
  Matrix        depth_factor(NLayers, N, 0.0);
  for (Index angle = 0; angle < N; ++angle) {
    const Numeric mu         = -downward_mu[angle];
    const bool    in_aureole = convention == ims_convention::pythonic_disort ||
                               std::abs(beam_theta - std::acos(downward_mu[angle])) <= Constant::pi / 18.0;
    if (!in_aureole) continue;
    for (Index l = 0; l < NLayers; ++l) {
      const Index   boundary = l + 1;
      const Numeric sign     = convention == ims_convention::disort ? -1.0 : 1.0;
      depth_factor[l, angle] = sign * IMS_scalar[boundary] * dc::ims_chi(tau_arr[l], mu, scaled_mu0[boundary]);
    }
  }
  for (Index j = 0; j < M; ++j) {
    calculate_nu(nu, downward_mu, phi[j], -mu0, phi0);
    for (Index angle = 0; angle < N; ++angle) {
      for (Index l = 0; l < NLayers; ++l) {
        const Index boundary = l + 1;
        ims[l, j, angle] = depth_factor[l, angle] * Legendre::legendre_sum(Leg_coeffs_residue_avg[boundary], nu[angle]);
      }
    }
  }
}

void main_data::gridded_u_corr(
    Tensor3View u_data, Tensor3View tms, Tensor3View ims, const Vector& phi, const ims_convention convention) const {
  ARTS_TIME_REPORT

  gridded_u(u_data, phi);

  if (has_beam_source) {
    gridded_TMS(tms, phi);
    gridded_IMS(ims, phi, convention);

    tms                         *= I0_orig;
    u_data                      += tms;
    ims                         *= I0_orig;
    u_data[joker, joker, rb(N)] += ims;
  }
}

flux_values main_data::flux(flux_data& data, const Numeric tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);
  const Index   l                = tau_index(tau);
  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau       = scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  const Numeric direct_beam        = has_beam_source ? dc::direct_beam_flux(I0, mu0, tau) : 0.0;
  const Numeric direct_beam_scaled = has_beam_source ? dc::direct_beam_flux(I0, mu0, scaled_tau) : 0.0;
  u0(data.u0, tau);

  const auto  diffuse = dc::integrate_diffuse(mu_arr[rf(N)], W, [&](const Index i) { return data.u0.u0[i]; });
  flux_values out;
  out.up           = Constant::two_pi * diffuse.upward;
  out.down_diffuse = Constant::two_pi * diffuse.downward - I0_orig * (direct_beam - direct_beam_scaled);
  out.down_direct  = I0_orig * direct_beam;

  Numeric mean_intensity = diffuse.mean_intensity;
  if (has_beam_source) mean_intensity += I0_orig * dc::direct_beam_radiance(I0, mu0, tau) / (4.0 * Constant::pi);

  const Numeric source =
      dc::horner_polynomial(Nscoeffs, tau, [&](const Index coefficient) { return source_poly_coeffs[l, coefficient]; });
  out.dfdt = (1.0 - omega_arr[l]) * 4.0 * Constant::pi * (mean_intensity - source);
  return out;
}

void main_data::gridded_flux(VectorView flux_up,
                             VectorView flux_down_diffuse,
                             VectorView flux_down_direct,
                             VectorView flux_dfdt) const {
  ARTS_TIME_REPORT

  assert(flux_up.size() == static_cast<Size>(NLayers));
  assert(flux_down_diffuse.size() == static_cast<Size>(NLayers));
  assert(flux_down_direct.size() == static_cast<Size>(NLayers));
  assert(flux_dfdt.size() == static_cast<Size>(NLayers));

  for (Index l = 0; l < NLayers; l++) {
    const auto&& u0 = um[l, 0];

    const Numeric direct_beam = has_beam_source ? dc::direct_beam_flux(I0, mu0, tau_arr[l]) : 0.0;
    const Numeric direct_beam_scaled =
        has_beam_source ? dc::direct_beam_flux(I0, mu0, scaled_tau_arr_with_0[l + 1]) : 0.0;
    const auto diffuse = dc::integrate_diffuse(mu_arr[rf(N)], W, [&](const Index i) { return u0[i]; });

    flux_up[l]           = Constant::two_pi * I0_orig * diffuse.upward;
    flux_down_diffuse[l] = I0_orig * (Constant::two_pi * diffuse.downward - direct_beam + direct_beam_scaled);
    flux_down_direct[l]  = I0_orig * direct_beam;

    Numeric mean_intensity = I0_orig * diffuse.mean_intensity;
    if (has_beam_source)
      mean_intensity += I0_orig * dc::direct_beam_radiance(I0, mu0, tau_arr[l]) / (4.0 * Constant::pi);
    const Numeric source = dc::horner_polynomial(
        Nscoeffs, tau_arr[l], [&](const Index coefficient) { return source_poly_coeffs[l, coefficient]; });
    flux_dfdt[l] = (1.0 - omega_arr[l]) * 4.0 * Constant::pi * (mean_intensity - source);
  }
}

ConstMatrixView main_data::layer_um(Size l) const {
  assert(l < static_cast<Size>(NLayers));

  return um[l];
}

void main_data::gridded_u(Tensor3View out, const Vector& phi) const {
  ARTS_TIME_REPORT

  for (Index l = 0; l < NLayers; l++) {
    static_assert(matpack::einsum_optpath<"pi", "im", "pm">(),
                  "On Failure, the einsum has been changed to not use optimal path");
    einsum<"pi", "im", "pm">(
        out[l],
        transpose(um[l]),
        eintra<Matrix, "pm", "p", "m">(
            [i0 = I0_orig, p0 = phi0](auto p, auto m) { return i0 * std::cos(static_cast<Numeric>(m) * (p0 - p)); },
            phi,
            stdv::iota(Index{0}, NFourier)));
  }
}

void main_data::ungridded_flux(VectorView           flux_up,
                               VectorView           flux_down_diffuse,
                               VectorView           flux_down_direct,
                               VectorView           flux_dfdt,
                               const AscendingGrid& tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(flux_up.size() != tau.size() or flux_down_diffuse.size() != tau.size() or
                         flux_down_direct.size() != tau.size() or flux_dfdt.size() != tau.size(),
                     "All ungridded flux outputs must have the same size as tau ({})",
                     tau.size());
  if (tau.empty()) return;

  ARTS_USER_ERROR_IF(tau.front() < 0, "the first tau ({}) must be positive", tau.front());
  ARTS_USER_ERROR_IF(tau.back() > tau_arr.back(),
                     "the last tau ({}) must be less than the last layer ({})",
                     tau.back(),
                     tau_arr.back());

  Vector u0(NQuad);

  Index l = tau_index(tau.front());
  for (Size il = 0; il < tau.size(); il++) {
    while (tau[il] > tau_arr[l]) l++;

    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau       = scaled_tau_arr_l - (tau_arr[l] - tau[il]) * scale_tau[l];

    homogeneous_field(u0, 0, l, scaled_tau);
    for (Index state = 0; state < NQuad; ++state) u0[state] += particular(0, l, state, scaled_tau);

    const Numeric direct_beam        = has_beam_source ? dc::direct_beam_flux(I0, mu0, tau[il]) : 0.0;
    const Numeric direct_beam_scaled = has_beam_source ? dc::direct_beam_flux(I0, mu0, scaled_tau) : 0.0;
    const auto    diffuse            = dc::integrate_diffuse(mu_arr[rf(N)], W, [&](const Index i) { return u0[i]; });
    flux_up[il]                      = Constant::two_pi * I0_orig * diffuse.upward;
    flux_down_diffuse[il] = I0_orig * (Constant::two_pi * diffuse.downward - direct_beam + direct_beam_scaled);
    flux_down_direct[il]  = I0_orig * direct_beam;

    Numeric mean_intensity = I0_orig * diffuse.mean_intensity;
    if (has_beam_source) mean_intensity += I0_orig * dc::direct_beam_radiance(I0, mu0, tau[il]) / (4.0 * Constant::pi);
    const Numeric source = dc::horner_polynomial(
        Nscoeffs, tau[il], [&](const Index coefficient) { return source_poly_coeffs[l, coefficient]; });
    flux_dfdt[il] = (1.0 - omega_arr[l]) * 4.0 * Constant::pi * (mean_intensity - source);
  }
}

void main_data::ungridded_u(Tensor3View out, const AscendingGrid& tau, const Vector& phi) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau.front() < 0, "the first tau ({}) must be positive", tau.front());
  ARTS_USER_ERROR_IF(tau.back() > tau_arr.back(),
                     "the last tau ({}) must be less than the last layer ({})",
                     tau.back(),
                     tau_arr.back());

  Matrix um(NFourier, NQuad);

  const Index Nphi = phi.size();
  Matrix      cp(Nphi, NFourier);
  for (Size p = 0; p < phi.size(); p++) {
    for (Index m = 0; m < NFourier; m++) { cp[p, m] = I0_orig * std::cos(static_cast<Numeric>(m) * (phi0 - phi[p])); }
  }

  Index l = tau_index(tau.front());
  for (Size il = 0; il < tau.size(); il++) {
    while (tau[il] > tau_arr[l]) l++;

    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau       = scaled_tau_arr_l - (tau_arr[l] - tau[il]) * scale_tau[l];

    for (Index m = 0; m < NFourier; ++m) {
      homogeneous_field(um[m], m, l, scaled_tau);
      for (Index state = 0; state < NQuad; ++state) um[m, state] += particular(m, l, state, scaled_tau);
    }

    static_assert(matpack::einsum_optpath<"pi", "im", "pm">(),
                  "On Failure, the einsum has been changed to not use optimal path");
    einsum<"pi", "im", "pm">(out[il], transpose(um), cp);
  }
}

ZenGriddedField1 main_data::gridded_weights() const {
  Vector mu = mu_arr;

  stdr::sort(mu);

  std::transform(mu.begin(), mu.end(), mu.begin(), [](const Numeric& m) { return 180.0 - Conversion::acosd(m); });

  ZenGriddedField1 disort_quadrature{
      .data_name = "Disort quadrature weights", .data = mu, .grid_names = {"Zenith grid"}, .grids = {std::move(mu)}};

  disort_quadrature[rf(N)] = W;
  disort_quadrature[rb(N)] = W;

  return disort_quadrature;
}

coupling_result couple(main_data&    atmosphere,
                       main_data&    subsurface,
                       const Numeric tolerance,
                       const Index   max_iterations,
                       const Numeric relaxation) {
  ARTS_USER_ERROR_IF(max_iterations < 1, "max_iterations must be at least 1, got {}", max_iterations);
  ARTS_USER_ERROR_IF(tolerance < 0, "tolerance must be non-negative, got {}", tolerance);
  ARTS_USER_ERROR_IF(relaxation <= 0 or relaxation > 1, "relaxation must be in (0, 1], got {}", relaxation);

  MatrixView  atm_up       = atmosphere.upward_boundary();
  MatrixView  subsurf_down = subsurface.downward_boundary();
  const Size  NLayer       = atmosphere.tau().size();
  const Size  Nquad        = atmosphere.mu().size();
  const Size  N            = Nquad / 2;
  const Size  Nfourier     = atm_up.nrows();
  const Size  nbrdfs       = atmosphere.brdf_modes().size();
  const Range front        = rf(N);
  const Range back         = rb(N);

  ARTS_USER_ERROR_IF(Nquad != subsurface.mu().size(),
                     "The coupled models must use the same quadrature dimension, got {} vs {}",
                     Nquad,
                     subsurface.mu().size());
  ARTS_USER_ERROR_IF(atm_up.shape() != subsurf_down.shape(),
                     "The coupled models must use identical boundary dimensions, got {:B,} vs. {:B,}",
                     atm_up.shape(),
                     subsurf_down.shape());

  coupling_result   out{};
  Matrix            res(Nfourier, Nquad), atm_sum(Nfourier, N), subsurf_sum(Nfourier, N);
  StridedMatrixView atm_res = res[joker, front], subsurf_res = res[joker, back];
  auto              eatm_sum     = matpack::eigen::as_eigen(atm_sum);
  auto              esubsurf_sum = matpack::eigen::as_eigen(subsurf_sum);

  const Tensor3 ImR = [nbrdfs, N, &atmosphere] -> Tensor3 {
    Tensor3 R(nbrdfs, N, N);
    Matrix  mathscr_D_neg(N, N);

    for (Size m = 0; m < nbrdfs; ++m) {
      const bool m_equals_0_bool = m == 0;
      atmosphere.brdf_modes()[m](mathscr_D_neg, atmosphere.mu()[rf(N)], atmosphere.mu()[rb(N)]),
          einsum<"ij", "", "ij", "j", "j">(
              R[m], -(1 + m_equals_0_bool), mathscr_D_neg, atmosphere.mu()[rf(N)], atmosphere.weights());
      diagonal(R[m]) += 1.0;
    }

    return R;
  }();

  for (; out.iterations < max_iterations and not out.converged; out.iterations++) {
    res = atmosphere.layer_um(NLayer - 1);
    clamp(res, 0.0, std::numeric_limits<Numeric>::max());
    eatm_sum.noalias() = relaxation * atm_res + (1 - relaxation) * atm_up;

    // The atmosphere-to-subsurface transmission factor is 1-R.
    for (Size m = 0; m < nbrdfs; ++m) { eatm_sum.row(m) = ImR[m] * eatm_sum.row(m); }

    res = subsurface.layer_um(0);
    clamp(res, 0.0, std::numeric_limits<Numeric>::max());
    esubsurf_sum.noalias() = relaxation * subsurf_res + (1 - relaxation) * subsurf_down;

    out.max_relative_change = std::max(copy_maxreliff(atm_sum, subsurf_down), copy_maxreliff(subsurf_sum, atm_up));
    out.converged           = out.max_relative_change <= tolerance;

    if (not out.converged) {
      if (out.iterations == 0) {
        // set_beam_source maintains the original and transport-scaled beam intensities.
        atmosphere.update_all(atmosphere.beam_source());
        subsurface.update_all(subsurface.beam_source());
      } else {
        atmosphere.solve_for_coefs();
        subsurface.solve_for_coefs();
        atmosphere.rad_field();
        subsurface.rad_field();
      }
    }
  }

  return out;
}
}  // namespace disort

void DisortSettings::resize(Index          quadrature_dimension_,
                            Index          legendre_polynomial_dimension_,
                            Index          fourier_mode_dimension_,
                            AscendingGrid  f_grid,
                            DescendingGrid alt_grid_) {
  quadrature_dimension          = quadrature_dimension_;
  legendre_polynomial_dimension = legendre_polynomial_dimension_;
  fourier_mode_dimension        = fourier_mode_dimension_;
  const Size nfreq              = f_grid.size();
  const Size nlay               = alt_grid_.size() - 1;

  freq_grid = std::move(f_grid);
  alt_grid  = std::move(alt_grid_);

  solar_source.resize(nfreq);
  solar_zenith_angle.resize(nfreq);
  solar_azimuth_angle.resize(nfreq);
  bidirectional_reflectance_distribution_functions.resize(nfreq, 0);
  optical_thicknesses.resize(nfreq, nlay);
  single_scattering_albedo.resize(nfreq, nlay);
  fractional_scattering.resize(nfreq, nlay);
  delta_m_peak_moments.resize(nfreq, nlay, legendre_polynomial_dimension);
  delta_m_peak_moments = 1.0;
  source_polynomial.resize(nfreq, nlay, 0);
  legendre_coefficients.resize(nfreq, nlay, legendre_polynomial_dimension);
  upward_boundary_condition.resize(nfreq, fourier_mode_dimension, quadrature_dimension / 2);
  downward_boundary_condition.resize(nfreq, fourier_mode_dimension, quadrature_dimension / 2);
}

void DisortSettings::check() const {
  const Index nfreq = freq_grid.size();
  const Index nlay  = alt_grid.size() - 1;

  ARTS_USER_ERROR_IF(
      solar_source.shape() != std::array{nfreq} or solar_zenith_angle.shape() != std::array{nfreq} or
          solar_azimuth_angle.shape() != std::array{nfreq} or
          (bidirectional_reflectance_distribution_functions.shape() !=
           std::array{nfreq, bidirectional_reflectance_distribution_functions.ncols()}) or
          (optical_thicknesses.shape() != std::array{nfreq, nlay}) or
          (single_scattering_albedo.shape() != std::array{nfreq, nlay}) or
          (fractional_scattering.shape() != std::array{nfreq, nlay}) or
          (delta_m_peak_moments.shape() != std::array{nfreq, nlay, legendre_polynomial_dimension}) or
          (source_polynomial.shape() != std::array{nfreq, nlay, source_polynomial.ncols()}) or
          (legendre_coefficients.shape() != std::array{nfreq, nlay, legendre_coefficients.ncols()}) or
          (upward_boundary_condition.shape() != std::array{nfreq, fourier_mode_dimension, quadrature_dimension / 2}) or
          (downward_boundary_condition.shape() !=
           std::array{nfreq, fourier_mode_dimension, quadrature_dimension / 2}) or
          legendre_polynomial_dimension > legendre_coefficients.ncols(),
      R"-x-(Input is of incorrect size.

{:Bs,}

Also note that the reduced Legendre polynomial dimension is {}.  It must be at most {}.
)-x-",
      *this,
      legendre_polynomial_dimension,
      legendre_coefficients.ncols());
}

disort::main_data DisortSettings::init() const {
  check();
  return disort::main_data(alt_grid.size() - 1,
                           quadrature_dimension,
                           legendre_polynomial_dimension,
                           fourier_mode_dimension,
                           source_polynomial.ncols(),
                           legendre_coefficients.ncols(),
                           bidirectional_reflectance_distribution_functions.ncols());
}

disort::main_data& DisortSettings::set(disort::main_data& dis, Index iv) const {
  using Conversion::cosd;
  using Conversion::deg2rad;

  for (Index i = 0; i < bidirectional_reflectance_distribution_functions.ncols(); i++) {
    dis.brdf_modes()[i] = bidirectional_reflectance_distribution_functions[iv, i];
  }

  dis.tau(optical_thicknesses[iv]);
  dis.solar_zenith()         = cosd(solar_zenith_angle[iv]);
  dis.beam_azimuth()         = deg2rad(solar_azimuth_angle[iv]);
  dis.omega()                = single_scattering_albedo[iv];
  dis.f()                    = fractional_scattering[iv];
  dis.delta_m_peak_moments() = delta_m_peak_moments[iv];
  dis.all_legendre_coeffs()  = legendre_coefficients[iv];
  dis.upward_boundary()      = upward_boundary_condition[iv];
  dis.downward_boundary()    = downward_boundary_condition[iv];
  dis.source_poly()          = source_polynomial[iv];

  dis.update_all(solar_source[iv]);

  return dis;
}

#ifdef ENABLE_CDISORT
disort::main_data& DisortSettings::set_cdisort(disort::main_data& dis, Index iv) const {
  using Conversion::cosd;
  using Conversion::deg2rad;

  for (Index i = 0; i < bidirectional_reflectance_distribution_functions.ncols(); i++) {
    dis.brdf_modes()[i] = bidirectional_reflectance_distribution_functions[iv, i];
  }

  dis.tau(optical_thicknesses[iv]);
  dis.solar_zenith()        = cosd(solar_zenith_angle[iv]);
  dis.beam_azimuth()        = deg2rad(solar_azimuth_angle[iv]);
  dis.omega()               = single_scattering_albedo[iv];
  dis.f()                   = fractional_scattering[iv];
  dis.all_legendre_coeffs() = legendre_coefficients[iv];
  dis.upward_boundary()     = upward_boundary_condition[iv];
  dis.downward_boundary()   = downward_boundary_condition[iv];

  return dis;
}
#endif

void xml_io_stream<DisortBDRF>::read(std::istream& is, DisortBDRF& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.f, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<DisortBDRF>::write(std::ostream& os, const DisortBDRF& x, bofstream* pbofs, std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.f, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<DisortSettings>::read(std::istream& is_xml, DisortSettings& v, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  xml_read_from_stream(is_xml, v.quadrature_dimension, pbifs);
  xml_read_from_stream(is_xml, v.legendre_polynomial_dimension, pbifs);
  xml_read_from_stream(is_xml, v.fourier_mode_dimension, pbifs);
  xml_read_from_stream(is_xml, v.freq_grid, pbifs);
  xml_read_from_stream(is_xml, v.alt_grid, pbifs);
  xml_read_from_stream(is_xml, v.solar_azimuth_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_zenith_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_source, pbifs);
  xml_read_from_stream(is_xml, v.bidirectional_reflectance_distribution_functions, pbifs);
  xml_read_from_stream(is_xml, v.optical_thicknesses, pbifs);
  xml_read_from_stream(is_xml, v.single_scattering_albedo, pbifs);
  xml_read_from_stream(is_xml, v.fractional_scattering, pbifs);
  xml_read_from_stream(is_xml, v.delta_m_peak_moments, pbifs);
  xml_read_from_stream(is_xml, v.source_polynomial, pbifs);
  xml_read_from_stream(is_xml, v.legendre_coefficients, pbifs);
  xml_read_from_stream(is_xml, v.upward_boundary_condition, pbifs);
  xml_read_from_stream(is_xml, v.downward_boundary_condition, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
}

void xml_io_stream<DisortSettings>::write(std::ostream&         os_xml,
                                          const DisortSettings& v,
                                          bofstream*            pbofs,
                                          std::string_view) {
  XMLTag open_tag;
  XMLTag close_tag;

  open_tag.name = type_name;
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(os_xml, v.quadrature_dimension, pbofs, "quadrature_dimension");
  xml_write_to_stream(os_xml, v.legendre_polynomial_dimension, pbofs, "legendre_polynomial_dimension");
  xml_write_to_stream(os_xml, v.fourier_mode_dimension, pbofs, "fourier_mode_dimension");
  xml_write_to_stream(os_xml, v.freq_grid, pbofs, "nfreq");
  xml_write_to_stream(os_xml, v.alt_grid, pbofs, "nlay");
  xml_write_to_stream(os_xml, v.solar_azimuth_angle, pbofs, "solazi");
  xml_write_to_stream(os_xml, v.solar_zenith_angle, pbofs, "solzen");
  xml_write_to_stream(os_xml, v.solar_source, pbofs, "solsrc");
  xml_write_to_stream(os_xml, v.bidirectional_reflectance_distribution_functions, pbofs, "BRDF");
  xml_write_to_stream(os_xml, v.optical_thicknesses, pbofs, "Tau");
  xml_write_to_stream(os_xml, v.single_scattering_albedo, pbofs, "albedo");
  xml_write_to_stream(os_xml, v.fractional_scattering, pbofs, "fractional_scattering");
  xml_write_to_stream(os_xml, v.delta_m_peak_moments, pbofs, "delta_m_peak_moments");
  xml_write_to_stream(os_xml, v.source_polynomial, pbofs, "source_polynomial");
  xml_write_to_stream(os_xml, v.legendre_coefficients, pbofs, "legendre_coefficients");
  xml_write_to_stream(os_xml, v.upward_boundary_condition, pbofs, "upward_boundary_condition");
  xml_write_to_stream(os_xml, v.downward_boundary_condition, pbofs, "downward_boundary_condition");

  close_tag.name = type_name;
  close_tag.write_to_end_stream(os_xml);
}
