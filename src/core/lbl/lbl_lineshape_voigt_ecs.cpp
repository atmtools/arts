#include "lbl_lineshape_voigt_ecs.h"

#include <arts_omp.h>
#include <atm.h>
#include <configtypes.h>
#include <debug.h>
#include <isotopologues.h>
#include <jacobian.h>
#include <partfun.h>
#include <physics_funcs.h>
#include <sorting.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <ranges>
#include <stdexcept>

#include "lbl_lineshape_linemixing.h"
#include "lbl_lineshape_model.h"
#include "lbl_lineshape_voigt_ecs_hartmann.h"
#include "lbl_lineshape_voigt_ecs_makarov.h"
#include "lbl_lineshape_voigt_lte.h"

#undef WIGNER3
#undef WIGNER6

namespace lbl::voigt::ecs {
ComputeData::ComputeData(const ConstVectorView&   f_grid,
                         const AtmPoint&          atm,
                         const Vector2&           los,
                         const ZeemanPolarization pol)
    : scl(f_grid.size()), shape(f_grid.size()) {
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 scl.begin(),
                 [N = number_density(atm.pressure, atm.temperature), T = atm.temperature](auto f) {
                   const Numeric r = (Constant::h * f) / (Constant::k * T);
                   return -N * f * std::expm1(-r);
                 });

  update_zeeman(los, atm.mag, pol);
}

void ComputeData::update_zeeman(const Vector2& los, const Vector3& mag, const ZeemanPolarization pol) {
  npm = zeeman::norm_view(pol, mag, los);
}

void ComputeData::core_calc_eqv() {
  // Ws is stored as the spectral operator itself: the same orientation the
  // resolvent and the eigendecomposition use, so nothing is transposed here.
  // The carrier frequency is subtracted from Ws and dW in place; both are
  // rebuilt by every adapt() call and are not read again after this one.
  const Size n = pop.size(), m = vmrs.size(), nt = dW.npages();
  ARTS_USER_ERROR_IF(dip.size() != n or Ws.npages() != static_cast<Index>(m) or Ws.nrows() != static_cast<Index>(n) or
                         Ws.ncols() != static_cast<Index>(n),
                     "Inconsistent ECS matrix and population dimensions")
  ARTS_USER_ERROR_IF(nt and (m != 1 or dW.nrows() != static_cast<Index>(n) or dW.ncols() != static_cast<Index>(n) or
                             dpop.nrows() != static_cast<Index>(nt) or dpop.ncols() != static_cast<Index>(n) or
                             ddip.shape() != dpop.shape()),
                     "Inconsistent ECS Jacobian dimensions")
  eqv_strs.resize(m, n);
  eqv_vals.resize(m, n);
  V.resize(n, n);
  deqv_strs.resize(nt, n);
  deqv_vals.resize(nt, n);
  eigenvector_rcond.resize(m);
  dVs.resize(nt, n, n);
  eqv_strs          = 0;
  deqv_strs         = 0;
  deqv_vals         = 0;
  eigenvector_rcond = 1;
  if (n == 0) return;

  ComplexVector rhs(n), coefficients(n);
  ComplexMatrix derivative_rhs(n, nt), derivative_coefficients(n, nt);
  Vector        dcenter(nt);
  for (Size j = 0; j < n; ++j) rhs[j] = pop[j] * dip[j];
  complex_diagonalize_workdata workspace(n);

  for (Size k = 0; k < m; ++k) {
    auto          W       = Ws[k];
    auto          eqv_str = eqv_strs[k];
    auto          eqv_val = eqv_vals[k];
    const Numeric center  = W[0, 0].real();
    for (Size j = 0; j < n; ++j) W[j, j] -= center;
    Numeric matrix_norm = 0;
    for (Size i = 0; i < n; ++i) {
      Numeric row_sum = 0;
      for (Size j = 0; j < n; ++j) row_sum += std::abs(W[i, j]);
      matrix_norm = std::max(matrix_norm, row_sum);
    }
    ARTS_USER_ERROR_IF(not std::isfinite(matrix_norm), "ECS centered matrix norm overflowed")

    bool changes_operator = false;
    if (nt) {
      for (Size t = 0; t < nt; ++t) {
        dcenter[t] = dW[t, 0, 0].real();
        for (Size j = 0; j < n; ++j) dW[t, j, j] -= dcenter[t];
        for (Size i = 0; i < n; ++i)
          for (Size j = 0; j < n; ++j) changes_operator |= dW[t, i, j] != Complex{};
      }
    }
    if (changes_operator) {
      diagonalize(V, eqv_val, dVs, deqv_vals, W, dW, workspace);
    } else {
      // Fixed operators (e.g. frequency or population-only targets) do not
      // require derivatives of an eigenbasis, even if modes are degenerate.
      diagonalize(V, eqv_val, W, workspace);
      dVs = 0;
    }
    eigenvector_rcond[k]            = solve(coefficients, V, rhs, 1e-12);
    const Numeric damping_tolerance = 64 * std::numeric_limits<Numeric>::epsilon() * matrix_norm;
    for (Size i = 0; i < n; ++i) {
      ARTS_USER_ERROR_IF(eqv_val[i].imag() < -damping_tolerance,
                         "ECS relaxation matrix has a negative damping eigenvalue: {} Hz",
                         eqv_val[i])
      if (eqv_val[i].imag() < 0) {
        eqv_val[i].imag(0);
        for (Size t = 0; t < nt; ++t) deqv_vals[t, i].imag(0);
      }
      eqv_val[i]         += center;
      Complex projection  = 0;
      for (Size j = 0; j < n; ++j) projection += dip[j] * V[j, i];
      eqv_str[i] = projection * coefficients[i];
    }
    for (Size t = 0; t < nt; ++t) {
      for (Size i = 0; i < n; ++i) {
        derivative_rhs[i, t] = dpop[t, i] * dip[i] + pop[i] * ddip[t, i];
        for (Size j = 0; j < n; ++j) derivative_rhs[i, t] -= dVs[t, i, j] * coefficients[j];
      }
    }
    if (nt) solve(derivative_coefficients, V, derivative_rhs, 1e-12);
    for (Size t = 0; t < nt; ++t) {
      for (Size i = 0; i < n; ++i) {
        deqv_vals[t, i]    += dcenter[t];
        Complex projection = 0, derivative_projection = 0;
        for (Size j = 0; j < n; ++j) {
          projection            += dip[j] * V[j, i];
          derivative_projection += ddip[t, j] * V[j, i] + dip[j] * dVs[t, j, i];
        }
        deqv_strs[t, i] = derivative_projection * coefficients[i] + projection * derivative_coefficients[i, t];
      }
    }
  }
}

void ComputeData::core_calc(const ConstVectorView& f_grid) try {
  ARTS_USER_ERROR_IF(dgd_fac.size() != static_cast<Size>(dW.npages()) or df.size() != dgd_fac.size(),
                     "Inconsistent ECS profile Jacobian dimensions")
  core_calc_eqv();
  const Size m = vmrs.size(), nf = f_grid.size(), nt = dW.npages();
  shape.resize(nf);
  shape = 0;
  dshape.resize(nt, nf);
  dshape = 0;
  for (Size k = 0; k < m; ++k) {
    if (vmrs[k] == 0) continue;
    for (Size i = 0; i < eqv_strs[k].size(); ++i) {
      const Numeric gamd = gd_fac * eqv_vals[k, i].real();
      ARTS_USER_ERROR_IF(not std::isfinite(gamd) or gamd <= 0,
                         "ECS Gaussian 1/e half-width must be positive and finite, got {} Hz",
                         gamd)
      const Numeric inv_gamd = 1 / gamd;
      for (Size iv = 0; iv < nf; ++iv) {
        const Complex z   = (eqv_vals[k, i] - f_grid[iv]) * inv_gamd;
        const Complex w   = lte::single_shape::F(z);
        shape[iv]        += vmrs[k] * eqv_strs[k, i] * w * inv_gamd;
        const Complex dw  = nt ? lte::single_shape::dF(z, w) : Complex{};
        for (Size t = 0; t < nt; ++t) {
          const Numeric dgamd  = dgd_fac[t] * eqv_vals[k, i].real() + gd_fac * deqv_vals[t, i].real();
          const Complex dz     = (deqv_vals[t, i] - df[t] - z * dgamd) * inv_gamd;
          dshape[t, iv]       += (deqv_strs[t, i] * w + eqv_strs[k, i] * (dw * dz - w * dgamd * inv_gamd)) * inv_gamd;
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
void get_vmrs(VectorView vmrs, const line_shape::model::map_t& mod, const AtmPoint& atm) {
  std::transform(mod.begin(), mod.end(), vmrs.begin(), [&atm](const auto& m) {
    const Numeric x = m.first == SpeciesEnum::Bath or not atm.has(m.first) ? 0.0 : atm[m.first];
    ARTS_USER_ERROR_IF(not std::isfinite(x) or x < 0, "Invalid ECS broadener VMR for {}: {}", m.first, x)
    return x;
  });

  const Numeric total = sum(vmrs);
  ARTS_USER_ERROR_IF(not std::isfinite(total), "Non-finite total ECS broadener VMR")
  const Size bath_spec = std::distance(mod.begin(), mod.find(SpeciesEnum::Bath));
  if (bath_spec != mod.size()) {
    ARTS_USER_ERROR_IF(total > 1 + 1e-12, "ECS explicit broadener VMRs exceed unity: {}", total)
    vmrs[bath_spec] = std::max(Numeric{0}, 1 - total);
  } else if (total > 0) {
    vmrs /= total;
  }
  // With no bath and no present collision partners the collision matrix is zero,
  // matching ordinary line-shape mixing, rather than producing 0/0.
}

// Composition derivatives include normalization in the same convention as get_vmrs.
void get_dvmrs(MatrixView                      out,
               ConstVectorView                 fractions,
               const line_shape::model::map_t& models,
               const AtmPoint&                 atm,
               const Jacobian::Targets&        targets) {
  out              = 0;
  const bool bath  = models.contains(SpeciesEnum::Bath);
  Numeric    total = 0;
  for (auto spec : models | stdv::keys)
    if (spec != SpeciesEnum::Bath and atm.has(spec)) total += atm[spec];
  for (const auto& target : targets.atm) {
    const auto* species = std::get_if<SpeciesEnum>(&target.type);
    if (not species or *species == SpeciesEnum::Bath or not models.contains(*species)) continue;
    Size i = 0;
    for (auto spec : models | stdv::keys) {
      if (bath) {
        out[target.target_pos, i] =
            spec == SpeciesEnum::Bath ? Numeric(total <= 1 ? -1 : 0) : Numeric(spec == *species);
      } else if (total > 0) {
        out[target.target_pos, i] = (Numeric(spec == *species) - fractions[i]) / total;
      } else {
        ARTS_USER_ERROR("ECS VMR derivatives require a nonzero total broadener VMR when no Bath is present")
      }
      ++i;
    }
  }
}

Numeric broadener_mass_derivative(const AtmPoint& atm, SpeciesEnum broadener, const AtmKeyVal& key) {
  if (const auto* species = std::get_if<SpeciesEnum>(&key)) {
    if (broadener != SpeciesEnum::Bath or *species == SpeciesEnum::Bath) return 0;
    Numeric total = 0;
    for (const auto& [spec, vmr] : atm.specs) total += vmr;
    return (atm.mean_mass(*species) - atm.mean_mass()) / total;
  }
  if (const auto* isot = std::get_if<SpeciesIsotope>(&key)) {
    if (isot->is_joker() or isot->is_predefined() or (broadener != SpeciesEnum::Bath and broadener != isot->spec))
      return 0;
    Numeric ratio = 0;
    for (const auto& [other, value] : atm.isots)
      if (other.spec == isot->spec and not(other.is_joker() or other.is_predefined())) ratio += value;
    Numeric derivative = (isot->mass - atm.mean_mass(isot->spec)) / ratio;
    if (broadener == SpeciesEnum::Bath) {
      Numeric total = 0;
      for (const auto& [spec, vmr] : atm.specs) total += vmr;
      derivative *= atm.has(isot->spec) ? atm[isot->spec] / total : 0;
    }
    return derivative;
  }
  return 0;
}

// Target dispatch is kept outside the angular loops. The kernels consume flat
// derivative arrays for every target, just as they consume Q and Omega.
void prepare_basis_jacobian(MatrixView                      dQ,
                            MatrixView                      dOmega,
                            const energy_data&              energies,
                            const linemixing::species_data& collision,
                            Numeric                         T0,
                            const SpeciesIsotope&           isot,
                            SpeciesEnum                     broadener,
                            const AtmPoint&                 atm,
                            const Jacobian::Targets&        targets) {
  dQ     = 0;
  dOmega = 0;
  if (dQ.nrows() == 0) return;
  const Numeric T            = atm.temperature;
  const Numeric mass         = broadener == SpeciesEnum::Bath ? atm.mean_mass() : atm.mean_mass(broadener);
  const Numeric inverse_mass = 1 / mass + 1 / isot.mass;
  const Numeric velocity2    = 8 * Constant::k * T * inverse_mass / (Constant::m_u * Constant::pi);
  const Numeric scaling = collision.scaling(T0, T), beta = collision.beta(T0, T), lambda = collision.lambda(T0, T),
                distance = collision.collisional_distance(T0, T);
  for (const auto& target : targets.atm) {
    Numeric dT = 0, ds = 0, db = 0, dl = 0, dd = 0;
    if (const auto* key = std::get_if<AtmKey>(&target.type)) {
      switch (*key) {
        case AtmKey::t:
          dT = 1;
          ds = collision.scaling.dT(T0, T);
          db = collision.beta.dT(T0, T);
          dl = collision.lambda.dT(T0, T);
          dd = collision.collisional_distance.dT(T0, T);
          break;
        default: break;
      }
    }
    const Numeric dm                 = broadener_mass_derivative(atm, broadener, target.type);
    const Numeric relative_velocity2 = dT / T - dm / (mass * mass * inverse_mass);
    for (Index L = 0; L < dQ.ncols(); ++L) {
      const Numeric gap = (energies.rotational[L] - energies.rotational_minus_two[L]) / Constant::h_bar;
      const Numeric x   = Math::pow2(gap * distance) / (24 * velocity2);
      const Numeric dx  = gap * gap * (2 * distance * dd - distance * distance * relative_velocity2) / (24 * velocity2);
      dOmega[target.target_pos, L] = -2 * dx / Math::pow3(1 + x);
      if (L == 0) continue;
      const Numeric e          = energies.rotational[L] / (Constant::k * T);
      const Numeric angular    = Numeric(L) * Numeric(L + 1);
      dQ[target.target_pos, L] = std::exp(-beta * e) / std::pow(angular, lambda) *
                                 (ds + scaling * (-db * e + beta * e * dT / T - dl * std::log(angular)));
    }
  }
}

using Offdiagonal = decltype(&hartmann::relaxation_matrix_offdiagonal);
Offdiagonal offdiagonal_kernel(LineByLineLineshape lineshape) {
  using enum LineByLineLineshape;
  switch (lineshape) {
    case VP_ECS_MAKAROV:  return makarov::relaxation_matrix_offdiagonal;
    case VP_ECS_HARTMANN: return hartmann::relaxation_matrix_offdiagonal;
    default:              ARTS_USER_ERROR("Unknown ECS line shape {}", lineshape)
  }
}

Numeric reduced_dipole(const QuantumIdentifier& qid, const band_data& bnd, const line& ln) {
  const auto& J = ln.qn.at(QuantumNumberType::J);
  using enum LineByLineLineshape;
  switch (bnd.lineshape) {
    case VP_ECS_MAKAROV:  return makarov::reduced_dipole(J.upper, J.lower, ln.qn.at(QuantumNumberType::N).upper);
    case VP_ECS_HARTMANN: {
      const auto& l = qid.state.at(QuantumNumberType::l2);
      return hartmann::reduced_dipole(J.upper, J.lower, l.upper, l.lower);
    }
    default: ARTS_USER_ERROR("Unknown ECS line shape {}", bnd.lineshape)
  }
}

Numeric closure_residual(ConstMatrixView W, ConstVectorView d) {
  Numeric result = 0;
  for (Index i = 0; i < W.nrows(); ++i) {
    Numeric residual = 0, scale = 0;
    for (Index j = 0; j < W.ncols(); ++j) {
      const Numeric term  = d[j] * W[i, j];
      residual           += term;
      scale              += std::abs(term);
    }
    if (scale > 0) result = std::max(result, std::abs(residual) / scale);
  }
  return result;
}
}  // namespace

void prepare_rotational_ladder(energy_data& energies, const int count, Numeric (*energy)(Rational)) {
  ARTS_USER_ERROR_IF(count < 0 or energy == nullptr, "Invalid ECS reference-rotor preparation")
  energies.rotational.resize(count);
  energies.rotational_minus_two.resize(count);
  for (Index L = 0; L < count; ++L) {
    energies.rotational[L]           = energy(Rational{L});
    energies.rotational_minus_two[L] = energy(Rational{L - 2});
  }
}

basis_data prepare_basis(const int                       count,
                         const energy_data&              energies,
                         const linemixing::species_data& collision,
                         const Numeric                   T0,
                         const SpeciesIsotope&           isot,
                         const SpeciesEnum               broadener,
                         const AtmPoint&                 atm) {
  ARTS_USER_ERROR_IF(count < 0 or energies.rotational.size() < static_cast<Size>(count) or
                         energies.rotational_minus_two.size() < static_cast<Size>(count),
                     "Incomplete ECS reference-rotor energies for {}",
                     isot)
  basis_data    out{.Q = Vector(count, 0.0), .Omega = Vector(count)};
  const Numeric mass = broadener == SpeciesEnum::Bath ? atm.mean_mass() : atm.mean_mass(broadener);
  for (Index L = 0; L < count; ++L) {
    out.Omega[L] =
        collision.Omega(atm.temperature, T0, mass, isot.mass, energies.rotational[L], energies.rotational_minus_two[L]);
  }
  // The power law is undefined at L=0. Keep positive channels from L=1,
  // including channels not selected by a particular species' angular sum.
  for (Index L = 1; L < count; ++L) {
    out.Q[L] = collision.Q(Rational{L}, atm.temperature, T0, energies.rotational[L]);
  }
  for (Index L = 0; L < count; ++L) {
    ARTS_USER_ERROR_IF(not std::isfinite(out.Omega[L]) or out.Omega[L] <= 0 or not std::isfinite(out.Q[L]),
                       "Invalid ECS basis rate or adiabaticity factor at L={} for {} and {}",
                       L,
                       isot,
                       broadener)
  }
  return out;
}

void apply_sum_rule(
    MatrixView W, ConstVectorView dipr, ConstVectorView e0, Numeric T, Tensor3View dW, ConstVectorView dT) {
  const Size  n  = dipr.size();
  const Index nq = dW.npages();
  ARTS_USER_ERROR_IF(W.nrows() != static_cast<Index>(n) or W.ncols() != static_cast<Index>(n) or e0.size() != n,
                     "Inconsistent ECS sum-rule dimensions")
  ARTS_USER_ERROR_IF((nq != 0 and (dW.nrows() != W.nrows() or dW.ncols() != W.ncols())) or
                         (not dT.empty() and dT.size() != static_cast<Size>(nq)),
                     "Inconsistent ECS sum-rule derivative dimensions")
  ARTS_USER_ERROR_IF(not std::isfinite(T) or T <= 0, "ECS sum-rule correction requires positive finite temperature")

  // The sequential correction retains the historical truncated-band closure.
  // In particular it cannot enforce the final row's sum rule. Overflow or
  // invalid rates leave a non-finite matrix behind rather than a silently
  // zeroed one; that is reported once for the propagation matrix instead of
  // being scanned for element by element here.
  Vector dsumlw(nq), dsumup(nq), dscale(nq);
  for (Size i = 0; i < n; ++i) {
    Numeric sumlw = 0.0;
    Numeric sumup = 0.0;
    dsumlw        = 0;
    dsumup        = 0;

    for (Size j = 0; j < n; ++j) {
      if (j > i) {
        sumlw += dipr[j] * W[i, j];
        for (Index q = 0; q < nq; ++q) dsumlw[q] += dipr[j] * dW[q, i, j];
      } else {
        sumup += dipr[j] * W[i, j];
        for (Index q = 0; q < nq; ++q) dsumup[q] += dipr[j] * dW[q, i, j];
      }
    }

    for (Index q = 0; q < nq; ++q) {
      dscale[q] = sumlw == 0 ? 0.0 : -(dsumup[q] + (-sumup / sumlw) * dsumlw[q]) / sumlw;
    }

    for (Size j = i + 1; j < n; ++j) {
      if (sumlw == 0) {
        W[i, j] = 0.0;
        W[j, i] = 0.0;
        for (Index q = 0; q < nq; ++q) {
          dW[q, i, j] = 0.0;
          dW[q, j, i] = 0.0;
        }
      } else {
        for (Index q = 0; q < nq; ++q) { dW[q, i, j] = dW[q, i, j] * (-sumup / sumlw) + W[i, j] * dscale[q]; }
        W[i, j] *= -sumup / sumlw;
        W[j, i]  = W[i, j] * std::exp((e0[i] - e0[j]) / (Constant::k * T));
        if (nq != 0) {
          const Numeric exponent = (e0[i] - e0[j]) / (Constant::k * T);
          const Numeric balance  = std::exp(exponent);
          for (Index q = 0; q < nq; ++q) {
            const Numeric dexponent = dT.empty() ? 0.0 : -exponent * dT[q] / T;
            dW[q, j, i]             = (dW[q, i, j] + W[i, j] * dexponent) * balance;
          }
        }
      }
    }
  }
}

void validate_band(const QuantumIdentifier& bnd_qid, const band_data& bnd) try {
  ARTS_USER_ERROR_IF(bnd.size() == 0, "Cannot use an empty ECS band")
  const auto& models = bnd.front().ls.single_models;
  ARTS_USER_ERROR_IF(models.empty(), "No broadening species in the ECS band")

  for (const auto& ln : bnd) {
    ARTS_USER_ERROR_IF(
        ln.ls.single_models.size() != models.size() or
            not stdr::all_of(models | stdv::keys, [&](auto spec) { return ln.ls.single_models.contains(spec); }),
        "All lines in an ECS band must have the same broadening species")
    ARTS_USER_ERROR_IF(not std::isfinite(ln.ls.T0) or ln.ls.T0 <= 0 or ln.ls.T0 != bnd.front().ls.T0,
                       "All lines in an ECS band must have the same positive reference temperature")
    ARTS_USER_ERROR_IF(not std::isfinite(ln.f0) or ln.f0 <= 0 or not std::isfinite(ln.a) or ln.a < 0 or
                           not std::isfinite(ln.gu) or ln.gu <= 0 or not std::isfinite(ln.e0),
                       "Invalid ECS line frequency, Einstein A, statistical weight, or lower-state energy")
  }

  using enum LineByLineLineshape;
  switch (bnd.lineshape) {
    case VP_ECS_MAKAROV:  makarov::validate_band(bnd_qid, bnd); break;
    case VP_ECS_HARTMANN: hartmann::validate_band(bnd_qid, bnd); break;
    default:              ARTS_USER_ERROR("Unknown ECS line shape {}", bnd.lineshape)
  }

  // The reduced dipoles are deliberately not evaluated here.  They need the
  // Wigner tables, which would make a pure catalogue check depend on
  // *WignerInit*, and a dipole that comes out non-finite is reported for the
  // propagation matrix anyway.
}
ARTS_METHOD_ERROR_CATCH

void ComputeData::adapt_multi(const QuantumIdentifier&        bnd_qid,
                              const band_data&                bnd,
                              const LinemixingSpeciesEcsData& rovib_data,
                              const AtmPoint&                 atm,
                              const bool                      presorted) {
  adapt(bnd_qid, bnd, rovib_data, atm, presorted, true, {});
}

void ComputeData::adapt_single(const QuantumIdentifier&        bnd_qid,
                               const band_data&                bnd,
                               const LinemixingSpeciesEcsData& rovib_data,
                               const AtmPoint&                 atm,
                               const bool                      presorted) {
  adapt(bnd_qid, bnd, rovib_data, atm, presorted, false, {});
}

void ComputeData::adapt_single(const QuantumIdentifier&        qid,
                               const band_data&                band,
                               const LinemixingSpeciesEcsData& data,
                               const AtmPoint&                 atm,
                               const Jacobian::Targets&        targets,
                               bool                            presorted) {
  adapt(qid, band, data, atm, presorted, false, targets);
}

void ComputeData::adapt(const QuantumIdentifier&        bnd_qid,
                        const band_data&                bnd,
                        const LinemixingSpeciesEcsData& rovib_data,
                        const AtmPoint&                 atm,
                        const bool                      presorted,
                        const bool                      per_broadener,
                        const Jacobian::Targets&        targets) try {
  const Size n = bnd.size();
  ARTS_USER_ERROR_IF(n == 0, "Cannot adapt an empty ECS band")
  ARTS_USER_ERROR_IF(
      not std::isfinite(atm.temperature) or atm.temperature <= 0 or not std::isfinite(atm.pressure) or atm.pressure < 0,
      "ECS requires positive finite temperature and nonnegative finite pressure")
  const auto& models          = bnd.front().ls.single_models;
  const Size  broadener_count = models.size();
  ARTS_USER_ERROR_IF(broadener_count == 0, "No broadening species in the ECS band")
  const Size m      = per_broadener ? broadener_count : 1;
  const auto kernel = offdiagonal_kernel(bnd.lineshape);
  ARTS_USER_ERROR_IF(bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN and bnd_qid.isot != "CO2-626"_isot,
                     "Hartmann ECS currently supports rotational energies only for CO2-626, got {}",
                     bnd_qid.isot)
  ARTS_USER_ERROR_IF(bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV and bnd_qid.isot != "O2-66"_isot,
                     "Makarov ECS currently supports only the O2-66 microwave band, got {}",
                     bnd_qid.isot)
  // Everything about the band itself is catalogue data that cannot change
  // between atmospheric points.  *validate_band* covers it, and is run by
  // *abs_bandsCheckEcs* rather than here.
  if (presorted) {
    ARTS_USER_ERROR_IF(sort.size() != n, "ECS presorting requires a previous adaptation of the same band size")
    auto indices = sort;
    stdr::sort(indices);
    for (Size i = 0; i < n; ++i) {
      ARTS_USER_ERROR_IF(indices[i] != static_cast<Index>(i), "Invalid ECS sorting permutation")
    }
  }

  const Size nt = targets.target_count();
  dW.resize(nt, n, n);
  dpop.resize(nt, n);
  ddip.resize(nt, n);
  dgd_fac.resize(nt);
  df.resize(nt);
  dW      = 0;
  dpop    = 0;
  ddip    = 0;
  dgd_fac = 0;
  df      = 0;
  Vector  dT(nt, 0.0);
  Tensor3 dWimag(nt, n, n);

  pop.resize(n);
  dip.resize(n);
  dipr.resize(n);
  sort.resize(n);
  Wimag.resize(n, n);
  vmrs.resize(m);
  eqv_strs.resize(m, n);
  eqv_vals.resize(m, n);
  Ws.resize(m, n, n);
  sum_rule_residual.resize(broadener_count);
  Ws                = 0;
  eqv_strs          = 0;
  eqv_vals          = 0;
  sum_rule_residual = 0;

  gd_fac = std::sqrt(Constant::doppler_broadening_const_squared * atm.temperature / bnd_qid.isot.mass);
  ARTS_USER_ERROR_IF(not std::isfinite(gd_fac) or gd_fac <= 0, "Invalid ECS Doppler width factor")
  const Numeric QT = PartitionFunctions::Q(atm.temperature, bnd_qid.isot);
  ARTS_USER_ERROR_IF(not std::isfinite(QT) or QT <= 0, "Invalid ECS partition function: {}", QT)
  for (Size i = 0; i < n; ++i) {
    const auto& ln = bnd.lines[i];
    pop[i]         = ln.gu * std::exp(-ln.e0 / (Constant::k * atm.temperature)) / QT;
    dipr[i]        = reduced_dipole(bnd_qid, bnd, ln);
    dip[i] = std::copysign(0.5 * Constant::c * std::sqrt(ln.a / (Math::pow3(ln.f0) * Constant::two_pi)), dipr[i]);
  }

  if (not presorted) {
    stdr::iota(sort, 0);
    stdr::sort(stdv::zip(sort, pop, dip, dipr), stdr::greater(), [&](const auto& v) {
      const auto& [i, pop_i, dip_i, dipr_i] = v;
      return bnd.lines[i].f0 * pop_i * dip_i * dip_i;
    });
  } else {
    const auto reorder = [this](const Vector& vec) {
      Vector out(vec.size());
      for (Size i = 0; i < sort.size(); ++i) { out[i] = vec[sort[i]]; }
      return out;
    };
    pop  = reorder(pop);
    dip  = reorder(dip);
    dipr = reorder(dipr);
  }

  const Numeric dQT = nt ? PartitionFunctions::dQdT(atm.temperature, bnd_qid.isot) : 0;
  for (const auto& target : targets.atm) {
    ARTS_USER_ERROR_IF(target.target_pos >= nt, "Invalid ECS Jacobian target position")
    const auto* key = std::get_if<AtmKey>(&target.type);
    if (not key) continue;
    switch (*key) {
      case AtmKey::t:
        dT[target.target_pos]      = 1;
        dgd_fac[target.target_pos] = gd_fac / (2 * atm.temperature);
        for (Size k = 0; k < n; ++k) {
          const auto& ln             = bnd.lines[sort[k]];
          dpop[target.target_pos, k] = pop[k] * (ln.e0 / (Constant::k * Math::pow2(atm.temperature)) - dQT / QT);
        }
        break;
      case AtmKey::wind_u:
      case AtmKey::wind_v:
      case AtmKey::wind_w:
        // Like the other LBL profiles, return the frequency derivative here;
        // propagation applies the wind projection and Doppler conversion.
        df[target.target_pos] = 1;
        break;
      default: break;
    }
  }
  for (const auto& target : targets.line) {
    ARTS_USER_ERROR_IF(target.target_pos >= nt, "Invalid ECS Jacobian target position")
    const auto& key = target.type;
    if (key.band != bnd_qid) continue;
    for (Size k = 0; k < n; ++k) {
      if (key.line != static_cast<Size>(sort[k])) continue;
      const auto& ln = bnd.lines[sort[k]];
      switch (key.var) {
        case LineByLineVariable::f0:
          ddip[target.target_pos, k] = -1.5 * dip[k] / ln.f0;
          dW[target.target_pos, k, k].real(1);
          break;
        case LineByLineVariable::e0: dpop[target.target_pos, k] = -pop[k] / (Constant::k * atm.temperature); break;
        case LineByLineVariable::a:
          ARTS_USER_ERROR_IF(ln.a == 0, "ECS Einstein-A derivatives require a positive Einstein A")
          ddip[target.target_pos, k] = dip[k] / (2 * ln.a);
          break;
        case LineByLineVariable::unused: break;
      }
    }
  }

  // Prepare kernel inputs once in matrix order. The kernels never access the
  // catalogue or its permutation; optical populations above retain catalogue e0.
  rotational_lines.resize(n);
  for (Size i = 0; i < n; ++i) {
    const auto& ln      = bnd.lines[sort[i]];
    const auto& J       = ln.qn.at(QuantumNumberType::J);
    rotational_lines[i] = {.Ju = J.upper, .Jl = J.lower};
    if (bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV) {
      const auto& N          = ln.qn.at(QuantumNumberType::N);
      rotational_lines[i].Nu = N.upper;
      rotational_lines[i].Nl = N.lower;
    }
  }
  if (bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN) {
    hartmann::prepare_energies(energies, bnd_qid, rotational_lines);
  } else {
    makarov::prepare_energies(energies, bnd_qid, rotational_lines);
  }

  Vector fractions(broadener_count);
  get_vmrs(fractions, models, atm);
  if (per_broadener)
    vmrs = fractions;
  else
    vmrs = 1;

  Matrix dfractions(nt, broadener_count);
  get_dvmrs(dfractions, fractions, models, atm, targets);
  Matrix dQ(nt, energies.rotational.size()), dOmega(nt, energies.rotational.size());
  Size   i = 0;
  for (auto spec : models | stdv::keys) {
    const Size    page   = per_broadener ? i : 0;
    const Numeric weight = per_broadener ? 1 : fractions[i];
    bool          active = weight != 0;
    for (Size t = 0; t < nt; ++t) active |= dfractions[t, i] != 0;
    if (not active) {
      ++i;
      continue;
    }
    const auto data = rovib_data.find(spec);
    ARTS_USER_ERROR_IF(data == rovib_data.end(), "No ECS collision data for species {}", spec)
    Wimag  = 0;
    dWimag = 0;
    for (Size k = 0; k < n; ++k) {
      const auto&   ln    = bnd.lines[sort[k]];
      const auto&   model = ln.ls.single_models.at(spec);
      const Numeric width = model.G0(ln.ls.T0, atm.temperature, atm.pressure);
      const Numeric shift = model.D0(ln.ls.T0, atm.temperature, atm.pressure);
      // A negative width is not detectable downstream: it stays finite and
      // turns the profile into gain.  Non-finite values are left to the
      // propagation-matrix check.
      ARTS_USER_ERROR_IF(width < 0, "Negative ECS pressure width for species {}", spec)
      Wimag[k, k]               = width;
      real_val(Ws[page, k, k]) += weight * shift;
      for (Size t = 0; t < nt; ++t) real_val(dW[t, k, k]) += dfractions[t, i] * shift;
      for (const auto& target : targets.atm) {
        const auto* key = std::get_if<AtmKey>(&target.type);
        if (not key) continue;
        Numeric dw = 0, ds = 0;
        switch (*key) {
          case AtmKey::t:
            dw = model.dG0_dT(ln.ls.T0, atm.temperature, atm.pressure);
            ds = model.dD0_dT(ln.ls.T0, atm.temperature, atm.pressure);
            break;
          case AtmKey::p:
            dw = model.dG0_dP(ln.ls.T0, atm.temperature, atm.pressure);
            ds = model.dD0_dP(ln.ls.T0, atm.temperature, atm.pressure);
            break;
          default: break;
        }
        dWimag[target.target_pos, k, k]        = dw;
        real_val(dW[target.target_pos, k, k]) += weight * ds;
      }
      for (const auto& target : targets.line) {
        const auto& key = target.type;
        if (key.band != bnd_qid or key.line != static_cast<Size>(sort[k]) or key.spec != spec) continue;
        switch (key.ls_var) {
          case LineShapeModelVariable::G0:
            dWimag[target.target_pos, k, k] = model.dG0_dX(ln.ls.T0, atm.temperature, atm.pressure, key.ls_coeff);
            break;
          case LineShapeModelVariable::D0:
            real_val(dW[target.target_pos, k, k]) +=
                weight * model.dD0_dX(ln.ls.T0, atm.temperature, atm.pressure, key.ls_coeff);
            break;
          default: break;
        }
      }
    }
    prepare_basis_jacobian(dQ, dOmega, energies, data->second, bnd.front().ls.T0, bnd_qid.isot, spec, atm, targets);
    kernel(Wimag,
           bnd_qid,
           rotational_lines,
           bnd.front().ls.T0,
           spec,
           data->second,
           dipr,
           energies,
           atm,
           dWimag,
           dT,
           dQ,
           dOmega);
    sum_rule_residual[i] = closure_residual(Wimag, dipr);
    for (Size r = 0; r < n; ++r) {
      for (Size c = 0; c < n; ++c) {
        imag_val(Ws[page, r, c]) += weight * Wimag[r, c];
        for (Size t = 0; t < nt; ++t)
          imag_val(dW[t, r, c]) += weight * dWimag[t, r, c] + dfractions[t, i] * Wimag[r, c];
      }
    }
    ++i;
  }
  for (Size i = 0; i < n; ++i) { Ws[joker, i, i] += bnd.lines[sort[i]].f0; }
}
ARTS_METHOD_ERROR_CATCH

void calculate(PropmatVectorView               pm_,
               PropmatMatrixView               dpm_,
               ComputeData&                    com_data,
               const ConstVectorView           f_grid_,
               const Range&                    f_range,
               const Jacobian::Targets&        jac_targets,
               const QuantumIdentifier&        bnd_qid,
               const band_data&                bnd,
               const LinemixingSpeciesEcsData& rovib_data,
               const AtmPoint&                 atm,
               const ZeemanPolarization        pol,
               const bool                      no_negative_absorption) try {
  if (pol != ZeemanPolarization::no) {
    ARTS_USER_ERROR_IF(stdr::any_of(
                           bnd, [](auto& zee) { return zee.on; }, &line::z),
                       "Zeeman effect and ECS in combination is not yet possible.")
    return;
  }

  PropmatVectorView     pm     = pm_[f_range];
  const ConstVectorView f_grid = f_grid_[f_range];

  ARTS_USER_ERROR_IF(dpm_.nrows() != static_cast<Index>(jac_targets.target_count()) or
                         dpm_.ncols() != static_cast<Index>(f_grid_.size()),
                     "Inconsistent ECS propagation Jacobian dimensions")

  if (bnd.size() == 0) return;

  com_data.adapt_single(bnd_qid, bnd, rovib_data, atm, jac_targets);

  com_data.core_calc(f_grid);

  const Size    nt        = jac_targets.target_count();
  const Numeric abundance = atm[bnd_qid.isot.spec], ratio = atm[bnd_qid.isot];
  const Numeric scale = Constant::inv_sqrt_pi * abundance * ratio;
  const Numeric N     = number_density(atm.pressure, atm.temperature);
  Matrix        dscale(nt, f_grid.size(), 0.0);
  for (const auto& target : jac_targets.atm) {
    if (const auto* key = std::get_if<AtmKey>(&target.type)) {
      for (Size i = 0; i < f_grid.size(); ++i) {
        const Numeric f = f_grid[i], r = Constant::h * f / (Constant::k * atm.temperature);
        const Numeric e = std::expm1(-r);
        switch (*key) {
          case AtmKey::t:      dscale[target.target_pos, i] = -scale * N * f * (r * (e + 1) - e) / atm.temperature; break;
          case AtmKey::p:      dscale[target.target_pos, i] = -scale * f * e / (Constant::k * atm.temperature); break;
          case AtmKey::wind_u:
          case AtmKey::wind_v:
          case AtmKey::wind_w: dscale[target.target_pos, i] = scale * N * (r * (e + 1) - e); break;
          default:             break;
        }
      }
    } else if (const auto* species = std::get_if<SpeciesEnum>(&target.type)) {
      if (*species == bnd_qid.isot.spec)
        for (Size i = 0; i < f_grid.size(); ++i)
          dscale[target.target_pos, i] = Constant::inv_sqrt_pi * ratio * com_data.scl[i];
    } else if (const auto* isot = std::get_if<SpeciesIsotope>(&target.type)) {
      if (*isot == bnd_qid.isot)
        for (Size i = 0; i < f_grid.size(); ++i)
          dscale[target.target_pos, i] = Constant::inv_sqrt_pi * abundance * com_data.scl[i];
    }
  }
  for (Size i = 0; i < f_grid.size(); ++i) {
    const Complex F = scale * com_data.scl[i] * com_data.shape[i];
    if (no_negative_absorption and F.real() < 0) continue;
    pm[i] += zeeman::scale(com_data.npm, F);
    for (Size t = 0; t < nt; ++t) {
      const Complex dF     = dscale[t, i] * com_data.shape[i] + scale * com_data.scl[i] * com_data.dshape[t, i];
      dpm_[t, f_range][i] += zeeman::scale(com_data.npm, dF);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void equivalent_values(ComplexTensor3View              eqv_str,
                       ComplexTensor3View              eqv_val,
                       ComputeData&                    com_data,
                       const QuantumIdentifier&        bnd_qid,
                       const band_data&                bnd,
                       const LinemixingSpeciesEcsData& rovib_data,
                       const AtmPoint&                 atm,
                       const Vector&                   T) try {
  const auto k = eqv_str.npages();
  const auto m = eqv_str.ncols();

  ARTS_USER_ERROR_IF(eqv_str.shape() != eqv_val.shape(), "eqv_str and eqv_val must have the same shape.")
  ARTS_USER_ERROR_IF(T.size() != static_cast<Size>(k), "T must have the same size as eqv_str pages.")
  ARTS_USER_ERROR_IF(bnd.size() != static_cast<Size>(m), "bnd must have the same size as eqv_str cols.")

  if (bnd.size() == 0) return;
  ARTS_USER_ERROR_IF(eqv_str.nrows() != static_cast<Index>(bnd.front().ls.single_models.size()),
                     "eqv_str rows must match the number of ECS broadening species")

  com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm, false);

  std::string err{};
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(com_data)
  for (Index i = 0; i < k; ++i) {
    try {
      AtmPoint atm_copy    = atm;
      atm_copy.temperature = T[i];
      com_data.adapt_multi(bnd_qid, bnd, rovib_data, atm_copy, true);
      com_data.core_calc_eqv();
      eqv_str[i] = com_data.eqv_strs;
      eqv_val[i] = com_data.eqv_vals;
    } catch (std::exception& e) {
#pragma omp critical
      err += std::format("{}\n", e.what());
    }
  }

  if (not err.empty()) throw std::runtime_error(err);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::voigt::ecs
