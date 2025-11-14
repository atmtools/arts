#include <workspace.h>

#include "matpack_mdspan_helpers_grid_t.h"

namespace {
void wind_shift(VectorView freq_grid,
                Vector3& freq_wind_shift_jac,
                const AtmPoint& atm_point,
                const PropagationPathPoint& ray_path_point) {
  constexpr Numeric c = Constant::speed_of_light;

  const auto& [u, v, w] = atm_point.wind;
  const auto [za, aa]   = path::mirror(ray_path_point.los);

  const Numeric u2v2 = u * u + v * v;
  const Numeric w2   = w * w;
  const Numeric f2   = u2v2 + w2;
  const Numeric f    = std::sqrt(f2);
  const Numeric za_f = f == w ? 0.0 : std::acos(w / f);
  const Numeric aa_f = std::atan2(u, v);
  const Numeric za_p = Conversion::deg2rad(za);
  const Numeric aa_p = Conversion::deg2rad(aa);
  const Numeric scl  = f2 * std::sqrt(f2 - w2);
  const Numeric czaf = std::cos(za_f);
  const Numeric szaf = std::sin(za_f);
  const Numeric czap = std::cos(za_p);
  const Numeric szap = std::sin(za_p);
  const Numeric caa  = std::cos(aa_f - aa_p);
  const Numeric saa  = std::sin(aa_p - aa_f);
  const Numeric dp   = czaf * czap + szaf * szap * caa;
  const Numeric fac  = 1.0 - (f * dp) / c;

  ARTS_USER_ERROR_IF(fac <= 0,
                     R"(Negative frequency scaling factor: {}

atm_point.wind: {:B,}
ray_path_point.los:     {:B,}
)",
                     fac,
                     atm_point.wind,
                     ray_path_point.los);

  //! Zero shift if nan
  if (std::isnan(fac)) {
    freq_wind_shift_jac = {0, 0, 0};
    return;
  }

  // shift the frequency grid
  {
    stdr::transform(freq_grid, freq_grid.begin(), [fac](const Numeric& f) {
      return fac * f;
    });
  }

  {
    const Numeric df_du    = (f == 0) ? 1.0 : u / f;
    const Numeric dczaf_du = (f2 == 0) ? 0.0 : (-w * df_du / f2);
    const Numeric dszaf_du = (f2 == w2) ? 0.0 : (w2 * df_du / scl);
    const Numeric dcaa_du  = (u2v2 == 0) ? 0.0 : (v * saa / u2v2);
    const Numeric ddp_du =
        czap * dczaf_du + szap * caa * dszaf_du + szap * szaf * dcaa_du;
    freq_wind_shift_jac[0] = -(dp * df_du + f * ddp_du) / c;
  }

  {
    const Numeric df_dv    = (f == 0) ? 1.0 : v / f;
    const Numeric dczaf_dv = (f2 == 0) ? 0.0 : (-w * df_dv / f2);
    const Numeric dszaf_dv = (f2 == w2) ? 0.0 : (w2 * df_dv / scl);
    const Numeric dcaa_dv  = (u2v2 == 0) ? 0.0 : (-u * saa / u2v2);
    const Numeric ddp_dv =
        czap * dczaf_dv + szap * caa * dszaf_dv + szap * szaf * dcaa_dv;
    freq_wind_shift_jac[1] = -(dp * df_dv + f * ddp_dv) / c;
  }

  {
    const Numeric df_dw    = (f == 0) ? 1.0 : w / f;
    const Numeric dczaf_dw = (f2 == 0) ? 0.0 : (-w * df_dw / f2 + 1.0 / f);
    const Numeric dszaf_dw = (scl == 0) ? 0.0 : ((w2 * df_dw - f * w) / scl);
    const Numeric ddp_dw   = czap * dczaf_dw + szap * caa * dszaf_dw;
    freq_wind_shift_jac[2] = -(dp * df_dw + f * ddp_dw) / c;
  }

  freq_wind_shift_jac /= fac;
}
}  // namespace

void freq_gridWindShift(AscendingGrid& freq_grid,
                        Vector3& freq_wind_shift_jac,
                        const AtmPoint& atm_point,
                        const PropagationPathPoint& ray_path_point) {
  ARTS_TIME_REPORT

  Vector tmp = std::move(freq_grid).rvec();
  wind_shift(tmp, freq_wind_shift_jac, atm_point, ray_path_point);
  freq_grid = std::move(tmp);
}

void freqWindShift(Numeric& frequency,
                   Vector3& freq_wind_shift_jac,
                   const AtmPoint& atm_point,
                   const PropagationPathPoint& ray_path_point) {
  ARTS_TIME_REPORT

  wind_shift(
      VectorView{frequency}, freq_wind_shift_jac, atm_point, ray_path_point);
}

void spectral_propmat_jacWindFix(PropmatMatrix& spectral_propmat_jac,
                                 StokvecMatrix& source_vector_nonlte_jacobian,
                                 const AscendingGrid& freq_grid,
                                 const JacobianTargets& jac_targets,
                                 const Vector3& freq_wind_shift_jac) {
  ARTS_TIME_REPORT

  using enum AtmKey;

  const auto& atm = jac_targets.atm;

  const auto [df_du, df_dv, df_dw] = freq_wind_shift_jac;

  if (auto ptr = std::ranges::find_if(
          atm, Cmp::eq(wind_u), &Jacobian::AtmTarget::type);
      ptr != atm.end()) {
    const Index i = ptr->target_pos;

    std::transform(
        freq_grid.begin(),
        freq_grid.end(),
        spectral_propmat_jac[i].begin(),
        spectral_propmat_jac[i].begin(),
        [df_du](const Numeric f, const Propmat& x) { return x * f * df_du; });
    std::transform(
        freq_grid.begin(),
        freq_grid.end(),
        source_vector_nonlte_jacobian[i].begin(),
        source_vector_nonlte_jacobian[i].begin(),
        [df_du](const Numeric f, const Stokvec& x) { return x * f * df_du; });
  }

  if (auto ptr = std::ranges::find_if(
          atm, Cmp::eq(wind_v), &Jacobian::AtmTarget::type);
      ptr != atm.end()) {
    const Index i = ptr->target_pos;

    std::transform(
        freq_grid.begin(),
        freq_grid.end(),
        spectral_propmat_jac[i].begin(),
        spectral_propmat_jac[i].begin(),
        [df_dv](const Numeric f, const Propmat& x) { return x * f * df_dv; });
    std::transform(
        freq_grid.begin(),
        freq_grid.end(),
        source_vector_nonlte_jacobian[i].begin(),
        source_vector_nonlte_jacobian[i].begin(),
        [df_dv](const Numeric f, const Stokvec& x) { return x * f * df_dv; });
  }

  if (auto ptr = std::ranges::find_if(
          atm, Cmp::eq(wind_w), &Jacobian::AtmTarget::type);
      ptr != atm.end()) {
    const Index i = ptr->target_pos;

    std::transform(
        freq_grid.begin(),
        freq_grid.end(),
        spectral_propmat_jac[i].begin(),
        spectral_propmat_jac[i].begin(),
        [df_dw](const Numeric f, const Propmat& x) { return x * f * df_dw; });
    std::transform(
        freq_grid.begin(),
        freq_grid.end(),
        source_vector_nonlte_jacobian[i].begin(),
        source_vector_nonlte_jacobian[i].begin(),
        [df_dw](const Numeric f, const Stokvec& x) { return x * f * df_dw; });
  }
}
