#include <atm.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>

#include "arts_constants.h"

void frequency_gridWindShift(AscendingGrid& frequency_grid,
                             Vector3& frequency_grid_wind_shift_jacobian,
                             const AtmPoint& atmospheric_point,
                             const PropagationPathPoint& ray_path_point) {
  constexpr Numeric c = Constant::speed_of_light;

  const auto& [u, v, w] = atmospheric_point.wind;
  const auto [za, aa]   = path::mirror(ray_path_point.los);

  const Numeric u2v2 = u * u + v * v;
  const Numeric w2   = w * w;
  const Numeric f2   = u2v2 + w2;
  const Numeric f    = std::sqrt(f2);
  const Numeric za_f = f == w ? (0.5 * Constant::pi) : std::acos(w / f);
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

  ARTS_USER_ERROR_IF(
      fac < 0 or nonstd::isnan(fac), "Bad frequency scaling factor: {}", fac)

  // shift the frequency grid
  std::transform(frequency_grid.begin(),
                 frequency_grid.end(),
                 frequency_grid.unsafe_begin(),
                 [fac](const Numeric& f) { return fac * f; });

  const Numeric df_du    = f == u ? 1 : u / f;
  const Numeric dczaf_du = f == u ? 0 : (-w * df_du / f2);
  const Numeric dszaf_du = f == w ? 0 : (w2 * df_du / scl);
  const Numeric dcaa_du  = u2v2 == 0 ? 0 : (v * saa / u2v2);
  const Numeric ddp_du =
      czap * dczaf_du + szap * caa * dszaf_du + szap * szaf * dcaa_du;
  frequency_grid_wind_shift_jacobian[0] = -(dp * df_du + f * ddp_du) / c;

  const Numeric df_dv    = f == v ? 1 : v / f;
  const Numeric dczaf_dv = f == v ? 0 : (-w * df_dv / f2);
  const Numeric dszaf_dv = f == v ? 1 : (w2 * df_dv / scl);
  const Numeric dcaa_dv  = u2v2 == 0 ? 0 : (-u * saa / u2v2);
  const Numeric ddp_dv =
      czap * dczaf_dv + szap * caa * dszaf_dv + szap * szaf * dcaa_dv;
  frequency_grid_wind_shift_jacobian[1] = -(dp * df_dv + f * ddp_dv) / c;

  const Numeric df_dw    = f == w ? 1 : w / f;
  const Numeric dczaf_dw = f == w ? 0 : (-w * df_dw / f2 + 1.0 / f);
  const Numeric dszaf_dw = f == w ? 0 : ((w2 * df_dw - f * w) / scl);
  const Numeric ddp_dw   = czap * dczaf_dw + szap * caa * dszaf_dw;
  frequency_grid_wind_shift_jacobian[2] = -(dp * df_dw + f * ddp_dw) / c;
}

void propagation_matrix_jacobianWindFix(
    PropmatMatrix& propagation_matrix_jacobian,
    StokvecMatrix& source_vector_nonlte_jacobian,
    const JacobianTargets& jacobian_targets,
    const Vector3& frequency_grid_wind_shift_jacobian) {
  using enum AtmKey;

  const auto& atm = jacobian_targets.atm();

  const auto [df_du, df_dv, df_dw] = frequency_grid_wind_shift_jacobian;

  if (auto ptr = std::ranges::find_if(
          atm, Cmp::eq(wind_u), &Jacobian::AtmTarget::type);
      ptr != atm.end()) {
    const Index i = ptr->target_pos;

    propagation_matrix_jacobian[i]   *= df_du;
    source_vector_nonlte_jacobian[i] *= df_du;
  }

  if (auto ptr = std::ranges::find_if(
          atm, Cmp::eq(wind_v), &Jacobian::AtmTarget::type);
      ptr != atm.end()) {
    const Index i = ptr->target_pos;

    propagation_matrix_jacobian[i]   *= df_dv;
    source_vector_nonlte_jacobian[i] *= df_dv;
  }

  if (auto ptr = std::ranges::find_if(
          atm, Cmp::eq(wind_w), &Jacobian::AtmTarget::type);
      ptr != atm.end()) {
    const Index i = ptr->target_pos;

    propagation_matrix_jacobian[i]   *= df_dw;
    source_vector_nonlte_jacobian[i] *= df_dw;
  }
}
