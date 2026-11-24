#include <workspace.h>

namespace {
void wind_shift(VectorView                  freq_grid,
                Vector3&                    freq_wind_shift_jac,
                const AtmPoint&             atm_point,
                const PropagationPathPoint& ray_point) {
  constexpr Numeric c = Constant::speed_of_light;

  // The line-of-sight of a path point is the direction of propagation, towards the observer.  Gas moving
  // along it approaches the observer, so a line is seen at a higher frequency than in the frame of the gas.
  // The gas moves with the wind and the rotation of the planet.
  const Vector3 wind  = atm_point.effective_wind();
  const auto [za, aa] = ray_point.los;

  // The direction of propagation in the local east, north and up, the frame of the wind
  const Vector3 los{Conversion::sind(za) * Conversion::sind(aa),
                    Conversion::sind(za) * Conversion::cosd(aa),
                    Conversion::cosd(za)};

  const Numeric fac = 1.0 - dot(wind, los) / c;

  if (fac <= 0)
    arts::user_error(R"(Negative frequency scaling factor: {}

atm_point.wind:  {:B,}
atm_point.rot_u: {}
ray_point.los:   {:B,}
)",
                     fac,
                     atm_point.wind,
                     atm_point.rot_u,
                     ray_point.los);

  //! Zero shift if nan
  if (std::isnan(fac)) {
    freq_wind_shift_jac = {0, 0, 0};
    return;
  }

  stdr::transform(freq_grid, freq_grid.begin(), [fac](const Numeric& f) { return fac * f; });

  // fac is linear in the wind, so its derivative is the same for any wind.  The rotation is not a variable
  // of the wind, so the derivatives with respect to wind_u, wind_v and wind_w are those of fac
  freq_wind_shift_jac = los;
  freq_wind_shift_jac /= -c * fac;
}
}  // namespace

void freq_gridWindShift(AscendingGrid&              freq_grid,
                        Vector3&                    freq_wind_shift_jac,
                        const AtmPoint&             atm_point,
                        const PropagationPathPoint& ray_point) {
  ARTS_TIME_REPORT

  Vector tmp = std::move(freq_grid).rvec();
  wind_shift(tmp, freq_wind_shift_jac, atm_point, ray_point);
  freq_grid = std::move(tmp);
}

void freqWindShift(Numeric&                    frequency,
                   Vector3&                    freq_wind_shift_jac,
                   const AtmPoint&             atm_point,
                   const PropagationPathPoint& ray_point) {
  ARTS_TIME_REPORT

  wind_shift(VectorView{frequency}, freq_wind_shift_jac, atm_point, ray_point);
}

void spectral_propmatWindFix(SpectralPropmat&       spectral_propmat,
                             const AscendingGrid&   freq_grid,
                             const JacobianTargets& jac_targets,
                             const Vector3&         freq_wind_shift_jac) {
  ARTS_TIME_REPORT

  PropmatMatrix& spectral_propmat_jac          = spectral_propmat.dK;
  StokvecMatrix& source_vector_nonlte_jacobian = spectral_propmat.dnlte;

  using enum AtmKey;

  constexpr auto is_wind_u = Cmp::eq<wind_u>();
  constexpr auto is_wind_v = Cmp::eq<wind_v>();
  constexpr auto is_wind_w = Cmp::eq<wind_w>();

  const auto& atm = jac_targets.atm;

  const auto [df_du, df_dv, df_dw] = freq_wind_shift_jac;

  if (auto ptr = stdr::find_if(atm, is_wind_u, &Jacobian::AtmTarget::type); ptr != atm.end()) {
    const Index i = ptr->target_pos;

    std::transform(freq_grid.begin(),
                   freq_grid.end(),
                   spectral_propmat_jac[i].begin(),
                   spectral_propmat_jac[i].begin(),
                   [df_du](const Numeric f, const Propmat& x) { return x * f * df_du; });
    std::transform(freq_grid.begin(),
                   freq_grid.end(),
                   source_vector_nonlte_jacobian[i].begin(),
                   source_vector_nonlte_jacobian[i].begin(),
                   [df_du](const Numeric f, const Stokvec& x) { return x * f * df_du; });
  }

  if (auto ptr = stdr::find_if(atm, is_wind_v, &Jacobian::AtmTarget::type); ptr != atm.end()) {
    const Index i = ptr->target_pos;

    std::transform(freq_grid.begin(),
                   freq_grid.end(),
                   spectral_propmat_jac[i].begin(),
                   spectral_propmat_jac[i].begin(),
                   [df_dv](const Numeric f, const Propmat& x) { return x * f * df_dv; });
    std::transform(freq_grid.begin(),
                   freq_grid.end(),
                   source_vector_nonlte_jacobian[i].begin(),
                   source_vector_nonlte_jacobian[i].begin(),
                   [df_dv](const Numeric f, const Stokvec& x) { return x * f * df_dv; });
  }

  if (auto ptr = stdr::find_if(atm, is_wind_w, &Jacobian::AtmTarget::type); ptr != atm.end()) {
    const Index i = ptr->target_pos;

    std::transform(freq_grid.begin(),
                   freq_grid.end(),
                   spectral_propmat_jac[i].begin(),
                   spectral_propmat_jac[i].begin(),
                   [df_dw](const Numeric f, const Propmat& x) { return x * f * df_dw; });
    std::transform(freq_grid.begin(),
                   freq_grid.end(),
                   source_vector_nonlte_jacobian[i].begin(),
                   source_vector_nonlte_jacobian[i].begin(),
                   [df_dw](const Numeric f, const Stokvec& x) { return x * f * df_dw; });
  }
}
