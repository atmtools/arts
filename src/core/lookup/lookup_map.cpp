#include "lookup_map.h"

#include <arts_omp.h>
#include <jacobian.h>

#include "matpack_view.h"

namespace lookup {
table::table(const SpeciesEnum& species,
             std::shared_ptr<const AscendingGrid> f_grid_,
             std::shared_ptr<const ArrayOfAtmPoint> atm_,
             const AbsorptionBands& absorption_bands,
             const LinemixingEcsData& ecs_data,
             AscendingGrid t_pert_,
             AscendingGrid water_pert_)
    : f_grid(std::move(f_grid_)),
      atmref(std::move(atm_)),
      log_p_grid(atmref->begin(),
                 atmref->end(),
                 [](const AtmPoint& x) { return std::log(x.pressure); }),
      t_pert(std::move(t_pert_)),
      water_pert(std::move(water_pert_)),
      water_atmref(atmref->size()),
      t_atmref(atmref->size()),
      xsec(std::max<Index>(1, t_pert.size()),
           std::max<Index>(1, water_pert.size()),
           f_grid->size(),
           atmref->size()) {
  const bool do_water = not water_pert.empty();

  std::transform(
      atmref->begin(),
      atmref->end(),
      water_atmref.begin(),
      [do_water](const auto& x) { return do_water ? x["H2O"_spec] : 0; });

  std::transform(
      atmref->begin(), atmref->end(), t_atmref.begin(), [](const auto& x) {
        return x.temperature;
      });

  String error;
  PropmatVector pm(f_grid->size());
  StokvecVector sv(f_grid->size());
  PropmatMatrix dpm(0, f_grid->size());
  StokvecMatrix dsv(0, f_grid->size());
  const JacobianTargets jacobian_targets = {};
  const Vector2 los                      = {180, 0};
  const bool no_negative_absorption      = true;

  const AscendingGrid water_vmr_local(water_pert.empty() ? AscendingGrid{1.0}
                                                         : water_pert);
  const AscendingGrid t_pert_local(t_pert.empty() ? AscendingGrid{0.0}
                                                  : t_pert);

#pragma omp parallel for collapse(3) if (not arts_omp_in_parallel()) \
    firstprivate(pm, sv, dpm, dsv)
  for (Index it = 0; it < t_pert_local.size(); ++it) {
    for (Index iw = 0; iw < water_vmr_local.size(); ++iw) {
      for (Size ip = 0; ip < atmref->size(); ++ip) {
        try {
          AtmPoint atm_point     = atmref->operator[](ip);
          atm_point.temperature += t_pert_local[it];
          if (do_water) atm_point["H2O"_spec] *= water_vmr_local[iw];

          pm = 0.0;
          lbl::calculate(pm,
                         sv,
                         dpm,
                         dsv,
                         *f_grid,
                         jacobian_targets,
                         species,
                         absorption_bands,
                         ecs_data,
                         atm_point,
                         los,
                         no_negative_absorption);

          const Numeric inv_nd = 1.0 / atm_point.number_density(species);
          for (Index ifreq = 0; ifreq < f_grid->size(); ++ifreq) {
            xsec(it, iw, ifreq, ip) = pm[ifreq].A() * inv_nd;
          }
        } catch (std::runtime_error& e) {
#pragma omp critical
          error += std::format("ERROR:\n{}\n\n", e.what());
        }
      }
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
}

void table::absorption(ExhaustiveVectorView absorption,
                       const SpeciesEnum& species,
                       const Index& p_interp_order,
                       const Index& t_interp_order,
                       const Index& water_interp_order,
                       const Index& f_interp_order,
                       const AtmPoint& atm_point,
                       const AscendingGrid& frequency_grid,
                       const Numeric& extpolfac) const {
  if (xsec.empty()) return;

  ARTS_USER_ERROR_IF(
      log_p_grid.nelem() != static_cast<Index>(atmref->size()),
      "The lookup table internal variable log_p_grid is incorrect.\n")

  // Frequency grid positions
  const Vector& frequency_grid_v(frequency_grid);
  const Vector& f_grid_v(*f_grid);
  const auto flag =
      my_interp::lagrange_interpolation_list<LagrangeInterpolation>(
          frequency_grid_v, f_grid_v, f_interp_order, extpolfac);

  // Pressure grid positions
  const auto plog_local = std::log(atm_point.pressure);
  const Vector& plog_v(log_p_grid);
  LagrangeInterpolation::check(plog_v, p_interp_order, plog_local, extpolfac);
  const Array<LagrangeInterpolation> plag{
      LagrangeInterpolation(0, plog_local, plog_v, p_interp_order)};

  const auto water_lag = [&]() {
    const auto x = atm_point["H2O"_spec] / interp(water_atmref, plag[0]);
    const Vector& xi(water_pert);
    LagrangeInterpolation::check(xi, water_interp_order, x, extpolfac);
    const Array<LagrangeInterpolation> wlag{
        LagrangeInterpolation(0, x, xi, water_interp_order)};
    return wlag;
  };

  const auto temperature_lag = [&]() {
    const auto x = atm_point.temperature - interp(t_atmref, plag[0]);
    const Vector& xi(t_pert);
    LagrangeInterpolation::check(xi, t_interp_order, x, extpolfac);
    const Array<LagrangeInterpolation> tlag{
        LagrangeInterpolation(0, x, xi, t_interp_order)};
    return tlag;
  };

  Vector xsec_local(frequency_grid.size());

  // Optional grid positions by switching
  if (water_pert.empty() and t_pert.empty()) {
    xsec_local =
        reinterp(xsec[0][0], flag, plag).reshape(frequency_grid.size());
  } else if (t_pert.empty()) {
    xsec_local = reinterp(xsec[0], water_lag(), flag, plag)
                     .reshape(frequency_grid.size());
  } else if (water_pert.empty()) {
    xsec_local =
        reinterp(xsec(joker, 0, joker, joker), temperature_lag(), flag, plag)
            .reshape(frequency_grid.size());
  } else {
    xsec_local = reinterp(xsec, temperature_lag(), water_lag(), flag, plag)
                     .reshape(frequency_grid.size());
  }

  const Numeric nd = atm_point.number_density(species);
  for (Index i = 0; i < frequency_grid.size(); ++i) {
    absorption[i] += xsec_local[i] * nd;
  }
}
}  // namespace lookup