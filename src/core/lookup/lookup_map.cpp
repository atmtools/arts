#include "lookup_map.h"

#include <arts_omp.h>
#include <jacobian.h>
#include <lagrange_interp.h>

namespace lookup {
table::table()                            = default;
table::table(const table&)                = default;
table::table(table&&) noexcept            = default;
table& table::operator=(const table&)     = default;
table& table::operator=(table&&) noexcept = default;

bool table::do_t() const { return t_pert and not t_pert->empty(); }

bool table::do_w() const { return w_pert and not w_pert->empty(); }

bool table::do_p() const { return log_p_grid and not log_p_grid->empty(); }

bool table::do_f() const { return f_grid and not f_grid->empty(); }

Index table::t_size() const {
  return do_t() ? static_cast<Index>(t_pert->size()) : 1;
}

Index table::w_size() const {
  return do_w() ? static_cast<Index>(w_pert->size()) : 1;
}

Index table::p_size() const {
  return do_p() ? static_cast<Index>(log_p_grid->size()) : 1;
}

Index table::f_size() const {
  return do_f() ? static_cast<Index>(f_grid->size()) : 1;
}

std::array<Index, 4> table::grid_shape() const {
  return {t_size(), w_size(), p_size(), f_size()};
}

table::table(const SpeciesEnum& species,
             const ArrayOfAtmPoint& atmref,
             std::shared_ptr<const AscendingGrid> f_grid_,
             const AbsorptionBands& absorption_bands,
             const LinemixingEcsData& ecs_data,
             std::shared_ptr<const AscendingGrid> t_pert_,
             std::shared_ptr<const AscendingGrid> w_pert_) try
    : f_grid(std::move(f_grid_)),
      t_pert(std::move(t_pert_)),
      w_pert(std::move(w_pert_)),
      water_atmref(atmref.size()),
      t_atmref(atmref.size()) {
  ARTS_USER_ERROR_IF(not do_f(), "Frequency grid is not set.")

  const bool do_water = do_w();

  log_p_grid = std::make_shared<const DescendingGrid>(
      atmref.begin(), atmref.end(), [](const AtmPoint& x) {
        return std::log(x.pressure);
      });

  xsec.resize(grid_shape());

  std::transform(
      atmref.begin(),
      atmref.end(),
      water_atmref.begin(),
      [do_water](const auto& x) { return do_water ? x["H2O"_spec] : NAN; });

  std::transform(
      atmref.begin(), atmref.end(), t_atmref.begin(), [](const auto& x) {
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

  const AscendingGrid empty_water({1.0});
  const AscendingGrid empty_t_pert({0.0});
  const AscendingGrid& water_vmr_local(do_water ? *w_pert : empty_water);
  const AscendingGrid& t_pert_local(do_t() ? *t_pert : empty_t_pert);

#pragma omp parallel for collapse(3) if (not arts_omp_in_parallel()) \
    firstprivate(pm, sv, dpm, dsv)
  for (Size it = 0; it < t_pert_local.size(); ++it) {
    for (Size iw = 0; iw < water_vmr_local.size(); ++iw) {
      for (Size ip = 0; ip < atmref.size(); ++ip) {
        try {
          AtmPoint atm_point     = atmref[ip];
          atm_point.temperature += t_pert_local[it];
          if (do_water) atm_point["H2O"_spec] *= water_vmr_local[iw];

          pm = 0.0;
          lbl::calculate(pm,
                         sv,
                         dpm,
                         dsv,
                         *f_grid,
                         Range(0, f_grid->size()),
                         jacobian_targets,
                         species,
                         absorption_bands,
                         ecs_data,
                         atm_point,
                         los,
                         no_negative_absorption);

          const Numeric inv_nd = 1.0 / atm_point.number_density(species);
          for (Size ifreq = 0; ifreq < f_grid->size(); ++ifreq) {
            xsec[it, iw, ip, ifreq] = pm[ifreq].A() * inv_nd;
          }
        } catch (std::runtime_error& e) {
#pragma omp critical
          error += std::format("ERROR:\n{}\n\n", e.what());
        }
      }
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

std::array<lagrange_interp::lag_t<-1>, 1> table::pressure_lagrange(
    const Numeric& pressure,
    const Index interpolation_order,
    const Numeric& extpolation_factor) const try {
  ARTS_USER_ERROR_IF(not do_p(), "No pressure grid set.");
  const std::array<Numeric, 1> plog_local = {std::log(pressure)};
  const Vector& plog_v(*log_p_grid);
  return lagrange_interp::make_lags(plog_v,
                                    plog_local,
                                    interpolation_order,
                                    extpolation_factor,
                                    "Log-Pressure");
}
ARTS_METHOD_ERROR_CATCH

std::vector<lagrange_interp::lag_t<-1>> table::frequency_lagrange(
    const Vector& frequency_grid,
    const Index interpolation_order,
    const Numeric& extpolation_factor) const try {
  ARTS_USER_ERROR_IF(not do_f(), "No frequency grid set.");
  const Vector& f_grid_v(*f_grid);
  return lagrange_interp::make_lags(f_grid_v,
                                    frequency_grid,
                                    interpolation_order,
                                    extpolation_factor,
                                    "Frequency");
}
ARTS_METHOD_ERROR_CATCH

std::array<lagrange_interp::lag_t<-1>, 1> table::water_lagrange(
    const Numeric& water_vmr,
    const std::array<lagrange_interp::lag_t<-1>, 1> & pressure_lagrange,
    const Index interpolation_order,
    const Numeric& extpolation_factor) const try {
  ARTS_USER_ERROR_IF(not do_w(), "No water grid set.");
  const std::array<Numeric, 1> x{water_vmr /
                                 interp(water_atmref, pressure_lagrange[0])};
  const Vector& xi(*w_pert);
  return lagrange_interp::make_lags(
      xi, x, interpolation_order, extpolation_factor, "Water VMR");
}
ARTS_METHOD_ERROR_CATCH

std::array<lagrange_interp::lag_t<-1>, 1> table::temperature_lagrange(
    const Numeric& temperature,
    const std::array<lagrange_interp::lag_t<-1>, 1> & pressure_lagrange,
    const Index interpolation_order,
    const Numeric& extpolation_factor) const try {
  ARTS_USER_ERROR_IF(not do_t(), "No temperature grid set.");
  const std::array<Numeric, 1> x{temperature -
                                 interp(t_atmref, pressure_lagrange[0])};
  const Vector& xi(*t_pert);
  return lagrange_interp::make_lags(
      xi, x, interpolation_order, extpolation_factor, "Temperature");
}
ARTS_METHOD_ERROR_CATCH

void table::absorption(VectorView absorption,
                       const SpeciesEnum& species,
                       const Index& p_interp_order,
                       const Index& t_interp_order,
                       const Index& water_interp_order,
                       const Index& f_interp_order,
                       const AtmPoint& atm_point,
                       const AscendingGrid& frequency_grid,
                       const Numeric& extpolfac) const try {
  check();

  if (xsec.empty()) return;

  // Frequency grid positions
  const auto flag =
      frequency_lagrange(frequency_grid, f_interp_order, extpolfac);

  // Pressure grid positions
  const auto plag =
      pressure_lagrange(atm_point.pressure, p_interp_order, extpolfac);

  Vector xsec_local(frequency_grid.size());

  // Optional grid positions by switching
  if (do_w() and do_t()) {
    const auto wlag = water_lagrange(
        atm_point["H2O"_spec], plag, water_interp_order, extpolfac);
    const auto tlag = temperature_lagrange(
        atm_point.temperature, plag, t_interp_order, extpolfac);
    xsec_local =
        reinterp(xsec, tlag, wlag, plag, flag).reshape(frequency_grid.size());
  } else if (do_w()) {
    const auto wlag = water_lagrange(
        atm_point["H2O"_spec], plag, water_interp_order, extpolfac);
    xsec_local =
        reinterp(xsec[0], wlag, plag, flag).reshape(frequency_grid.size());
  } else if (do_t()) {
    const auto tlag = temperature_lagrange(
        atm_point.temperature, plag, t_interp_order, extpolfac);
    xsec_local = reinterp(xsec[joker, 0, joker, joker], tlag, plag, flag)
                     .reshape(frequency_grid.size());
  } else {
    xsec_local =
        reinterp(xsec[0][0], plag, flag).reshape(frequency_grid.size());
  }

  const Numeric nd = atm_point.number_density(species);
  for (Size i = 0; i < frequency_grid.size(); ++i) {
    absorption[i] += xsec_local[i] * nd;
  }
}
ARTS_METHOD_ERROR_CATCH

void table::check() const {
  ARTS_USER_ERROR_IF(not do_f() or not do_p(),
                     R"(Must have frequency and pressure grids.
  Frequency grid: {}
  Pressure grid:  {}
)",
                     do_f(),
                     do_p());

  const auto [t_size, w_size, p_size, f_size] = grid_shape();
  ARTS_USER_ERROR_IF(
      (xsec.shape() != std::array{t_size, w_size, p_size, f_size}),
      R"(The shape of the absorption cross section table is incorrect.

  Found:    {4:B,},
  Expected: [{0}, {1}, {2}, {3}]
    where
      {0} = t_size
      {1} = w_size
      {2} = p_size
      {3} = f_size
)",
      t_size,
      w_size,
      p_size,
      f_size,
      xsec.shape());

  ARTS_USER_ERROR_IF(water_atmref.size() != static_cast<Size>(p_size),
                     R"(Bad size of water_atmref
  Expected: {}
  Has size: {}
)",
                     p_size,
                     water_atmref.size());

  ARTS_USER_ERROR_IF(t_atmref.size() != static_cast<Size>(p_size),
                     R"(Bad size of t_atmref
  Expected: {}
  Has size: {}
)",
                     p_size,
                     t_atmref.size());
}

void extend_atmosphere(ArrayOfAtmPoint& atm,
                       const InterpolationExtrapolation extrapolation_type,
                       const Numeric new_max_pressure,
                       const Numeric new_min_pressure) {
  if (not std::isnan(new_max_pressure)) {
    ARTS_USER_ERROR_IF(new_max_pressure <= atm.back().pressure,
                       "The new maximum pressure must be greater than the "
                       "current maximum pressure: {} <= {}",
                       new_max_pressure,
                       atm.back().pressure)
    Atm::extend_in_pressure(atm, new_max_pressure, extrapolation_type);
  }

  if (not std::isnan(new_min_pressure)) {
    ARTS_USER_ERROR_IF(new_max_pressure >= atm.front().pressure,
                       "The new minimum pressure must be smaller than the "
                       "current minimum pressure: {} >= {}",
                       new_min_pressure,
                       atm.back().pressure)
    Atm::extend_in_pressure(atm, new_min_pressure, extrapolation_type);
  }
}
}  // namespace lookup

void xml_io_stream<AbsorptionLookupTable>::write(std::ostream& os,
                                                 const AbsorptionLookupTable& x,
                                                 bofstream* pbofs,
                                                 std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.f_grid, pbofs, "frequency"sv);
  xml_write_to_stream(os, x.log_p_grid, pbofs, "log-p"sv);
  xml_write_to_stream(os, x.t_pert, pbofs, "t-pert"sv);
  xml_write_to_stream(os, x.w_pert, pbofs, "water-pert"sv);
  xml_write_to_stream(os, x.water_atmref, pbofs, "water-ref"sv);
  xml_write_to_stream(os, x.t_atmref, pbofs, "t-ref"sv);
  xml_write_to_stream(os, x.xsec, pbofs, "xsec");

  tag.write_to_end_stream(os);
}

void xml_io_stream<AbsorptionLookupTable>::read(std::istream& is,
                                                AbsorptionLookupTable& x,
                                                bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.f_grid, pbifs);
  xml_read_from_stream(is, x.log_p_grid, pbifs);
  xml_read_from_stream(is, x.t_pert, pbifs);
  xml_read_from_stream(is, x.w_pert, pbifs);
  xml_read_from_stream(is, x.water_atmref, pbifs);
  xml_read_from_stream(is, x.t_atmref, pbifs);
  xml_read_from_stream(is, x.xsec, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);

  x.check();
}
