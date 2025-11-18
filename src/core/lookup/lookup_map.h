#pragma once

#include <atm.h>
#include <isotopologues.h>
#include <lagrange_interp.h>
#include <lbl.h>
#include <matpack.h>
#include <xml.h>

#include <unordered_map>

namespace lookup {
struct table {
  //! The frequency grid in Hz
  std::shared_ptr<const AscendingGrid> f_grid;

  //! The pressure grid in log(P) [any number of elements]
  std::shared_ptr<const DescendingGrid> log_p_grid;

  //! The temperautre perturbation grid in K [any number of elements or empty for nothing]
  std::shared_ptr<const AscendingGrid> t_pert;

  //! The humidity perturbation grid in fractional units [any number of elements or empty for nothing]
  std::shared_ptr<const AscendingGrid> w_pert;

  //! Local grids so that pressure interpolation may work [all in log_p_grid size]
  Vector water_atmref;
  Vector t_atmref;

  /*! The absorption cross section table

      Dimensions: t_pert.size() x w_pert.size() x log_p_grid -> size() x f_grid -> size()
  */
  Tensor4 xsec;

  table();
  table(const table&);
  table(table&&) noexcept;
  table& operator=(const table&);
  table& operator=(table&&) noexcept;

  table(const SpeciesEnum& species,
        const ArrayOfAtmPoint& atmref,
        std::shared_ptr<const AscendingGrid> f_grid,
        const AbsorptionBands& abs_bands,
        const LinemixingEcsData& abs_ecs_data,
        std::shared_ptr<const AscendingGrid> t_pert = nullptr,
        std::shared_ptr<const AscendingGrid> w_pert = nullptr);

  void absorption(VectorView absorption,
                  const SpeciesEnum& species,
                  const Index& p_interp_order,
                  const Index& t_interp_order,
                  const Index& water_interp_order,
                  const Index& f_interp_order,
                  const AtmPoint& atm_point,
                  const AscendingGrid& freq_grid,
                  const Numeric& extpolfac) const;

  [[nodiscard]] bool do_t() const;
  [[nodiscard]] bool do_w() const;
  [[nodiscard]] bool do_p() const;
  [[nodiscard]] bool do_f() const;

  [[nodiscard]] Index t_size() const;
  [[nodiscard]] Index w_size() const;
  [[nodiscard]] Index p_size() const;
  [[nodiscard]] Index f_size() const;

  [[nodiscard]] std::array<Index, 4> grid_shape() const;
  void check() const;

  [[nodiscard]] std::array<lagrange_interp::lag_t<-1>, 1> pressure_lagrange(
      const Numeric& pressure,
      const Index interpolation_order,
      const Numeric& extpolation_factor) const;

  [[nodiscard]] std::vector<lagrange_interp::lag_t<-1>> frequency_lagrange(
      const Vector& freq_grid,
      const Index interpolation_order,
      const Numeric& extpolation_factor) const;

  [[nodiscard]] std::array<lagrange_interp::lag_t<-1>, 1> water_lagrange(
      const Numeric& water_vmr,
      const std::array<lagrange_interp::lag_t<-1>, 1>& pressure_lagrange,
      const Index interpolation_order,
      const Numeric& extpolation_factor) const;

  [[nodiscard]] std::array<lagrange_interp::lag_t<-1>, 1> temperature_lagrange(
      const Numeric& temperature,
      const std::array<lagrange_interp::lag_t<-1>, 1>& pressure_lagrange,
      const Index interpolation_order,
      const Numeric& extpolation_factor) const;
};

/** Wraps calling Atm::extend_in_pressure but for an atmosphere fitting a lookup table.
 * 
 * Additional checks are performed to ensure that the input fits the ideas of the lookup table.
 * 
 * @param atm The atmospheric grid that goes into creating a lookup table
 * @param extrapolation_type See Atm::extend_in_pressure
 * @param new_max_pressure The maximum pressure to extend to, must be greater than the last pressure in atm
 * @param new_min_pressure The minimum pressure to extend to, must be less than the first pressure in atm
 */
void extend_atmosphere(ArrayOfAtmPoint& atm,
                       const InterpolationExtrapolation extrapolation_type =
                           InterpolationExtrapolation::Nearest,
                       const Numeric new_max_pressure = NAN,
                       const Numeric new_min_pressure = NAN);
}  // namespace lookup

using AbsorptionLookupTable = lookup::table;
using AbsorptionLookupTables =
    std::unordered_map<SpeciesEnum, AbsorptionLookupTable>;

template <>
struct std::formatter<AbsorptionLookupTable> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const AbsorptionLookupTable& v,
                              FmtContext& ctx) const {
    tags.format(ctx, "f_grid: "sv);
    if (v.f_grid) tags.format(ctx, *v.f_grid);
    tags.format(ctx, "\nlog_p_grid: "sv);
    if (v.log_p_grid) tags.format(ctx, *v.log_p_grid);
    tags.format(ctx, "\nt_pert: "sv);
    if (v.t_pert) tags.format(ctx, *v.t_pert);
    tags.format(ctx, "\nw_pert: "sv);
    if (v.w_pert) tags.format(ctx, *v.w_pert);
    return tags.format(ctx,

                       "\nwater_atmref: "sv,
                       v.water_atmref,
                       "\nt_atmref: "sv,
                       v.t_atmref,
                       "\nxsec:\n"sv,
                       v.xsec);
  }
};

template <>
struct xml_io_stream<AbsorptionLookupTable> {
  static constexpr std::string_view type_name = "AbsorptionLookupTable"sv;

  static void write(std::ostream& os,
                    const AbsorptionLookupTable& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   AbsorptionLookupTable& x,
                   bifstream* pbifs = nullptr);
};
