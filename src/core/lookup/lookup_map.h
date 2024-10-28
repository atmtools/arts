#pragma once

#include <atm.h>
#include <isotopologues.h>
#include <lbl.h>
#include <matpack.h>

#include <unordered_map>

namespace lookup {
struct table {
  //! The frequency grid in Hz
  std::shared_ptr<const AscendingGrid> f_grid{std::make_shared<AscendingGrid>()};

  //! The atmospheric states (Pressure must be in decreasing order)
  std::shared_ptr<const ArrayOfAtmPoint> atmref{std::make_shared<ArrayOfAtmPoint>()};

  //! The pressure grid in std::log Pa [same dimension as atm]
  DescendingGrid log_p_grid;

  //! The temperautre perturbation grid in K [any number of elements or empty for nothing]
  AscendingGrid t_pert;

  //! The humidity perturbation grid in fractional units [any number of elements or empty for nothing]
  AscendingGrid water_pert;

  //! Local grids so that pressure interpolation may work
  Vector water_atmref;
  Vector t_atmref;

  /*! The absorption cross section table

      Dimensions: t_pert.size() x water_pert.size() x f_grid -> size() x atm -> size()
  */
  Tensor4 xsec;

  table()                        = default;
  table(const table&)            = default;
  table(table&&)                 = default;
  table& operator=(const table&) = default;
  table& operator=(table&&)      = default;

  table(const SpeciesEnum& species,
        std::shared_ptr<const AscendingGrid> f_grid,
        std::shared_ptr<const ArrayOfAtmPoint> atm,
        const AbsorptionBands& absorption_bands,
        const LinemixingEcsData& ecs_data,
        AscendingGrid t_pert     = {},
        AscendingGrid water_pert = {});

  void absorption(ExhaustiveVectorView absorption,
                  const SpeciesEnum& species,
                  const Index& p_interp_order,
                  const Index& t_interp_order,
                  const Index& water_interp_order,
                  const Index& f_interp_order,
                  const AtmPoint& atm_point,
                  const AscendingGrid& frequency_grid,
                  const Numeric& extpolfac) const;
};
}  // namespace lookup

using AbsorptionLookupTable  = lookup::table;
using AbsorptionLookupTables = std::unordered_map<SpeciesEnum, AbsorptionLookupTable>;

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
  FmtContext::iterator format(const AbsorptionLookupTable& v, FmtContext& ctx) const {
    return tags.format(ctx,
                       "f_grid: "sv,
                       *v.f_grid,
                       "\natmref: "sv,
                       *v.atmref,
                       "\nlog_p_grid: "sv,
                       v.log_p_grid,
                       "\nt_pert: "sv,
                       v.t_pert,
                       "\nwater_pert: "sv,
                       v.water_pert,
                       "\nwater_atmref: "sv,
                       v.water_atmref,
                       "\nt_atmref: "sv,
                       v.t_atmref,
                       "\nxsec: "sv,
                       v.xsec);
  }
};
