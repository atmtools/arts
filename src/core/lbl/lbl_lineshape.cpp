#include "lbl_lineshape.h"

#include <algorithm>
#include <memory>
#include <ranges>

#include "artstime.h"
#include "debug.h"
#include "lbl_data.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_lineshape_voigt_nlte.h"
#include "lbl_zeeman.h"
#include "quantum_numbers.h"

namespace lbl {
std::unique_ptr<voigt::lte::ComputeData> init_voigt_lte_data(
    const ExhaustiveConstVectorView& f_grid,
    const std::span<const lbl::band>& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (std::ranges::any_of(
          bnds,
          [](auto& bnd) { return bnd.lineshape == Lineshape::VP_LTE; },
          &band::data))
    return std::make_unique<voigt::lte::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}
std::unique_ptr<voigt::nlte::ComputeData> init_voigt_line_nlte_data(
    const ExhaustiveConstVectorView& f_grid,
    const std::span<const lbl::band>& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (std::ranges::any_of(
          bnds,
          [](auto& bnd) { return bnd.lineshape == Lineshape::VP_LINE_NLTE; },
          &band::data))
    return std::make_unique<voigt::nlte::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}

void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               matpack::matpack_view<Stokvec, 2, false, true> dsv,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const std::span<const lbl::band>& bnds,
               const AtmPoint& atm,
               const Vector2 los) {
  auto voigt_lte_data = init_voigt_lte_data(f_grid, bnds, atm, los);
  auto voigt_line_nlte_data = init_voigt_line_nlte_data(f_grid, bnds, atm, los);

  const auto calc_voigt_lte = [&](const QuantumIdentifier& bnd_key,
                                  const band_data& bnd,
                                  const zeeman::pol pol) {
    voigt::lte::calculate(pm,
                          dpm,
                          *voigt_lte_data,
                          f_grid,
                          jacobian_targets,
                          bnd_key,
                          bnd,
                          atm,
                          pol);
  };

  const auto calc_voigt_line_nlte = [&](const QuantumIdentifier& bnd_key,
                                        const band_data& bnd,
                                        const zeeman::pol pol) {
    voigt::nlte::calculate(pm,
                           sv,
                           dpm,
                           dsv,
                           *voigt_line_nlte_data,
                           f_grid,
                           jacobian_targets,
                           bnd_key,
                           bnd,
                           atm,
                           pol);
  };

  const auto calc_switch = [&](const QuantumIdentifier& bnd_key,
                               const band_data& bnd,
                               const zeeman::pol pol) {
    switch (bnd.lineshape) {
      case Lineshape::VP_LTE:
        calc_voigt_lte(bnd_key, bnd, pol);
        break;
      case Lineshape::VP_LINE_NLTE:
        calc_voigt_line_nlte(bnd_key, bnd, pol);
        break;
      case Lineshape::FINAL:
        ARTS_USER_ERROR("Bad line shape state");
    }
  };

  for (auto& [bnd_key, bnd] : bnds) {
    calc_switch(bnd_key, bnd, zeeman::pol::no);
  }

  for (auto pol : {zeeman::pol::pi, zeeman::pol::sm, zeeman::pol::sp}) {
    if (voigt_lte_data) voigt_lte_data->update_zeeman(los, atm.mag, pol);

    for (auto& [bnd_key, bnd] : bnds) {
      calc_switch(bnd_key, bnd, pol);
    }
  }
}
}  // namespace lbl
