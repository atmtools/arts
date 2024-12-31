#include "lbl_lineshape.h"

#include <algorithm>
#include <memory>
#include <ranges>

#include "debug.h"
#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"
#include "lbl_lineshape_voigt_ecs.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_lineshape_voigt_lte_mirrored.h"
#include "lbl_lineshape_voigt_nlte.h"
#include "lbl_zeeman.h"
#include "quantum_numbers.h"
#include "species.h"

namespace lbl {
std::unique_ptr<voigt::lte::ComputeData> init_voigt_lte_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (std::ranges::any_of(
          bnds |std::ranges::views::values,
          [](auto& bnd) { return bnd.lineshape == LineByLineLineshape::VP_LTE; }))
    return std::make_unique<voigt::lte::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}
std::unique_ptr<voigt::lte_mirror::ComputeData> init_voigt_lte_mirrored_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (std::ranges::any_of(
          bnds|std::ranges::views::values,
          [](auto& bnd) { return bnd.lineshape == LineByLineLineshape::VP_LTE_MIRROR; }))
    return std::make_unique<voigt::lte_mirror::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}

std::unique_ptr<voigt::nlte::ComputeData> init_voigt_line_nlte_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (std::ranges::any_of(
          bnds|std::ranges::views::values,
          [](auto& bnd) { return bnd.lineshape == LineByLineLineshape::VP_LINE_NLTE; }))
    return std::make_unique<voigt::nlte::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}

std::unique_ptr<voigt::ecs::ComputeData> init_voigt_ecs_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (std::ranges::any_of(
          bnds|std::ranges::views::values,
          [](auto& bnd) {
            return bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV or
                   bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN;
          }))
    return std::make_unique<voigt::ecs::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}

void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               PropmatMatrixView dpm,
               StokvecMatrixView dsv,
               const ConstVectorView f_grid,
               const Range& f_range,
               const Jacobian::Targets& jacobian_targets,
               const SpeciesEnum species,
               const AbsorptionBands& bnds,
               const linemixing::isot_map& ecs_data,
               const AtmPoint& atm,
               const Vector2 los,
               const bool no_negative_absorption) {
  auto voigt_lte_data = init_voigt_lte_data(f_grid[f_range], bnds, atm, los);
  auto voigt_lte_mirror_data =
      init_voigt_lte_mirrored_data(f_grid[f_range], bnds, atm, los);
  auto voigt_line_nlte_data = init_voigt_line_nlte_data(f_grid[f_range], bnds, atm, los);
  auto voigt_ecs_data = init_voigt_ecs_data(f_grid[f_range], bnds, atm, los);

  const auto calc_voigt_lte = [&](const QuantumIdentifier& bnd_key,
                                  const band_data& bnd,
                                  const zeeman::pol pol) {
    voigt::lte::calculate(pm,
                          dpm,
                          *voigt_lte_data,
                          f_grid,
                          f_range,
                          jacobian_targets,
                          bnd_key,
                          bnd,
                          atm,
                          pol,
                          no_negative_absorption);
  };

  const auto calc_voigt_lte_mirrored = [&](const QuantumIdentifier& bnd_key,
                                           const band_data& bnd,
                                           const zeeman::pol pol) {
    voigt::lte_mirror::calculate(pm,
                                 dpm,
                                 *voigt_lte_mirror_data,
                                 f_grid,
                                 f_range,
                                 jacobian_targets,
                                 bnd_key,
                                 bnd,
                                 atm,
                                 pol,
                                 no_negative_absorption);
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
                           f_range,
                           jacobian_targets,
                           bnd_key,
                           bnd,
                           atm,
                           pol,
                           no_negative_absorption);
  };

  const auto calc_voigt_ecs_linemixing = [&](const QuantumIdentifier& bnd_key,
                                             const band_data& bnd,
                                             const zeeman::pol pol) {
    auto it = ecs_data.find(bnd_key.Isotopologue());
    if (it == ecs_data.end()) {
      ARTS_USER_ERROR("No ECS data for isotopologue {}", bnd_key.Isotopologue());
    }

    voigt::ecs::calculate(pm,
                          dpm,
                          *voigt_ecs_data,
                          f_grid,
                          f_range,
                          jacobian_targets,
                          bnd_key,
                          bnd,
                          it->second,
                          atm,
                          pol,
                          no_negative_absorption);
  };

  const auto calc_switch = [&](const QuantumIdentifier& bnd_key,
                               const band_data& bnd,
                               const zeeman::pol pol) {
    switch (bnd.lineshape) {
      case LineByLineLineshape::VP_LTE:
        calc_voigt_lte(bnd_key, bnd, pol);
        break;
      case LineByLineLineshape::VP_LTE_MIRROR:
        calc_voigt_lte_mirrored(bnd_key, bnd, pol);
        break;
      case LineByLineLineshape::VP_LINE_NLTE:
        calc_voigt_line_nlte(bnd_key, bnd, pol);
        break;
      case LineByLineLineshape::VP_ECS_MAKAROV:
        [[fallthrough]];
      case LineByLineLineshape::VP_ECS_HARTMANN:
        calc_voigt_ecs_linemixing(bnd_key, bnd, pol);
        break;
    }
  };

  for (auto& [bnd_key, bnd] : bnds) {
    if (species == bnd_key.Species() or species == SpeciesEnum::Bath) {
      calc_switch(bnd_key, bnd, zeeman::pol::no);
    }
  }

  for (auto pol : {zeeman::pol::pi, zeeman::pol::sm, zeeman::pol::sp}) {
    if (voigt_lte_data) voigt_lte_data->update_zeeman(los, atm.mag, pol);

    for (auto& [bnd_key, bnd] : bnds) {
      if (species == bnd_key.Species() or species == SpeciesEnum::Bath) {
        calc_switch(bnd_key, bnd, pol);
      }
    }
  }
}
}  // namespace lbl
