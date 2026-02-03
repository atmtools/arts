#include "lbl_lineshape.h"

#include <debug.h>
#include <quantum.h>

#include <algorithm>
#include <memory>
#include <ranges>

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"
#include "lbl_lineshape_voigt_ecs.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_lineshape_voigt_lte_mirrored.h"
#include "lbl_lineshape_voigt_nlte.h"

namespace lbl {
namespace {
std::unique_ptr<voigt::lte::ComputeData> init_voigt_lte_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (stdr::any_of(bnds | stdv::values, [](auto& bnd) {
        return bnd.lineshape == LineByLineLineshape::VP_LTE;
      }))
    return std::make_unique<voigt::lte::ComputeData>(
        f_grid, atm, los, ZeemanPolarization::no);
  return nullptr;
}

std::unique_ptr<voigt::lte_mirror::ComputeData> init_voigt_lte_mirrored_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (stdr::any_of(bnds | stdv::values, [](auto& bnd) {
        return bnd.lineshape == LineByLineLineshape::VP_LTE_MIRROR;
      }))
    return std::make_unique<voigt::lte_mirror::ComputeData>(
        f_grid, atm, los, ZeemanPolarization::no);
  return nullptr;
}

std::unique_ptr<voigt::nlte::ComputeData> init_voigt_line_nlte_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (stdr::any_of(bnds | stdv::values, [](auto& bnd) {
        return bnd.lineshape == LineByLineLineshape::VP_LINE_NLTE;
      }))
    return std::make_unique<voigt::nlte::ComputeData>(
        f_grid, atm, los, ZeemanPolarization::no);
  return nullptr;
}

std::unique_ptr<voigt::ecs::ComputeData> init_voigt_abs_ecs_data(
    const ConstVectorView& f_grid,
    const AbsorptionBands& bnds,
    const AtmPoint& atm,
    const Vector2 los) {
  if (stdr::any_of(bnds | stdv::values, [](auto& bnd) {
        return bnd.lineshape == LineByLineLineshape::VP_ECS_MAKAROV or
               bnd.lineshape == LineByLineLineshape::VP_ECS_HARTMANN;
      }))
    return std::make_unique<voigt::ecs::ComputeData>(
        f_grid, atm, los, ZeemanPolarization::no);
  return nullptr;
}
}  // namespace

void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               PropmatMatrixView dpm,
               StokvecMatrixView dsv,
               const ConstVectorView f_grid,
               const Range& f_range,
               const Jacobian::Targets& jac_targets,
               const SpeciesEnum species,
               const AbsorptionBands& bnds,
               const LinemixingEcsData& abs_ecs_data,
               const AtmPoint& atm,
               const Vector2 los,
               const bool no_negative_absorption) {
  auto voigt_lte_data = init_voigt_lte_data(f_grid[f_range], bnds, atm, los);
  auto voigt_lte_mirror_data =
      init_voigt_lte_mirrored_data(f_grid[f_range], bnds, atm, los);
  auto voigt_line_nlte_data =
      init_voigt_line_nlte_data(f_grid[f_range], bnds, atm, los);
  auto voigt_abs_ecs_data =
      init_voigt_abs_ecs_data(f_grid[f_range], bnds, atm, los);

  const auto calc_voigt_lte = [&](const QuantumIdentifier& bnd_key,
                                  const band_data& bnd,
                                  const ZeemanPolarization pol) {
    voigt::lte::calculate(pm,
                          dpm,
                          *voigt_lte_data,
                          f_grid,
                          f_range,
                          jac_targets,
                          bnd_key,
                          bnd,
                          atm,
                          pol,
                          no_negative_absorption);
  };

  const auto calc_voigt_lte_mirrored = [&](const QuantumIdentifier& bnd_key,
                                           const band_data& bnd,
                                           const ZeemanPolarization pol) {
    voigt::lte_mirror::calculate(pm,
                                 dpm,
                                 *voigt_lte_mirror_data,
                                 f_grid,
                                 f_range,
                                 jac_targets,
                                 bnd_key,
                                 bnd,
                                 atm,
                                 pol,
                                 no_negative_absorption);
  };

  const auto calc_voigt_line_nlte = [&](const QuantumIdentifier& bnd_key,
                                        const band_data& bnd,
                                        const ZeemanPolarization pol) {
    voigt::nlte::calculate(pm,
                           sv,
                           dpm,
                           dsv,
                           *voigt_line_nlte_data,
                           f_grid,
                           f_range,
                           jac_targets,
                           bnd_key,
                           bnd,
                           atm,
                           pol,
                           no_negative_absorption);
  };

  const auto calc_voigt_ecs_linemixing = [&](const QuantumIdentifier& bnd_key,
                                             const band_data& bnd,
                                             const ZeemanPolarization pol) {
    auto it = abs_ecs_data.find(bnd_key.isot);
    if (it == abs_ecs_data.end()) {
      ARTS_USER_ERROR("No ECS data for isotopologue {}", bnd_key.isot);
    }

    voigt::ecs::calculate(pm,
                          dpm,
                          *voigt_abs_ecs_data,
                          f_grid,
                          f_range,
                          jac_targets,
                          bnd_key,
                          bnd,
                          it->second,
                          atm,
                          pol,
                          no_negative_absorption);
  };

  const auto calc_switch = [&](const QuantumIdentifier& bnd_key,
                               const band_data& bnd,
                               const ZeemanPolarization pol) {
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
      case LineByLineLineshape::VP_ECS_MAKAROV: [[fallthrough]];
      case LineByLineLineshape::VP_ECS_HARTMANN:
        calc_voigt_ecs_linemixing(bnd_key, bnd, pol);
        break;
    }
  };

  for (auto& [bnd_key, bnd] : bnds) {
    if (species == bnd_key.isot.spec or species == SpeciesEnum::Bath) {
      calc_switch(bnd_key, bnd, ZeemanPolarization::no);
    }
  }

  for (auto pol : {ZeemanPolarization::pi,
                   ZeemanPolarization::sm,
                   ZeemanPolarization::sp}) {
    if (voigt_lte_data) voigt_lte_data->update_zeeman(los, atm.mag, pol);
    if (voigt_lte_mirror_data) voigt_lte_mirror_data->update_zeeman(los, atm.mag, pol);
    if (voigt_line_nlte_data) voigt_line_nlte_data->update_zeeman(los, atm.mag, pol);

    for (auto& [bnd_key, bnd] : bnds) {
      if (species == bnd_key.isot.spec or species == SpeciesEnum::Bath) {
        calc_switch(bnd_key, bnd, pol);
      }
    }
  }
}
}  // namespace lbl
