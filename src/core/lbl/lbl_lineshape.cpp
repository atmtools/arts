#include "lbl_lineshape.h"

#include <algorithm>
#include <memory>
#include <ranges>

#include "artstime.h"
#include "debug.h"
#include "lbl_data.h"
#include "lbl_lineshape_voigt.h"
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
          [](auto& bnd) {
            return bnd.lineshape == Lineshape::VP and
                   bnd.linestrength == Linestrength::LTE;
          },
          &band::data))
    return std::make_unique<voigt::lte::ComputeData>(
        f_grid, atm, los, zeeman::pol::no);
  return nullptr;
}

void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const std::span<const lbl::band>& bnds,
               const AtmPoint& atm,
               const Vector2 los) {
  auto voigt_lte_data = init_voigt_lte_data(f_grid, bnds, atm, los);

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

  const auto calc_switch = [&](const QuantumIdentifier& bnd_key,
                               const band_data& bnd,
                               const zeeman::pol pol) {
    switch (bnd.lineshape) {
      case Lineshape::VP:
        switch (bnd.linestrength) {
          case Linestrength::LTE:
            calc_voigt_lte(bnd_key, bnd, pol);
            break;
          case Linestrength::FINAL:
            ARTS_USER_ERROR("Bad line strength state");
        }
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
