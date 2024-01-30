#include "lbl_fwd.h"

#include <physics_funcs.h>

#include <limits>
#include <ranges>
#include <type_traits>

#include "atm.h"
#include "configtypes.h"
#include "debug.h"
#include "lbl_data.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_zeeman.h"

namespace lbl::fwd {
void line_storage::adapt() {
  ARTS_USER_ERROR_IF(not bands, "No bands set")
  ARTS_USER_ERROR_IF(not atm, "No atm set")
  ARTS_USER_ERROR_IF(not good_enum(pol), "Bad polarization set")

  lte_shapes.resize(0);
  cutoff_lte_shapes.resize(0);
  lte_mirror_shapes.resize(0);
  cutoff_lte_mirror_shapes.resize(0);
  nlte_shapes.resize(0);
  cutoff_nlte_shapes.resize(0);

  std::vector<line_pos> tmp_pos;
  std::vector<voigt::lte::single_shape> tmp_lte_shapes;
  std::vector<voigt::lte_mirror::single_shape> tmp_lte_mirror_shapes;
  std::vector<voigt::nlte::single_shape> tmp_nlte_shapes;

  const auto pushback =
      [](const auto& inlines, auto& out, auto& cutoff_out, const auto& band) {
        switch (band.cutoff) {
          case CutoffType::None:
            for (auto& line : inlines) out.push_back(line);
            break;
          case CutoffType::ByLine:
            ARTS_USER_ERROR_IF(band.cutoff_value != 750e9,
                               "Only 750 GHz cutoff is supported")
            for (auto& line : inlines) cutoff_out.push_back(line);
            break;
          case CutoffType::FINAL:
            ARTS_USER_ERROR("Undefined cutoff type")
        }
      };

  for (auto& [qid, band] : *bands) {
    switch (band.lineshape) {
      case Lineshape::VP_LTE:
        band_shape_helper(tmp_lte_shapes,
                          tmp_pos,
                          qid.Isotopologue(),
                          band,
                          *atm,
                          std::numeric_limits<Numeric>::lowest(),
                          std::numeric_limits<Numeric>::max(),
                          pol);
        pushback(tmp_lte_shapes, lte_shapes, cutoff_lte_shapes, band);
        break;
      case Lineshape::VP_LTE_MIRROR:
        band_shape_helper(tmp_lte_mirror_shapes,
                          tmp_pos,
                          qid.Isotopologue(),
                          band,
                          *atm,
                          std::numeric_limits<Numeric>::lowest(),
                          std::numeric_limits<Numeric>::max(),
                          pol);
        pushback(tmp_lte_mirror_shapes,
                 lte_mirror_shapes,
                 cutoff_lte_mirror_shapes,
                 band);
        break;
      case Lineshape::VP_LINE_NLTE:
        band_shape_helper(tmp_nlte_shapes,
                          tmp_pos,
                          qid,
                          band,
                          *atm,
                          std::numeric_limits<Numeric>::lowest(),
                          std::numeric_limits<Numeric>::max(),
                          pol);
        pushback(tmp_nlte_shapes, nlte_shapes, cutoff_nlte_shapes, band);
        break;
      case Lineshape::VP_ECS_HARTMANN:
        [[fallthrough]];
      case Lineshape::VP_ECS_MAKAROV:
        ARTS_USER_ERROR("Cannot use ECS models in forward model")
      case Lineshape::FINAL: {
        ARTS_USER_ERROR("Undefined lineshape")
      }
    }

    lte_shapes.shrink_to_fit();
    nlte_shapes.shrink_to_fit();
    lte_mirror_shapes.shrink_to_fit();

    cutoff_lte_shapes.shrink_to_fit();
    cutoff_nlte_shapes.shrink_to_fit();
    cutoff_lte_mirror_shapes.shrink_to_fit();

    const auto by_freq = [](auto& l1, auto& l2) { return l1.f0 < l2.f0; };
    std::ranges::sort(cutoff_lte_shapes, by_freq);
    std::ranges::sort(cutoff_nlte_shapes, by_freq);
    std::ranges::sort(cutoff_lte_mirror_shapes, by_freq);
  }
}

line_storage::line_storage(std::shared_ptr<AbsorptionBands> bands_,
                           std::shared_ptr<AtmPoint> atm_,
                           const zeeman::pol pol_)
    : bands(std::move(bands_)), atm(std::move(atm_)), pol(pol_) {
  adapt();
}

void line_storage::set_lines(std::shared_ptr<AbsorptionBands> bands_) {
  bands = std::move(bands_);
  adapt();
}

void line_storage::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}

void line_storage::set_pol(zeeman::pol pol_) {
  pol = pol_;
  adapt();
}

struct NlteSumup {
  Complex abs{};
  Complex src{};
  NlteSumup operator+(const NlteSumup& other) const {
    return {abs + other.abs, src + other.src};
  }
  friend NlteSumup operator*(const Numeric& other, const NlteSumup& x) {
    return {other * x.abs, other * x.src};
  }
};

std::pair<Complex, Complex> line_storage::operator()(const Numeric f) const {
  static constexpr Numeric cutoff_frequency = 750e9;
  constexpr auto constant = Constant::c * Constant::c / (8 * Constant::pi);
  const auto nlte_fac =
      number_density(atm->pressure, atm->temperature) * constant * f;
  const auto lte_fac = -nlte_fac * std::expm1(-(Constant::h * f) /
                                              (Constant::k * atm->temperature));

  const auto valid_span = [fmin = f - cutoff_frequency,
                           fmax = f + cutoff_frequency](const auto& lines) {
    using T = std::decay_t<decltype(lines.front())>;
    const auto start = std::ranges::lower_bound(lines, fmin, {}, &T::f0);
    const auto end =
        std::ranges::upper_bound(start, lines.end(), fmax, {}, &T::f0);
    return std::span{start, end};
  };

  const auto lte_sumup = [f](auto&& lines, const bool cutoff) {
    Complex out{0.0, 0.0};
    if (cutoff) {
      for (auto& line : lines) {
        out += line(f) - line(line.f0 + cutoff_frequency);
      }
    } else {
      for (auto& line : lines) {
        out += line(f);
      }
    }
    return out;
  };

  const auto nlte_sumup = [f](auto&& lines, const bool cutoff) -> NlteSumup {
    Complex aout{0.0, 0.0}, sout{0.0, 0.0};
    if (cutoff) {
      for (auto& line : lines) {
        const auto [a, s] = line(f);
        const auto [ac, sc] = line(line.f0 + cutoff_frequency);
        aout += a - ac;
        sout += s - sc;
      }
    } else {
      for (auto& line : lines) {
        const auto [a, s] = line(f);
        aout += a;
        sout += s;
      }
    }
    return {aout, sout};
  };

  const Complex lte_abs =
      lte_fac * (lte_sumup(lte_shapes, false) +
                 lte_sumup(valid_span(cutoff_lte_shapes), true) +
                 lte_sumup(lte_mirror_shapes, false) +
                 lte_sumup(valid_span(cutoff_lte_mirror_shapes), true));

  const NlteSumup nlte =
      nlte_fac * (nlte_sumup(nlte_shapes, false) +
                  nlte_sumup(valid_span(cutoff_nlte_shapes), true));

  return {lte_abs + nlte.abs, nlte.src};
}

Propmat line_storage::polarization(Vector2 los) const {
  return zeeman::norm_view(pol, atm->mag, los);
}
}  // namespace lbl::fwd