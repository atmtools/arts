#include "lbl_fwd.h"

#include <physics_funcs.h>

#include <iomanip>
#include <limits>

#include "configtypes.h"
#include "debug.h"
#include "lbl_data.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_zeeman.h"

namespace lbl::fwd {
namespace models {
void lte::adapt() try {
  lines.lines.resize(0);
  cutoff_lines.lines.resize(0);
  cutoff.resize(0);

  if (not bands) {
    return;
  }

  if (bands->empty()) {
    return;
  }

  ARTS_USER_ERROR_IF(not atm, "Must have an atmosphere")

  std::vector<voigt::lte::single_shape> shapes;
  std::vector<line_pos> shapes_pos;
  decltype(cutoff) cutoff_this;

  for (auto& [qid, band] : *bands) {
    if (band.lineshape != LineByLineLineshape::VP_LTE) continue;

    band_shape_helper(shapes,
                      shapes_pos,
                      qid.Isotopologue(),
                      band,
                      *atm,
                      std::numeric_limits<Numeric>::lowest(),
                      std::numeric_limits<Numeric>::max(),
                      pol);

    if (shapes.size() == 0) continue;

    voigt::lte::band_shape b{std::move(shapes), band.cutoff_value};
    switch (band.cutoff) {
      case LineByLineCutoffType::ByLine:
        cutoff_this.resize(b.lines.size());
        b(cutoff_this);

        for (auto& line : b.lines) {
          cutoff_lines.lines.push_back(line);
        }

        for (auto& c : cutoff_this) {
          cutoff.push_back(c);
        }
        break;
      case LineByLineCutoffType::None:
        for (auto& line : b.lines) {
          lines.lines.push_back(line);
        }
        break;
    }

    shapes = std::move(b.lines);
  }
}
ARTS_METHOD_ERROR_CATCH

void lte_mirror::adapt() {
  lines.lines.resize(0);
  cutoff_lines.lines.resize(0);
  cutoff.resize(0);

  if (not bands) {
    return;
  }

  ARTS_USER_ERROR_IF(not atm, "Must have an atmosphere")

  std::vector<voigt::lte_mirror::single_shape> shapes;
  std::vector<line_pos> shapes_pos;
  decltype(cutoff) cutoff_this;

  for (auto& [qid, band] : *bands) {
    if (band.lineshape != LineByLineLineshape::VP_LTE_MIRROR) continue;

    band_shape_helper(shapes,
                      shapes_pos,
                      qid.Isotopologue(),
                      band,
                      *atm,
                      std::numeric_limits<Numeric>::lowest(),
                      std::numeric_limits<Numeric>::max(),
                      pol);

    if (shapes.size() == 0) continue;

    voigt::lte_mirror::band_shape b{std::move(shapes), band.cutoff_value};
    switch (band.cutoff) {
      case LineByLineCutoffType::ByLine:
        cutoff_this.resize(b.lines.size());
        b(cutoff_this);

        for (auto& line : b.lines) {
          cutoff_lines.lines.push_back(line);
        }

        for (auto& c : cutoff_this) {
          cutoff.push_back(c);
        }
        break;
      case LineByLineCutoffType::None:
        for (auto& line : b.lines) {
          lines.lines.push_back(line);
        }
        break;
    }

    shapes = std::move(b.lines);
  }
}

void nlte::adapt() {
  lines.lines.resize(0);
  cutoff_lines.lines.resize(0);
  cutoff.resize(0);

  if (not bands) {
    return;
  }

  ARTS_USER_ERROR_IF(not atm, "Must have an atmosphere")

  std::vector<voigt::nlte::single_shape> shapes;
  std::vector<line_pos> shapes_pos;
  decltype(cutoff) cutoff_this;

  for (auto& [qid, band] : *bands) {
    if (band.lineshape != LineByLineLineshape::VP_LINE_NLTE) continue;

    band_shape_helper(shapes,
                      shapes_pos,
                      qid,
                      band,
                      *atm,
                      std::numeric_limits<Numeric>::lowest(),
                      std::numeric_limits<Numeric>::max(),
                      pol);

    if (shapes.size() == 0) continue;

    voigt::nlte::band_shape b{std::move(shapes), band.cutoff_value};
    switch (band.cutoff) {
      case LineByLineCutoffType::ByLine:
        cutoff_this.resize(b.lines.size());
        b(cutoff_this);

        for (auto& line : b.lines) {
          cutoff_lines.lines.push_back(line);
        }

        for (auto& c : cutoff_this) {
          cutoff.push_back(c);
        }
        break;
      case LineByLineCutoffType::None:
        for (auto& line : b.lines) {
          lines.lines.push_back(line);
        }
        break;
    }

    shapes = std::move(b.lines);
  }
}

std::pair<Complex, Complex> lte::operator()(const Numeric frequency) const {
  const auto scl = [N = number_density(atm->pressure, atm->temperature),
                    T = atm->temperature](auto f) {
    constexpr Numeric c = Constant::c * Constant::c / (8 * Constant::pi);
    const Numeric r = (Constant::h * f) / (Constant::k * T);
    return -N * f * std::expm1(-r) * c;
  }(frequency);

  return {scl * (lines(frequency) + cutoff_lines(cutoff, frequency)),
          Complex{0.0, 0.0}};
}

std::pair<Complex, Complex> lte_mirror::operator()(
    const Numeric frequency) const {
  const auto scl = [N = number_density(atm->pressure, atm->temperature),
                    T = atm->temperature](auto f) {
    constexpr Numeric c = Constant::c * Constant::c / (8 * Constant::pi);
    const Numeric r = (Constant::h * f) / (Constant::k * T);
    return -N * f * std::expm1(-r) * c;
  }(frequency);

  ARTS_ASSERT(lines(frequency) == cutoff_lines(cutoff, frequency) and
              (cutoff_lines(cutoff, frequency) == Complex{0.0, 0.0}))

  return {scl * (lines(frequency) + cutoff_lines(cutoff, frequency)),
          Complex{0.0, 0.0}};
}

std::pair<Complex, Complex> nlte::operator()(const Numeric frequency) const {
  const auto scl =
      [N = number_density(atm->pressure, atm->temperature)](auto f) {
        constexpr Numeric c = Constant::c * Constant::c / (8 * Constant::pi);
        return N * f * c;  // * std::expm1(-r); ??? It feels like LTE term...
      }(frequency);

  auto [a, s] = lines(frequency);
  auto [ac, sc] = cutoff_lines(cutoff, frequency);

  ARTS_ASSERT(a == s and ac == sc and (a == Complex{0.0, 0.0}))
  return {scl * (a + ac), scl * (s + sc)};
}

void lte::set_model(std::shared_ptr<ArrayOfAbsorptionBand> bands_) {
  bands = std::move(bands_);
  adapt();
}

void lte_mirror::set_model(std::shared_ptr<ArrayOfAbsorptionBand> bands_) {
  bands = std::move(bands_);
  adapt();
}

void nlte::set_model(std::shared_ptr<ArrayOfAbsorptionBand> bands_) {
  bands = std::move(bands_);
  adapt();
}

void lte::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}

void lte_mirror::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}

void nlte::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}

void lte::set_pol(zeeman::pol pol_) {
  pol = pol_;
  adapt();
}

void lte_mirror::set_pol(zeeman::pol pol_) {
  pol = pol_;
  adapt();
}

void nlte::set_pol(zeeman::pol pol_) {
  pol = pol_;
  adapt();
}

void lte::set(std::shared_ptr<ArrayOfAbsorptionBand> bands_,
              std::shared_ptr<AtmPoint> atm_,
              zeeman::pol pol_) {
  bands = std::move(bands_);
  atm = std::move(atm_);
  pol = pol_;
  adapt();
}

void lte_mirror::set(std::shared_ptr<ArrayOfAbsorptionBand> bands_,
                     std::shared_ptr<AtmPoint> atm_,
                     zeeman::pol pol_) {
  bands = std::move(bands_);
  atm = std::move(atm_);
  pol = pol_;
  adapt();
}

void nlte::set(std::shared_ptr<ArrayOfAbsorptionBand> bands_,
               std::shared_ptr<AtmPoint> atm_,
               zeeman::pol pol_) {
  bands = std::move(bands_);
  atm = std::move(atm_);
  pol = pol_;
  adapt();
}
}  // namespace models

line_storage::line_storage(std::shared_ptr<AtmPoint> atm_,
                           std::shared_ptr<ArrayOfAbsorptionBand> bands_)
    : atm(std::move(atm_)), bands(std::move(bands_)) {
  for (auto& [qid, band] : *bands) {
    ARTS_USER_ERROR_IF(
        band.lineshape != LineByLineLineshape::VP_LTE and
            band.lineshape != LineByLineLineshape::VP_LTE_MIRROR and
            band.lineshape != LineByLineLineshape::VP_LINE_NLTE,
        "Lineshape not supported ",
        '"', band.lineshape, '"',
        " for band: ",
        qid)
  }

  lte[static_cast<Size>(zeeman::pol::sm)].set(bands, atm, zeeman::pol::sm);
  lte[static_cast<Size>(zeeman::pol::pi)].set(bands, atm, zeeman::pol::pi);
  lte[static_cast<Size>(zeeman::pol::sp)].set(bands, atm, zeeman::pol::sp);
  lte[static_cast<Size>(zeeman::pol::no)].set(bands, atm, zeeman::pol::no);

  lte_mirror[static_cast<Size>(zeeman::pol::sm)].set(
      bands, atm, zeeman::pol::sm);
  lte_mirror[static_cast<Size>(zeeman::pol::pi)].set(
      bands, atm, zeeman::pol::pi);
  lte_mirror[static_cast<Size>(zeeman::pol::sp)].set(
      bands, atm, zeeman::pol::sp);
  lte_mirror[static_cast<Size>(zeeman::pol::no)].set(
      bands, atm, zeeman::pol::no);

  nlte[static_cast<Size>(zeeman::pol::sm)].set(bands, atm, zeeman::pol::sm);
  nlte[static_cast<Size>(zeeman::pol::pi)].set(bands, atm, zeeman::pol::pi);
  nlte[static_cast<Size>(zeeman::pol::sp)].set(bands, atm, zeeman::pol::sp);
  nlte[static_cast<Size>(zeeman::pol::no)].set(bands, atm, zeeman::pol::no);
}

void line_storage::set_model(std::shared_ptr<ArrayOfAbsorptionBand> bands_) {
  for (auto& m : lte) m.set_model(bands_);
  for (auto& m : lte_mirror) m.set_model(bands_);
  for (auto& m : nlte) m.set_model(bands_);
  bands = std::move(bands_);
}

void line_storage::set_atm(std::shared_ptr<AtmPoint> atm_) {
  for (auto& m : lte) m.set_atm(atm_);
  for (auto& m : lte_mirror) m.set_atm(atm_);
  for (auto& m : nlte) m.set_atm(atm_);
  atm = std::move(atm_);
}

std::pair<Complex, Complex> line_storage::operator()(
    const Numeric f, const zeeman::pol pol) const {
  std::array res{lte[static_cast<Size>(pol)](f),
                 lte_mirror[static_cast<Size>(pol)](f),
                 nlte[static_cast<Size>(pol)](f)};

  return std::reduce(res.begin(),
                     res.end(),
                     std::pair<Complex, Complex>{0.0, 0.0},
                     [](const auto& a, const auto& b) {
                       return std::pair{a.first + b.first, a.second + b.second};
                     });
}
}  // namespace lbl::fwd