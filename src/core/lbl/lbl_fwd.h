#pragma once

#include <memory>

#include "atm.h"
#include "lbl_data.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_lineshape_voigt_lte_mirrored.h"
#include "lbl_lineshape_voigt_nlte.h"

namespace lbl::fwd {
namespace models {
class lte {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};
  zeeman::pol pol{};

  voigt::lte::band_shape lines{};

  voigt::lte::band_shape cutoff_lines{};
  ComplexVector cutoff;

  void adapt();

 public:
  std::pair<Complex, Complex> operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(zeeman::pol pol);
  void set(std::shared_ptr<AbsorptionBands> bands, std::shared_ptr<AtmPoint> atm,
  zeeman::pol pol);
};

class lte_mirror {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};
  zeeman::pol pol{};

  voigt::lte_mirror::band_shape lines{};

  voigt::lte_mirror::band_shape cutoff_lines{};
  ComplexVector cutoff;

  void adapt();

 public:
  std::pair<Complex, Complex> operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(zeeman::pol pol);
  void set(std::shared_ptr<AbsorptionBands> bands, std::shared_ptr<AtmPoint> atm,
  zeeman::pol pol);
};

class nlte {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};
  zeeman::pol pol{};

  voigt::nlte::band_shape lines{};

  voigt::nlte::band_shape cutoff_lines{};
  matpack::matpack_data<std::pair<Complex, Complex>, 1> cutoff;

  void adapt();

 public:
  std::pair<Complex, Complex> operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(zeeman::pol pol);
  void set(std::shared_ptr<AbsorptionBands> bands, std::shared_ptr<AtmPoint> atm,
  zeeman::pol pol);
};
}  // namespace models

class line_storage {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};

  std::array<models::lte, 4> lte{};
  std::array<models::lte_mirror, 4> lte_mirror{};
  std::array<models::nlte, 4> nlte{};

 public:
  line_storage() = default;
  line_storage(const line_storage&) = default;
  line_storage(line_storage&&) = default;
  line_storage& operator=(const line_storage&) = default;
  line_storage& operator=(line_storage&&) = default;

  line_storage(std::shared_ptr<AtmPoint> atm,
               std::shared_ptr<AbsorptionBands> bands);

  std::pair<Complex, Complex> operator()(const Numeric frequency, const zeeman::pol pol) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};  // struct frequency
}  // namespace lbl::fwd
