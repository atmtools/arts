#pragma once

#include <atm.h>

#include <memory>

#include "lbl_data.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_lineshape_voigt_lte_mirrored.h"
#include "lbl_lineshape_voigt_nlte.h"

namespace lbl::fwd {
namespace models {
class lte {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};
  ZeemanPolarization pol{};

  voigt::lte::band_shape lines{};

  voigt::lte::band_shape cutoff_lines{};
  ComplexVector cutoff;

  void adapt();

 public:
  lte();
  lte(const lte&);
  lte(lte&&) noexcept;
  lte& operator=(const lte&);
  lte& operator=(lte&&) noexcept;

  std::pair<Complex, Complex> operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(ZeemanPolarization pol);
  void set(std::shared_ptr<AbsorptionBands> bands,
           std::shared_ptr<AtmPoint> atm,
           ZeemanPolarization pol);
};

class lte_mirror {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};
  ZeemanPolarization pol{};

  voigt::lte_mirror::band_shape lines{};

  voigt::lte_mirror::band_shape cutoff_lines{};
  ComplexVector cutoff;

  void adapt();

 public:
  lte_mirror();
  lte_mirror(const lte_mirror&);
  lte_mirror(lte_mirror&&) noexcept;
  lte_mirror& operator=(const lte_mirror&);
  lte_mirror& operator=(lte_mirror&&) noexcept;

  std::pair<Complex, Complex> operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(ZeemanPolarization pol);
  void set(std::shared_ptr<AbsorptionBands> bands,
           std::shared_ptr<AtmPoint> atm,
           ZeemanPolarization pol);
};

class nlte {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};
  ZeemanPolarization pol{};

  voigt::nlte::band_shape lines{};

  voigt::nlte::band_shape cutoff_lines{};
  matpack::data_t<std::pair<Complex, Complex>, 1> cutoff;

  void adapt();

 public:
  nlte();
  nlte(const nlte&);
  nlte(nlte&&) noexcept;
  nlte& operator=(const nlte&);
  nlte& operator=(nlte&&) noexcept;

  std::pair<Complex, Complex> operator()(const Numeric frequency) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(ZeemanPolarization pol);
  void set(std::shared_ptr<AbsorptionBands> bands,
           std::shared_ptr<AtmPoint> atm,
           ZeemanPolarization pol);
};
}  // namespace models

class line_storage {
  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<AbsorptionBands> bands{};

  std::array<models::lte, 4> lte{};
  std::array<models::lte_mirror, 4> lte_mirror{};
  std::array<models::nlte, 4> nlte{};

 public:
  line_storage();
  line_storage(const line_storage&);
  line_storage(line_storage&&) noexcept;
  line_storage& operator=(const line_storage&);
  line_storage& operator=(line_storage&&) noexcept;

  line_storage(std::shared_ptr<AtmPoint> atm,
               std::shared_ptr<AbsorptionBands> bands);

  std::pair<Complex, Complex> operator()(const Numeric frequency,
                                         const ZeemanPolarization pol) const;

  void set_model(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};  // struct frequency
}  // namespace lbl::fwd
