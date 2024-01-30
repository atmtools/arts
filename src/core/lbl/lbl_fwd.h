#pragma once

#include <memory>

#include "atm.h"
#include "lbl_data.h"
#include "lbl_lineshape_voigt_lte.h"
#include "lbl_lineshape_voigt_lte_mirrored.h"
#include "lbl_lineshape_voigt_nlte.h"

namespace lbl::fwd {
class line_storage {
  std::shared_ptr<AbsorptionBands> bands{nullptr};
  std::shared_ptr<AtmPoint> atm{nullptr};
  zeeman::pol pol{};

  void adapt();

public:
  std::vector<voigt::lte::single_shape> lte_shapes{};
  std::vector<voigt::lte::single_shape> cutoff_lte_shapes{};

  std::vector<voigt::lte_mirror::single_shape> lte_mirror_shapes{};
  std::vector<voigt::lte_mirror::single_shape> cutoff_lte_mirror_shapes{};

  std::vector<voigt::nlte::single_shape> nlte_shapes{};
  std::vector<voigt::nlte::single_shape> cutoff_nlte_shapes{};

  line_storage() = default;
  line_storage(const line_storage&) = default;
  line_storage(line_storage&&) = default;
  line_storage& operator=(const line_storage&) = default;
  line_storage& operator=(line_storage&&) = default;

  line_storage(std::shared_ptr<AbsorptionBands> bands,
               std::shared_ptr<AtmPoint> atm,
               const zeeman::pol pol=zeeman::pol::no);

  std::pair<Complex, Complex> operator()(const Numeric frequency) const;
  
  void set_lines(std::shared_ptr<AbsorptionBands> bands);
  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(zeeman::pol pol);
  [[nodiscard]] Propmat polarization(Vector2 los) const;
};  // struct frequency
}  // namespace lbl::fwd