#include "spectral_radiance_transform_operator.h"

#include <debug.h>
#include <physics_funcs.h>

#include "arts_constexpr_math.h"

namespace {
struct spectral_unit_op {
  void operator()(StokvecVector& iy,
                         StokvecMatrix& diy,
                         const AscendingGrid&,
                         const PropagationPathPoint& point) {
    if (point.nreal != 1.0) {
      diy *= Math::pow2(point.nreal);
      iy  *= Math::pow2(point.nreal);
    }
  }
};

struct spectral_rjbt_op {
  void operator()(StokvecVector& iy,
                         StokvecMatrix& diy,
                         const AscendingGrid& freqs,
                         const PropagationPathPoint&) {
    ARTS_USER_ERROR_IF(
        iy.size() != freqs.size() or
            static_cast<Size>(diy.ncols()) != freqs.size(),
        R"(Mismatch in size of spectral radiance, spectral radiance jacobian, and frequency grid

frequency_grid.size() = {}
iy.size()             = {}
diy.shape()           = {:B} [column size must match frequency grid size]
)",
        freqs.size(),
        iy.size(),
        diy.shape());

    for (Size j = 0; j < freqs.size(); j++) {
      const Numeric df  = invrayjean(1.0, freqs[j]);
      iy[j]            *= df;
      diy[joker, j]    *= df;
    }
  }
};

struct spectral_planck_op {
  void operator()(StokvecVector& iy,
                         StokvecMatrix& diy,
                         const AscendingGrid& freqs,
                         const PropagationPathPoint&) {
    ARTS_USER_ERROR_IF(
        iy.size() != freqs.size() or
            static_cast<Size>(diy.ncols()) != freqs.size(),
        R"(Mismatch in size of spectral radiance, spectral radiance jacobian, and frequency grid

frequency_grid.size() = {}
iy.size()             = {}
diy.shape()           = {:B} [column size must match frequency grid size]
)",
        freqs.size(),
        iy.size(),
        diy.shape());

    for (Size j = 0; j < freqs.size(); j++) {
      auto& v  = iy[j];
      auto&& f = freqs[j];
      const Stokvec deriv{dinvplanckdI(v.I(), f),
                          dinvplanckdI(0.5 * (v.I() + v.Q()), f) -
                              dinvplanckdI(0.5 * (v.I() - v.Q()), f),
                          dinvplanckdI(0.5 * (v.I() + v.U()), f) -
                              dinvplanckdI(0.5 * (v.I() - v.U()), f),
                          dinvplanckdI(0.5 * (v.I() + v.V()), f) -
                              dinvplanckdI(0.5 * (v.I() - v.V()), f)};
      const Stokvec normal{invplanck(v.I(), f),
                           invplanck(0.5 * (v.I() + v.Q()), f) -
                               invplanck(0.5 * (v.I() - v.Q()), f),
                           invplanck(0.5 * (v.I() + v.U()), f) -
                               invplanck(0.5 * (v.I() - v.U()), f),
                           invplanck(0.5 * (v.I() + v.V()), f) -
                               invplanck(0.5 * (v.I() - v.V()), f)};

      v = normal;
      for (Index i = 0; i < diy.nrows(); i++) diy[i, j] *= deriv;
    }
  }
};

struct spectral_W_m2_m_sr_op {
  void operator()(StokvecVector& iy,
                         StokvecMatrix& diy,
                         const AscendingGrid& freqs,
                         const PropagationPathPoint& point) {
    ARTS_USER_ERROR_IF(
        iy.size() != freqs.size() or
            static_cast<Size>(diy.ncols()) != freqs.size(),
        R"(Mismatch in size of spectral radiance, spectral radiance jacobian, and frequency grid

frequency_grid.size() = {}
iy.size()             = {}
diy.shape()           = {:B} [column size must match frequency grid size]
)",
        freqs.size(),
        iy.size(),
        diy.shape());

    for (Size j = 0; j < freqs.size(); j++) {
      const Numeric df  = (freqs[j] * (freqs[j] / Constant::c));
      diy[joker, j]    *= df * Math::pow2(point.nreal);
      iy[j]            *= df * Math::pow2(point.nreal);
    }
  }
};

struct spectral_W_m2_m1_sr_op {
  void operator()(StokvecVector& iy,
                         StokvecMatrix& diy,
                         const AscendingGrid&,
                         const PropagationPathPoint& point) {
    diy *= Math::pow2(point.nreal) * Constant::c;
    iy  *= Math::pow2(point.nreal) * Constant::c;
  }
};
}  // namespace

SpectralRadianceTransformOperator::SpectralRadianceTransformOperator() =
    default;
SpectralRadianceTransformOperator::SpectralRadianceTransformOperator(
    const SpectralRadianceTransformOperator& x) = default;

SpectralRadianceTransformOperator::SpectralRadianceTransformOperator(
    SpectralRadianceTransformOperator&& x) noexcept = default;

SpectralRadianceTransformOperator& SpectralRadianceTransformOperator::operator=(
    const SpectralRadianceTransformOperator& x) = default;

SpectralRadianceTransformOperator& SpectralRadianceTransformOperator::operator=(
    SpectralRadianceTransformOperator&& x) noexcept = default;

SpectralRadianceTransformOperator::SpectralRadianceTransformOperator(
    SpectralRadianceTransformOperator::Op&& x) noexcept
    : f(std::move(x)) {}

SpectralRadianceTransformOperator::SpectralRadianceTransformOperator(
    const SpectralRadianceTransformOperator::Op& x)
    : f(x) {}

SpectralRadianceTransformOperator::Op from(SpectralRadianceUnitType x) {
  switch (x) {
    case SpectralRadianceUnitType::unit:
      return SpectralRadianceTransformOperator::Op{spectral_unit_op{}};
    case SpectralRadianceUnitType::RJBT:
      return SpectralRadianceTransformOperator::Op{spectral_rjbt_op{}};
    case SpectralRadianceUnitType::PlanckBT:
      return SpectralRadianceTransformOperator::Op{spectral_planck_op{}};
    case SpectralRadianceUnitType::W_m2_m_sr:
      return SpectralRadianceTransformOperator::Op{spectral_W_m2_m_sr_op{}};
    case SpectralRadianceUnitType::W_m2_m1_sr:
      return SpectralRadianceTransformOperator::Op{spectral_W_m2_m1_sr_op{}};
  }

  ARTS_USER_ERROR("Unsupported SpectralRadianceUnitType: {}", x);
  std::unreachable();
}

SpectralRadianceTransformOperator::SpectralRadianceTransformOperator(
    SpectralRadianceUnitType x)
    : f(from(x)), type(toString(x)) {}

SpectralRadianceTransformOperator::SpectralRadianceTransformOperator(
    const std::string_view& x)
    : SpectralRadianceTransformOperator(to<SpectralRadianceUnitType>(x)) {}

void SpectralRadianceTransformOperator::operator()(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point) const {
  f(spectral_radiance,
    spectral_radiance_jacobian,
    frequency_grid,
    ray_path_point);
}