#include "spectral_radiance_transform_operator.h"

#include <arts_constexpr_math.h>
#include <debug.h>
#include <physics_funcs.h>

namespace {
struct spectral_unit_op {
  void operator()(StokvecVector& iy,
                  StokvecMatrix& diy,
                  const AscendingGrid&,
                  const PropagationPathPoint& point) const {
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
                  const PropagationPathPoint&) const {
    ARTS_USER_ERROR_IF(
        iy.size() != freqs.size() or
            static_cast<Size>(diy.ncols()) != freqs.size(),
        R"(Mismatch in size of spectral radiance, spectral radiance jacobian, and frequency grid

freq_grid.size() = {}
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
                  const PropagationPathPoint&) const {
    ARTS_USER_ERROR_IF(
        iy.size() != freqs.size() or
            static_cast<Size>(diy.ncols()) != freqs.size(),
        R"(Mismatch in size of spectral radiance, spectral radiance jacobian, and frequency grid

freq_grid.size() = {}
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
                  const PropagationPathPoint& point) const {
    ARTS_USER_ERROR_IF(
        iy.size() != freqs.size() or
            static_cast<Size>(diy.ncols()) != freqs.size(),
        R"(Mismatch in size of spectral radiance, spectral radiance jacobian, and frequency grid

freq_grid.size() = {}
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
                  const PropagationPathPoint& point) const {
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
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const PropagationPathPoint& ray_path_point) const {
  f(spectral_rad, spectral_rad_jac, freq_grid, ray_path_point);
}

void xml_io_stream<SpectralRadianceTransformOperator>::write(
    std::ostream& os,
    const SpectralRadianceTransformOperator& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.type, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<SpectralRadianceTransformOperator>::read(
    std::istream& is, SpectralRadianceTransformOperator& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  x = SpectralRadianceTransformOperator(x.type);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
