#include "phase_matrix.h"

#include <arts_constants.h>

#include <algorithm>
#include <map>
#include <ranges>

namespace scattering {

namespace detail {
namespace {
bool equal(Numeric a, Numeric b, Numeric epsilon = 1e-6) {
  return std::abs(a - b) <= ((std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b)) * epsilon);
}

bool small(Numeric a, Numeric epsilon = 1e-6) { return std::abs(a) < epsilon; }

/** acos of a cosine that rounding may have put just outside [-1, 1]
 *
 * Only clamps: snapping cosines within a relative 1e-6 of +-1 to 0 or pi, as
 * this did, moves angles by up to 1.4e-3 rad and makes the laboratory-frame
 * phase matrix jump there.
 */
Numeric save_acos(Numeric a) { return std::acos(std::clamp<Numeric>(a, -1.0, 1.0)); }
}  // namespace

std::array<Numeric, 5> rotation_coefficients(Numeric aa_inc_d, Numeric za_inc_d, Numeric aa_scat_d, Numeric za_scat_d) {
  Numeric aa_inc  = Conversion::deg2rad(aa_inc_d);
  Numeric za_inc  = Conversion::deg2rad(za_inc_d);
  Numeric aa_scat = Conversion::deg2rad(aa_scat_d);
  Numeric za_scat = Conversion::deg2rad(za_scat_d);

  Numeric cos_theta = cos(za_inc) * cos(za_scat) + sin(za_inc) * sin(za_scat) * cos(aa_scat - aa_inc);
  Numeric theta     = save_acos(cos_theta);
  if ((small(std::abs(aa_scat - aa_inc))) || (equal(std::abs(aa_scat - aa_inc), Constant::two_pi))) {
    theta = std::abs(za_inc - za_scat);
  } else if ((equal(aa_scat - aa_inc, Constant::two_pi))) {
    theta = za_scat + za_inc;
    if (theta > Constant::pi) { theta = Constant::two_pi - theta; }
  }
  const Numeric delta_aa = std::abs(std::remainder(static_cast<Numeric>(aa_scat - aa_inc), Constant::two_pi));
  const auto    ray      = [](Numeric za, Numeric aa) {
    const Numeric s = std::sin(za);
    return rtepack::meridional_basis(Vector3{s * std::sin(aa), s * std::cos(aa), std::cos(za)},
                                     Vector3{std::cos(aa), -std::sin(aa), 0.0});
  };
  const auto [c_1, s_1, c_2, s_2] = rtepack::scattering_plane_rotation(ray(static_cast<Numeric>(za_inc), 0.0),
                                                                       ray(static_cast<Numeric>(za_scat), delta_aa));

  return {theta,
          static_cast<Numeric>(c_1),
          static_cast<Numeric>(c_2),
          static_cast<Numeric>(s_1),
          static_cast<Numeric>(s_2)};
}
}  // namespace detail

ScatteringDataGrids::ScatteringDataGrids(std::shared_ptr<const Vector> t_grid_, std::shared_ptr<const Vector> f_grid_)
    : t_grid(std::move(t_grid_)),
      f_grid(std::move(f_grid_)),
      aa_inc_grid(nullptr),
      za_inc_grid(nullptr),
      aa_scat_grid(nullptr),
      za_scat_grid(nullptr) {}

ScatteringDataGrids::ScatteringDataGrids(std::shared_ptr<const Vector>          t_grid_,
                                         std::shared_ptr<const Vector>          f_grid_,
                                         std::shared_ptr<const ZenithAngleGrid> za_scat_grid_)
    : t_grid(std::move(t_grid_)),
      f_grid(std::move(f_grid_)),
      aa_inc_grid(nullptr),
      za_inc_grid(nullptr),
      aa_scat_grid(nullptr),
      za_scat_grid(std::move(za_scat_grid_)) {}

ScatteringDataGrids::ScatteringDataGrids(std::shared_ptr<const Vector>          t_grid_,
                                         std::shared_ptr<const Vector>          f_grid_,
                                         std::shared_ptr<const Vector>          za_inc_grid_,
                                         std::shared_ptr<const Vector>          delta_aa_grid_,
                                         std::shared_ptr<const ZenithAngleGrid> za_scat_grid_)
    : t_grid(std::move(t_grid_)),
      f_grid(std::move(f_grid_)),
      aa_inc_grid(nullptr),
      za_inc_grid(std::move(za_inc_grid_)),
      aa_scat_grid(std::move(delta_aa_grid_)),
      za_scat_grid(std::move(za_scat_grid_)) {}

RegridWeights calc_regrid_weights(std::shared_ptr<const Vector>          t_grid,
                                  std::shared_ptr<const Vector>          f_grid,
                                  std::shared_ptr<const Vector>          aa_inc_grid,
                                  std::shared_ptr<const Vector>          za_inc_grid,
                                  std::shared_ptr<const Vector>          aa_scat_grid,
                                  std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                  ScatteringDataGrids                    new_grids) {
  RegridWeights res{};

  if (!t_grid) { ARTS_USER_ERROR("The old t_grid must be provided for calculating regridding weights."); }
  if (!new_grids.t_grid) { ARTS_USER_ERROR("The new t_grid must be provided for calculating regridding weights."); }
  if (!f_grid) { ARTS_USER_ERROR("The old f_grid must be provided for calculating regridding weights."); }
  if (!new_grids.f_grid) { ARTS_USER_ERROR("The new f_grid must be provided for calculating regridding weights."); }

  const auto positions = [](ArrayOfGridPos &out, const Vector &old_grid, const Vector &new_grid) {
    if (old_grid.size() == 1) {
      stdr::fill(out, GridPos{.idx = 0, .fd = {0.0, 1.0}});
    } else {
      gridpos(out, old_grid, new_grid, 1e99);
    }
  };

  res.t_grid_weights = ArrayOfGridPos(new_grids.t_grid->size());
  positions(res.t_grid_weights, *t_grid, *new_grids.t_grid);
  res.f_grid_weights = ArrayOfGridPos(new_grids.f_grid->size());
  positions(res.f_grid_weights, *f_grid, *new_grids.f_grid);

  if ((aa_inc_grid) && (new_grids.aa_inc_grid)) {
    res.aa_inc_grid_weights = ArrayOfGridPos(new_grids.aa_inc_grid->size());
    positions(res.aa_inc_grid_weights, *aa_inc_grid, *new_grids.aa_inc_grid);
  }
  if ((za_inc_grid) && (new_grids.za_inc_grid)) {
    res.za_inc_grid_weights = ArrayOfGridPos(new_grids.za_inc_grid->size());
    positions(res.za_inc_grid_weights, *za_inc_grid, *new_grids.za_inc_grid);
  }
  if ((aa_scat_grid) && (new_grids.aa_scat_grid)) {
    res.aa_scat_grid_weights = ArrayOfGridPos(new_grids.aa_scat_grid->size());
    positions(res.aa_scat_grid_weights, *aa_scat_grid, *new_grids.aa_scat_grid);
  }
  if ((za_scat_grid) && (new_grids.za_scat_grid)) {
    res.za_scat_grid_weights =
        ArrayOfGridPos(std::visit([](const auto &grd) { return grd.angles.size(); }, *new_grids.za_scat_grid));
    positions(res.za_scat_grid_weights,
              std::visit([](const auto &grd) { return grd.angles.vec(); }, *za_scat_grid),
              std::visit([](const auto &grd) { return grd.angles.vec(); }, *new_grids.za_scat_grid));
  }
  return res;
}

BackscatterMatrixData<Format::TRO> &BackscatterMatrixData<Format::TRO>::operator=(const TensorType &data) {
  ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                     "Provided backscatter coefficient data do not match temperature grid.");
  ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided backscatter coefficient data do not match frequency grid.");
  TensorType::operator=(data);
  return *this;
}

BackscatterMatrixData<Format::TRO> BackscatterMatrixData<Format::TRO>::regrid(const ScatteringDataGrids &grids,
                                                                              const RegridWeights &weights) const {
  BackscatterMatrixData result(grids.t_grid, grids.f_grid);
  const TensorType     &coeffs_this = *this;
  TensorType           &coeffs_res  = result;
  for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
    GridPos    gp_t  = weights.t_grid_weights[i_t];
    Numeric    w_t_l = gp_t.fd[1];
    Numeric    w_t_r = gp_t.fd[0];
    const auto i_t_l = std::clamp<Index>(gp_t.idx, 0, n_temps_ - 1);
    const auto i_t_r = std::min(i_t_l + 1, n_temps_ - 1);
    for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
      GridPos    gp_f      = weights.f_grid_weights[i_f];
      Numeric    w_f_l     = gp_f.fd[1];
      Numeric    w_f_r     = gp_f.fd[0];
      const auto i_f_l     = std::clamp<Index>(gp_f.idx, 0, n_freqs_ - 1);
      const auto i_f_r     = std::min(i_f_l + 1, n_freqs_ - 1);
      coeffs_res[i_t, i_f] = (w_t_l * w_f_l * coeffs_this[i_t_l, i_f_l] + w_t_l * w_f_r * coeffs_this[i_t_l, i_f_r] +
                              w_t_r * w_f_l * coeffs_this[i_t_r, i_f_l] + w_t_r * w_f_r * coeffs_this[i_t_r, i_f_r]);
    }
  }
  return result;
}

BackscatterMatrixData<Format::TRO> BackscatterMatrixData<Format::TRO>::regrid(const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
  return regrid(grids, weights);
}

BackscatterMatrixData<Format::ARO>::BackscatterMatrixData(const BackscatterMatrixData<Format::TRO> &bsmat,
                                                          std::shared_ptr<const Vector>             za_inc_grid)
    : BackscatterMatrixData(bsmat.get_t_grid_ptr(), bsmat.get_f_grid_ptr(), za_inc_grid) {
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) this->operator[](i_t, i_f, joker) = bsmat[i_t, i_f].expand();
}

BackscatterMatrixData<Format::ARO> &BackscatterMatrixData<Format::ARO>::operator=(const TensorType &data) {
  ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                     "Provided backscatter coefficient data do not match temperature grid.");
  ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided backscatter coefficient data do not match frequency grid.");
  ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                     "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
  TensorType::operator=(data);
  return *this;
}

BackscatterMatrixData<Format::ARO> BackscatterMatrixData<Format::ARO>::regrid(const ScatteringDataGrids &grids,
                                                                              const RegridWeights &weights) const {
  BackscatterMatrixData result(grids.t_grid, grids.f_grid, grids.za_inc_grid);
  const TensorType     &coeffs_this = *this;
  TensorType           &coeffs_res  = result;
  for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
    GridPos     gp_t    = weights.t_grid_weights[i_t];
    const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
    Numeric     w_t_l   = gp_t.fd[1];
    Numeric     w_t_r   = gp_t.fd[0];
    for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
      GridPos     gp_f    = weights.f_grid_weights[i_f];
      const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
      Numeric     w_f_l   = gp_f.fd[1];
      Numeric     w_f_r   = gp_f.fd[0];
      for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size(); ++i_za_inc) {
        GridPos     gp_za_inc          = weights.za_inc_grid_weights[i_za_inc];
        const Index za_inc_upper       = std::min<Index>(gp_za_inc.idx + 1, coeffs_this.extent(2) - 1);
        Numeric     w_za_inc_l         = gp_za_inc.fd[1];
        Numeric     w_za_inc_r         = gp_za_inc.fd[0];
        coeffs_res[i_t, i_f, i_za_inc] = (w_t_l * w_f_l * w_za_inc_l * coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx] +
                                          w_t_l * w_f_l * w_za_inc_r * coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper] +
                                          w_t_l * w_f_r * w_za_inc_l * coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx] +
                                          w_t_l * w_f_r * w_za_inc_r * coeffs_this[gp_t.idx, f_upper, za_inc_upper] +

                                          w_t_r * w_f_l * w_za_inc_l * coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx] +
                                          w_t_r * w_f_l * w_za_inc_r * coeffs_this[t_upper, gp_f.idx, za_inc_upper] +
                                          w_t_r * w_f_r * w_za_inc_l * coeffs_this[t_upper, f_upper, gp_za_inc.idx] +
                                          w_t_r * w_f_r * w_za_inc_r * coeffs_this[t_upper, f_upper, za_inc_upper]);
      }
    }
  }
  return result;
}

BackscatterMatrixData<Format::ARO> BackscatterMatrixData<Format::ARO>::regrid(const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, nullptr, nullptr, grids);
  return regrid(grids, weights);
}

PhaseMatrixData<Format::TRO, Representation::Gridded> &PhaseMatrixData<Format::TRO, Representation::Gridded>::operator=(
    const TensorType &data) {
  ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                     "Provided backscatter coefficient data do not match temperature grid.");
  ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided backscatter coefficient data do not match frequency grid.");
  ARTS_USER_ERROR_IF(data.shape()[2] != n_za_scat_,
                     "Provided backscatter coefficient data do not match expected number of scattering zenith angles.");
  TensorType::operator=(data);
  return *this;
}

BackscatterMatrixData<Format::TRO> PhaseMatrixData<Format::TRO, Representation::Gridded>::extract_backscatter_matrix() {
  BackscatterMatrixData<Format::TRO> result(t_grid_, f_grid_);
  GridPos                            interp = find_interp_weights(grid_vector(*za_scat_grid_), 180.0);
  //gridpos(interp, grid_vector(*za_scat_grid_), 180.0);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      result[i_t, i_f] =
          interp.fd[1] * (*this)[i_t, i_f, interp.idx] + interp.fd[0] * (*this)[i_t, i_f, interp.idx + 1];
    }
  }
  return result;
}

ForwardscatterMatrixData<Format::TRO>
PhaseMatrixData<Format::TRO, Representation::Gridded>::extract_forwardscatter_matrix() {
  ForwardscatterMatrixData<Format::TRO> result(t_grid_, f_grid_);
  GridPos                               interp;
  gridpos(interp, grid_vector(*za_scat_grid_), 0.0, 1e99);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      result[i_t, i_f] =
          interp.fd[1] * (*this)[i_t, i_f, interp.idx] + interp.fd[0] * (*this)[i_t, i_f, interp.idx + 1];
    }
  }
  return result;
}

rtepack::compact_planar_muelmat_matrix PhaseMatrixData<Format::TRO, Representation::Gridded>::integrate_phase_matrix() {
  rtepack::compact_planar_muelmat_matrix results(n_temps_, n_freqs_);
  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      results[i_t, i_f] = Constant::two_pi * integrate_zenith_angle(this->operator[](i_t, i_f, joker), *za_scat_grid_);
    }
  }
  return results;
}

PhaseMatrixData<Format::TRO, Representation::Gridded> PhaseMatrixData<Format::TRO, Representation::Gridded>::regrid(
    const ScatteringDataGrids &grids, const RegridWeights &weights) const {
  PhaseMatrixData   result(grids.t_grid, grids.f_grid, grids.za_scat_grid);
  const TensorType &coeffs_this = *this;
  TensorType       &coeffs_res  = result;
  for (Index i_t = 0; i_t < static_cast<Index>(weights.t_grid_weights.size()); ++i_t) {
    GridPos     gp_t  = weights.t_grid_weights[i_t];
    Numeric     w_t_l = gp_t.fd[1];
    Numeric     w_t_r = gp_t.fd[0];
    const Index it0   = std::clamp<Index>(gp_t.idx, 0, n_temps_ - 1);
    const Index it1   = std::min(it0 + 1, n_temps_ - 1);
    for (Index i_f = 0; i_f < static_cast<Index>(weights.f_grid_weights.size()); ++i_f) {
      GridPos     gp_f  = weights.f_grid_weights[i_f];
      Numeric     w_f_l = gp_f.fd[1];
      Numeric     w_f_r = gp_f.fd[0];
      const Index if0   = std::clamp<Index>(gp_f.idx, 0, n_freqs_ - 1);
      const Index if1   = std::min(if0 + 1, n_freqs_ - 1);
      for (Index i_za_scat = 0; i_za_scat < static_cast<Index>(weights.za_scat_grid_weights.size()); ++i_za_scat) {
        GridPos     gp_za_scat  = weights.za_scat_grid_weights[i_za_scat];
        Numeric     w_za_scat_l = gp_za_scat.fd[1];
        Numeric     w_za_scat_r = gp_za_scat.fd[0];
        const Index iza0        = std::clamp<Index>(gp_za_scat.idx, 0, n_za_scat_ - 1);
        const Index iza1        = std::min(iza0 + 1, n_za_scat_ - 1);

        coeffs_res[i_t, i_f, i_za_scat] = rtepack::compact_planar_muelmat{};

        if (w_t_l > 0.0) {
          if (w_f_l > 0.0) {
            coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_l * w_za_scat_l * coeffs_this[it0, if0, iza0];
            coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_l * w_za_scat_r * coeffs_this[it0, if0, iza1];
          }
          if (w_f_r > 0.0) {
            coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_r * w_za_scat_l * coeffs_this[it0, if1, iza0];
            coeffs_res[i_t, i_f, i_za_scat] += w_t_l * w_f_r * w_za_scat_r * coeffs_this[it0, if1, iza1];
          }
        }
        if (w_t_r > 0.0) {
          if (w_f_l > 0.0) {
            coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_l * w_za_scat_l * coeffs_this[it1, if0, iza0];
            coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_l * w_za_scat_r * coeffs_this[it1, if0, iza1];
          }
          if (w_f_r > 0.0) {
            coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_r * w_za_scat_l * coeffs_this[it1, if1, iza0];
            coeffs_res[i_t, i_f, i_za_scat] += w_t_r * w_f_r * w_za_scat_r * coeffs_this[it1, if1, iza1];
          }
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::TRO, Representation::Gridded> PhaseMatrixData<Format::TRO, Representation::Gridded>::regrid(
    const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, za_scat_grid_, grids);
  return regrid(grids, weights);
}

auto PhaseMatrixData<Format::TRO, Representation::Gridded>::to_spectral(Index degree, Index order) const
    -> PhaseMatrixDataSpectral {
  ARTS_USER_ERROR_IF(order != 0,
                     "A TRO phase matrix depends on the scattering angle only, so its spectral form has order 0; "
                     "got order {}",
                     order)
  PhaseMatrixDataSpectral result(t_grid_, f_grid_, degree);
  const Vector            angles{grid_vector(*za_scat_grid_)};
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f)
      result[i_t, i_f, joker] = tro_legendre::project(angles, (*this)[i_t, i_f, joker], degree);
  return result;
}

LegendreReport PhaseMatrixData<Format::TRO, Representation::Gridded>::legendre_report(
    const PhaseMatrixDataSpectral &spectral) const {
  ARTS_USER_ERROR_IF(spectral.extent(0) != n_temps_ or spectral.extent(1) != n_freqs_,
                     "The Legendre series has {} temperatures and {} frequencies, but the gridded phase matrix has "
                     "{} and {}",
                     spectral.extent(0),
                     spectral.extent(1),
                     n_temps_,
                     n_freqs_)
  LegendreReport report(n_temps_, n_freqs_);
  const Vector   angles{grid_vector(*za_scat_grid_)};
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f)
      report.set(i_t, i_f, tro_legendre::assess(spectral[i_t, i_f, joker], angles, (*this)[i_t, i_f, joker]));
  return report;
}

auto PhaseMatrixData<Format::TRO, Representation::Gridded>::to_lab_frame(
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const Vector>          delta_aa_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) const -> PhaseMatrixDataLabFrame {
  const auto source_angles = grid_vector(*za_scat_grid_);
  const auto interpolate   = [&](Numeric theta, rtepack::compact_planar_muelmat_matrix &scattering_matrix) {
    const Numeric scat_angle = Conversion::rad2deg(theta);
    Index         angle0;
    Index         angle1;
    Numeric       weight0;
    Numeric       weight1;
    if (n_za_scat_ == 1 or scat_angle <= source_angles.front()) {
      angle0 = angle1 = 0;
      weight0         = 1.0;
      weight1         = 0.0;
    } else if (scat_angle >= source_angles.back()) {
      angle0 = angle1 = n_za_scat_ - 1;
      weight0         = 1.0;
      weight1         = 0.0;
    } else {
      const auto upper = std::ranges::upper_bound(source_angles, scat_angle);
      angle1           = static_cast<Index>(upper - source_angles.begin());
      angle0           = angle1 - 1;
      weight1          = (scat_angle - source_angles[angle0]) / (source_angles[angle1] - source_angles[angle0]);
      weight0          = 1.0 - weight1;
    }

    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        rtepack::compact_planar_muelmat value{};
        if (weight0 > 0.0) value += weight0 * (*this)[i_t, i_f, angle0];
        if (weight1 > 0.0) value += weight1 * (*this)[i_t, i_f, angle1];
        scattering_matrix[i_t, i_f] = value;
      }
    }
  };
  return tro_lab_frame(t_grid_, f_grid_, za_inc_grid, delta_aa_grid, za_scat_grid_new, interpolate);
}

PhaseMatrixData<Format::TRO, Representation::Spectral> &
PhaseMatrixData<Format::TRO, Representation::Spectral>::operator=(const TensorType &data) {
  ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_, "Provided phase matrix coefficients do not match temperature grid.");
  ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided phase matrix coefficients do not match frequency grid.");
  ARTS_USER_ERROR_IF(data.shape()[2] != degree_ + 1,
                     "Provided phase matrix coefficients do not match the Legendre degree {}.",
                     degree_);
  TensorType::operator=(data);
  return *this;
}

auto PhaseMatrixData<Format::TRO, Representation::Spectral>::to_gridded(
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid) const -> PhaseMatrixDataGridded {
  PhaseMatrixDataGridded result(t_grid_, f_grid_, za_scat_grid);
  const Vector           angles{grid_vector(*za_scat_grid)};
  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      result[i_t, i_f, joker] = tro_legendre::evaluate((*this)[i_t, i_f, joker], angles);
    }
  }
  return result;
}

ScatteringMatrixFunction PhaseMatrixData<Format::TRO, Representation::Spectral>::scattering_matrix_function() const {
  Vector norm(degree_ + 1);
  for (Index l = 0; l <= degree_; ++l) norm[l] = std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * Constant::pi));
  return [this, p = Vector(degree_ + 1), norm = std::move(norm)](
             Numeric theta, rtepack::compact_planar_muelmat_matrix &scattering_matrix) mutable {
    Legendre::legendre_polynomials(p, std::clamp<Numeric>(std::cos(theta), -1.0, 1.0));
    for (Index l = 0; l <= degree_; ++l) p[l] *= norm[l];
    for (Index i_t = 0; i_t < n_temps_; ++i_t) {
      for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
        rtepack::compact_planar_muelmat value{};
        for (Index l = 0; l <= degree_; ++l) value += (*this)[i_t, i_f, l] * p[l];
        scattering_matrix[i_t, i_f] = value;
      }
    }
  };
}

rtepack::compact_planar_muelmat_matrix PhaseMatrixData<Format::TRO, Representation::Spectral>::integrate_phase_matrix()
    const {
  rtepack::compact_planar_muelmat_matrix results(n_temps_, n_freqs_);
  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      results[i_t, i_f] = this->operator[](i_t, i_f, 0) * sqrt(4.0 * Constant::pi);
    }
  }
  return results;
}

PhaseMatrixData<Format::TRO, Representation::Spectral>
PhaseMatrixData<Format::TRO, Representation::Spectral>::to_spectral(Index degree, Index order) const {
  ARTS_USER_ERROR_IF(order != 0,
                     "A TRO phase matrix depends on the scattering angle only, so its spectral form has order 0; "
                     "got order {}",
                     order)
  ARTS_USER_ERROR_IF(degree > degree_,
                     "The Legendre series has degree {}, so it cannot give the coefficients to degree {}; convert "
                     "the gridded data to that degree instead",
                     degree_,
                     degree)
  PhaseMatrixData out(t_grid_, f_grid_, degree);
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f)
      for (Index l = 0; l <= degree; ++l) out[i_t, i_f, l] = this->operator[](i_t, i_f, l);
  return out;
}

PhaseMatrixData<Format::TRO, Representation::Spectral> PhaseMatrixData<Format::TRO, Representation::Spectral>::regrid(
    const ScatteringDataGrids &grids, const RegridWeights weights) const {
  PhaseMatrixData   result(grids.t_grid, grids.f_grid, degree_);
  const TensorType &coeffs_this = *this;
  TensorType       &coeffs_res  = result;
  for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
    GridPos     gp_t    = weights.t_grid_weights[i_t];
    const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
    Numeric     w_t_l   = gp_t.fd[1];
    Numeric     w_t_r   = gp_t.fd[0];
    for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
      GridPos     gp_f    = weights.f_grid_weights[i_f];
      const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
      Numeric     w_f_l   = gp_f.fd[1];
      Numeric     w_f_r   = gp_f.fd[0];

      for (Index l = 0; l <= degree_; ++l) {
        coeffs_res[i_t, i_f, l] = rtepack::compact_planar_muelmat{};
        if (w_t_l > 0.0) {
          if (w_f_l > 0.0) coeffs_res[i_t, i_f, l] += w_t_l * w_f_l * coeffs_this[gp_t.idx, gp_f.idx, l];
          if (w_f_r > 0.0) coeffs_res[i_t, i_f, l] += w_t_l * w_f_r * coeffs_this[gp_t.idx, f_upper, l];
        }
        if (w_t_r > 0.0) {
          if (w_f_l > 0.0) coeffs_res[i_t, i_f, l] += w_t_r * w_f_l * coeffs_this[t_upper, gp_f.idx, l];
          if (w_f_r > 0.0) coeffs_res[i_t, i_f, l] += w_t_r * w_f_r * coeffs_this[t_upper, f_upper, l];
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::TRO, Representation::Spectral> PhaseMatrixData<Format::TRO, Representation::Spectral>::regrid(
    const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
  return regrid(grids, weights);
}

Index PhaseMatrixData<Format::TRO, Representation::Spectral>::checked_degree(Index degree) {
  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree)
  return degree;
}

BackscatterMatrixData<Format::TRO> PhaseMatrixData<Format::TRO, Representation::Spectral>::extract_at(
    Numeric angle) const {
  BackscatterMatrixData<Format::TRO> result(t_grid_, f_grid_);
  const Vector                       angles{angle};
  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      result[i_t, i_f] = tro_legendre::evaluate((*this)[i_t, i_f, joker], angles)[0];
    }
  }
  return result;
}

auto PhaseMatrixData<Format::TRO, Representation::Spectral>::to_lab_frame(
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const Vector>          delta_aa_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid_new) const -> PhaseMatrixDataLabFrame {
  return tro_lab_frame(t_grid_, f_grid_, za_inc_grid, delta_aa_grid, za_scat_grid_new, scattering_matrix_function());
}

PhaseMatrixData<Format::ARO, Representation::Fourier>
PhaseMatrixData<Format::TRO, Representation::Spectral>::to_lab_frame_fourier_modes(
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    Index                                  max_mode) const {
  const auto integral = integrate_phase_matrix();
  Matrix     f11(n_temps_, n_freqs_);
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) f11[i_t, i_f] = integral[i_t, i_f].F11();
  return tro_lab_frame_fourier_modes(
      t_grid_, f_grid_, std::move(za_inc_grid), std::move(za_scat_grid), max_mode, f11, scattering_matrix_function());
}

PhaseMatrixData<Format::ARO, Representation::Gridded> &PhaseMatrixData<Format::ARO, Representation::Gridded>::operator=(
    const TensorType &data) {
  ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                     "Provided backscatter coefficient data do not match temperature grid.");
  ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided backscatter coefficient data do not match frequency grid.");
  ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                     "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
  ARTS_USER_ERROR_IF(
      data.shape()[3] != n_delta_aa_,
      "Provided backscatter coefficient data do not match expected number of scattering azimuth angles.");
  ARTS_USER_ERROR_IF(data.shape()[4] != n_za_scat_,
                     "Provided backscatter coefficient data do not match expected number of scattering zenith angles.");
  TensorType::operator=(data);
  return *this;
}

BackscatterMatrixData<Format::ARO> PhaseMatrixData<Format::ARO, Representation::Gridded>::extract_backscatter_matrix() {
  BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
  GridPos                            za_scat_interp, delta_aa_interp;
  gridpos(delta_aa_interp, *delta_aa_grid_, 180.0, 1e99);
  Vector     weights(4);
  const auto elements = rtepack::elements(*this);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        auto za_inc = (*za_inc_grid_)[i_za_inc];
        gridpos(za_scat_interp, grid_vector(*za_scat_grid_), 180.0 - za_inc, 1e99);
        interpweights(weights, delta_aa_interp, za_scat_interp);

        for (Index i_s = 0; i_s < 16; ++i_s) {
          auto mat                             = elements[i_t, i_f, i_za_inc, joker, joker, i_s];
          result[i_t, i_f, i_za_inc].data[i_s] = interp(weights, mat, delta_aa_interp, za_scat_interp);
        }
      }
    }
  }
  return result;
}

ForwardscatterMatrixData<Format::ARO>
PhaseMatrixData<Format::ARO, Representation::Gridded>::extract_forwardscatter_matrix() {
  BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
  GridPos                            za_scat_interp, delta_aa_interp;
  gridpos(delta_aa_interp, *delta_aa_grid_, 0.0, 1e99);
  Vector     weights(4);
  const auto elements = rtepack::elements(*this);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        auto za_inc = (*za_inc_grid_)[i_za_inc];
        gridpos(za_scat_interp, grid_vector(*za_scat_grid_), za_inc, 1e99);
        interpweights(weights, delta_aa_interp, za_scat_interp);

        for (Index i_s = 0; i_s < 16; ++i_s) {
          auto mat                             = elements[i_t, i_f, i_za_inc, joker, joker, i_s];
          result[i_t, i_f, i_za_inc].data[i_s] = interp(weights, mat, delta_aa_interp, za_scat_interp);
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Gridded> PhaseMatrixData<Format::ARO, Representation::Gridded>::regrid(
    const ScatteringDataGrids &grids, const RegridWeights &weights) const {
  PhaseMatrixData   result(grids.t_grid, grids.f_grid, grids.za_inc_grid, grids.aa_scat_grid, grids.za_scat_grid);
  const TensorType &coeffs_this = *this;
  TensorType       &coeffs_res  = result;

  for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
    GridPos     gp_t    = weights.t_grid_weights[i_t];
    const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
    Numeric     w_t_l   = gp_t.fd[1];
    Numeric     w_t_r   = gp_t.fd[0];
    for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
      GridPos     gp_f    = weights.f_grid_weights[i_f];
      const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
      Numeric     w_f_l   = gp_f.fd[1];
      Numeric     w_f_r   = gp_f.fd[0];
      for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size(); ++i_za_inc) {
        GridPos     gp_za_inc    = weights.za_inc_grid_weights[i_za_inc];
        const Index za_inc_upper = std::min<Index>(gp_za_inc.idx + 1, coeffs_this.extent(2) - 1);
        Numeric     w_za_inc_l   = gp_za_inc.fd[1];
        Numeric     w_za_inc_r   = gp_za_inc.fd[0];
        for (Size i_aa_scat = 0; i_aa_scat < weights.aa_scat_grid_weights.size(); ++i_aa_scat) {
          GridPos     gp_aa_scat    = weights.aa_scat_grid_weights[i_aa_scat];
          const Index aa_scat_upper = std::min<Index>(gp_aa_scat.idx + 1, coeffs_this.extent(3) - 1);
          Numeric     w_aa_scat_l   = gp_aa_scat.fd[1];
          Numeric     w_aa_scat_r   = gp_aa_scat.fd[0];
          for (Size i_za_scat = 0; i_za_scat < weights.za_scat_grid_weights.size(); ++i_za_scat) {
            GridPos     gp_za_scat    = weights.za_scat_grid_weights[i_za_scat];
            const Index za_scat_upper = std::min<Index>(gp_za_scat.idx + 1, coeffs_this.extent(4) - 1);
            Numeric     w_za_scat_l   = gp_za_scat.fd[1];
            Numeric     w_za_scat_r   = gp_za_scat.fd[0];

            coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] = rtepack::muelmat{0.0};

            if (w_t_l > 0.0) {
              if (w_f_l > 0.0) {
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, aa_scat_upper, za_scat_upper];
              }
              if (w_f_r > 0.0) {
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[gp_t.idx, f_upper, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[gp_t.idx, f_upper, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[gp_t.idx, f_upper, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_l * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[gp_t.idx, f_upper, za_inc_upper, aa_scat_upper, za_scat_upper];
              }
            }
            if (w_t_r > 0.0) {
              if (w_f_l > 0.0) {
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[t_upper, gp_f.idx, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[t_upper, gp_f.idx, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[t_upper, gp_f.idx, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_l * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[t_upper, gp_f.idx, za_inc_upper, aa_scat_upper, za_scat_upper];
              }
              if (w_f_r > 0.0) {
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[t_upper, f_upper, gp_za_inc.idx, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_l * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[t_upper, f_upper, gp_za_inc.idx, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[t_upper, f_upper, gp_za_inc.idx, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_l * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[t_upper, f_upper, gp_za_inc.idx, aa_scat_upper, za_scat_upper];

                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_l *
                    coeffs_this[t_upper, f_upper, za_inc_upper, gp_aa_scat.idx, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_r * w_aa_scat_l * w_za_scat_r *
                    coeffs_this[t_upper, f_upper, za_inc_upper, gp_aa_scat.idx, za_scat_upper];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_l *
                    coeffs_this[t_upper, f_upper, za_inc_upper, aa_scat_upper, gp_za_scat.idx];
                coeffs_res[i_t, i_f, i_za_inc, i_aa_scat, i_za_scat] +=
                    w_t_r * w_f_r * w_za_inc_r * w_aa_scat_r * w_za_scat_r *
                    coeffs_this[t_upper, f_upper, za_inc_upper, aa_scat_upper, za_scat_upper];
              }
            }
          }
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Gridded> PhaseMatrixData<Format::ARO, Representation::Gridded>::regrid(
    const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, delta_aa_grid_, za_scat_grid_, grids);
  return regrid(grids, weights);
}

auto PhaseMatrixData<Format::ARO, Representation::Gridded>::to_spectral(std::shared_ptr<SHT> sht) const
    -> PhaseMatrixDataSpectral {
  assert(sht->get_n_azimuth_angles() == n_delta_aa_);
  assert(sht->get_n_zenith_angles() == n_za_scat_);

  PhaseMatrixDataSpectral result(t_grid_, f_grid_, za_inc_grid_, sht);
  const auto              elements = rtepack::elements(*this);
  auto                    coeffs   = rtepack::elements(result);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        for (Index i_s = 0; i_s < 16; ++i_s) {
          coeffs[i_t, i_f, i_za_inc, joker, i_s] = sht->transform(elements[i_t, i_f, i_za_inc, joker, joker, i_s]);
        }
      }
    }
  }
  return result;
}

auto PhaseMatrixData<Format::ARO, Representation::Gridded>::to_spectral(Index degree, Index order) const
    -> PhaseMatrixDataSpectral {
  auto sht_ptr = sht::provider.get_instance_lm(degree, order);
  return to_spectral(sht_ptr);
}

auto PhaseMatrixData<Format::ARO, Representation::Gridded>::to_spectral() const -> PhaseMatrixDataSpectral {
  return to_spectral(sht::provider.get_instance(n_delta_aa_, n_za_scat_));
}

auto PhaseMatrixData<Format::ARO, Representation::Gridded>::to_fourier(Index max_mode) const -> PhaseMatrixDataFourier {
  ARTS_USER_ERROR_IF(max_mode < 0, "The highest Fourier mode must be >= 0, got {}", max_mode)
  ARTS_USER_ERROR_IF(
      n_delta_aa_ < 2 or std::abs((*delta_aa_grid_)[n_delta_aa_ - 1] - (*delta_aa_grid_)[0] - 360.0) > 1e-9,
      "Azimuthal Fourier modes need the ARO data over one full period of azimuth differences (the "
      "last angle of the grid 360 deg above the first), but the grid spans [{}, {}] deg",
      n_delta_aa_ > 0 ? (*delta_aa_grid_)[0] : 0.0,
      n_delta_aa_ > 0 ? (*delta_aa_grid_)[n_delta_aa_ - 1] : 0.0)
  for (Index k = 0; k + 1 < n_delta_aa_; ++k)
    ARTS_USER_ERROR_IF(not((*delta_aa_grid_)[k] < (*delta_aa_grid_)[k + 1]),
                       "The azimuth-difference grid must ascend strictly")

  PhaseMatrixDataFourier result(t_grid_, f_grid_, za_inc_grid_, za_scat_grid_, max_mode);
  // w[k, m, cs]: the weight of node k in C_m (cs = 0) and S_m (cs = 1) of the piecewise-linear function
  Tensor3 w(n_delta_aa_, max_mode + 1, 2, 0.0);
  for (Index k = 0; k + 1 < n_delta_aa_; ++k) {
    const Numeric a = Conversion::deg2rad((*delta_aa_grid_)[k]), b = Conversion::deg2rad((*delta_aa_grid_)[k + 1]);
    const Numeric h  = b - a;
    w[k, 0, 0]      += 0.5 * h / Constant::two_pi;
    w[k + 1, 0, 0]  += 0.5 * h / Constant::two_pi;
    for (Index m = 1; m <= max_mode; ++m) {
      // int f cos(m x) = [f sin(m x) / m] + s [cos(m x)] / m^2, int f sin(m x) = [-f cos(m x) / m] + s [sin(m x)] / m^2,
      // with f linear from f_a at a to f_b at b and s = (f_b - f_a) / h
      const Numeric dm = static_cast<Numeric>(m);
      const Numeric ca = std::cos(dm * a), cb = std::cos(dm * b), sa = std::sin(dm * a), sb = std::sin(dm * b);
      const Numeric dc = (cb - ca) / (dm * dm * h), ds = (sb - sa) / (dm * dm * h);
      w[k, m, 0]     += (-sa / dm - dc) / Constant::pi;
      w[k + 1, m, 0] += (sb / dm + dc) / Constant::pi;
      w[k, m, 1]     += (ca / dm - ds) / Constant::pi;
      w[k + 1, m, 1] += (-cb / dm + ds) / Constant::pi;
    }
  }
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f)
      for (Index i_i = 0; i_i < n_za_inc_; ++i_i)
        for (Index i_s = 0; i_s < n_za_scat_; ++i_s)
          for (Index k = 0; k < n_delta_aa_; ++k)
            for (Index m = 0; m <= max_mode; ++m)
              for (Index cs = 0; cs < 2; ++cs)
                for (Index e = 0; e < 16; ++e)
                  result[i_t, i_f, i_i, i_s, m, cs].data[e] +=
                      w[k, m, cs] * this->operator[](i_t, i_f, i_i, k, i_s).data[e];

  // The phase integral 2 pi int C_0,11 sin(za) dza of the data, linear in the scattering zenith angle between
  // the nodes, exactly; unknown (NaN) unless the nodes span all scattering zenith angles
  const Vector za{grid_vector(*za_scat_grid_)};
  Tensor3      integral(n_temps_, n_freqs_, n_za_inc_, std::numeric_limits<Numeric>::quiet_NaN());
  if (n_za_scat_ >= 2 and std::abs(za[0]) <= 1e-9 and std::abs(za[n_za_scat_ - 1] - 180.0) <= 1e-9) {
    integral = 0.0;
    for (Index k = 0; k + 1 < n_za_scat_; ++k) {
      // int sin = cos(a) - cos(b), int x sin(x) = [sin(x) - x cos(x)], for the two linear weights of the segment
      const Numeric a = Conversion::deg2rad(za[k]), b = Conversion::deg2rad(za[k + 1]), h = b - a;
      const Numeric i0 = std::cos(a) - std::cos(b);
      const Numeric i1 = (std::sin(b) - b * std::cos(b)) - (std::sin(a) - a * std::cos(a));
      const Numeric wa = Constant::two_pi * (b * i0 - i1) / h, wb = Constant::two_pi * (i1 - a * i0) / h;
      for (Index i_t = 0; i_t < n_temps_; ++i_t)
        for (Index i_f = 0; i_f < n_freqs_; ++i_f)
          for (Index i_i = 0; i_i < n_za_inc_; ++i_i)
            integral[i_t, i_f, i_i] +=
                wa * result[i_t, i_f, i_i, k, 0, 0][0, 0] + wb * result[i_t, i_f, i_i, k + 1, 0, 0][0, 0];
    }
  }
  result.set_phase_integral(std::move(integral));
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Gridded> tro_lab_frame(std::shared_ptr<const Vector> t_grid,
                                                                    std::shared_ptr<const Vector> f_grid,
                                                                    std::shared_ptr<const Vector> za_inc_grid,
                                                                    std::shared_ptr<const Vector> delta_aa_grid,
                                                                    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                                                    const ScatteringMatrixFunction &scattering_matrix) {
  PhaseMatrixData<Format::ARO, Representation::Gridded> result(
      t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);

  const auto                             za_scat = grid_vector(*za_scat_grid);
  rtepack::compact_planar_muelmat_matrix tro(t_grid->size(), f_grid->size());
  for (Size i_za_inc = 0; i_za_inc < za_inc_grid->size(); ++i_za_inc) {
    for (Size i_delta_aa = 0; i_delta_aa < delta_aa_grid->size(); ++i_delta_aa) {
      // In [0, 360) deg, so that delta_aa > 180 tells the side of the principal plane for any grid, e.g. [-180, 180]
      Numeric delta_aa = std::fmod((*delta_aa_grid)[i_delta_aa], Numeric{360});
      if (delta_aa < 0) delta_aa += 360;
      for (Size i_za_scat = 0; i_za_scat < za_scat.size(); ++i_za_scat) {
        const auto [theta, c_1, c_2, s_1, s_2] =
            detail::rotation_coefficients(0.0, (*za_inc_grid)[i_za_inc], delta_aa, za_scat[i_za_scat]);
        scattering_matrix(theta, tro);
        for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
          for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
            // Beyond 180 deg, the mirror image of the other side of the principal plane
            const rtepack::muelmat Z                          = rtepack::rotated(tro[i_t, i_f], c_1, s_1, c_2, s_2);
            result[i_t, i_f, i_za_inc, i_delta_aa, i_za_scat] = delta_aa > 180.0 ? rtepack::mirror(Z) : Z;
          }
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Spectral> &
PhaseMatrixData<Format::ARO, Representation::Spectral>::operator=(const TensorType &data) {
  ARTS_USER_ERROR_IF(data.shape()[0] != n_temps_,
                     "Provided backscatter coefficient data do not match temperature grid.");
  ARTS_USER_ERROR_IF(data.shape()[1] != n_freqs_, "Provided backscatter coefficient data do not match frequency grid.");
  ARTS_USER_ERROR_IF(data.shape()[2] != n_za_inc_,
                     "Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
  ARTS_USER_ERROR_IF(data.shape()[3] != n_spectral_coeffs_,
                     "Provided backscatter coefficient data do not match expected number of SHT coefficients.");
  TensorType::operator=(data);
  return *this;
}

auto PhaseMatrixData<Format::ARO, Representation::Spectral>::to_gridded() const -> PhaseMatrixDataGridded {
  PhaseMatrixDataGridded result(t_grid_, f_grid_, za_inc_grid_, sht_->get_aa_grid_ptr(), sht_->get_za_grid_ptr());
  auto                   elements = rtepack::elements(result);
  const auto             coeffs   = rtepack::elements(*this);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        for (Index i_s = 0; i_s < 16; ++i_s) {
          elements[i_t, i_f, i_za_inc, joker, joker, i_s] = sht_->synthesize(coeffs[i_t, i_f, i_za_inc, joker, i_s]);
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Spectral>
PhaseMatrixData<Format::ARO, Representation::Spectral>::to_spectral(Index l_new, Index m_new) const {
  auto            sht_new = sht::provider.get_instance_lm(l_new, m_new);
  PhaseMatrixData pm_new(t_grid_, f_grid_, za_inc_grid_, sht_new);
  auto            coeffs_new = rtepack::elements(pm_new);
  const auto      coeffs     = rtepack::elements(*this);
  for (Size f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
    for (Size t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
      for (Size za_inc_ind = 0; za_inc_ind < za_inc_grid_->size(); ++za_inc_ind) {
        for (Index i_s = 0; i_s < 16; ++i_s) {
          coeffs_new[t_ind, f_ind, za_inc_ind, joker, i_s] =
              sht::add_coeffs(*sht_new,
                              coeffs_new[t_ind, f_ind, za_inc_ind, joker, i_s],
                              *sht_,
                              coeffs[t_ind, f_ind, za_inc_ind, joker, i_s]);
        }
      }
    }
  }
  return pm_new;
}

BackscatterMatrixData<Format::ARO> PhaseMatrixData<Format::ARO, Representation::Spectral>::extract_backscatter_matrix()
    const {
  BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
  const auto                         elements = rtepack::elements(*this);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        auto    za_inc  = (*za_inc_grid_)[i_za_inc];
        Numeric za_scat = 180.0 - za_inc;
        for (Index i_s = 0; i_s < 16; ++i_s) {
          auto coeffs = elements[i_t, i_f, i_za_inc, joker, i_s];
          result[i_t, i_f, i_za_inc].data[i_s] =
              sht_->evaluate(coeffs, Conversion::deg2rad(180.0), Conversion::deg2rad(za_scat));
        }
      }
    }
  }
  return result;
}

ForwardscatterMatrixData<Format::ARO>
PhaseMatrixData<Format::ARO, Representation::Spectral>::extract_forwardscatter_matrix() {
  BackscatterMatrixData<Format::ARO> result(t_grid_, f_grid_, za_inc_grid_);
  const auto                         elements = rtepack::elements(*this);

  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        auto za_inc = (*za_inc_grid_)[i_za_inc];
        for (Index i_s = 0; i_s < 16; ++i_s) {
          auto coeffs = elements[i_t, i_f, i_za_inc, joker, i_s];
          result[i_t, i_f, i_za_inc].data[i_s] =
              sht_->evaluate(coeffs, Conversion::deg2rad(0.0), Conversion::deg2rad(za_inc));
        }
      }
    }
  }
  return result;
}

rtepack::muelmat_tensor3 PhaseMatrixData<Format::ARO, Representation::Spectral>::integrate_phase_matrix() const {
  const Numeric            factor = sht_->get_l_max() == 0 ? 4.0 * Constant::pi : std::sqrt(4.0 * Constant::pi);
  rtepack::muelmat_tensor3 results(n_temps_, n_freqs_, n_za_inc_);
  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_za_inc = 0; i_za_inc < n_za_inc_; ++i_za_inc) {
        const rtepack::specmat &a_0 = this->operator[](i_t, i_f, i_za_inc, 0);
        for (Index i_s = 0; i_s < 16; ++i_s) { results[i_t, i_f, i_za_inc].data[i_s] = a_0.data[i_s].real() * factor; }
      }
    }
  }
  return results;
}

PhaseMatrixData<Format::ARO, Representation::Spectral> PhaseMatrixData<Format::ARO, Representation::Spectral>::regrid(
    const ScatteringDataGrids &grids, const RegridWeights weights) const {
  PhaseMatrixData   result(grids.t_grid, grids.f_grid, grids.za_inc_grid, sht_);
  const TensorType &coeffs_this = *this;
  TensorType       &coeffs_res  = result;

  for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
    GridPos     gp_t    = weights.t_grid_weights[i_t];
    const Index t_upper = std::min<Index>(gp_t.idx + 1, coeffs_this.extent(0) - 1);
    Numeric     w_t_l   = gp_t.fd[1];
    Numeric     w_t_r   = gp_t.fd[0];
    for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
      GridPos     gp_f    = weights.f_grid_weights[i_f];
      const Index f_upper = std::min<Index>(gp_f.idx + 1, coeffs_this.extent(1) - 1);
      Numeric     w_f_l   = gp_f.fd[1];
      Numeric     w_f_r   = gp_f.fd[0];
      for (Size i_za_inc = 0; i_za_inc < weights.za_inc_grid_weights.size(); ++i_za_inc) {
        GridPos     gp_za_inc    = weights.za_inc_grid_weights[i_za_inc];
        const Index za_inc_upper = std::min<Index>(gp_za_inc.idx + 1, coeffs_this.extent(2) - 1);
        Numeric     w_za_inc_l   = gp_za_inc.fd[1];
        Numeric     w_za_inc_r   = gp_za_inc.fd[0];

        for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
          coeffs_res[i_t, i_f, i_za_inc, i_sht] = rtepack::specmat{0.0};
        }

        if (w_t_l > 0.0) {
          if (w_f_l > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_l * w_f_l * w_za_inc_l * coeffs_this[gp_t.idx, gp_f.idx, gp_za_inc.idx, i_sht];
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_l * w_f_l * w_za_inc_r * coeffs_this[gp_t.idx, gp_f.idx, za_inc_upper, i_sht];
            }
          }
          if (w_f_r > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_l * w_f_r * w_za_inc_l * coeffs_this[gp_t.idx, f_upper, gp_za_inc.idx, i_sht];
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_l * w_f_r * w_za_inc_r * coeffs_this[gp_t.idx, f_upper, za_inc_upper, i_sht];
            }
          }
        }
        if (w_t_r > 0.0) {
          if (w_f_l > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_r * w_f_l * w_za_inc_l * coeffs_this[t_upper, gp_f.idx, gp_za_inc.idx, i_sht];
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_r * w_f_l * w_za_inc_r * coeffs_this[t_upper, gp_f.idx, za_inc_upper, i_sht];
            }
          }
          if (w_f_r > 0.0) {
            for (Index i_sht = 0; i_sht < n_spectral_coeffs_; ++i_sht) {
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_r * w_f_r * w_za_inc_l * coeffs_this[t_upper, f_upper, gp_za_inc.idx, i_sht];
              coeffs_res[i_t, i_f, i_za_inc, i_sht] +=
                  w_t_r * w_f_r * w_za_inc_r * coeffs_this[t_upper, f_upper, za_inc_upper, i_sht];
            }
          }
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Spectral> PhaseMatrixData<Format::ARO, Representation::Spectral>::regrid(
    const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, nullptr, nullptr, grids);
  return regrid(grids, weights);
}

PhaseMatrixData<Format::ARO, Representation::Fourier>
PhaseMatrixData<Format::ARO, Representation::Spectral>::to_fourier(std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
                                                                   Index max_mode) const {
  PhaseMatrixData<Format::ARO, Representation::Fourier> out(t_grid_, f_grid_, za_inc_grid_, za_scat_grid, max_mode);
  const Vector                                          za{grid_vector(*za_scat_grid)};
  const Index                                           N        = 2 * (sht_->get_m_max() + max_mode) + 2;
  const auto                                            elements = rtepack::elements(*this);
  for (Index i_t = 0; i_t < n_temps_; ++i_t) {
    for (Index i_f = 0; i_f < n_freqs_; ++i_f) {
      for (Index i_i = 0; i_i < n_za_inc_; ++i_i) {
        for (Index e = 0; e < 16; ++e) {
          const auto coeffs = elements[i_t, i_f, i_i, joker, e];
          for (Size i_s = 0; i_s < za.size(); ++i_s) {
            const Numeric theta = Conversion::deg2rad(za[i_s]);
            for (Index k = 0; k < N; ++k) {
              const Numeric phi = Constant::two_pi * static_cast<Numeric>(k) / static_cast<Numeric>(N);
              const Numeric v   = sht_->evaluate(coeffs, phi, theta) / static_cast<Numeric>(N);
              for (Index m = 0; m <= max_mode; ++m) {
                const Numeric w                        = m == 0 ? 1.0 : 2.0;
                out[i_t, i_f, i_i, i_s, m, 0].data[e] += w * v * std::cos(static_cast<Numeric>(m) * phi);
                out[i_t, i_f, i_i, i_s, m, 1].data[e] += w * v * std::sin(static_cast<Numeric>(m) * phi);
              }
            }
          }
        }
      }
    }
  }
  const auto integral = integrate_phase_matrix();
  Tensor3    sigma(n_temps_, n_freqs_, n_za_inc_);
  for (Index i_t = 0; i_t < n_temps_; ++i_t)
    for (Index i_f = 0; i_f < n_freqs_; ++i_f)
      for (Index i_i = 0; i_i < n_za_inc_; ++i_i) sigma[i_t, i_f, i_i] = integral[i_t, i_f, i_i][0, 0];
  out.set_phase_integral(std::move(sigma));
  return out;
}

PhaseMatrixData<Format::ARO, Representation::Fourier> &PhaseMatrixData<Format::ARO, Representation::Fourier>::operator=(
    const rtepack::muelmat_tensor6 &data) {
  ARTS_USER_ERROR_IF(data.shape() != this->shape(),
                     "The Fourier modes have shape {:B,}, but these grids need {:B,}",
                     data.shape(),
                     this->shape())
  TensorType::operator=(data);
  return *this;
}

void PhaseMatrixData<Format::ARO, Representation::Fourier>::set_phase_integral(Tensor3 phase_integral) {
  ARTS_USER_ERROR_IF(phase_integral.shape() != phase_integral_.shape(),
                     "The phase integral has shape {:B,}, but these grids need {:B,}",
                     phase_integral.shape(),
                     phase_integral_.shape())
  phase_integral_ = std::move(phase_integral);
}

PhaseMatrixData<Format::ARO, Representation::Fourier> &
PhaseMatrixData<Format::ARO, Representation::Fourier>::operator+=(const PhaseMatrixData &other) {
  ARTS_USER_ERROR_IF(other.shape() != this->shape(),
                     "Cannot add Fourier modes of shape {:B,} to modes of shape {:B,}",
                     other.shape(),
                     this->shape())
  static_cast<TensorType &>(*this) += static_cast<const TensorType &>(other);
  phase_integral_                  += other.phase_integral_;
  return *this;
}

PhaseMatrixData<Format::ARO, Representation::Fourier> &
PhaseMatrixData<Format::ARO, Representation::Fourier>::operator*=(Numeric x) {
  static_cast<TensorType &>(*this) *= x;
  phase_integral_                  *= x;
  return *this;
}

auto PhaseMatrixData<Format::ARO, Representation::Fourier>::to_gridded(
    std::shared_ptr<const Vector> delta_aa_grid) const -> PhaseMatrixDataGridded {
  PhaseMatrixDataGridded result(t_grid_, f_grid_, za_inc_grid_, delta_aa_grid, za_scat_grid_);
  const auto [nt, nf, ni, ns, nm, ncs] = this->shape();
  for (Size k = 0; k < delta_aa_grid->size(); ++k) {
    const Numeric delta = Conversion::deg2rad((*delta_aa_grid)[k]);
    for (Index m = 0; m < nm; ++m) {
      const Numeric c = std::cos(static_cast<Numeric>(m) * delta), s = std::sin(static_cast<Numeric>(m) * delta);
      for (Index i_t = 0; i_t < nt; ++i_t)
        for (Index i_f = 0; i_f < nf; ++i_f)
          for (Index i_i = 0; i_i < ni; ++i_i)
            for (Index i_s = 0; i_s < ns; ++i_s) {
              const rtepack::muelmat &C = this->operator[](i_t, i_f, i_i, i_s, m, 0);
              const rtepack::muelmat &S = this->operator[](i_t, i_f, i_i, i_s, m, 1);
              for (Index e = 0; e < 16; ++e) result[i_t, i_f, i_i, k, i_s].data[e] += c * C.data[e] + s * S.data[e];
            }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Fourier> PhaseMatrixData<Format::ARO, Representation::Fourier>::to_fourier(
    Index max_mode) const {
  ARTS_USER_ERROR_IF(max_mode > max_mode_,
                     "The Fourier modes go to m = {}, so they cannot give the modes to m = {}; convert the gridded "
                     "data to that mode instead",
                     max_mode_,
                     max_mode)
  PhaseMatrixData out(t_grid_, f_grid_, za_inc_grid_, za_scat_grid_, max_mode);
  const auto [nt, nf, ni, ns, nm, ncs] = out.shape();
  for (Index i_t = 0; i_t < nt; ++i_t)
    for (Index i_f = 0; i_f < nf; ++i_f)
      for (Index i_i = 0; i_i < ni; ++i_i)
        for (Index i_s = 0; i_s < ns; ++i_s)
          for (Index m = 0; m < nm; ++m) out[i_t, i_f, i_i, i_s, m] = this->operator[](i_t, i_f, i_i, i_s, m);
  out.phase_integral_ = phase_integral_;
  return out;
}

PhaseMatrixData<Format::ARO, Representation::Fourier> PhaseMatrixData<Format::ARO, Representation::Fourier>::regrid(
    const ScatteringDataGrids &grids, const RegridWeights &weights) const {
  auto            new_inc  = grids.za_inc_grid ? grids.za_inc_grid : za_inc_grid_;
  auto            new_scat = grids.za_scat_grid ? grids.za_scat_grid : za_scat_grid_;
  const auto      inc      = node_indices(*new_inc, *za_inc_grid_, "incidence");
  const auto      scat     = node_indices(grid_vector(*new_scat), grid_vector(*za_scat_grid_), "scattering");
  PhaseMatrixData result(grids.t_grid, grids.f_grid, new_inc, new_scat, max_mode_);
  const Index     nt = this->extent(0), nf = this->extent(1);
  for (Size i_t = 0; i_t < weights.t_grid_weights.size(); ++i_t) {
    const GridPos gp_t = weights.t_grid_weights[i_t];
    const Index   t0 = std::clamp<Index>(gp_t.idx, 0, nt - 1), t1 = std::min<Index>(t0 + 1, nt - 1);
    for (Size i_f = 0; i_f < weights.f_grid_weights.size(); ++i_f) {
      const GridPos gp_f = weights.f_grid_weights[i_f];
      const Index   f0 = std::clamp<Index>(gp_f.idx, 0, nf - 1), f1 = std::min<Index>(f0 + 1, nf - 1);
      const std::array<std::tuple<Numeric, Index, Index>, 4> corners{{{gp_t.fd[1] * gp_f.fd[1], t0, f0},
                                                                      {gp_t.fd[1] * gp_f.fd[0], t0, f1},
                                                                      {gp_t.fd[0] * gp_f.fd[1], t1, f0},
                                                                      {gp_t.fd[0] * gp_f.fd[0], t1, f1}}};
      for (const auto &[w, it, jf] : corners) {
        if (w == 0.0) continue;
        for (Size i_i = 0; i_i < inc.size(); ++i_i) {
          result.phase_integral_[i_t, i_f, i_i] += w * phase_integral_[it, jf, inc[i_i]];
          for (Size i_s = 0; i_s < scat.size(); ++i_s) {
            auto       out = result[i_t, i_f, i_i, i_s];
            const auto in  = this->operator[](it, jf, inc[i_i], scat[i_s]);
            for (auto [o, x] : std::views::zip(out | by_elem, in | by_elem))
              for (auto [oe, xe] : std::views::zip(o.data, x.data)) oe += w * xe;
          }
        }
      }
    }
  }
  return result;
}

PhaseMatrixData<Format::ARO, Representation::Fourier> PhaseMatrixData<Format::ARO, Representation::Fourier>::regrid(
    const ScatteringDataGrids &grids) const {
  auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
  return regrid(grids, weights);
}

Index PhaseMatrixData<Format::ARO, Representation::Fourier>::checked_mode(Index max_mode) {
  ARTS_USER_ERROR_IF(max_mode < 0, "The highest Fourier mode must be >= 0, got {}", max_mode)
  return max_mode;
}

std::vector<Index> PhaseMatrixData<Format::ARO, Representation::Fourier>::node_indices(
    const StridedConstVectorView &wanted, const StridedConstVectorView &nodes, const char *kind) {
  std::vector<Index> out(wanted.size());
  for (Size i = 0; i < wanted.size(); ++i) {
    const auto it = stdr::find_if(nodes, [x = wanted[i]](Numeric y) { return std::abs(x - y) <= 1e-9; });
    ARTS_USER_ERROR_IF(it == nodes.end(),
                       "ARO Fourier modes are known on their own {} zenith angles only, {:B,} deg, which do not "
                       "include {} deg; give the data on the wanted angles",
                       kind,
                       nodes,
                       wanted[i])
    out[i] = static_cast<Index>(it - nodes.begin());
  }
  return out;
}

PhaseMatrixData<Format::ARO, Representation::Fourier> tro_lab_frame_fourier_modes(
    std::shared_ptr<const Vector>          t_grid,
    std::shared_ptr<const Vector>          f_grid,
    std::shared_ptr<const Vector>          za_inc_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid,
    Index                                  max_mode,
    const ConstMatrixView                 &phase_integral,
    const ScatteringMatrixFunction        &scattering_matrix) {
  PhaseMatrixData<Format::ARO, Representation::Fourier> out(t_grid, f_grid, za_inc_grid, za_scat_grid, max_mode);
  ARTS_USER_ERROR_IF(phase_integral.nrows() != static_cast<Index>(t_grid->size()) or
                         phase_integral.ncols() != static_cast<Index>(f_grid->size()),
                     "The phase integral has shape [{}, {}], but there are {} temperatures and {} frequencies",
                     phase_integral.nrows(),
                     phase_integral.ncols(),
                     t_grid->size(),
                     f_grid->size()) {
    Tensor3 integral(t_grid->size(), f_grid->size(), za_inc_grid->size());
    for (Index i_t = 0; i_t < integral.extent(0); ++i_t)
      for (Index i_f = 0; i_f < integral.extent(1); ++i_f) integral[i_t, i_f, joker] = phase_integral[i_t, i_f];
    out.set_phase_integral(std::move(integral));
  }

  const Vector za_scat{grid_vector(*za_scat_grid)};
  for (const Numeric za : *za_inc_grid)
    ARTS_USER_ERROR_IF(not(za >= 0.0 and za <= 180.0), "Incidence zenith angles must be in [0, 180] deg, got {}", za)
  for (const Numeric za : za_scat)
    ARTS_USER_ERROR_IF(not(za >= 0.0 and za <= 180.0), "Scattering zenith angles must be in [0, 180] deg, got {}", za)

  constexpr Index   max_nodes = 4096;
  constexpr Numeric tolerance = 1e-13;
  const Index       nt = t_grid->size(), nf = f_grid->size(), M = max_mode, nset = nt * nf;

  // Z(2 pi - delta) = D Z(delta) D with D = diag(1, 1, -1, -1), the side flip rtepack::mirror: the
  // elements with exactly one index in {U, V} are odd in delta and have sine modes only, the others cosine modes
  // only.  So [0, pi] gives the modes: C_m = (2 - delta_m0) / 2 sum_k w_k Z(phi_k) cos(m phi_k) and S_m = sum_k
  // w_k Z(phi_k) sin(m phi_k) for the n-point Gauss-Legendre rule (w_k, phi_k) on [0, pi], whose ends are the
  // only directions where Z can have a kink (exact forward and backward scattering).
  constexpr std::array<bool, 16> odd{
      false, false, true, true, false, false, true, true, true, true, false, false, true, true, false, false};

  //! The rule of n nodes on [0, pi], and its weights of every mode, [n, M + 1]
  struct rule {
    Vector phi;
    Matrix cos_weight, sin_weight;
  };

  Index start = 16;
  while (start < M + 8) start *= 2;

  // [set, m, (cos, sin), element]
  const Index nmode = (M + 1) * 2 * 16;

  // The direction pairs are independent; each thread has its own copy of the scattering matrix (which may keep
  // scratch), its own rules and its own work arrays, and every pair writes its own part of out
  std::string error;
#pragma omp parallel if (not arts_omp_in_parallel())
  {
    auto                                   local_matrix = scattering_matrix;
    std::map<Index, rule>                  rules;
    rtepack::compact_planar_muelmat_matrix F(nt, nf);
    std::vector<Numeric>                   coarse(nset * nmode), fine(nset * nmode);

    const auto get_rule = [&](Index n) -> const rule & {
      auto [it, inserted] = rules.try_emplace(n);
      if (inserted) {
        Vector x(n), w(n);
        Legendre::GaussLegendre(x, w);
        rule &r = it->second;
        r.phi   = Vector(n);
        r.cos_weight.resize(n, M + 1);
        r.sin_weight.resize(n, M + 1);
        for (Index k = 0; k < n; k++) {
          r.phi[k] = Constant::pi * 0.5 * (x[k] + 1.0);
          for (Index m = 0; m <= M; m++) {
            const Numeric dm   = static_cast<Numeric>(m);
            r.cos_weight[k, m] = (m == 0 ? 0.5 : 1.0) * w[k] * std::cos(dm * r.phi[k]);
            r.sin_weight[k, m] = w[k] * std::sin(dm * r.phi[k]);
          }
        }
      }
      return it->second;
    };

#pragma omp for collapse(2) schedule(dynamic)
    for (Size ii = 0; ii < za_inc_grid->size(); ii++) {
      for (Size is = 0; is < za_scat.size(); is++) {
        try {
          Numeric    scale = 0.0;
          const auto modes = [&](Index n, std::vector<Numeric> &c) {
            const rule &r = get_rule(n);
            stdr::fill(c, 0.0);
            for (Index k = 0; k < n; k++) {
              const Numeric delta = Conversion::rad2deg(r.phi[k]);
              const auto [theta, c_1, c_2, s_1, s_2] =
                  detail::rotation_coefficients(0.0, (*za_inc_grid)[ii], delta, za_scat[is]);
              local_matrix(theta, F);
              for (Index set = 0; set < nset; set++) {
                // Its 16 elements, row-major
                const auto Z = rtepack::rotated(F[set / nf, set % nf], c_1, s_1, c_2, s_2).data;
                for (Index e = 0; e < 16; e++) scale = std::max<Numeric>(scale, std::abs(Z[e]));
                Numeric *cs = c.data() + set * nmode;
                for (Index m = 0; m <= M; m++, cs += 32) {
                  const Numeric wc = r.cos_weight[k, m], ws = r.sin_weight[k, m];
                  for (Index e = 0; e < 16; e++) {
                    cs[e]      += wc * Z[e];
                    cs[16 + e] += ws * Z[e];
                  }
                }
              }
            }
            // The parts that vanish by the symmetry
            for (Index j = 0; j < nset * (M + 1); j++)
              for (Index e = 0; e < 16; e++) c[j * 32 + (odd[e] ? 0 : 16) + e] = 0.0;
          };

          Index n = start;
          modes(n, coarse);
          while (true) {
            n *= 2;
            modes(n, fine);
            Numeric change = 0.0;
            for (Size j = 0; j < fine.size(); j++) change = std::max(change, std::abs(fine[j] - coarse[j]));
            if (change <= tolerance * scale) break;
            ARTS_USER_ERROR_IF(2 * n > max_nodes,
                               "The azimuthal Fourier modes of the laboratory-frame phase matrix at za_inc = {} deg "
                               "and za_scat = {} deg did not converge with {} azimuths (they still change by {:.3e} "
                               "of the largest phase-matrix element); the scattering matrix is too sharp in the "
                               "scattering angle to resolve",
                               (*za_inc_grid)[ii],
                               za_scat[is],
                               2 * n,
                               change / scale)
            std::swap(coarse, fine);
          }

          for (Index i_t = 0; i_t < nt; i_t++)
            for (Index i_f = 0; i_f < nf; i_f++)
              for (Index m = 0; m <= M; m++)
                for (Index cs = 0; cs < 2; cs++)
                  for (Index e = 0; e < 16; e++)
                    out[i_t, i_f, ii, is, m, cs].data[e] = fine[(i_t * nf + i_f) * nmode + (m * 2 + cs) * 16 + e];
        } catch (const std::exception &e) {
#pragma omp critical(tro_lab_frame_fourier_modes)
          if (error.empty()) error = e.what();
        }
      }
    }
  }
  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  return out;
}

std::ostream &operator<<(std::ostream &out, Format format) {
  switch (format) {
    case Format::TRO:     out << "TRO"; break;
    case Format::ARO:     out << "ARO"; break;
    case Format::General: out << "General"; break;
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, Representation repr) {
  switch (repr) {
    case Representation::Gridded:        out << "gridded"; break;
    case Representation::Spectral:       out << "spectral"; break;
    case Representation::DoublySpectral: out << "doubly-spectral"; break;
    case Representation::Fourier:        out << "fourier"; break;
  }
  return out;
}

}  // namespace scattering
