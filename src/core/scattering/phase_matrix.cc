#include "phase_matrix.h"

#include <ranges>

namespace scattering {
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
    w[k, 0, 0]      += 0.5 * h / (2.0 * pi_v<Numeric>);
    w[k + 1, 0, 0]  += 0.5 * h / (2.0 * pi_v<Numeric>);
    for (Index m = 1; m <= max_mode; ++m) {
      // int f cos(m x) = [f sin(m x) / m] + s [cos(m x)] / m^2, int f sin(m x) = [-f cos(m x) / m] + s [sin(m x)] / m^2,
      // with f linear from f_a at a to f_b at b and s = (f_b - f_a) / h
      const Numeric dm = static_cast<Numeric>(m);
      const Numeric ca = std::cos(dm * a), cb = std::cos(dm * b), sa = std::sin(dm * a), sb = std::sin(dm * b);
      const Numeric dc = (cb - ca) / (dm * dm * h), ds = (sb - sa) / (dm * dm * h);
      w[k, m, 0]     += (-sa / dm - dc) / pi_v<Numeric>;
      w[k + 1, m, 0] += (sb / dm + dc) / pi_v<Numeric>;
      w[k, m, 1]     += (ca / dm - ds) / pi_v<Numeric>;
      w[k + 1, m, 1] += (-cb / dm + ds) / pi_v<Numeric>;
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
      const Numeric wa = 2.0 * pi_v<Numeric> * (b * i0 - i1) / h, wb = 2.0 * pi_v<Numeric> * (i1 - a * i0) / h;
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
              const Numeric phi = 2.0 * pi_v<Numeric> * static_cast<Numeric>(k) / static_cast<Numeric>(N);
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
