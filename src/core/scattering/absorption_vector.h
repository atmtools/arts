#ifndef SCATTERING_ABSORPTION_VECTOR_H_
#define SCATTERING_ABSORPTION_VECTOR_H_

#include <matpack.h>
#include <rtepack.h>

#include <memory>

#include "phase_matrix.h"

namespace scattering {

template <Format format, Representation representation> class AbsorptionVectorData;

/** The absorption vectors [t, f] of totally randomly oriented particles
 *
 * Their absorption vector is [a1, 0, 0, 0]: I of each
 * compact_planar_stokvec is a1, and Q = 0.
 */
template <Representation repr> class AbsorptionVectorData<Format::TRO, repr>
    : public rtepack::compact_planar_stokvec_matrix {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::compact_planar_stokvec_matrix::reshape;
  using rtepack::compact_planar_stokvec_matrix::resize;

 public:
  using AbsorptionVectorDataLabFrame = AbsorptionVectorData<Format::ARO, Representation::Gridded>;
  /// The tensor type used to store the absorption vector data.
  using TensorType = rtepack::compact_planar_stokvec_matrix;

  using TensorType::operator[];

  AbsorptionVectorData() = default;

  /** Create a new AbsorptionVectorData container.
   *
   * Creates a container to hold extinction matrix data for the
   * provided grids. The extinction matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   */
  AbsorptionVectorData(std::shared_ptr<const Vector> t_grid, std::shared_ptr<const Vector> f_grid)
      : TensorType(t_grid->size(), f_grid->size(), rtepack::compact_planar_stokvec{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid) {}

  AbsorptionVectorData& operator=(const TensorType& data) {
    if (data.shape()[0] != n_temps_)
      arts::user_error("Provided backscatter coefficient data do not match temperature grid.");
    if (data.shape()[1] != n_freqs_)
      arts::user_error("Provided backscatter coefficient data do not match frequency grid.");
    TensorType::operator=(data);
    return *this;
  }

  [[nodiscard]] AbsorptionVectorData<Format::ARO, repr> to_lab_frame(std::shared_ptr<const Vector> za_inc_grid) const {
    AbsorptionVectorData<Format::ARO, repr> av_new{t_grid_, f_grid_, za_inc_grid};
    for (Size t_ind = 0; t_ind < t_grid_->size(); ++t_ind) {
      for (Size f_ind = 0; f_ind < f_grid_->size(); ++f_ind) {
        for (Size za_inc_ind = 0; za_inc_ind < za_inc_grid->size(); ++za_inc_ind) {
          av_new[t_ind, f_ind, za_inc_ind] = this->operator[](t_ind, f_ind);
        }
      }
    }
    return av_new;
  }

  [[nodiscard]] AbsorptionVectorData<Format::TRO, Representation::Spectral> to_spectral() const {
    AbsorptionVectorData<Format::TRO, Representation::Spectral> avd_new(t_grid_, f_grid_);
    avd_new = *this;
    return avd_new;
  }

  [[nodiscard]] AbsorptionVectorData regrid(const ScatteringDataGrids grids, const RegridWeights weights) const {
    AbsorptionVectorData result(grids.t_grid, grids.f_grid);
    regrid_with(
        static_cast<TensorType&>(result), static_cast<const TensorType&>(*this), weights.t_lags, weights.f_lags);
    return result;
  }

  [[nodiscard]] AbsorptionVectorData regrid(const ScatteringDataGrids& grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, nullptr, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  [[nodiscard]] AbsorptionVectorData<Format::TRO, Representation::Gridded> to_gridded() const {
    auto abs_new = AbsorptionVectorData<Format::TRO, Representation::Gridded>(t_grid_, f_grid_);
    abs_new      = *this;
    return abs_new;
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;

  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;
};

/** The absorption vectors [t, f, za_inc] of azimuthally randomly oriented particles
 *
 * In the laboratory (meridional) basis I and Q of each
 * compact_planar_stokvec are a1 and a2.
 */
template <Representation repr> class AbsorptionVectorData<Format::ARO, repr>
    : public rtepack::compact_planar_stokvec_tensor3 {
 private:
  // Hiding resize and reshape functions to avoid inconsistencies.
  // between grids and data.
  using rtepack::compact_planar_stokvec_tensor3::reshape;
  using rtepack::compact_planar_stokvec_tensor3::resize;

 public:
  /// The tensor type used to store the absorption vector data.
  using TensorType = rtepack::compact_planar_stokvec_tensor3;

  AbsorptionVectorData() = default;
  /** Create a new AbsorptionVectorData container.
   *
   * Creates a container to hold extinction matrix data for the
   * provided grids. The extinction matrix data in the container is
   * initialized to 0.
   *
   * @param t_grid: A pointer to the temperature grid over which the
   * data is defined.
   * @param f_grid: A pointer to the frequency grid over which the
   * data is defined.
   */
  AbsorptionVectorData(std::shared_ptr<const Vector> t_grid,
                       std::shared_ptr<const Vector> f_grid,
                       std::shared_ptr<const Vector> za_inc_grid)
      : TensorType(t_grid->size(), f_grid->size(), za_inc_grid->size(), rtepack::compact_planar_stokvec{}),
        n_temps_(t_grid->size()),
        t_grid_(t_grid),
        n_freqs_(f_grid->size()),
        f_grid_(f_grid),
        n_za_inc_(za_inc_grid->size()),
        za_inc_grid_(za_inc_grid) {}

  AbsorptionVectorData& operator=(const TensorType& data) {
    if (data.shape()[0] != n_temps_)
      arts::user_error("Provided backscatter coefficient data do not match temperature grid.");
    if (data.shape()[1] != n_freqs_)
      arts::user_error("Provided backscatter coefficient data do not match frequency grid.");
    if (data.shape()[2] != n_za_inc_)
      arts::user_error("Provided backscatter coefficient data do not match expected number of incoming zenith angles.");
    TensorType::operator=(data);
    return *this;
  }

  [[nodiscard]] AbsorptionVectorData<Format::ARO, Representation::Spectral> to_spectral() const {
    AbsorptionVectorData<Format::ARO, Representation::Spectral> avd_new(t_grid_, f_grid_, za_inc_grid_);
    avd_new = *this;
    return avd_new;
  }

  //! The same data, as the absorption vector of azimuthal Fourier-mode optics (it has no azimuth dependence)
  [[nodiscard]] AbsorptionVectorData<Format::ARO, Representation::Fourier> to_fourier() const {
    AbsorptionVectorData<Format::ARO, Representation::Fourier> avd_new(t_grid_, f_grid_, za_inc_grid_);
    avd_new = *this;
    return avd_new;
  }

  [[nodiscard]] AbsorptionVectorData regrid(const ScatteringDataGrids& grids, const RegridWeights& weights) const {
    AbsorptionVectorData result(grids.t_grid, grids.f_grid, std::make_shared<Vector>(grid_vector(*grids.za_inc_grid)));
    regrid_with(static_cast<TensorType&>(result),
                static_cast<const TensorType&>(*this),
                weights.t_lags,
                weights.f_lags,
                weights.za_inc_lags);
    return result;
  }

  [[nodiscard]] AbsorptionVectorData regrid(const ScatteringDataGrids& grids) const {
    auto weights = calc_regrid_weights(t_grid_, f_grid_, nullptr, za_inc_grid_, nullptr, nullptr, grids);
    return regrid(grids, weights);
  }

  [[nodiscard]] AbsorptionVectorData<Format::ARO, Representation::Gridded> to_gridded() const {
    auto abs_new = AbsorptionVectorData<Format::ARO, Representation::Gridded>(t_grid_, f_grid_, za_inc_grid_);
    abs_new      = *this;
    return abs_new;
  }

 protected:
  /// The size of the temperature grid.
  Index n_temps_;
  /// The temperature grid.
  std::shared_ptr<const Vector> t_grid_;

  /// The size of the frequency grid.
  Index n_freqs_;
  /// The frequency grid.
  std::shared_ptr<const Vector> f_grid_;

  /// The size of the incoming zenith angle grid.
  Index n_za_inc_;
  /// The incoming zenith angle grid.
  std::shared_ptr<const Vector> za_inc_grid_;
};

}  // namespace scattering

#endif  // ABSORPTION_VECTOR_H_
